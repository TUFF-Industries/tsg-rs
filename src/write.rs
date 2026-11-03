//! CSV and `.npy` output. Finite floats are written in their shortest form
//! that parses back to the same float32, so CSV output is lossless for them;
//! NaN and infinities are written as `NaN`, `inf` and `-inf`.

use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::{ScalarColumn, ScalarValues, SpectralFile};

/// Shortest round-trip text for an `f32`: plain notation in `[1e-4, 1e7)`,
/// scientific outside it (Rust prints both forms shortest-round-trip).
pub fn fmt_f32(v: f32) -> String {
    let a = v.abs();
    if v == 0.0 || !v.is_finite() || (1e-4..1e7).contains(&a) {
        format!("{v}")
    } else {
        format!("{v:e}")
    }
}

/// Quote a CSV field when it contains a comma, quote or line break.
pub fn csv_field(s: &str) -> String {
    if s.contains([',', '"', '\n', '\r']) {
        format!("\"{}\"", s.replace('"', "\"\""))
    } else {
        s.to_string()
    }
}

/// Write to `<path>.partial` and rename into place, so a failed write never
/// leaves a truncated file under the requested name.
fn atomically(
    path: &Path,
    body: impl FnOnce(&mut BufWriter<File>) -> io::Result<()>,
) -> io::Result<()> {
    let mut tmp = path.as_os_str().to_owned();
    tmp.push(".partial");
    let tmp = PathBuf::from(tmp);
    let result = File::create(&tmp).and_then(|f| {
        let mut w = BufWriter::new(f);
        body(&mut w)?;
        w.flush()
    });
    let result = result.and_then(|()| std::fs::rename(&tmp, path));
    if result.is_err() {
        let _ = std::fs::remove_file(&tmp);
    }
    result
}

/// `.npy` header (format 1.0) for a C-order little-endian float32 array.
fn npy_header(shape: &[usize]) -> Vec<u8> {
    let dims: Vec<String> = shape.iter().map(ToString::to_string).collect();
    let shape_txt = match dims.as_slice() {
        [one] => format!("({one},)"),
        _ => format!("({})", dims.join(", ")),
    };
    let mut dict = format!("{{'descr': '<f4', 'fortran_order': False, 'shape': {shape_txt}, }}");
    // Magic (6) + version (2) + length (2) + dict + '\n' is a multiple of 64.
    while (10 + dict.len() + 1) % 64 != 0 {
        dict.push(' ');
    }
    dict.push('\n');
    let len = u16::try_from(dict.len()).expect("npy header for a 2-D shape fits in u16");
    let mut out = b"\x93NUMPY\x01\x00".to_vec();
    out.extend_from_slice(&len.to_le_bytes());
    out.extend_from_slice(dict.as_bytes());
    out
}

/// A depth as CSV text: shortest round-trip form, empty when absent.
fn depth_text(d: Option<f64>) -> String {
    d.map(|v| v.to_string()).unwrap_or_default()
}

/// Spectra as an `(n_samples, n_bands)` float32 `.npy`, values bit-for-bit.
pub fn spectra_npy(f: &SpectralFile, path: &Path) -> io::Result<()> {
    atomically(path, |w| {
        w.write_all(&npy_header(&[f.n_samples(), f.n_bands()]))?;
        f.spectra()
            .iter()
            .try_for_each(|v| w.write_all(&v.to_le_bytes()))
    })
}

/// Spectra as CSV: `sample,depth_m,<wavelength>...`, one row per sample.
pub fn spectra_csv(f: &SpectralFile, path: &Path) -> io::Result<()> {
    atomically(path, |w| {
        write!(w, "sample,depth_m")?;
        for wl in f.wavelengths() {
            write!(w, ",{wl}")?;
        }
        writeln!(w)?;
        for s in f.samples() {
            write!(w, "{},{}", s.index, depth_text(s.depth_m))?;
            for v in f.spectrum(s.index) {
                write!(w, ",{}", fmt_f32(*v))?;
            }
            writeln!(w)?;
        }
        Ok(())
    })
}

/// Columns written before the scalars in `scalars` CSV output.
const LEADING_COLUMNS: [&str; 2] = ["sample", "depth_m"];

/// CSV column names: the band-header names, with ` #<column>` appended to a
/// name that repeats or clashes with a leading column, until every header in
/// the file is unique.
pub fn scalar_column_names(cols: &[ScalarColumn<'_>]) -> Vec<String> {
    let mut count: BTreeMap<&str, usize> = LEADING_COLUMNS.iter().map(|n| (*n, 1)).collect();
    for c in cols {
        *count.entry(c.header.name.as_str()).or_default() += 1;
    }
    let mut taken: BTreeSet<String> = LEADING_COLUMNS.iter().map(ToString::to_string).collect();
    let mut names = Vec::with_capacity(cols.len());
    for c in cols {
        let mut name = c.header.name.clone();
        if count[name.as_str()] > 1 {
            name = format!("{name} #{}", c.header.band);
        }
        while taken.contains(&name) {
            name = format!("{name} #{}", c.header.band);
        }
        taken.insert(name.clone());
        names.push(name);
    }
    names
}

/// Scalars as CSV: `sample,depth_m,<scalar>...`. Class-coded columns hold
/// their labels; nulls are empty fields.
pub fn scalars_csv(f: &SpectralFile, cols: &[ScalarColumn<'_>], path: &Path) -> io::Result<()> {
    let names = scalar_column_names(cols);
    atomically(path, |w| {
        write!(w, "sample,depth_m")?;
        for n in &names {
            write!(w, ",{}", csv_field(n))?;
        }
        writeln!(w)?;
        for s in f.samples() {
            write!(w, "{},{}", s.index, depth_text(s.depth_m))?;
            for c in cols {
                let cell = match &c.values {
                    ScalarValues::Numeric(v) => v[s.index].map(fmt_f32).unwrap_or_default(),
                    ScalarValues::Class { labels, .. } => {
                        labels[s.index].map(csv_field).unwrap_or_default()
                    }
                };
                write!(w, ",{cell}")?;
            }
            writeln!(w)?;
        }
        Ok(())
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn float_text_round_trips_exactly() {
        for v in [
            0.0f32,
            -0.0,
            0.1,
            1.0 / 3.0,
            0.000_123,
            1e-30,
            123_456_780.0,
            crate::NULL_VALUE,
            f32::MIN_POSITIVE,
            7.0,
        ] {
            let s = fmt_f32(v);
            assert_eq!(
                s.parse::<f32>().unwrap().to_bits(),
                v.to_bits(),
                "{v} -> {s}"
            );
        }
        assert_eq!(fmt_f32(0.25), "0.25");
        assert_eq!(fmt_f32(crate::NULL_VALUE), "-3.4028233e38");
    }

    #[test]
    fn csv_fields_are_quoted_when_needed() {
        assert_eq!(csv_field("Depth (m)"), "Depth (m)");
        assert_eq!(csv_field("a,b"), "\"a,b\"");
        assert_eq!(csv_field("say \"x\""), "\"say \"\"x\"\"\"");
    }

    #[test]
    fn npy_header_is_aligned() {
        let h = npy_header(&[22356, 531]);
        assert_eq!(h.len() % 64, 0);
        assert!(std::str::from_utf8(&h[10..])
            .unwrap()
            .contains("'shape': (22356, 531)"));
    }
}
