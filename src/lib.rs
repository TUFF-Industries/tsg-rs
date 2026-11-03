//! Reader for the NIR/SWIR part of a TSG (The Spectral Geologist) package:
//! a `.tsg` text header plus a `.bip` binary file holding the reflectance
//! spectra and the scalar table, including any TSA results stored in it.
//!
//! The layout was learned from pytsg (Geological Survey of Western Australia,
//! MIT) and checked against real NVCL HyLogger packages; see
//! `docs/tsg-format.md`. Files that fail a consistency check are refused with
//! a specific message rather than read on a guess.
//!
//! ```no_run
//! let f = tsg_rs::read_pair("hole_tsg.tsg".as_ref(), "hole_tsg.bip".as_ref())?;
//! println!("{} spectra of {} bands", f.n_samples(), f.n_bands());
//! for col in f.scalar_columns()? {
//!     println!("{}", col.header.name);
//! }
//! # Ok::<(), tsg_rs::TsgError>(())
//! ```

#![forbid(unsafe_code)]

mod header;
mod scalars;
pub mod write;

use std::collections::{BTreeMap, BTreeSet};
use std::fmt;
use std::path::{Path, PathBuf};

pub use header::{decode_cp1252, Header};
pub use scalars::{is_null, BandHeader, ClassTable, ScalarColumn, ScalarValues, NULL_VALUE};

/// Backward depth steps up to this size (metres) below the running maximum
/// are ignored by [`SpectralFile::depth_backsteps`]: real files jitter by a
/// few micrometres within a scan line.
pub const DEPTH_TOLERANCE_M: f64 = 0.001;

/// Wavelength grids outside this range (nm) are refused: only the VNIR/SWIR
/// pair has been checked, and TIR files are not supported.
pub const WAVELENGTH_RANGE_NM: (f64, f64) = (300.0, 2600.0);

#[derive(Debug)]
pub enum TsgError {
    Io(PathBuf, std::io::Error),
    Format(String),
}

impl fmt::Display for TsgError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(p, e) => write!(f, "{}: {e}", p.display()),
            Self::Format(m) => write!(f, "refusing TSG file: {m}"),
        }
    }
}

impl std::error::Error for TsgError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(_, e) => Some(e),
            Self::Format(_) => None,
        }
    }
}

pub(crate) fn fmt_err(msg: impl Into<String>) -> TsgError {
    TsgError::Format(msg.into())
}

/// One `[sample headers]` record, e.g.
/// `0:05KCD001_0001_1  T=0001 L=1 P=1 D=59.500008 X=6.263977 H=05KCD001`.
#[derive(Debug, Clone, PartialEq)]
pub struct SampleHeader {
    pub index: usize,
    pub name: String,
    /// The `D=` field; `None` when the file has no depths at all.
    pub depth_m: Option<f64>,
    pub fields: BTreeMap<String, String>,
}

fn parse_sample_header(line: &str) -> Result<SampleHeader, TsgError> {
    let (idx, rest) = line
        .split_once(':')
        .ok_or_else(|| fmt_err(format!("sample header without `:`: {line}")))?;
    let index = idx
        .trim()
        .parse()
        .map_err(|_| fmt_err(format!("bad sample index in sample header: {line}")))?;
    // The first token is always part of the name (real names can contain `=`);
    // later tokens up to the first `key=value` extend it. From there on every
    // token must be a `key=value` field whose key has not been seen before.
    let rest = rest.trim();
    let tokens: Vec<(usize, &str)> = rest
        .split_whitespace()
        .map(|t| (t.as_ptr() as usize - rest.as_ptr() as usize, t))
        .collect();
    let split = tokens
        .iter()
        .skip(1)
        .position(|(_, t)| t.contains('='))
        .map_or(tokens.len(), |p| p + 1);
    let name = match tokens.get(split) {
        Some((offset, _)) => rest[..*offset].trim_end(),
        None => rest,
    }
    .to_string();
    let mut fields = BTreeMap::new();
    for (_, t) in &tokens[split..] {
        let (k, v) = t.split_once('=').ok_or_else(|| {
            fmt_err(format!(
                "sample {index}: `{t}` follows the key=value fields but is not one"
            ))
        })?;
        if fields.insert(k.to_string(), v.to_string()).is_some() {
            return Err(fmt_err(format!("sample {index} has `{k}=` twice")));
        }
    }
    let depth_m = match fields.get("D") {
        None => None,
        Some(d) => match d.parse::<f64>() {
            Ok(v) if v.is_finite() => Some(v),
            _ => {
                return Err(fmt_err(format!(
                    "sample {index} depth D={d} is not a finite number"
                )))
            }
        },
    };
    Ok(SampleHeader {
        index,
        name,
        depth_m,
        fields,
    })
}

/// Evenly spaced band centres `start + i * (end - start) / (n - 1)` from the
/// single `start end unit` line of `[wavelength specs]`.
fn wavelengths_from_specs(lines: &[&str], n_bands: usize) -> Result<(Vec<f64>, String), TsgError> {
    let [line] = lines else {
        return Err(fmt_err(format!(
            "[wavelength specs] should be one line, found {}",
            lines.len()
        )));
    };
    let parts: Vec<&str> = line.split_whitespace().collect();
    let [s, e, unit] = parts.as_slice() else {
        return Err(fmt_err(format!(
            "[wavelength specs] should be `start end unit`, found `{line}`"
        )));
    };
    let parse = |v: &str, what: &str| {
        v.parse::<f64>()
            .map_err(|_| fmt_err(format!("[wavelength specs] {what} `{v}` is not a number")))
    };
    let (start, end) = (parse(s, "start")?, parse(e, "end")?);
    if n_bands < 2 {
        return Err(fmt_err(format!(
            "{n_bands} band(s); at least 2 are needed to place a wavelength grid"
        )));
    }
    if !start.is_finite() || !end.is_finite() || end <= start {
        return Err(fmt_err(format!(
            "wavelengths not increasing: start={start} end={end}"
        )));
    }
    if *unit != "nm" {
        return Err(fmt_err(format!("wavelength unit `{unit}` is not nm")));
    }
    let (lo, hi) = WAVELENGTH_RANGE_NM;
    if start < lo || end > hi {
        return Err(fmt_err(format!(
            "wavelength grid {start}-{end} nm is outside {lo}-{hi} nm; \
             only the VNIR/SWIR pair is supported (TIR is not)"
        )));
    }
    let step = (end - start) / (n_bands - 1) as f64;
    let wl = (0..n_bands).map(|i| start + i as f64 * step).collect();
    Ok((wl, (*unit).to_string()))
}

/// Either every sample header has a depth or none does; a mix is refused.
fn check_depths(samples: &[SampleHeader]) -> Result<(), TsgError> {
    let with = samples.iter().filter(|s| s.depth_m.is_some()).count();
    if with != 0 && with != samples.len() {
        let missing = samples
            .iter()
            .find(|s| s.depth_m.is_none())
            .map_or(0, |s| s.index);
        return Err(fmt_err(format!(
            "{with} of {} sample headers have a depth (D=) but sample {missing} has none",
            samples.len()
        )));
    }
    Ok(())
}

/// Everything the `.tsg` header says, checked, before any `.bip` data is read.
#[derive(Debug, Clone)]
struct Layout {
    header: Header,
    n_samples: usize,
    n_bands: usize,
    wavelengths: Vec<f64>,
    wavelength_unit: String,
    samples: Vec<SampleHeader>,
    band_headers: Vec<BandHeader>,
    classes: BTreeMap<u32, ClassTable>,
    bip_len: usize,
}

/// A decoded `.tsg`/`.bip` pair.
#[derive(Debug, Clone)]
pub struct SpectralFile {
    layout: Layout,
    spectra: Vec<f32>,
    scalars: Vec<f32>,
}

impl SpectralFile {
    /// The parsed `.tsg` header, for sections this crate does not interpret.
    pub fn header(&self) -> &Header {
        &self.layout.header
    }

    pub fn n_samples(&self) -> usize {
        self.layout.n_samples
    }

    pub fn n_bands(&self) -> usize {
        self.layout.n_bands
    }

    /// Band centres, one per band.
    pub fn wavelengths(&self) -> &[f64] {
        &self.layout.wavelengths
    }

    pub fn wavelength_unit(&self) -> &str {
        &self.layout.wavelength_unit
    }

    pub fn samples(&self) -> &[SampleHeader] {
        &self.layout.samples
    }

    /// `[band headers]` records in file order.
    pub fn band_headers(&self) -> &[BandHeader] {
        &self.layout.band_headers
    }

    /// `[class N]` tables by `N`.
    pub fn classes(&self) -> &BTreeMap<u32, ClassTable> {
        &self.layout.classes
    }

    /// Reflectance spectrum of sample `i` (length `n_bands`).
    pub fn spectrum(&self, i: usize) -> &[f32] {
        &self.spectra[i * self.layout.n_bands..(i + 1) * self.layout.n_bands]
    }

    /// All spectra, row-major `[sample][band]`.
    pub fn spectra(&self) -> &[f32] {
        &self.spectra
    }

    /// Raw values of scalar-plane column `band`, nulls included.
    pub fn scalar_column(&self, band: usize) -> Option<Vec<f32>> {
        (band < self.layout.n_bands).then(|| {
            (0..self.layout.n_samples)
                .map(|i| self.scalars[i * self.layout.n_bands + band])
                .collect()
        })
    }

    /// Depths in sample order; all `None` when the file has no `D=` fields.
    pub fn depths(&self) -> Vec<Option<f64>> {
        self.layout.samples.iter().map(|s| s.depth_m).collect()
    }

    /// Samples whose depth is more than [`DEPTH_TOLERANCE_M`] below the
    /// running maximum, i.e. where depth goes backwards. Reported, not refused:
    /// the values are read as stored.
    pub fn depth_backsteps(&self) -> Vec<usize> {
        let mut max_so_far = f64::NEG_INFINITY;
        let mut out = Vec::new();
        for s in &self.layout.samples {
            let Some(d) = s.depth_m else { continue };
            if max_so_far - d > DEPTH_TOLERANCE_M {
                out.push(s.index);
            }
            max_so_far = max_so_far.max(d);
        }
        out
    }

    /// Every column named in `[band headers]`, in column order, with class
    /// codes replaced by their labels. Refuses the file if any class value is
    /// not a whole-number index into its label table.
    pub fn scalar_columns(&self) -> Result<Vec<ScalarColumn<'_>>, TsgError> {
        let mut headers: Vec<&BandHeader> = self.layout.band_headers.iter().collect();
        headers.sort_by_key(|h| h.band);
        headers
            .into_iter()
            .map(|h| {
                let raw = self.scalar_column(h.band).unwrap_or_default();
                let values = match h.class_section().and_then(|c| self.layout.classes.get(&c)) {
                    Some(table) => ScalarValues::Class {
                        table,
                        labels: scalars::decode_class(&raw, h, table)?,
                    },
                    None => ScalarValues::Numeric(
                        raw.into_iter()
                            .map(|v| (!is_null(v)).then_some(v))
                            .collect(),
                    ),
                };
                Ok(ScalarColumn { header: h, values })
            })
            .collect()
    }

    /// Class tables that hold TSA results (name contains `TSA`), with the
    /// scalar columns coded by each.
    pub fn tsa_tables(&self) -> Vec<(&ClassTable, Vec<&BandHeader>)> {
        self.layout
            .classes
            .values()
            .filter(|t| t.name.contains("TSA"))
            .map(|t| {
                let mut cols: Vec<&BandHeader> = self
                    .layout
                    .band_headers
                    .iter()
                    .filter(|h| h.class_section() == Some(t.number))
                    .collect();
                cols.sort_by_key(|h| h.band);
                (t, cols)
            })
            .collect()
    }
}

fn parse_usize(h: &Header, section: &str, key: &str) -> Result<usize, TsgError> {
    let v = h.require(section, key)?;
    v.parse().map_err(|_| {
        fmt_err(format!(
            "[{section}] {key} = `{v}` is not a non-negative integer"
        ))
    })
}

fn parse_samples(header: &Header, n_samples: usize) -> Result<Vec<SampleHeader>, TsgError> {
    let samples = header
        .section("sample headers")?
        .ok_or_else(|| fmt_err("missing [sample headers]"))?
        .into_iter()
        .map(parse_sample_header)
        .collect::<Result<Vec<_>, _>>()?;
    if samples.len() != n_samples {
        return Err(fmt_err(format!(
            "{} sample headers but [coordinates] lastsample = {n_samples}",
            samples.len()
        )));
    }
    if let Some(s) = samples.iter().enumerate().find(|(i, s)| s.index != *i) {
        return Err(fmt_err(format!(
            "sample header {} is at position {}; indices must be 0..n in order",
            s.1.index, s.0
        )));
    }
    check_depths(&samples)?;
    Ok(samples)
}

/// A `key = value` in `section` that, when present, must equal `want`.
fn check_optional(h: &Header, section: &str, key: &str, want: &str) -> Result<(), TsgError> {
    match h.value(section, key)? {
        Some(v) if v != want => Err(fmt_err(format!(
            "[{section}] {key} = {v}; only {want} is supported"
        ))),
        _ => Ok(()),
    }
}

/// Parse and check the `.tsg` header on its own.
fn parse_layout(header_bytes: &[u8]) -> Result<Layout, TsgError> {
    if header_bytes.starts_with(b"\xEF\xBB\xBF") {
        return Err(fmt_err(
            "header starts with a UTF-8 byte-order mark; TSG headers are cp1252 text",
        ));
    }
    let header = Header::parse(&decode_cp1252(header_bytes)?);
    let datatype = header.require("structure", "datatype")?;
    if datatype != "float" {
        return Err(fmt_err(format!(
            "unsupported [structure] datatype `{datatype}`"
        )));
    }
    check_optional(&header, "structure", "dataheader", "0")?;
    check_optional(&header, "platform", "byteorder", "3")?;
    check_optional(&header, "coordinates", "firstband", "1")?;
    check_optional(&header, "coordinates", "firstsample", "1")?;
    let n_bands = parse_usize(&header, "coordinates", "lastband")?;
    let n_samples = parse_usize(&header, "coordinates", "lastsample")?;
    check_optional(&header, "dimensions", "bands", &n_bands.to_string())?;
    check_optional(&header, "dimensions", "samples", &n_samples.to_string())?;
    if n_bands == 0 || n_samples == 0 {
        return Err(fmt_err("zero bands or zero samples"));
    }
    let bip_len = n_samples
        .checked_mul(n_bands)
        .and_then(|v| v.checked_mul(2 * 4))
        .ok_or_else(|| fmt_err("sample x band count overflows"))?;

    let specs = header
        .section("wavelength specs")?
        .ok_or_else(|| fmt_err("missing [wavelength specs]"))?;
    let (wavelengths, wavelength_unit) = wavelengths_from_specs(&specs, n_bands)?;
    let samples = parse_samples(&header, n_samples)?;

    let classes = scalars::parse_class_tables(header.sections())?;
    let band_headers = header
        .section("band headers")?
        .unwrap_or_default()
        .into_iter()
        .map(scalars::parse_band_header)
        .collect::<Result<Vec<_>, _>>()?;
    scalars::check_band_headers(&band_headers, n_bands, &classes)?;
    Ok(Layout {
        header,
        n_samples,
        n_bands,
        wavelengths,
        wavelength_unit,
        samples,
        band_headers,
        classes,
        bip_len,
    })
}

impl Layout {
    fn check_bip_len(&self, len: u64) -> Result<(), TsgError> {
        if u64::try_from(self.bip_len).ok() != Some(len) {
            return Err(fmt_err(format!(
                ".bip is {len} bytes, expected {} (2 planes x {} samples x {} bands x 4)",
                self.bip_len, self.n_samples, self.n_bands
            )));
        }
        Ok(())
    }

    /// Split checked `.bip` bytes into the spectra and scalar planes.
    fn with_data(self, bip: &[u8]) -> Result<SpectralFile, TsgError> {
        self.check_bip_len(bip.len() as u64)?;
        let to_f32 = |b: &[u8]| -> Vec<f32> {
            b.chunks_exact(4)
                .map(|c| f32::from_le_bytes([c[0], c[1], c[2], c[3]]))
                .collect()
        };
        let (spectra, scalars) = bip.split_at(self.bip_len / 2);
        Ok(SpectralFile {
            spectra: to_f32(spectra),
            scalars: to_f32(scalars),
            layout: self,
        })
    }
}

/// Parse a header and `.bip` bytes that are already in memory.
pub fn decode(header_bytes: &[u8], bip: &[u8]) -> Result<SpectralFile, TsgError> {
    parse_layout(header_bytes)?.with_data(bip)
}

/// Read a `.tsg`/`.bip` pair from disk. The header is checked, and the `.bip`
/// size compared with it, before the `.bip` is read.
pub fn read_pair(tsg: &Path, bip: &Path) -> Result<SpectralFile, TsgError> {
    let layout = parse_layout(&std::fs::read(tsg).map_err(io_err(tsg))?)?;
    layout.check_bip_len(std::fs::metadata(bip).map_err(io_err(bip))?.len())?;
    layout.with_data(&std::fs::read(bip).map_err(io_err(bip))?)
}

fn io_err(p: &Path) -> impl FnOnce(std::io::Error) -> TsgError + '_ {
    move |e| TsgError::Io(p.to_path_buf(), e)
}

fn has_ext(p: &Path, ext: &str) -> bool {
    p.extension()
        .and_then(|e| e.to_str())
        .is_some_and(|e| e.eq_ignore_ascii_case(ext))
}

/// The partner file with extension `ext` (either case), which must exist.
fn sibling(p: &Path, ext: &str) -> Result<PathBuf, TsgError> {
    [ext.to_lowercase(), ext.to_uppercase()]
        .iter()
        .map(|e| p.with_extension(e))
        .find(|c| c.is_file())
        .ok_or_else(|| fmt_err(format!("no .{ext} file next to {}", p.display())))
}

/// Resolve the `.tsg`/`.bip` pair from either file of the pair, or from a
/// directory holding exactly one `*tsg.tsg` (the NIR/SWIR header).
pub fn resolve_pair(path: &Path) -> Result<(PathBuf, PathBuf), TsgError> {
    if path.is_dir() {
        return find_nir_pair(path);
    }
    if !path.is_file() {
        return Err(TsgError::Io(
            path.to_path_buf(),
            std::io::Error::new(std::io::ErrorKind::NotFound, "no such file"),
        ));
    }
    if has_ext(path, "tsg") {
        Ok((path.to_path_buf(), sibling(path, "bip")?))
    } else if has_ext(path, "bip") {
        Ok((sibling(path, "tsg")?, path.to_path_buf()))
    } else {
        Err(fmt_err(format!(
            "{} is neither a .tsg nor a .bip file",
            path.display()
        )))
    }
}

/// Find the single NIR/SWIR pair (`*tsg.tsg` + `*tsg.bip`) in a directory.
pub fn find_nir_pair(dir: &Path) -> Result<(PathBuf, PathBuf), TsgError> {
    let entries = std::fs::read_dir(dir).map_err(|e| TsgError::Io(dir.to_path_buf(), e))?;
    let mut tsg = BTreeSet::new();
    for entry in entries {
        let path = entry
            .map_err(|e| TsgError::Io(dir.to_path_buf(), e))?
            .path();
        let name = path
            .file_name()
            .and_then(|n| n.to_str())
            .unwrap_or_default()
            .to_lowercase();
        if name.ends_with("tsg.tsg") {
            tsg.insert(path);
        }
    }
    let mut it = tsg.into_iter();
    match (it.next(), it.next()) {
        (Some(one), None) => Ok((one.clone(), sibling(&one, "bip")?)),
        (None, _) => Err(fmt_err(format!("no *tsg.tsg in {}", dir.display()))),
        (Some(_), Some(_)) => Err(fmt_err(format!(
            "more than one *tsg.tsg in {}; name the file",
            dir.display()
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn header(samples: &[(usize, f64)], bands: usize) -> String {
        let sh: String = samples
            .iter()
            .map(|(i, d)| format!("{i}:H_{i}  T=1 L=1 P={i} D={d} H=H\n"))
            .collect();
        format!(
            "[structure]\ndatatype = float\ndataheader = 0\n\n[coordinates]\nlastsample = {}\n\
             lastband = {bands}\n\n[band headers]\n0:Flags\n1:Grp;0;1;2;3;x\n\n\
             [sample headers]\n{sh}\n[class 3]\nname = TSA Groups\nmax = 2\n0:A\n1:B\n\n\
             [wavelength specs]\n1000.0 1003.0 nm\n",
            samples.len()
        )
    }

    fn bip(values: &[f32]) -> Vec<u8> {
        values.iter().flat_map(|v| v.to_le_bytes()).collect()
    }

    #[test]
    fn decodes_two_planes_wavelengths_and_classes() {
        let h = header(&[(0, 1.0), (1, 1.01)], 4);
        let mut vals: Vec<f32> = (0..16).map(|v| v as f32).collect();
        vals[9] = 1.0;
        vals[13] = NULL_VALUE;
        let f = decode(h.as_bytes(), &bip(&vals)).unwrap();
        assert_eq!(f.wavelengths(), [1000.0, 1001.0, 1002.0, 1003.0]);
        assert_eq!(f.spectrum(1), [4.0, 5.0, 6.0, 7.0]);
        assert_eq!(f.depths(), [Some(1.0), Some(1.01)]);
        let cols = f.scalar_columns().unwrap();
        assert_eq!(
            cols[0].values,
            ScalarValues::Numeric(vec![Some(8.0), Some(12.0)])
        );
        match &cols[1].values {
            ScalarValues::Class { labels, table } => {
                assert_eq!(labels, &[Some("B"), None]);
                assert_eq!(table.name, "TSA Groups");
            }
            other => panic!("expected class column, got {other:?}"),
        }
        assert_eq!(f.tsa_tables()[0].1[0].name, "Grp");
    }

    #[test]
    fn refuses_bip_with_wrong_size() {
        let h = header(&[(0, 1.0)], 4);
        let err = decode(h.as_bytes(), &bip(&[0.0; 7])).unwrap_err();
        assert!(err.to_string().contains("expected 32"), "{err}");
    }

    #[test]
    fn reports_backward_depths_beyond_tolerance_and_allows_missing_depths() {
        let h = header(&[(0, 5.0), (1, 4.9999), (2, 4.9), (3, 6.0)], 4);
        let f = decode(h.as_bytes(), &bip(&[0.0; 32])).unwrap();
        assert_eq!(f.depth_backsteps(), [2]);
        let none = header(&[(0, 5.0), (1, 6.0)], 4)
            .replace(" D=5", "")
            .replace(" D=6", "");
        let f = decode(none.as_bytes(), &bip(&[0.0; 16])).unwrap();
        assert_eq!(f.depths(), [None, None]);
        let mixed = header(&[(0, 5.0), (1, 6.0)], 4).replace(" D=6", "");
        let err = decode(mixed.as_bytes(), &bip(&[0.0; 16])).unwrap_err();
        assert!(err.to_string().contains("but sample 1 has none"), "{err}");
        let bad = header(&[(0, 5.0)], 4).replace("D=5", "D=five");
        let err = decode(bad.as_bytes(), &bip(&[0.0; 8])).unwrap_err();
        assert!(
            err.to_string().contains("D=five is not a finite number"),
            "{err}"
        );
    }

    #[test]
    fn refuses_bad_grids_and_counts() {
        let h = header(&[(0, 1.0)], 4).replace("1000.0 1003.0", "1003.0 1000.0");
        assert!(decode(h.as_bytes(), &bip(&[0.0; 8])).is_err());
        let h = header(&[(0, 1.0)], 4).replace("1000.0 1003.0 nm", "6000.0 14500.0 nm");
        let err = decode(h.as_bytes(), &bip(&[0.0; 8])).unwrap_err();
        assert!(err.to_string().contains("TIR is not"), "{err}");
        let h = header(&[(0, 1.0)], 4).replace("lastsample = 1", "lastsample = 2");
        let err = decode(h.as_bytes(), &bip(&[0.0; 16])).unwrap_err();
        assert!(err.to_string().contains("1 sample headers but"), "{err}");
        let h = header(&[(0, 1.0)], 4).replace("datatype = float", "datatype = int");
        assert!(decode(h.as_bytes(), &bip(&[0.0; 8])).is_err());
    }
}
