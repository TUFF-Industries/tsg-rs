//! The scalar table: `[band headers]` describe the columns of the second
//! `.bip` plane, and `[class N]` sections hold the labels of class-coded
//! columns (TSA mineral and group results among them).

use std::collections::{BTreeMap, BTreeSet};

use crate::{fmt_err, TsgError};

/// TSG's missing-value marker in the scalar plane: the float32 bit pattern
/// `0xFF7FFFFE`, one step above `f32::MIN`.
pub const NULL_VALUE: f32 = f32::from_bits(0xFF7F_FFFE);

/// True for TSG's missing-value marker. Uses the same test as pytsg
/// (`numpy.isclose(value, finfo(float32).min)`, default tolerances), which
/// accepts both `NULL_VALUE` and `f32::MIN`.
pub fn is_null(v: f32) -> bool {
    let min = f64::from(f32::MIN);
    (f64::from(v) - min).abs() <= 1e-8 + 1e-5 * min.abs()
}

/// One `[band headers]` record, e.g. `91:Grp1 sTSAS;1.000000;16.000000;2;1;...`.
#[derive(Debug, Clone, PartialEq)]
pub struct BandHeader {
    /// Column of the scalar plane.
    pub band: usize,
    pub name: String,
    /// Field 3 (`2` = plain or class-coded scalar); `None` for name-only records.
    pub flag: Option<i64>,
    /// Field 4 when `flag <= 2`: the `[class N]` section number, `-1` for none.
    pub class_number: Option<i64>,
}

impl BandHeader {
    /// The `[class N]` number whose labels this column's values index, if any.
    /// As in pytsg, only `flag == 2` columns with a class number are coded.
    pub fn class_section(&self) -> Option<u32> {
        match (self.flag, self.class_number) {
            (Some(2), Some(n)) => u32::try_from(n).ok(),
            _ => None,
        }
    }
}

pub(crate) fn parse_band_header(line: &str) -> Result<BandHeader, TsgError> {
    let (idx, rest) = line
        .split_once(':')
        .ok_or_else(|| fmt_err(format!("band header without `:`: {line}")))?;
    let band = idx
        .trim()
        .parse()
        .map_err(|_| fmt_err(format!("bad band index in band header: {line}")))?;
    let parts: Vec<&str> = rest.split(';').collect();
    let name = parts[0].to_string();
    if parts.len() == 1 {
        return Ok(BandHeader {
            band,
            name,
            flag: None,
            class_number: None,
        });
    }
    let field = |i: usize, what: &str| -> Result<i64, TsgError> {
        parts
            .get(i)
            .and_then(|s| s.trim().parse::<i64>().ok())
            .ok_or_else(|| {
                fmt_err(format!(
                    "band header {band} (`{name}`) has no integer {what}"
                ))
            })
    };
    let flag = field(3, "flag in field 3")?;
    let class_number = if flag <= 2 {
        Some(field(4, "class number in field 4")?)
    } else {
        None
    };
    if let (2, Some(c)) = (flag, class_number) {
        if c != -1 && u32::try_from(c).is_err() {
            return Err(fmt_err(format!(
                "band header {band} (`{name}`) has class number {c}; expected -1 or a [class N] number"
            )));
        }
    }
    Ok(BandHeader {
        band,
        name,
        flag: Some(flag),
        class_number,
    })
}

/// A `[class N]` section: a name and the labels `0:LABEL`, `1:LABEL`, ...
#[derive(Debug, Clone, PartialEq)]
pub struct ClassTable {
    pub number: u32,
    pub name: String,
    pub labels: Vec<String>,
}

fn is_label_line(line: &str) -> bool {
    line.split_once(':')
        .is_some_and(|(i, _)| !i.is_empty() && i.bytes().all(|b| b.is_ascii_digit()))
}

pub(crate) fn parse_class_table(number: u32, lines: &[&str]) -> Result<ClassTable, TsgError> {
    let sec = format!("[class {number}]");
    let mut name = None;
    let mut max = None;
    let mut labels = Vec::new();
    for line in lines {
        if is_label_line(line) {
            let (i, label) = line.split_once(':').unwrap_or_default();
            if i.parse::<usize>().ok() != Some(labels.len()) {
                return Err(fmt_err(format!(
                    "{sec} label `{line}` is out of order (expected index {})",
                    labels.len()
                )));
            }
            if label.is_empty() {
                return Err(fmt_err(format!(
                    "{sec} label {} is empty, which CSV output could not tell from a null",
                    labels.len()
                )));
            }
            labels.push(label.to_string());
        } else if let Some((k, v)) = line.split_once('=') {
            match k.trim() {
                "name" => name = Some(v.trim().to_string()),
                "max" => max = Some(v.trim().to_string()),
                _ => {}
            }
        } else {
            return Err(fmt_err(format!("{sec} has an unrecognised line: {line}")));
        }
    }
    let name = name.ok_or_else(|| fmt_err(format!("{sec} has no `name =`")))?;
    let max = max.ok_or_else(|| fmt_err(format!("{sec} has no `max =`")))?;
    if max.parse::<usize>().ok() != Some(labels.len()) {
        return Err(fmt_err(format!(
            "{sec} says max = {max} but lists {} labels",
            labels.len()
        )));
    }
    Ok(ClassTable {
        number,
        name,
        labels,
    })
}

/// Parse every `[class N]` section; a malformed `N` or a repeated one is refused.
pub(crate) fn parse_class_tables<'a>(
    sections: impl Iterator<Item = (&'a str, Vec<&'a str>)>,
) -> Result<BTreeMap<u32, ClassTable>, TsgError> {
    let mut out = BTreeMap::new();
    for (name, lines) in sections {
        let Some(num) = name.strip_prefix("class ") else {
            continue;
        };
        let number: u32 = num
            .trim()
            .parse()
            .map_err(|_| fmt_err(format!("section [{name}] has a non-numeric class number")))?;
        if out
            .insert(number, parse_class_table(number, &lines)?)
            .is_some()
        {
            return Err(fmt_err(format!("[class {number}] appears more than once")));
        }
    }
    Ok(out)
}

/// Band headers must name distinct columns inside the scalar plane and every
/// class-coded column must point at a class section that exists.
pub(crate) fn check_band_headers(
    headers: &[BandHeader],
    n_bands: usize,
    classes: &BTreeMap<u32, ClassTable>,
) -> Result<(), TsgError> {
    let mut seen = BTreeSet::new();
    for h in headers {
        if h.band >= n_bands {
            return Err(fmt_err(format!(
                "band header {} (`{}`) is outside the {n_bands} scalar columns",
                h.band, h.name
            )));
        }
        if !seen.insert(h.band) {
            return Err(fmt_err(format!(
                "band header {} appears more than once",
                h.band
            )));
        }
        if let Some(c) = h.class_section() {
            if !classes.contains_key(&c) {
                return Err(fmt_err(format!(
                    "scalar `{}` (column {}) refers to [class {c}], which is missing",
                    h.name, h.band
                )));
            }
        }
    }
    Ok(())
}

/// Values of one scalar column, nulls as `None`.
#[derive(Debug, Clone, PartialEq)]
pub enum ScalarValues<'a> {
    Numeric(Vec<Option<f32>>),
    /// Class-coded: label strings borrowed from the class table.
    Class {
        table: &'a ClassTable,
        labels: Vec<Option<&'a str>>,
    },
}

/// A decoded scalar column.
#[derive(Debug, Clone, PartialEq)]
pub struct ScalarColumn<'a> {
    pub header: &'a BandHeader,
    pub values: ScalarValues<'a>,
}

/// A class value must be a whole number indexing the table's labels.
pub(crate) fn decode_class<'a>(
    raw: &[f32],
    header: &BandHeader,
    table: &'a ClassTable,
) -> Result<Vec<Option<&'a str>>, TsgError> {
    raw.iter()
        .enumerate()
        .map(|(sample, &v)| {
            if is_null(v) {
                return Ok(None);
            }
            let valid = v.fract() == 0.0 && v >= 0.0 && (v as usize) < table.labels.len();
            if !valid {
                return Err(fmt_err(format!(
                    "scalar `{}` (column {}) sample {sample} holds {v}, which is not a label \
                     index of [class {}] `{}` ({} labels)",
                    header.name,
                    header.band,
                    table.number,
                    table.name,
                    table.labels.len()
                )));
            }
            Ok(Some(table.labels[v as usize].as_str()))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn null_marker_matches_pytsg_tolerance() {
        assert!(is_null(NULL_VALUE));
        assert!(is_null(f32::MIN));
        assert!(!is_null(-3.0e38));
        assert!(!is_null(0.0));
        assert!(!is_null(f32::NEG_INFINITY));
    }

    #[test]
    fn parses_band_header_kinds() {
        let coded = parse_band_header("91:Grp1 sTSAS;1.0;16.0;2;1;1;1;0").unwrap();
        assert_eq!(coded.class_section(), Some(1));
        let plain = parse_band_header("80:Depth (m);59.5;201.5;2;-1;0").unwrap();
        assert_eq!(plain.class_section(), None);
        let other = parse_band_header("82:sec_end_mask;0;1;9;5;-1").unwrap();
        assert_eq!((other.flag, other.class_number), (Some(9), None));
        let bare = parse_band_header("0:Flags").unwrap();
        assert_eq!((bare.flag, bare.class_section()), (None, None));
        for bad in ["4:G;0;1;2;-5", "4:G;0;1;2;4294967296"] {
            let err = parse_band_header(bad).unwrap_err().to_string();
            assert!(
                err.contains("band header 4 (`G`) has class number"),
                "{err}"
            );
        }
        let err = parse_band_header("5:X;0;1;two;1").unwrap_err().to_string();
        assert!(
            err.contains("band header 5 (`X`) has no integer flag"),
            "{err}"
        );
    }

    #[test]
    fn class_table_checks_order_and_max() {
        let t = parse_class_table(1, &["name = G", "max = 2", "colours = 1 2", "0:A", "1:B:C"])
            .unwrap();
        assert_eq!(t.labels, ["A", "B:C"]);
        let err = parse_class_table(1, &["name = G", "max = 3", "0:A", "1:B"]).unwrap_err();
        assert!(
            err.to_string().contains("max = 3 but lists 2 labels"),
            "{err}"
        );
        let err = parse_class_table(1, &["name = G", "max = 2", "1:A", "0:B"]).unwrap_err();
        assert!(err.to_string().contains("out of order"), "{err}");
        assert!(parse_class_table(1, &["max = 0"]).is_err());
        let err = parse_class_table(1, &["name = G", "max = 1", "0:"]).unwrap_err();
        assert!(err.to_string().contains("label 0 is empty"), "{err}");
    }

    #[test]
    fn decodes_class_values_and_refuses_bad_ones() {
        let t = parse_class_table(1, &["name = G", "max = 2", "0:A", "1:B"]).unwrap();
        let h = parse_band_header("3:Grp1 sTSAS;0;1;2;1").unwrap();
        let got = decode_class(&[1.0, NULL_VALUE, 0.0], &h, &t).unwrap();
        assert_eq!(got, [Some("B"), None, Some("A")]);
        for bad in [2.0, 0.5, -1.0, f32::NAN] {
            let err = decode_class(&[bad], &h, &t).unwrap_err().to_string();
            assert!(err.contains("`Grp1 sTSAS` (column 3) sample 0"), "{err}");
        }
    }

    #[test]
    fn band_header_checks() {
        let classes =
            parse_class_tables([("class 1", vec!["name = G", "max = 1", "0:A"])].into_iter())
                .unwrap();
        let h = |l: &str| parse_band_header(l).unwrap();
        assert!(check_band_headers(&[h("0:A;0;1;2;1")], 2, &classes).is_ok());
        let err = check_band_headers(&[h("2:A")], 2, &classes).unwrap_err();
        assert!(
            err.to_string().contains("outside the 2 scalar columns"),
            "{err}"
        );
        let err = check_band_headers(&[h("0:A"), h("0:B")], 2, &classes).unwrap_err();
        assert!(err.to_string().contains("appears more than once"), "{err}");
        let err = check_band_headers(&[h("0:A;0;1;2;7")], 2, &classes).unwrap_err();
        assert!(
            err.to_string().contains("[class 7], which is missing"),
            "{err}"
        );
    }
}
