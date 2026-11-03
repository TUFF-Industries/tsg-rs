//! Parsing of the `.tsg` text header: cp1252 text split into `[section]` blocks.

use crate::{fmt_err, TsgError};

/// Bytes 0x80..=0x9F of Windows-1252. `None` marks the five undefined code
/// points, which Python's `cp1252` codec (and therefore pytsg) also rejects.
const CP1252_HIGH: [Option<char>; 32] = [
    Some('\u{20AC}'),
    None,
    Some('\u{201A}'),
    Some('\u{0192}'),
    Some('\u{201E}'),
    Some('\u{2026}'),
    Some('\u{2020}'),
    Some('\u{2021}'),
    Some('\u{02C6}'),
    Some('\u{2030}'),
    Some('\u{0160}'),
    Some('\u{2039}'),
    Some('\u{0152}'),
    None,
    Some('\u{017D}'),
    None,
    None,
    Some('\u{2018}'),
    Some('\u{2019}'),
    Some('\u{201C}'),
    Some('\u{201D}'),
    Some('\u{2022}'),
    Some('\u{2013}'),
    Some('\u{2014}'),
    Some('\u{02DC}'),
    Some('\u{2122}'),
    Some('\u{0161}'),
    Some('\u{203A}'),
    Some('\u{0153}'),
    None,
    Some('\u{017E}'),
    Some('\u{0178}'),
];

/// Decode Windows-1252 bytes, refusing undefined code points.
pub fn decode_cp1252(bytes: &[u8]) -> Result<String, TsgError> {
    bytes
        .iter()
        .enumerate()
        .map(|(pos, &b)| match b {
            0x80..=0x9F => CP1252_HIGH[usize::from(b - 0x80)].ok_or_else(|| {
                fmt_err(format!(
                    "header byte 0x{b:02X} at offset {pos} is undefined in cp1252"
                ))
            }),
            _ => Ok(char::from(b)),
        })
        .collect()
}

/// A `.tsg` header: ordered `[name]` sections, each a list of lines with
/// trailing whitespace removed. Lines before the first section are ignored.
#[derive(Debug, Clone)]
pub struct Header {
    sections: Vec<(String, Vec<String>)>,
}

fn section_name(line: &str) -> Option<&str> {
    let rest = line.strip_prefix('[')?;
    let end = rest.find(']')?;
    let name = &rest[..end];
    let valid = !name.is_empty() && name.chars().all(|c| c.is_ascii_alphanumeric() || c == ' ');
    valid.then_some(name)
}

fn non_blank(body: &[String]) -> Vec<&str> {
    body.iter()
        .map(String::as_str)
        .filter(|l| !l.trim().is_empty())
        .collect()
}

impl Header {
    /// Split decoded header text into sections.
    pub fn parse(text: &str) -> Self {
        let mut sections: Vec<(String, Vec<String>)> = Vec::new();
        for raw in text.lines() {
            let line = raw.trim_end();
            if let Some(name) = section_name(line) {
                sections.push((name.to_string(), Vec::new()));
            } else if let Some((_, body)) = sections.last_mut() {
                body.push(line.to_string());
            }
        }
        Self { sections }
    }

    /// Non-blank lines of the section called `name`; `Ok(None)` if absent.
    /// A name that appears twice is refused: which copy is meant would be a guess.
    pub fn section(&self, name: &str) -> Result<Option<Vec<&str>>, TsgError> {
        let mut found = self.sections.iter().filter(|(n, _)| n == name);
        let first = found.next();
        if found.next().is_some() {
            return Err(fmt_err(format!("section [{name}] appears more than once")));
        }
        Ok(first.map(|(_, body)| non_blank(body)))
    }

    /// All sections in file order with their non-blank lines.
    pub fn sections(&self) -> impl Iterator<Item = (&str, Vec<&str>)> {
        self.sections
            .iter()
            .map(|(n, body)| (n.as_str(), non_blank(body)))
    }

    /// Value of `key = value` inside `section`. A key given twice is refused.
    pub fn value(&self, section: &str, key: &str) -> Result<Option<&str>, TsgError> {
        let Some(lines) = self.section(section)? else {
            return Ok(None);
        };
        let mut found = lines.into_iter().filter_map(|l| {
            let (k, v) = l.split_once('=')?;
            (k.trim() == key).then(|| v.trim())
        });
        let first = found.next();
        if found.next().is_some() {
            return Err(fmt_err(format!("`{key}` appears twice in [{section}]")));
        }
        Ok(first)
    }

    /// Like [`Header::value`] but an absent key is a format error.
    pub fn require(&self, section: &str, key: &str) -> Result<&str, TsgError> {
        self.value(section, key)?
            .ok_or_else(|| fmt_err(format!("missing `{key}` in [{section}]")))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn decodes_cp1252_specials_and_rejects_undefined() {
        assert_eq!(decode_cp1252(b"a\x80\xB5").unwrap(), "a\u{20AC}\u{00B5}");
        let err = decode_cp1252(b"ab\x81").unwrap_err().to_string();
        assert!(err.contains("0x81 at offset 2"), "{err}");
    }

    #[test]
    fn splits_sections_and_reads_values() {
        let h = Header::parse(
            "junk\n[coordinates]\nlastband = 531 \r\n\n[wavelength specs]\n380 2500 nm\n",
        );
        let names: Vec<_> = h.sections().map(|(n, _)| n).collect();
        assert_eq!(names, ["coordinates", "wavelength specs"]);
        assert_eq!(h.value("coordinates", "lastband").unwrap(), Some("531"));
        assert_eq!(
            h.section("wavelength specs").unwrap().unwrap(),
            ["380 2500 nm"]
        );
        assert!(h.require("coordinates", "lastsample").is_err());
        assert!(h.section("absent").unwrap().is_none());
    }

    #[test]
    fn keeps_last_line_when_no_blank_line_precedes_next_section() {
        let h = Header::parse("[a]\n1:x\n2:y\n[b]\nz\n");
        assert_eq!(h.section("a").unwrap().unwrap(), ["1:x", "2:y"]);
    }

    #[test]
    fn refuses_duplicate_section() {
        let h = Header::parse("[a]\nx = 1\n[a]\nx = 2\n");
        let err = h.value("a", "x").unwrap_err().to_string();
        assert!(err.contains("[a] appears more than once"), "{err}");
    }

    #[test]
    fn refuses_duplicate_key() {
        let h = Header::parse("[a]\nx = 1\nx = 2\n");
        let err = h.value("a", "x").unwrap_err().to_string();
        assert!(err.contains("`x` appears twice in [a]"), "{err}");
    }

    #[test]
    fn bracket_lines_with_symbols_are_not_sections() {
        let h = Header::parse("[a]\n[x=1]\n");
        assert_eq!(h.section("a").unwrap().unwrap(), ["[x=1]"]);
    }
}
