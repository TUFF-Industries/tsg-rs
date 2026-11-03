//! End-to-end tests of the `tsg` binary on synthetic TSG files.

mod common;

use std::path::Path;
use std::process::{Command, Output};

use common::{temp_dir, Synth, BANDS, NULL, SAMPLES};

/// Mutates a synthetic file into an inconsistent one.
type Corrupt = Box<dyn Fn(&mut Synth)>;

fn tsg(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_tsg"))
        .args(args)
        .output()
        .expect("run tsg")
}

fn ok(args: &[&str]) -> String {
    let out = tsg(args);
    assert!(
        out.status.success(),
        "tsg {args:?} failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).unwrap()
}

fn refused(args: &[&str]) -> String {
    let out = tsg(args);
    assert!(!out.status.success(), "tsg {args:?} should have failed");
    String::from_utf8(out.stderr).unwrap()
}

fn s(p: &Path) -> &str {
    p.to_str().unwrap()
}

#[test]
fn info_reports_dimensions_grid_and_tsa_libraries() {
    let dir = temp_dir("info");
    let tsg_path = Synth::new().write(&dir);
    let out = ok(&["info", s(&tsg_path)]);
    assert!(out.contains("samples:     3"), "{out}");
    assert!(out.contains("bands:       5"), "{out}");
    assert!(
        out.contains("wavelengths: 1300 to 2500 nm, step 300 nm"),
        "{out}"
    );
    assert!(out.contains("depth:       10 to 10.2 m"), "{out}");
    assert!(out.contains("4 named columns, 3 class-coded"), "{out}");
    assert!(out.contains("reflibs:     system = SYS_TEST_SWIR"), "{out}");
    assert!(
        out.contains("[class 1] S_SWIR_TSA705 Groups (3 labels): Grp1 sTSAS"),
        "{out}"
    );
    assert!(
        out.contains("[class 2] S_SWIR_TSA705 Minerals (2 labels): Min1 sTSAS"),
        "{out}"
    );
    assert!(!out.contains("Domain"), "non-TSA class listed: {out}");
    // The .bip and the directory resolve to the same pair.
    let bip = dir.join("SYN_tsg.bip");
    assert_eq!(ok(&["info", s(&bip)]), out);
    assert_eq!(ok(&["info", s(&dir)]), out);
}

#[test]
fn spectra_npy_is_bit_identical() {
    let dir = temp_dir("npy");
    let synth = Synth::new();
    let tsg_path = synth.write(&dir);
    let out = dir.join("spectra.npy");
    ok(&["spectra", s(&tsg_path), "-o", s(&out)]);
    let bytes = std::fs::read(&out).unwrap();
    assert_eq!(&bytes[..6], b"\x93NUMPY");
    let hlen = usize::from(u16::from_le_bytes([bytes[8], bytes[9]]));
    let header = std::str::from_utf8(&bytes[10..10 + hlen]).unwrap();
    assert!(
        header.contains("'descr': '<f4'") && header.contains("'shape': (3, 5)"),
        "{header}"
    );
    let body = &bytes[10 + hlen..];
    let expected: Vec<u8> = synth.spectra.iter().flat_map(|v| v.to_le_bytes()).collect();
    assert_eq!(body, expected.as_slice());
}

#[test]
fn spectra_csv_round_trips_every_value() {
    let dir = temp_dir("csv");
    let synth = Synth::new();
    let tsg_path = synth.write(&dir);
    let out = dir.join("spectra.csv");
    ok(&["spectra", s(&tsg_path), "-o", s(&out)]);
    let text = std::fs::read_to_string(&out).unwrap();
    let mut lines = text.lines();
    assert_eq!(
        lines.next().unwrap(),
        "sample,depth_m,1300,1600,1900,2200,2500"
    );
    let rows: Vec<&str> = lines.collect();
    assert_eq!(rows.len(), SAMPLES);
    for (i, row) in rows.iter().enumerate() {
        let f: Vec<&str> = row.split(',').collect();
        assert_eq!(f[0], i.to_string());
        for b in 0..BANDS {
            let v: f32 = f[2 + b].parse().unwrap();
            assert_eq!(v.to_bits(), synth.spectra[i * BANDS + b].to_bits());
        }
    }
}

#[test]
fn scalars_decode_class_labels_and_nulls() {
    let dir = temp_dir("scalars");
    let tsg_path = Synth::new().write(&dir);
    let out = dir.join("scalars.csv");
    ok(&["scalars", s(&tsg_path), "-o", s(&out)]);
    let text = std::fs::read_to_string(&out).unwrap();
    let expected = "\
sample,depth_m,Depth (m),Grp1 sTSAS,Min1 sTSAS,Grp1 sTSAV
0,10,10,KAOLIN,Muscovite,OXIDE
1,10.1,10.1,WHITE-MICA,Muscovite,INVALID
2,10.2,10.2,,,
";
    assert_eq!(text, expected);
}

#[test]
fn duplicate_scalar_names_get_column_suffixes() {
    let dir = temp_dir("dup");
    let mut synth = Synth::new();
    synth.header = synth.header.replace(
        "3:Grp1 sTSAV;1.000000;3.000000;2;0;2;1;0",
        "3:Depth (m);0;1;2;-1",
    );
    let tsg_path = synth.write(&dir);
    let out = dir.join("scalars.csv");
    ok(&["scalars", s(&tsg_path), "-o", s(&out)]);
    let text = std::fs::read_to_string(&out).unwrap();
    assert!(
        text.starts_with("sample,depth_m,Depth (m) #0,Grp1 sTSAS,Min1 sTSAS,Depth (m) #3\n"),
        "{text}"
    );
}

/// Each corruption must be refused with a message naming the problem, and
/// must not leave an output file behind.
#[test]
fn inconsistent_files_are_refused_with_specific_messages() {
    let base = Synth::new();
    let cases: Vec<(&str, Corrupt, &str)> = vec![
        (
            "short bip",
            Box::new(|s| s.scalars.pop().map(|_| ()).unwrap()),
            "expected 120",
        ),
        (
            "sample count",
            Box::new(|s| s.header = s.header.replace("lastsample = 3", "lastsample = 4")),
            "3 sample headers but [coordinates] lastsample = 4",
        ),
        (
            "missing sample header",
            Box::new(|s| {
                s.header = s
                    .header
                    .replace("2:SYN_0001_3  T=0001 L=1 P=3 D=10.200000 X=3.0 H=SYN\n", "")
            }),
            "2 sample headers but [coordinates] lastsample = 3",
        ),
        (
            "depth on some samples only",
            Box::new(|s| s.header = s.header.replace(" D=10.100000", "")),
            "2 of 3 sample headers have a depth (D=) but sample 1 has none",
        ),
        (
            "class value out of range",
            Box::new(|s| s.scalars[1] = 3.0),
            "scalar `Grp1 sTSAS` (column 1) sample 0 holds 3",
        ),
        (
            "fractional class value",
            Box::new(|s| s.scalars[2] = 0.5),
            "scalar `Min1 sTSAS` (column 2) sample 0 holds 0.5",
        ),
        (
            "class max disagrees",
            Box::new(|s| s.header = s.header.replace("max = 3\n", "max = 4\n")),
            "[class 1] says max = 4 but lists 3 labels",
        ),
        (
            "missing class section",
            Box::new(|s| s.header = s.header.replace("[class 2]", "[class 9]")),
            "refers to [class 2], which is missing",
        ),
        (
            "duplicate class section",
            Box::new(|s| s.header = s.header.replace("[class 3]", "[class 2]")),
            "[class 2] appears more than once",
        ),
        (
            "band header outside the plane",
            Box::new(|s| s.header = s.header.replace("0:Depth (m)", "5:Depth (m)")),
            "band header 5 (`Depth (m)`) is outside the 5 scalar columns",
        ),
        (
            "TIR grid",
            Box::new(|s| {
                s.header = s
                    .header
                    .replace("1300.000000 2500.000000", "6000.0 14500.0")
            }),
            "TIR is not",
        ),
        (
            "bip one sample too long",
            Box::new(|s| s.scalars.extend([0.0; 10])),
            ".bip is 160 bytes, expected 120",
        ),
        (
            "byte order",
            Box::new(|s| {
                s.header = s
                    .header
                    .replace("[reflibs]", "[platform]\nbyteorder = 1\n\n[reflibs]")
            }),
            "[platform] byteorder = 1; only 3 is supported",
        ),
        (
            "dimensions disagree",
            Box::new(|s| {
                s.header = s
                    .header
                    .replace("[reflibs]", "[dimensions]\nbands = 6\n\n[reflibs]")
            }),
            "[dimensions] bands = 6; only 5 is supported",
        ),
        (
            "duplicate key",
            Box::new(|s| {
                s.header = s
                    .header
                    .replace("lastband = 5", "lastband = 5\nlastband = 9")
            }),
            "`lastband` appears twice in [coordinates]",
        ),
        (
            "data header bytes",
            Box::new(|s| s.header = s.header.replace("dataheader = 0", "dataheader = 64")),
            "[structure] dataheader = 64; only 0 is supported",
        ),
        (
            "first band",
            Box::new(|s| s.header = s.header.replace("firstband = 1", "firstband = 0")),
            "[coordinates] firstband = 0; only 1 is supported",
        ),
        (
            "two wavelength lines",
            Box::new(|s| {
                s.header = s
                    .header
                    .replace("2500.000000 nm\n", "2500.000000 nm\n1 2 nm\n")
            }),
            "[wavelength specs] should be one line, found 2",
        ),
        (
            "extra wavelength token",
            Box::new(|s| s.header = s.header.replace("2500.000000 nm", "2500.000000 300 nm")),
            "[wavelength specs] should be `start end unit`",
        ),
        (
            "negative class number",
            Box::new(|s| {
                s.header = s.header.replace(
                    "2:Min1 sTSAS;1.000000;3.000000;2;2",
                    "2:Min1 sTSAS;1.000000;3.000000;2;-5",
                )
            }),
            "band header 2 (`Min1 sTSAS`) has class number -5",
        ),
        (
            "empty class label",
            Box::new(|s| s.header = s.header.replace("0:Default", "0:")),
            "[class 3] label 0 is empty",
        ),
        (
            "duplicate band header",
            Box::new(|s| s.header = s.header.replace("3:Grp1 sTSAV", "2:Grp1 sTSAV")),
            "band header 2 appears more than once",
        ),
        (
            "stray sample token",
            Box::new(|s| {
                s.header = s
                    .header
                    .replace("D=10.100000 X=2.0", "D=10.100000 junk X=2.0")
            }),
            "sample 1: `junk` follows the key=value fields but is not one",
        ),
        (
            "byte-order mark",
            Box::new(|s| s.header = format!("\u{EF}\u{BB}\u{BF}{}", s.header)),
            "UTF-8 byte-order mark",
        ),
        (
            "undefined cp1252 byte",
            Box::new(|s| s.header = s.header.replace("SYNTHETIC", "SYN\u{81}")),
            "undefined in cp1252",
        ),
    ];
    for (what, corrupt, message) in cases {
        let dir = temp_dir("bad");
        let mut synth = Synth::new();
        corrupt(&mut synth);
        let tsg_path = dir.join("SYN_tsg.tsg");
        // Header bytes are written as cp1252: chars below U+0100 map to one byte.
        let bytes: Vec<u8> = synth.header.chars().map(|c| c as u32 as u8).collect();
        std::fs::write(&tsg_path, bytes).unwrap();
        std::fs::write(dir.join("SYN_tsg.bip"), synth.bip()).unwrap();
        let out = dir.join("out.csv");
        let npy = dir.join("out.npy");
        for args in [
            vec!["scalars", s(&tsg_path), "-o", s(&out)],
            vec!["spectra", s(&tsg_path), "-o", s(&npy)],
            vec!["info", s(&tsg_path)],
        ] {
            let err = refused(&args);
            assert!(
                err.contains(message),
                "{what} ({}): expected `{message}` in: {err}",
                args[0]
            );
        }
        let leftovers: Vec<_> = std::fs::read_dir(&dir)
            .unwrap()
            .map(|e| e.unwrap().file_name())
            .collect();
        assert_eq!(
            leftovers.len(),
            2,
            "{what}: output written despite refusal: {leftovers:?}"
        );
    }
    // The untouched synthetic file is accepted.
    let dir = temp_dir("good");
    let tsg_path = base.write(&dir);
    ok(&["scalars", s(&tsg_path), "-o", s(&dir.join("out.csv"))]);
}

#[test]
fn keeps_last_line_of_sections_without_trailing_blank_line() {
    let dir = temp_dir("noblank");
    let mut synth = Synth::new();
    synth.header = synth.header.replace("\n\n[", "\n[");
    let tsg_path = synth.write(&dir);
    let out = dir.join("scalars.csv");
    ok(&["scalars", s(&tsg_path), "-o", s(&out)]);
    let text = std::fs::read_to_string(&out).unwrap();
    assert_eq!(text.lines().count(), 1 + SAMPLES, "{text}");
    assert!(
        text.lines().next().unwrap().ends_with(",Grp1 sTSAV"),
        "{text}"
    );
}

#[test]
fn backward_depths_are_reported_and_missing_depths_left_empty() {
    let dir = temp_dir("depth");
    let mut synth = Synth::new();
    synth.header = synth.header.replace("D=10.200000", "D=9.000000");
    let tsg_path = synth.write(&dir);
    let out = ok(&["info", s(&tsg_path)]);
    assert!(
        out.contains(
            "depth falls more than 0.001 m below an earlier sample at 1 samples, first at sample 2"
        ),
        "{out}"
    );
    let csv = dir.join("spectra.csv");
    ok(&["spectra", s(&tsg_path), "-o", s(&csv)]);
    assert!(std::fs::read_to_string(&csv).unwrap().contains("\n2,9,"));

    let dir = temp_dir("nodepth");
    let mut synth = Synth::new();
    for d in ["10.000000", "10.100000", "10.200000"] {
        synth.header = synth.header.replace(&format!(" D={d}"), "");
    }
    let tsg_path = synth.write(&dir);
    assert!(ok(&["info", s(&tsg_path)]).contains("depth:       none (no D= fields)"));
    let csv = dir.join("scalars.csv");
    ok(&["scalars", s(&tsg_path), "-o", s(&csv)]);
    let text = std::fs::read_to_string(&csv).unwrap();
    assert!(text.contains("\n0,,10,KAOLIN,"), "{text}");
}

#[test]
fn bad_usage_is_reported() {
    let dir = temp_dir("usage");
    let tsg_path = Synth::new().write(&dir);
    let p = s(&tsg_path);
    assert!(refused(&["spectra", p]).contains("spectra needs -o"));
    assert!(refused(&["spectra", p, "-o", s(&dir.join("x.txt"))]).contains(".csv or .npy"));
    assert!(refused(&["scalars", p, "-o", s(&dir.join("x.npy"))]).contains("must end in .csv"));
    assert!(refused(&["frobnicate", p]).contains("unknown command"));
    assert!(refused(&["info", s(&dir.join("absent.tsg"))]).contains("no such file"));
    assert!(ok(&["--help"]).contains("tsg info <file>"));
}

#[test]
fn null_marker_constant_matches_library() {
    assert_eq!(NULL.to_bits(), tsg_rs::NULL_VALUE.to_bits());
}

#[test]
fn labels_keep_colons_and_cp1252_text_and_names_never_collide() {
    let dir = temp_dir("labels");
    let mut synth = Synth::new();
    synth.header = synth
        .header
        .replace("1:Muscovite", "1:Musc:ovite \u{B5}m")
        .replace("0:Depth (m);", "0:sample;")
        .replace("3:Grp1 sTSAV;", "3:depth_m;");
    let bytes: Vec<u8> = synth.header.chars().map(|c| c as u32 as u8).collect();
    let tsg_path = dir.join("SYN_tsg.tsg");
    std::fs::write(&tsg_path, bytes).unwrap();
    std::fs::write(dir.join("SYN_tsg.bip"), synth.bip()).unwrap();
    let out = dir.join("scalars.csv");
    ok(&["scalars", s(&tsg_path), "-o", s(&out)]);
    let text = std::fs::read_to_string(&out).unwrap();
    let mut lines = text.lines();
    assert_eq!(
        lines.next().unwrap(),
        "sample,depth_m,sample #0,Grp1 sTSAS,Min1 sTSAS,depth_m #3"
    );
    assert_eq!(
        lines.next().unwrap(),
        "0,10,10,KAOLIN,Musc:ovite \u{B5}m,OXIDE"
    );
}

#[test]
fn special_floats_round_trip_through_csv() {
    let dir = temp_dir("floats");
    let mut synth = Synth::new();
    let special = [
        -0.0,
        1e-45,
        f32::MIN_POSITIVE,
        f32::INFINITY,
        f32::NEG_INFINITY,
        f32::NAN,
        3.0e38,
    ];
    synth.spectra[..special.len()].copy_from_slice(&special);
    let tsg_path = synth.write(&dir);
    let out = dir.join("spectra.csv");
    ok(&["spectra", s(&tsg_path), "-o", s(&out)]);
    let text = std::fs::read_to_string(&out).unwrap();
    let values: Vec<f32> = text
        .lines()
        .skip(1)
        .flat_map(|l| {
            l.split(',')
                .skip(2)
                .map(|v| v.parse::<f32>().unwrap())
                .collect::<Vec<_>>()
        })
        .collect();
    assert_eq!(values.len(), synth.spectra.len());
    for (got, want) in values.iter().zip(&synth.spectra) {
        if want.is_nan() {
            assert!(got.is_nan());
        } else {
            assert_eq!(got.to_bits(), want.to_bits(), "{want}");
        }
    }
}

#[test]
fn pair_resolution_failures_are_specific() {
    let dir = temp_dir("pairs");
    assert!(refused(&["info", s(&dir)]).contains("no *tsg.tsg in"));
    let tsg_path = Synth::new().write(&dir);
    std::fs::copy(&tsg_path, dir.join("OTHER_tsg.tsg")).unwrap();
    assert!(refused(&["info", s(&dir)]).contains("more than one *tsg.tsg"));
    let lone = temp_dir("lone");
    std::fs::copy(&tsg_path, lone.join("X_tsg.tsg")).unwrap();
    assert!(refused(&["info", s(&lone.join("X_tsg.tsg"))]).contains("no .bip file next to"));
    let txt = lone.join("notes.txt");
    std::fs::write(&txt, "x").unwrap();
    assert!(refused(&["info", s(&txt)]).contains("neither a .tsg nor a .bip"));
}

#[test]
fn failed_rename_leaves_no_partial_file() {
    let dir = temp_dir("rename");
    let tsg_path = Synth::new().write(&dir);
    let target = dir.join("taken.csv");
    std::fs::create_dir(&target).unwrap();
    refused(&["spectra", s(&tsg_path), "-o", s(&target)]);
    assert!(!dir.join("taken.csv.partial").exists());
}

#[test]
fn double_dash_allows_names_starting_with_a_dash() {
    let dir = temp_dir("dash");
    let synth = Synth::new();
    std::fs::write(dir.join("-x_tsg.tsg"), &synth.header).unwrap();
    std::fs::write(dir.join("-x_tsg.bip"), synth.bip()).unwrap();
    let out = Command::new(env!("CARGO_BIN_EXE_tsg"))
        .current_dir(&dir)
        .args(["info", "--", "-x_tsg.tsg"])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
}

#[test]
fn closed_stdout_is_not_a_crash() {
    let dir = temp_dir("pipe");
    let tsg_path = Synth::new().write(&dir);
    let mut child = Command::new(env!("CARGO_BIN_EXE_tsg"))
        .args(["info", s(&tsg_path)])
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::piped())
        .spawn()
        .unwrap();
    drop(child.stdout.take());
    let out = child.wait_with_output().unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
}

#[cfg(unix)]
#[test]
fn non_utf8_paths_are_refused_not_panicked() {
    use std::os::unix::ffi::OsStrExt;
    let dir = temp_dir("utf8");
    let bad = dir.join(std::ffi::OsStr::from_bytes(b"\xff_tsg.tsg"));
    let out = Command::new(env!("CARGO_BIN_EXE_tsg"))
        .arg("info")
        .arg(&bad)
        .output()
        .unwrap();
    assert_eq!(out.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&out.stderr).contains("no such file"));
}
