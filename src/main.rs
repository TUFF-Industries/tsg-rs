//! `tsg`: inspect and export the NIR/SWIR spectra and scalars of a TSG file.

use std::ffi::OsString;
use std::fmt::Write as _;
use std::io::{self, Write as _};
use std::path::{Path, PathBuf};
use std::process::ExitCode;

use tsg_rs::{resolve_pair, write, SpectralFile, TsgError};

const USAGE: &str = "\
tsg: read the NIR/SWIR spectra and scalars of a TSG (.tsg + .bip) file

Usage:
  tsg info <file>
  tsg spectra <file> -o <out.csv|out.npy>
  tsg scalars <file> -o <out.csv>

<file> is the .tsg or the .bip of a pair, or a directory holding one *tsg.tsg.
Use `--` before a file name that starts with `-`.
Files that fail a consistency check are refused with a reason.";

enum Command {
    Info(PathBuf),
    Spectra(PathBuf, PathBuf),
    Scalars(PathBuf, PathBuf),
}

fn parse_args(args: &[OsString]) -> Result<Command, String> {
    let (cmd, rest) = args.split_first().ok_or("missing command")?;
    let cmd = cmd.to_str().ok_or("the command is not valid UTF-8")?;
    if !matches!(cmd, "info" | "spectra" | "scalars") {
        return Err(format!("unknown command `{cmd}`"));
    }
    let mut input = None;
    let mut output = None;
    let mut options_done = false;
    let mut it = rest.iter();
    while let Some(a) = it.next() {
        match a.to_str() {
            Some("--") if !options_done => options_done = true,
            Some("-o" | "--output") if !options_done => {
                let o = it.next().ok_or("-o needs a path")?;
                if output.replace(PathBuf::from(o)).is_some() {
                    return Err("-o given twice".into());
                }
            }
            Some(s) if !options_done && s.starts_with('-') && s.len() > 1 => {
                return Err(format!("unknown option {s}"));
            }
            _ => {
                if input.replace(PathBuf::from(a)).is_some() {
                    return Err("more than one input file".into());
                }
            }
        }
    }
    let input = input.ok_or("missing <file>")?;
    match (cmd, output) {
        ("info", None) => Ok(Command::Info(input)),
        ("info", Some(_)) => Err("info takes no -o".into()),
        ("spectra", Some(o)) => Ok(Command::Spectra(input, o)),
        ("scalars", Some(o)) => Ok(Command::Scalars(input, o)),
        _ => Err(format!("{cmd} needs -o <out>")),
    }
}

fn ext(p: &Path) -> String {
    p.extension()
        .and_then(|e| e.to_str())
        .unwrap_or_default()
        .to_ascii_lowercase()
}

/// Read and fully check a pair, class values included, whatever the command.
fn load(input: &Path) -> Result<(PathBuf, PathBuf, SpectralFile), String> {
    let (tsg, bip) = resolve_pair(input).map_err(|e| e.to_string())?;
    let f = tsg_rs::read_pair(&tsg, &bip).map_err(|e| e.to_string())?;
    f.scalar_columns().map_err(|e| e.to_string())?;
    Ok((tsg, bip, f))
}

fn io_err(path: &Path) -> impl FnOnce(io::Error) -> String + '_ {
    move |e| format!("writing {}: {e}", path.display())
}

fn info_text(tsg: &Path, bip: &Path, f: &SpectralFile) -> Result<String, TsgError> {
    let wl = f.wavelengths();
    let (first, last) = (wl[0], wl[wl.len() - 1]);
    let step = (last - first) / (wl.len() - 1) as f64;
    let unit = f.wavelength_unit();
    let n_class = f
        .band_headers()
        .iter()
        .filter(|h| h.class_section().is_some())
        .count();
    let mut o = String::new();
    // Writing to a String cannot fail.
    let _ = writeln!(o, "header:      {}", tsg.display());
    let _ = writeln!(o, "data:        {}", bip.display());
    let _ = writeln!(o, "samples:     {}", f.n_samples());
    let _ = writeln!(o, "bands:       {}", f.n_bands());
    let _ = writeln!(
        o,
        "wavelengths: {first} to {last} {unit}, step {step} {unit}"
    );
    let known: Vec<f64> = f.depths().into_iter().flatten().collect();
    if known.is_empty() {
        let _ = writeln!(o, "depth:       none (no D= fields)");
    } else {
        let lo = known.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = known.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let _ = writeln!(o, "depth:       {lo} to {hi} m");
    }
    let back = f.depth_backsteps();
    if let Some(first) = back.first() {
        let _ = writeln!(
            o,
            "             depth falls more than {} m below an earlier sample at {} samples, \
             first at sample {first}",
            tsg_rs::DEPTH_TOLERANCE_M,
            back.len()
        );
    }
    let _ = writeln!(
        o,
        "scalars:     {} named columns, {n_class} class-coded",
        f.band_headers().len()
    );
    for key in ["system", "system2"] {
        if let Some(v) = f.header().value("reflibs", key)? {
            let _ = writeln!(o, "reflibs:     {key} = {v}");
        }
    }
    let tsa = f.tsa_tables();
    let _ = writeln!(
        o,
        "TSA libraries:{}",
        if tsa.is_empty() { " none" } else { "" }
    );
    for (t, cols) in tsa {
        let mut names: Vec<&str> = cols.iter().map(|h| h.name.as_str()).collect();
        names.sort_unstable();
        let used = if names.is_empty() {
            "no scalar uses it".to_string()
        } else {
            names.join(", ")
        };
        let _ = writeln!(
            o,
            "  [class {}] {} ({} labels): {used}",
            t.number,
            t.name,
            t.labels.len()
        );
    }
    Ok(o)
}

/// Print to stdout; a reader that has gone away (`tsg info x | head`) is not an error.
fn print(text: &str) -> Result<(), String> {
    match io::stdout().lock().write_all(text.as_bytes()) {
        Err(e) if e.kind() != io::ErrorKind::BrokenPipe => Err(format!("stdout: {e}")),
        _ => Ok(()),
    }
}

fn info(input: &Path) -> Result<(), String> {
    let (tsg, bip, f) = load(input)?;
    print(&info_text(&tsg, &bip, &f).map_err(|e| e.to_string())?)
}

fn spectra(input: &Path, out: &Path) -> Result<(), String> {
    let kind = ext(out);
    if kind != "csv" && kind != "npy" {
        return Err(format!(
            "{}: output must end in .csv or .npy",
            out.display()
        ));
    }
    let (_, _, f) = load(input)?;
    if kind == "npy" {
        write::spectra_npy(&f, out).map_err(io_err(out))?;
    } else {
        write::spectra_csv(&f, out).map_err(io_err(out))?;
    }
    eprintln!(
        "wrote {} spectra x {} bands to {}",
        f.n_samples(),
        f.n_bands(),
        out.display()
    );
    Ok(())
}

fn scalars(input: &Path, out: &Path) -> Result<(), String> {
    if ext(out) != "csv" {
        return Err(format!("{}: output must end in .csv", out.display()));
    }
    let (_, _, f) = load(input)?;
    let cols = f.scalar_columns().map_err(|e| e.to_string())?;
    write::scalars_csv(&f, &cols, out).map_err(io_err(out))?;
    eprintln!(
        "wrote {} scalars for {} samples to {}",
        cols.len(),
        f.n_samples(),
        out.display()
    );
    Ok(())
}

fn main() -> ExitCode {
    let args: Vec<OsString> = std::env::args_os().skip(1).collect();
    let result = match args.first().and_then(|a| a.to_str()) {
        None if args.is_empty() => print(&format!("{USAGE}\n")),
        Some("-h" | "--help" | "help") => print(&format!("{USAGE}\n")),
        Some("-V" | "--version") => print(&format!("tsg {}\n", env!("CARGO_PKG_VERSION"))),
        _ => parse_args(&args)
            .map_err(|e| format!("{e}\n\n{USAGE}"))
            .and_then(|cmd| match cmd {
                Command::Info(i) => info(&i),
                Command::Spectra(i, o) => spectra(&i, &o),
                Command::Scalars(i, o) => scalars(&i, &o),
            }),
    };
    match result {
        Ok(()) => ExitCode::SUCCESS,
        Err(e) => {
            eprintln!("tsg: {e}");
            ExitCode::FAILURE
        }
    }
}
