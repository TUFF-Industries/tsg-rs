# tsg-rs

A reader and command-line tool for the NIR/SWIR spectra in TSG (The Spectral
Geologist) files, as published by Australia's National Virtual Core Library
(NVCL) for HyLogger drill-core scans.

A TSG NIR/SWIR dataset is a pair of files:
- a `.tsg` text header;
- a `.bip` binary file holding two planes: the reflectance spectra, and a
  table of per-sample scalars.

The scalar table includes the TSA mineral and group results that TSG
computed, stored as class codes. tsg-rs reads both planes, decodes the class
codes to their labels, and writes CSV or NumPy `.npy`. The TSA results come
out as ordinary columns (`Grp1 sTSAS`, `Min1 sTSAS`, `Grp1 uTSAS`, ...).

Every command (`info`, `spectra` and `scalars`) refuses a file that fails a
consistency check, with a message saying what is wrong, rather than guessing. Examples: a `.bip` of the wrong size, a sample
count that disagrees with the header, or a class value that is not a valid
label.

## Install

Rust 1.83 or newer. The crate has no dependencies.

```bash
git clone https://github.com/TUFF-Industries/tsg-rs.git && cd tsg-rs
cargo install --path .
```

## Use

`<file>` can be the `.tsg` or the `.bip` of a pair, or a directory holding
exactly one `*tsg.tsg`.

```text
$ tsg info 05KCD001_tsg.tsg
header:      05KCD001_tsg.tsg
data:        05KCD001_tsg.bip
samples:     22356
bands:       531
wavelengths: 380 to 2500 nm, step 4 nm
depth:       59.5 to 201.500009 m
scalars:     151 named columns, 21 class-coded
reflibs:     system = SYS_TSA63_SWIR
reflibs:     system2 = SYS_TSA62_VNIR
TSA libraries:
  [class 0] S_VNIR_TSA704 Groups (8 labels): Grp1 sTSAV, Grp2 sTSAV
  [class 1] S_SWIR_TSA705 Groups (17 labels): Grp1 sTSAS, Grp2 sTSAS, Grp3 sTSAS
  [class 2] U_SWIR_TSA705 Minerals (61 labels): Min1 uTSAS, Min2 uTSAS, Min3 uTSAS
  [class 3] S_VNIR_TSA704 Minerals (21 labels): Min1 sTSAV, Min2 sTSAV
  [class 4] U_SWIR_TSA705 Groups (17 labels): Grp1 uTSAS, Grp2 uTSAS, Grp3 uTSAS
  [class 7] S_SWIR_TSA705 Minerals (61 labels): Min1 sTSAS, Min2 sTSAS, Min3 sTSAS

$ tsg spectra 05KCD001_tsg.tsg -o spectra.npy     # float32, shape (samples, bands)
$ tsg spectra 05KCD001_tsg.tsg -o spectra.csv     # sample,depth_m,<one column per wavelength>
$ tsg scalars 05KCD001_tsg.tsg -o scalars.csv     # sample,depth_m,<one column per scalar>
```

Here are a few columns of `scalars.csv`:

```text
sample,depth_m,Grp1 sTSAS,Min1 sTSAS,Grp1 uTSAS
8999,116.72377,CHLORITE,Chlorite-Mg,CHLORITE
9000,116.730707,CHLORITE,Chlorite-Mg,CHLORITE
9001,116.737634,CHLORITE,Chlorite-Mg,CHLORITE
```

Output conventions:
- `depth_m` is the `D=` field of each sample header. It is empty if the file
  has no depths.
- Null scalars (TSG's missing-value marker) are written as empty fields.
- Spectra are written exactly as stored, including any marker values.
- Finite floats in CSV are written in the shortest form that parses back to
  the same float32, so CSV output is lossless for them. NaN and infinities are
  written as `NaN`, `inf` and `-inf`; a NaN's sign and payload are not kept.
- A scalar name that occurs more than once, or is `sample` or `depth_m`, gets
  ` #<column>` appended, so every CSV header is unique.
- The `.npy` holds spectra only. The wavelengths are evenly spaced, so the
  `info` line gives the whole grid.
- Output is written to `<out>.partial` and renamed when complete. A refused
  file leaves no output behind.

As a library:

```rust
let f = tsg_rs::read_pair("x_tsg.tsg".as_ref(), "x_tsg.bip".as_ref())?;
let first = f.spectrum(0); // &[f32], one value per band in f.wavelengths
for col in f.scalar_columns()? {
    if let tsg_rs::ScalarValues::Class { labels, .. } = &col.values {
        println!("{}: {:?}", col.header.name, labels.iter().take(3).collect::<Vec<_>>());
    }
}
```

## Verified against

[VERIFY.md](VERIFY.md) covers five real NVCL holes, one each from WA, QLD,
SA, NT and NSW. On all five:
- spectra are bit-identical to pytsg 1.0.0, in both `.npy` and CSV output;
- wavelengths and depths are equal to pytsg's;
- every scalar column matches pytsg's, except one class label that pytsg
  truncates at a colon;
- the decoded `Grp1 sTSAS` is identical, row for row, to the survey's
  published `Grp1_sTSAS.csv`.

A sixth hole, WA BDRCD019, was checked as a known exception. Its spectra and
scalars match pytsg, but its embedded `Grp1 sTSAS` differs from the survey CSV
on 5,043 of 40,200 rows. An earlier scan of 197 files found one more such
file. tsg-rs reports what the file stores; it cannot tell you which TSA run is
the right one.

As a smoke run, tsg-rs read 196 of 197 real NVCL TSG files. The one it
refused was a TIR file.

[docs/tsg-format.md](docs/tsg-format.md) lists every format assumption, how
each was checked, and what the reader refuses.

## Speed compared with pytsg

On an Apple M3 Pro with a warm cache, a 95 MB `.bip` (22,356 spectra × 531
bands) took:

| Step | Time |
|---|---|
| `tsg info` (whole process: parse, read, check) | 0.11 s |
| `tsg spectra` to `.npy` | 0.19 s |
| `tsg scalars` to CSV | 0.53 s |
| `tsg spectra` to CSV | 1.7 s |
| pytsg `read_spectra`, in a running Python process | 0.16 s |
| Python start-up plus `import pytsg` | about 0.35 s more |

Parsing is somewhat faster than pytsg. Exporting to `.npy` takes about as
long as pytsg's read, and CSV export is slower. Use tsg-rs if you want a
single binary with no Python, strict checks, or decoded class columns, rather
than for speed. The full table is in VERIFY.md.

## Limitations

- **NIR/SWIR pair only** (`*_tsg.tsg` + `*_tsg.bip`). TIR is refused by its
  wavelength grid, whatever the file is called: any grid outside 300–2600 nm
  is rejected.
- **Not read:** CRAS core imagery (`*_cras.bip`), profilometer data
  (`*_hires.dat`), and the `.cal`/`.ini` side files.
- **Only `datatype = float`, `dataheader = 0` and `byteorder = 3` files.**
  The `.bip` is read as little-endian float32. Byte order 3 is the only value
  seen on real files, and the bit-identical match with pytsg confirms it.
  Other values are refused, not converted.
- **Class codes are taken as written.** tsg-rs does not run TSA, and it does
  not judge whether the embedded results are right. It reports what TSG
  stored.
- **The format is undocumented.** The layout comes from pytsg's source and was
  checked on real files (docs/tsg-format.md). Files written by other TSG
  versions or tools may break an assumption. If they do, tsg-rs refuses them
  rather than reading them on a guess.
- **Metadata sections are kept in `SpectralFile::header()` but not
  interpreted.** Apart from the `[reflibs]` library names, the CLI prints none
  of them. `[events]`, for example, can hold user names and licence strings.

## Licence and credits

Licensed under either the MIT licence ([LICENSE-MIT](LICENSE-MIT)) or the
Apache License 2.0 ([LICENSE-APACHE](LICENSE-APACHE)), at your option.

The file layout was learned from
[pytsg](https://github.com/Geological-Survey-of-Western-Australia/pytsg)
(Geological Survey of Western Australia, MIT). See [NOTICE](NOTICE).

The examples and the verification use public HyLogger data from the
Geological Survey of Western Australia and the other state and territory
geological surveys, accessed through AuScope's National Virtual Core Library
(NVCL). No data are included in this repository.
