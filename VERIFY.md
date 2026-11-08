# Verification on real NVCL files

Run on 2026-11-03 with `verify/run.sh`, which calls `verify/compare.py` once
per hole. The raw output is in `verify/results-2026-11-03.jsonl`. The data are public NVCL HyLogger downloads that were already on
disk; they are not in this repository.

Setup:
- tsg-rs 0.1.0, release build.
- Reference reader: pytsg 1.0.0, Python 3.12.14, numpy 2.5.4, pandas 3.1.0, run with `python -I`.
- Machine: Apple M3 Pro, warm file cache.

## What was checked

For each `.tsg`/`.bip` pair:
1. **Spectra:** `tsg spectra -o x.npy` and `tsg spectra -o x.csv`, each compared
   with pytsg's `read_spectra(...).spectra` as float32 bit patterns. The CSV
   values were parsed back to float32 before comparing.
2. **Wavelengths and depths:** compared with pytsg's `wavelength` and its
   sample-header `D` column.
3. **Every scalar column**, from `tsg scalars -o x.csv`, against pytsg's
   `scalars` DataFrame:
   - numeric columns: bit patterns, with the null marker treated as empty;
   - class-coded columns: label strings.
   pytsg leaves `[class 0]` columns undecoded, so for those the script maps the
   raw values through pytsg's own class table.
4. **`Grp1 sTSAS` against the survey:** the label decoded from the scalar plane
   was compared with the survey's published `Grp1_sTSAS.csv` for that hole, row
   by row over every sample. The survey's `null` matches an empty field.
   The survey CSV is the `Grp1 sTSAS` scalar log published by each state's
   NVCL data service (`downloadscalars`), one row per TSG sample.

The five holes are the three used in the earliest pytsg cross-check (one each
from WA, QLD and SA) plus the alphabetically first NT and NSW holes that had a
survey `Grp1_sTSAS.csv` on disk. A sixth hole, WA BDRCD019, is listed
separately as a known exception: an earlier scan of 197 NVCL TSG files, made
outside this repository, had found its survey CSV disagreeing with its
embedded TSA.

## Results

| Hole | Samples × bands | Spectra .npy = pytsg | Spectra .csv = pytsg | Wavelengths, depths | Scalar columns = pytsg | Grp1 sTSAS = survey CSV |
|---|---|---|---|---|---|---|
| WA 05KCD001 | 22,356 × 531 | bit-identical | bit-identical | equal | 151 / 151 | 22,356 / 22,356 rows (14,499 non-null) |
| QLD ChartersTowers1 | 25,544 × 531 | bit-identical | bit-identical | equal | 154 / 154 | 25,544 / 25,544 (19,519 non-null) |
| SA 103960 (ANDDH1) | 27,776 × 531 | bit-identical | bit-identical | equal | 156 / 157 (see note) | 27,776 / 27,776 (14,452 non-null) |
| NT 1113632 (DD85GL5) | 16,375 × 531 | bit-identical | bit-identical | equal | 125 / 125 | 16,375 / 16,375 (15,530 non-null) |
| NSW MIN_000610 (BU-1) | 244 × 531 | bit-identical | bit-identical | equal | 152 / 152 | 244 / 244 (142 non-null) |
| *Exception:* WA BDRCD019 | 40,200 × 531 | bit-identical | bit-identical | equal | 142 / 142 | 35,157 / 40,200 |

Survey rows and TSG samples line up one-to-one: on every hole the survey's
`StartDepth` and the TSG `D=` depth differ by at most 0.0001 m.

**Note on SA 103960.** The one scalar column that differs is `RockMarks`. Its
`[class 11]` has a single label, `Operator|tray0032: E.O.H`. pytsg splits label
lines on every `:`, so it returns `Operator|tray0032`; tsg-rs returns the whole
label. The values agree everywhere else in that column: 27,775 nulls and one
label. I checked this by hand; it is not in the results file.

**Exception: WA BDRCD019.** The decoded embedded `Grp1 sTSAS` disagrees with
the survey CSV on 5,043 of 40,200 rows. The most common pairs (survey vs
embedded) are:

| Survey | Embedded | Rows |
|---|---|---|
| CHLORITE | AMPHIBOLE | 1,526 |
| AMPHIBOLE | CHLORITE | 483 |
| DARK-MICA | CHLORITE | 420 |
| CHLORITE | EPIDOTE | 326 |
| WHITE-MICA | CHLORITE | 260 |

The spectra and all scalars match pytsg on this hole, and its depths line up
with the survey, so the reader is not the cause. The likely explanation is
that the published CSV and the TSG file carry different TSA runs. That has not
been confirmed.

The same scan over 197 `*_tsg.tsg` files found:
- 2 files whose embedded and survey `Grp1 sTSAS` disagree: this hole and WA
  ARMCRCD002A;
- 177 files where they agree on every sample where both have a valid label.
  That scan skipped rows where the survey is null, so it is a looser test
  than the row-by-row check above;
- 18 files with no survey CSV or no `Grp1 sTSAS` to compare.

## Whole-collection smoke run

`verify/smoke.sh` ran `tsg scalars` and `tsg info` on all 197 `*_tsg.tsg`
files on disk. The per-file log is `verify/smoke-2026-11-03.tsv`.
- **196 were read.** Every class-coded value on them was a valid label index.
- **1 was refused:** a TIR file, `MS9-0 May2019_tsg.tsg`, whose 6000–14500 nm
  grid is outside the supported range.
- **3 files have depths more than 1 mm below an earlier sample**, at 4,399,
  2,006 and 2 samples. Earlier builds refused these.
- **8 files have no `D=` depths at all**: 7 from SA and 1 from WA. Earlier
  builds refused these too.

Both behaviours were changed on purpose; see docs/tsg-format.md claims 9–10.
Spectra and labels were compared only on the six holes above, not on all 197.

## Speed

Median of 3 wall-clock runs, in seconds:
- The `tsg` columns time the whole process, including writing the output file.
- `pytsg read` times `read_spectra` inside an already-running Python process,
  with nothing written.
- A fresh Python process spends a further ~0.35 s on start-up and
  `import pytsg`.

| Hole | .bip MB | `tsg info` (parse, read, check) | `tsg spectra` → .npy | `tsg spectra` → .csv | `tsg scalars` → .csv | pytsg read |
|---|---|---|---|---|---|---|
| WA 05KCD001 | 95.0 | 0.105 | 0.186 | 1.726 | 0.527 | 0.160 |
| QLD ChartersTowers1 | 108.5 | 0.113 | 0.174 | 1.948 | 0.613 | 0.177 |
| SA 103960 | 118.0 | 0.123 | 0.196 | 2.126 | 0.648 | 0.212 |
| NT DD85GL5 | 69.6 | 0.073 | 0.109 | 1.264 | 0.338 | 0.097 |
| NSW BU-1 | 1.0 | 0.006 | 0.006 | 0.024 | 0.010 | 0.008 |
| WA BDRCD019 | 170.8 | 0.178 | 0.274 | 3.047 | 0.925 | 0.295 |

Reading and checking a file (`tsg info`) takes 58–75% of pytsg's in-process
read time, with no Python start-up to add. Exporting to `.npy` takes about as
long as pytsg's read alone. CSV export is dominated by float-to-text
formatting: on the five large files it is 9–12 times slower than `.npy`.
tsg-rs is somewhat faster than pytsg at parsing but not meaningfully faster
for bulk export. Its advantages are a dependency-free binary, refusals on
inconsistent files, and decoded class columns.

## Reproduce

```bash
verify/run.sh <raw-data-root> <python-with-pytsg-1.0.0> <workdir>
verify/smoke.sh <raw-data-root> <workdir>
```

`<raw-data-root>/<hole>/tsg/*_tsg.{tsg,bip}` and
`<raw-data-root>/<hole>/scalars/Grp1_sTSAS.csv` must already exist; the hole
list is in `verify/run.sh`. Treat the downloaded files as untrusted, which is
why the script runs Python with `-I`.
