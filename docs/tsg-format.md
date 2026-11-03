# TSG NIR/SWIR file: what tsg-rs assumes

Scope: the NIR/SWIR pair `<name>_tsg.tsg` + `<name>_tsg.bip` of a TSG package
(HyLogger VNIR+SWIR). The reader does not read TIR spectra (usually
`_tsg_tir.*`, but recognised by their wavelength grid, not their name), CRAS core
imagery (`_tsg_cras.bip`), profilometer data (`_tsg_hires.dat`), or the
`.cal`/`.ini` side files.

The layout comes from the source of pytsg 1.0.0 (GSWA, MIT),
`pytsg/readers/spectra.py` and `models.py`, and was then checked on real NVCL
packages. Each claim below says how it was checked:

- **"Real files"** means the five verification holes in [VERIFY.md](../VERIFY.md):
  - WA 05KCD001
  - QLD ChartersTowers1
  - SA 103960
  - NT DD85GL5
  - NSW BU-1
- **"All 197"** means a header and data scan of 197 `*_tsg.tsg` files from
  193 NVCL holes in WA, QLD, SA, NT and NSW. The files themselves are not
  redistributed.

Where tsg-rs refuses a file, the claim's last column gives the refusal.

| # | Claim | How checked | tsg-rs |
|---|---|---|---|
| 1 | The `.tsg` file is cp1252 text; pytsg opens it with `encoding="cp1252"`. The bytes 0x81, 0x8D, 0x8F, 0x90 and 0x9D are undefined in cp1252. | pytsg source. Every real header decoded. | Refuses an undefined byte and reports its offset. Also refuses a UTF-8 byte-order mark. |
| 2 | The header is split into sections by lines matching `^\[[A-Za-z0-9 ]+\]`. Section bodies hold `key = value` or `index:payload` lines. | pytsg regex. WA has 31 sections, among them `[coordinates]`, `[band headers]`, `[sample headers]`, `[class N]` and `[wavelength specs]`. | Refuses a duplicate of any section it reads, and a key given twice in a section it reads. |
| 3 | `[coordinates] lastsample` and `lastband` are the sample and band **counts**: pytsg reshapes with them, and `firstsample = firstband = 1`. They match `[dimensions] samples` and `bands`. | WA: lastsample 22356 = samples 22356, and lastband 531 = bands 531. They are equal on all 197 (header scan). | Refuses a missing or non-integer count, any `firstsample` or `firstband` other than 1, and a `[dimensions]` count that disagrees. |
| 4 | The `.bip` is headerless little-endian float32 in 2 planes × samples × bands. Plane 0 holds the reflectance spectra. Plane 1 is the scalar table, where column k is described by the `[band headers]` record with index k. The file size is therefore exactly 2·S·B·4. | WA: 2·22356·531·4 = 94,968,288 bytes, the file size. Spectra are bit-identical to pytsg on the real files. | Refuses any other size and states the expected size. |
| 5 | `[structure] datatype = float` and `dataheader = 0`, so the `.bip` has no header bytes. | All 197 have both keys (header scan). | Refuses any other datatype, or a `dataheader` other than 0. An absent `dataheader` is accepted. |
| 6 | `[platform] fpcheck = -971227136` is 0xC61C4000, the float32 bit pattern of −10000.0, which looks like a byte-order canary. | Arithmetic on the WA value; this is a hypothesis only. `byteorder = 3` on all 197 (header scan). | `fpcheck` is not used. A `byteorder` other than 3 is refused. The data are read as little-endian, which claim 4's bit-identity supports. |
| 7 | Wavelengths are not stored per band. `[wavelength specs]` gives `start end ... unit`, and band i is at `start + i·(end−start)/(B−1)`. pytsg uses `np.arange(start, end+step, step)`, which can yield B+1 values when the floats round badly. tsg-rs computes each band by index. | The real files have `380.000000 2500.000000 nm` and 531 bands, so 4 nm spacing. Equal to pytsg on all five. All 197 have exactly one line of three tokens (header scan). | Refuses anything but one `start end unit` line, non-increasing or non-numeric bounds, fewer than 2 bands, a unit other than `nm`, and any grid outside 300–2600 nm (TIR is 6000–14500). |
| 8 | `[sample headers]` lines read `i:name  T=tray L=line P=pos D=depth_m X=... H=hole`, and depth is the `D=` field. The count must equal `lastsample`, and indices must run 0..S−1 in order. | Real files. Depths are equal to pytsg's `D` column on all five. Names can contain spaces (WA GUD175: `Mean VNIR-SWIR magnetite_0001  GUD174_tsg.A.1  (3) I=-3`) or `=` (SA 131270: `131270_0001_370__T=0001_L=2_..._D=28.428935_...`, one token with no `D=` field). | The first token is always part of the name, and so are later tokens up to the first `key=value`. Refuses a count mismatch, an index out of order, a non-numeric `D=`, a repeated key, or a non-`key=value` token after the fields. |
| 9 | Some files have no `D=` field at all, for example spectra-library datasets. | 8 of 197 files, 7 from SA and 1 from WA (verify/smoke-2026-11-03.tsv). | All-or-none: no depths gives an empty `depth_m` column. Depths on only some samples is refused. |
| 10 | Depths are not strictly ordered. Consecutive `D=` values within a scan line step back by a few µm, for example WA 59.500008 → 59.500005. | In 3 of 197 files some depths are more than 1 mm below an earlier sample, at 4,399, 2,006 and 2 samples (verify/smoke-2026-11-03.tsv). | Read as stored. `info` reports samples more than 1 mm below the running maximum. |
| 11 | `[band headers]` lines read `k:Name;min;max;flag;class;...`. When `flag ≤ 2`, field 4 is an integer class number (−1 for none). A scalar whose `flag = 2` and class number is ≥ 0 holds a class index into `[class N]`. A name-only line such as `0:Flags` is a built-in column. Band numbers can skip; WA has no 78 or 79. | pytsg `_parse_bandheaders` and `_parse_scalars`. Every class-coded column decoded on 196 of 197 files; the one file not read is TIR. | Refuses a band index outside the plane, a band index used twice, or a missing integer flag or class field. For `flag = 2` it also refuses a class number other than −1 or a `[class N]` number, and a class number ≥ 0 with no section. |
| 12 | `[class N]` has `name =`, `max =` (the label count) and optionally `colours =`, followed by `0:LABEL`, `1:LABEL` and so on. A label may itself contain `:`; SA 103960 has `Operator\|tray0032: E.O.H`. | On all 197, `max` equals the label count and the labels are contiguous (header scan); no label is empty. | Refuses a missing name or max, `max` ≠ the label count, labels out of order, an empty label, or a line that is neither `index:label` nor `key = value` (other keys are ignored). The label is everything after the first `:`. A line is read as a label if it starts with digits and `:`; pytsg checks for `=` first, so a label containing `=` would be read differently. |
| 13 | TSA results are ordinary class-coded scalars. The WA example is `91:Grp1 sTSAS;...;2;1;...` with `[class 1] name = S_SWIR_TSA705 Groups`. Library names contain `TSA`, for example `S_SWIR_TSA705 Groups` and `U_SWIR_TSA705 Minerals`. | Decoded `Grp1 sTSAS` is identical to the survey's published `Grp1_sTSAS.csv` on all five real files (see VERIFY.md). | `info` lists class tables whose name contains `TSA`, with their scalars. |
| 14 | The missing value in the scalar plane is float32 0xFF7FFFFE (−3.4028233e38), one step above `f32::MIN`. pytsg tests for it with `np.isclose(v, finfo(float32).min)`, but only in class-coded columns; its numeric columns keep the raw value. | WA: the only non-label value in the TSA columns. | tsg-rs uses the same `isclose` test in every scalar column. A null becomes an empty CSV field. The spectra plane is written as stored. |
| 15 | A class value must be a whole number in 0..max−1. pytsg truncates with `astype(int)`. | No real file breaks this. | Every command refuses any other value and names the column and the sample. |
| 16 | Scalar names are not unique. WA YDH12 repeats 15 system scalars (`Domain`, `Final Mask`, `HoleID`, ...). | 1 of 197 files. | CSV headers get ` #<column>` appended to every repeated name and to a scalar called `sample` or `depth_m`, until all headers are unique. |
| 17 | pytsg 1.0.0 slices each section as `lines[start+1 : next_header-1]`, which drops the last line of every section and the last line of the file. Real TSG-written files put a blank line there, so nothing is lost. A header without those blank lines loses its last sample header, band header and class label. | Section bounds in the pytsg source; synthetic repro. pytsg reads every real file correctly. | Keeps every non-blank line; tested on a synthetic file without blank lines. |

Known differences from pytsg 1.0.0. Neither affects the spectra.
- pytsg leaves `[class 0]` scalars as raw floats, because its check is
  `class_number > 0`. tsg-rs decodes them.
- pytsg cuts a class label at its second `:`. tsg-rs keeps the whole label.
