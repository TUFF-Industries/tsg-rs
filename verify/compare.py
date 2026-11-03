"""Check tsg-rs output against pytsg and a survey's published scalar CSV.

For one real .tsg/.bip pair this runs `tsg spectra` (.npy and .csv) and
`tsg scalars`, reads the same pair with pytsg's `read_spectra`, and reports:
spectra bit-identity (npy and csv), wavelength and depth agreement, every
scalar column against pytsg, and the decoded `Grp1 sTSAS` against the
survey CSV. Prints one JSON object. TSG files are untrusted input, so run
with `python -I`:

    python -I verify/compare.py <tsg binary> <x_tsg.tsg> <survey Grp1_sTSAS.csv> <workdir>
"""
import csv
import json
import statistics
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np
from pytsg import read_spectra

REPEATS = 3


def timed(fn) -> tuple[float, object]:
    """Median wall time of REPEATS calls, and the last result."""
    times, out = [], None
    for _ in range(REPEATS):
        t0 = time.perf_counter()
        out = fn()
        times.append(time.perf_counter() - t0)
    return statistics.median(times), out


def run(*args: str) -> None:
    subprocess.run(args, check=True, capture_output=True)


def bits_equal(a: np.ndarray, b: np.ndarray) -> bool:
    a32 = np.ascontiguousarray(a, dtype=np.float32)
    b32 = np.ascontiguousarray(b, dtype=np.float32)
    return a32.shape == b32.shape and np.array_equal(a32.view(np.uint32), b32.view(np.uint32))


def is_null(v: np.ndarray) -> np.ndarray:
    return np.isclose(v, np.finfo("float32").min)


def read_csv(path: Path) -> tuple[list[str], list[list[str]]]:
    with open(path, newline="") as f:
        rows = list(csv.reader(f))
    return rows[0], rows[1:]


def compare_scalars(sp, header: list[str], rows: list[list[str]]) -> dict:
    """Every band-header column: ours (CSV text) against pytsg's decode."""
    ours = {name: [r[j] for r in rows] for j, name in enumerate(header)}
    checked, bad = 0, []
    names = Counter(bh.name for bh in sp.bandheaders)
    for bh in sp.bandheaders:
        if names[bh.name] > 1:
            continue  # pytsg's DataFrame cannot tell duplicate names apart
        mine = ours[bh.name]
        ref = sp.scalars[bh.name]
        coded = bh.flag == 2 and isinstance(bh.class_number, int) and bh.class_number >= 0
        if coded and bh.class_number == 0:
            # pytsg leaves [class 0] columns as raw floats; decode with its own table.
            raw = np.asarray(ref, dtype=np.float32)
            ints = np.where(is_null(raw), -1, raw).astype(int)
            ref = sp.classes[0].map_ints(ints)
        if coded:
            ok = list(ref) == mine
        else:
            vals = np.asarray(ref, dtype=np.float32)
            null = is_null(vals)
            mine_null = np.array([m == "" for m in mine])
            parsed = np.array([np.float32(m) if m else np.float32(0) for m in mine], dtype=np.float32)
            ok = bool(np.array_equal(null, mine_null)) and bits_equal(parsed[~null], vals[~null])
        checked += 1
        if not ok:
            bad.append(bh.name)
    return {"scalar_columns_checked": checked, "scalar_columns_differing": bad,
            "duplicate_names_skipped": sorted(n for n, c in names.items() if c > 1)}


def compare_survey(header: list[str], rows: list[list[str]], survey: Path) -> dict:
    s_head, s_rows = read_csv(survey)
    mine = [r[header.index("Grp1 sTSAS")] for r in rows]
    if len(s_rows) != len(mine):
        return {"survey_rows": len(s_rows), "tsg_rows": len(mine), "survey_equal": None}
    theirs = ["" if r[2] == "null" else r[2] for r in s_rows]
    diff = Counter((t, m) for t, m in zip(theirs, mine) if t != m)
    depth = np.array([float(r[header.index("depth_m")]) for r in rows])
    start = np.array([float(r[0]) for r in s_rows])
    return {
        "survey_rows": len(s_rows),
        "survey_equal": sum(t == m for t, m in zip(theirs, mine)),
        "survey_non_null": sum(t != "" for t in theirs),
        "survey_differences": [f"{t or 'null'} vs {m or 'null'}: {n}" for (t, m), n in diff.most_common(5)],
        "max_abs_depth_diff_vs_survey_m": float(np.max(np.abs(depth - start))),
    }


def main() -> None:
    tsg_bin, tsg, survey, work = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4])
    work.mkdir(parents=True, exist_ok=True)
    npy, sp_csv, sc_csv = work / "spectra.npy", work / "spectra.csv", work / "scalars.csv"
    t_info, _ = timed(lambda: run(tsg_bin, "info", str(tsg)))
    t_npy, _ = timed(lambda: run(tsg_bin, "spectra", str(tsg), "-o", str(npy)))
    t_csv, _ = timed(lambda: run(tsg_bin, "spectra", str(tsg), "-o", str(sp_csv)))
    t_sc, _ = timed(lambda: run(tsg_bin, "scalars", str(tsg), "-o", str(sc_csv)))
    t_py, sp = timed(lambda: read_spectra(tsg, tsg.with_suffix(".bip")))

    ref = np.asarray(sp.spectra)
    ours = np.load(npy)
    c_head, c_rows = read_csv(sp_csv)
    csv_vals = np.array([[np.float32(v) for v in r[2:]] for r in c_rows], dtype=np.float32)
    wl_ours = np.array([float(w) for w in c_head[2:]])
    wl_ref = np.asarray(sp.wavelength, dtype=np.float64)
    s_head, s_rows = read_csv(sc_csv)
    depth_ref = np.array([float(d) for d in sp.sampleheaders["D"]])
    depth_ours = np.array([float(r[1]) for r in c_rows])

    out = {
        "file": tsg.name,
        "samples": int(ref.shape[0]),
        "bands": int(ref.shape[1]),
        "spectra_npy_bit_identical": bits_equal(ours, ref),
        "spectra_csv_bit_identical": bits_equal(csv_vals, ref),
        "wavelengths": "equal" if wl_ours.shape == wl_ref.shape and np.array_equal(wl_ours, wl_ref)
        else f"pytsg {wl_ref.shape[0]} values, max diff "
             f"{np.max(np.abs(wl_ours - wl_ref[:len(wl_ours)])) if len(wl_ref) >= len(wl_ours) else 'n/a'}",
        "depths_equal": bool(np.array_equal(depth_ours, depth_ref)),
        **compare_scalars(sp, s_head, s_rows),
        **compare_survey(s_head, s_rows, survey),
        "seconds_tsg_info": round(t_info, 3),
        "seconds_tsg_spectra_npy": round(t_npy, 3),
        "seconds_tsg_spectra_csv": round(t_csv, 3),
        "seconds_tsg_scalars_csv": round(t_sc, 3),
        "seconds_pytsg_read_spectra": round(t_py, 3),
        "bip_mb": round(tsg.with_suffix(".bip").stat().st_size / 1e6, 1),
    }
    print(json.dumps(out))


if __name__ == "__main__":
    main()
