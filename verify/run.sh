#!/usr/bin/env bash
# Run verify/compare.py on the holes listed in VERIFY.md.
# Usage: verify/run.sh <raw-data-root> <python-with-pytsg> <workdir>
# <raw-data-root>/<hole>/tsg/*_tsg.{tsg,bip} and <hole>/scalars/Grp1_sTSAS.csv
# must already be on disk (downloaded from the NVCL services; not in this repo).
set -euo pipefail
root=$1 python=$2 work=$3
here=$(cd "$(dirname "$0")/.." && pwd)
cargo build --release --manifest-path "$here/Cargo.toml" >&2
holes=(
  WA_05KCD001
  QLD_ChartersTowers1
  SA_103960
  r5/NT_1113632_DD85GL5
  r5/NSW_MIN_000610
  r2/WA_BDRCD019
)
for h in "${holes[@]}"; do
  tsg=$(find "$root/$h/tsg" -name '*_tsg.tsg' | head -1)
  "$python" -I "$here/verify/compare.py" "$here/target/release/tsg" "$tsg" \
    "$root/$h/scalars/Grp1_sTSAS.csv" "$work/${h//\//_}" | sed "s|^{|{\"hole\": \"$h\", |"
done
