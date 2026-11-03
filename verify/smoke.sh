#!/usr/bin/env bash
# Run `tsg info` and `tsg scalars` on every *_tsg.tsg under a directory and
# print one tab-separated line per file:
#   file  status  depth  backward_samples  message
# Usage: verify/smoke.sh <raw-data-root> <workdir>
set -euo pipefail
root=$1 work=$2
here=$(cd "$(dirname "$0")/.." && pwd)
cargo build --release --manifest-path "$here/Cargo.toml" >&2
tsg="$here/target/release/tsg"
mkdir -p "$work"
printf 'file\tstatus\tdepth\tbackward_samples\tmessage\n'
while IFS= read -r -d '' f; do
  rel=${f#"$root"/}
  if "$tsg" scalars "$f" -o "$work/smoke.csv" 2>"$work/err.txt" && "$tsg" info "$f" >"$work/info.txt"; then
    depth=$(grep '^depth:' "$work/info.txt" | sed 's/^depth: *//')
    back=$(grep -o 'at [0-9]* samples' "$work/info.txt" | grep -o '[0-9]*' || echo 0)
    printf '%s\tread\t%s\t%s\t\n' "$rel" "$depth" "$back"
  else
    printf '%s\trefused\t\t\t%s\n' "$rel" "$(tr '\n' ' ' <"$work/err.txt")"
  fi
done < <(find "$root" -name '*_tsg.tsg' -print0 | sort -z)
rm -f "$work/smoke.csv" "$work/err.txt" "$work/info.txt"
