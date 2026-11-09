#!/bin/bash
#
# Byte regression on a few small examples: for plain knor -a and for
# --auto-encode, the exit code, the hash of the AIGER output and the SAT calls
# of the minimiser (see fingerprint.sh) must be those in fingerprints.txt. The
# examples are those whose results are the same on macOS (arm64) and Linux
# (x86_64); in general the results are deterministic per platform only (the
# order of BDD nodes in memory can break ties), so check the controllers
# (test controllers) before suspecting a regression on another platform.
# Update the file only for a deliberate change, explained in the commit:
# test/fingerprints.sh path/to/knor --update
#
# Usage: test/fingerprints.sh path/to/knor [--update]
# Needs python3.

dir=$(cd "$(dirname "$0")" && pwd)
knor=${1:-$dir/../build/knor}
expected=$dir/fingerprints.txt

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

: > "$tmp/actual"
for mode in plain auto-encode; do
    awk -F'\t' -v m=$mode '$1 == m {print $2}' "$expected" > "$tmp/examples"
    opts=()
    [ $mode = auto-encode ] && opts=(--auto-encode)
    EXAMPLES="$tmp/examples" "$dir/fingerprint.sh" "$knor" "$tmp/out" "${opts[@]}" || exit 2
    cut -f1-4 "$tmp/out" | sed "s/^/$mode	/" >> "$tmp/actual"
done

if [ "$2" = --update ]; then
    cp "$tmp/actual" "$expected"
    echo "updated $expected"
    exit 0
fi
if ! diff "$expected" "$tmp/actual"; then
    echo "the fingerprints differ (< expected, > actual)"
    exit 1
fi
echo "all fingerprints as expected"
