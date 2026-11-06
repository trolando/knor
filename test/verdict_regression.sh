#!/usr/bin/env bash
# Verdict regression: every examples/*.ehoa must produce the realizability
# exit code recorded in test/expected_verdicts.txt (10 REALIZABLE / 20
# UNREALIZABLE). This is the Phase-0 invariant — the backend must stay
# byte-for-byte identical on the shipped automata regardless of front-end work.
#
# Usage: test/verdict_regression.sh [KNOR_BIN]   (default: ./build/knor)
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
KNOR="${1:-./build/knor}"
EXPECTED="$HERE/expected_verdicts.txt"

if [[ ! -x "$KNOR" ]]; then echo "error: knor binary not found at $KNOR" >&2; exit 1; fi
if [[ ! -f "$EXPECTED" ]]; then echo "error: $EXPECTED not found" >&2; exit 1; fi

fail=0
while read -r name want; do
    [[ -z "$name" ]] && continue
    f="$HERE/../examples/$name"
    set +e
    "$KNOR" "$f" --real >/dev/null 2>&1
    got=$?
    set -e
    if [[ "$got" != "$want" ]]; then
        echo "MISMATCH $name: expected $want got $got"
        fail=$((fail+1))
    fi
done < "$EXPECTED"

total=$(grep -c . "$EXPECTED")
if [[ "$fail" -eq 0 ]]; then
    echo "verdict regression: all $total examples match"
else
    echo "verdict regression: $fail / $total MISMATCH"
fi
[[ "$fail" -eq 0 ]]
