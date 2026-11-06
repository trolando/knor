#!/bin/bash
#
# Regression test for the eHOA parser.
# Every file in test/parser/accept must be parsed; every file in
# test/parser/reject must be rejected with a parse error (exit code 1).
#
# Usage: test/parser.sh [path/to/knor]   (default: build/knor)

dir=$(cd "$(dirname "$0")" && pwd)
knor=${1:-$dir/../build/knor}

if ! [ -x "$knor" ]; then
    echo "knor not found at $knor"
    exit 2
fi

failures=0

for f in "$dir"/parser/accept/*.hoa; do
    err=$("$knor" --real "$f" 2>&1 >/dev/null)
    res=$?
    if [ $res -eq 1 ] || [[ $err == *"Parsing error"* ]]; then
        echo "FAIL (rejected): $f"
        echo "$err" | head -3
        failures=$((failures+1))
    fi
done

for f in "$dir"/parser/reject/*.hoa; do
    err=$("$knor" --real "$f" 2>&1 >/dev/null)
    res=$?
    if [ $res -ne 1 ] || [[ $err != *"Parsing error"* ]]; then
        echo "FAIL (accepted, exit code $res): $f"
        failures=$((failures+1))
    fi
done

if [ $failures -ne 0 ]; then
    echo "$failures parser test(s) failed"
    exit 1
fi
echo "all parser tests passed"
