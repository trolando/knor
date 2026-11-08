#!/bin/bash
#
# Regression test for the eHOA parser.
# Every file in test/parser/accept must be parsed; every file in
# test/parser/reject must be rejected with a parse error (exit code 1).
# Inputs nested as deep as the parser allows (MAX_NESTING in
# src/simplehoa.cpp) must be parsed, one level deeper must be rejected; these
# are generated.
#
# Usage: test/parser.sh [path/to/knor]   (default: build/knor)

dir=$(cd "$(dirname "$0")" && pwd)
knor=${1:-$dir/../build/knor}

if ! [ -x "$knor" ]; then
    echo "knor not found at $knor"
    exit 2
fi

failures=0

accepted() {
    local err res
    err=$("$knor" --real "$1" 2>&1 >/dev/null)
    res=$?
    if [ $res -eq 1 ] || [[ $err == *"Parsing error"* ]]; then
        echo "FAIL (rejected): $1"
        echo "$err" | head -3
        failures=$((failures+1))
    fi
}

rejected() {
    local err res
    err=$("$knor" --real "$1" 2>&1 >/dev/null)
    res=$?
    if [ $res -ne 1 ] || [[ $err != *"Parsing error"* ]]; then
        echo "FAIL (accepted, exit code $res): $1"
        failures=$((failures+1))
    fi
}

for f in "$dir"/parser/accept/*.hoa; do accepted "$f"; done
for f in "$dir"/parser/reject/*.hoa; do rejected "$f"; done

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
repeat() { printf "%$2s" "" | tr ' ' "$1"; }
# deep D: an acceptance condition, a label with negations, and a label with
# parentheses, each nested D deep
deep() {
    printf 'HOA: v1\nStates: 1\nStart: 0\nAcceptance: 1 %sInf(0)%s\n--BODY--\nState: 0\n[t] 0\n--END--\n' \
        "$(repeat '(' "$1")" "$(repeat ')' "$1")" > "$tmp/acc_$1.hoa"
    printf 'HOA: v1\nStates: 2\nStart: 0\nAP: 2 "a" "b"\nAcceptance: 1 Inf(0)\n--BODY--\nState: 0\n[%s0] 1\nState: 1\n[t] 1\n--END--\n' \
        "$(repeat '!' "$1")" > "$tmp/not_$1.hoa"
    printf 'HOA: v1\nStates: 2\nStart: 0\nAP: 2 "a" "b"\nAcceptance: 1 Inf(0)\n--BODY--\nState: 0\n[%s0%s] 1\nState: 1\n[t] 1\n--END--\n' \
        "$(repeat '(' "$1")" "$(repeat ')' "$1")" > "$tmp/parens_$1.hoa"
}
limit=100000
deep $limit
deep $((limit + 1))
for f in "$tmp"/*_$limit.hoa; do accepted "$f"; done
for f in "$tmp"/*_$((limit + 1)).hoa; do rejected "$f"; done

if [ $failures -ne 0 ]; then
    echo "$failures parser test(s) failed"
    exit 1
fi
echo "all parser tests passed"
