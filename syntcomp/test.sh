#!/bin/bash
#
# Test of the SYNTCOMP image: every configuration in /home/solvers (or in
# SOLVERS) on every eHOA file in DIR, through runlim as the competition runs
# it. The first line of the output must be REALIZABLE or UNREALIZABLE,
# matching the exit code (10 or 20), the configurations must agree, and with
# EXPECTED (lines "FILE CODE", as test/expected_verdicts.txt) the verdict
# must be the expected one.
#
# Usage: syntcomp-test DIR [EXPECTED]
# The limits are TIME_LIMIT seconds (default 60) and SPACE_LIMIT MB (default
# 1300, the memory of a single-core job in the competition).

dir=${1:?usage: syntcomp-test DIR [EXPECTED]}
expected=$2
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

failures=0
fail() {
    echo "FAIL: $1"
    failures=$((failures+1))
}

count=0
for f in "$dir"/*.ehoa; do
    name=$(basename "$f")
    want=
    [ -n "$expected" ] && want=$(awk -v n="$name" '$1 == n {print $2}' "$expected")
    first=
    for solver in "${SOLVERS:-/home/solvers}"/*; do
        config=$(basename "$solver")
        runlim -o "$tmp/info" --real-time-limit="${TIME_LIMIT:-60}" --space-limit="${SPACE_LIMIT:-1300}" \
            "$solver" "$f" > "$tmp/out" 2> /dev/null
        status=$(grep '^\[runlim\] status:' "$tmp/info" | cut -d: -f2- | xargs)
        rc=$(grep '^\[runlim\] result:' "$tmp/info" | awk '{print $3}')
        line=$(head -n 1 "$tmp/out")
        case "$rc:$line" in
            10:REALIZABLE | 20:UNREALIZABLE) ;;
            *) fail "$config $name: status \"$status\", exit code $rc, output \"$line\""; continue ;;
        esac
        [ -n "$want" ] && [ "$rc" != "$want" ] && fail "$config $name: exit code $rc, expected $want"
        [ -z "$first" ] && first=$rc
        [ "$rc" != "$first" ] && fail "$config $name: exit code $rc, another configuration gave $first"
    done
    count=$((count+1))
done

if [ $failures -ne 0 ]; then
    echo "$failures failure(s) on $count file(s)"
    exit 1
fi
echo "all configurations agree on $count file(s)"
