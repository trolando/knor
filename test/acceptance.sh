#!/bin/bash
#
# Regression test for the conversion of acceptance conditions to parity.
# Every file in test/acceptance starts with a comment "/* expect: CODE OPTIONS */":
# knor --real OPTIONS on the file must exit with CODE (10 realizable,
# 20 unrealizable, 104 needs acceptance memory forbidden by --no-paritize).
# The opt-out must reject conditions that need memory and still accept
# parity conditions and conversions that only relabel the existing edges.
#
# Usage: test/acceptance.sh [path/to/knor]   (default: build/knor)

dir=$(cd "$(dirname "$0")" && pwd)
knor=${1:-$dir/../build/knor}

if ! [ -x "$knor" ]; then
    echo "knor not found at $knor"
    exit 2
fi

failures=0
check() {
    local code=$1 f=$2 res
    shift 2
    "$knor" --real "$@" "$f" >/dev/null 2>&1
    res=$?
    if [ "$res" != "$code" ]; then
        echo "FAIL: $f $* (exit code $res, expected $code)"
        failures=$((failures+1))
    fi
}

for f in "$dir"/acceptance/*.hoa; do
    expect=$(sed -n '1s|^/\* expect: \(.*\) \*/$|\1|p' "$f")
    code=${expect%% *}
    options=""
    [ "$expect" != "$code" ] && options=${expect#* }
    # shellcheck disable=SC2086
    check "$code" "$f" $options
done

for name in genbuchi_alternate genbuchi_env rabin2 rabin2_env; do
    check 104 "$dir/acceptance/$name.hoa" --no-paritize
done
for name in buchi buchi_state cobuchi_env parity_permuted true false; do
    f="$dir/acceptance/$name.hoa"
    code=$(sed -n '1s|^/\* expect: \([0-9]*\).*|\1|p' "$f")
    check "$code" "$f" --no-paritize
done

if [ $failures -ne 0 ]; then
    echo "$failures acceptance test(s) failed"
    exit 1
fi
echo "all acceptance tests passed"
