#!/bin/bash
#
# Regression test for the conversion of acceptance conditions to parity.
# Every file in test/acceptance starts with a comment "/* expect: CODE OPTIONS */":
# knor --real OPTIONS on the file must exit with CODE (10 realizable,
# 20 unrealizable, 104 needs --paritize).
#
# Usage: test/acceptance.sh [path/to/knor]   (default: build/knor)

dir=$(cd "$(dirname "$0")" && pwd)
knor=${1:-$dir/../build/knor}

if ! [ -x "$knor" ]; then
    echo "knor not found at $knor"
    exit 2
fi

failures=0
for f in "$dir"/acceptance/*.hoa; do
    expect=$(sed -n '1s|^/\* expect: \(.*\) \*/$|\1|p' "$f")
    code=${expect%% *}
    options=""
    [ "$expect" != "$code" ] && options=${expect#* }
    # shellcheck disable=SC2086
    "$knor" --real $options "$f" >/dev/null 2>&1
    res=$?
    if [ "$res" != "$code" ]; then
        echo "FAIL: $f (exit code $res, expected $code)"
        failures=$((failures+1))
    fi
done

if [ $failures -ne 0 ]; then
    echo "$failures acceptance test(s) failed"
    exit 1
fi
echo "all acceptance tests passed"
