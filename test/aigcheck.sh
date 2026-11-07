#!/bin/bash
#
# Correctness test of the synthesized controllers: for every realizable
# example (see expected_verdicts.txt) and every realizable automaton in
# test/acceptance and test/controllers, the controller of knor -a OPTIONS must
# be accepted by aigcheck, the independent checker of controllers against
# their automaton. A controller with a negated output must be rejected, so that
# the checker is known to find errors. Checking takes time proportional to
# the size of the circuit times the number of transitions of its product with
# the automaton (about 6 minutes for the 131,805 AND gates of full_arbiter_8),
# so controllers with more than AIGCHECK_MAX_ANDS AND gates (default 20000)
# are skipped.
#
# Usage: test/aigcheck.sh [path/to/knor [path/to/aigcheck [OPTIONS...]]]
#        (default: build/knor, build/aigcheck, no options)

dir=$(cd "$(dirname "$0")" && pwd)
knor=${1:-$dir/../build/knor}
aigcheck=${2:-$dir/../build/aigcheck}
shift $(($# < 2 ? $# : 2))
options=("$@")
max_ands=${AIGCHECK_MAX_ANDS:-20000}

for tool in "$knor" "$aigcheck"; do
    if ! [ -x "$tool" ]; then
        echo "not found: $tool"
        exit 2
    fi
done

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

failures=0
checked=0
skipped=0
fail() {
    echo "FAIL: $1"
    failures=$((failures+1))
}

# check FILE [KNOR OPTIONS...]: the controller must be accepted by aigcheck
check() {
    local f=$1
    shift
    "$knor" -a "$@" "${options[@]}" "$f" >"$tmp/c.aag" 2>/dev/null
    local res=$?
    if [ "$res" != 10 ]; then
        fail "$f: knor -a $* ${options[*]} exited with $res, expected 10"
        return
    fi
    read -r _ _ _ _ _ ands < "$tmp/c.aag"
    if [ "$ands" -gt "$max_ands" ]; then
        skipped=$((skipped+1))
        return
    fi
    if ! out=$("$aigcheck" "$f" "$tmp/c.aag" 2>&1); then
        fail "$f: knor -a $* ${options[*]}: $out"
        return
    fi
    checked=$((checked+1))
}

while read -r name verdict; do
    [ "$verdict" = 10 ] && check "$dir/../examples/$name"
done < "$dir/expected_verdicts.txt"

for f in "$dir"/acceptance/*.hoa; do
    expect=$(sed -n '1s|^/\* expect: \(.*\) \*/$|\1|p' "$f")
    code=${expect%% *}
    acceptance_options=()
    [ "$expect" != "$code" ] && read -ra acceptance_options <<< "${expect#* }"
    [ "$code" = 10 ] && check "$f" "${acceptance_options[@]}"
done

check "$dir/controllers/no_outputs.ehoa"

# the checker must find the error in a controller with a negated output
f="$dir/../examples/full_arbiter.tlsf.ehoa"
"$knor" -a "$f" >"$tmp/c.aag" 2>/dev/null
read -r _ _ i l _ _ < "$tmp/c.aag"
awk -v line=$((i + l + 2)) 'NR == line { $1 = $1 % 2 == 0 ? $1 + 1 : $1 - 1 } { print }' "$tmp/c.aag" >"$tmp/negated.aag"
"$aigcheck" "$f" "$tmp/negated.aag" >/dev/null 2>&1
res=$?
[ "$res" = 1 ] || fail "$f with a negated output: aigcheck exited with $res, expected 1"

if [ $failures -ne 0 ]; then
    echo "$failures controller check(s) failed"
    exit 1
fi
echo "all $checked controllers passed aigcheck ($skipped with more than $max_ands AND gates skipped)"
