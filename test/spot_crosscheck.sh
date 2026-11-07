#!/bin/bash
#
# Cross-check with Spot, which shares no code with knor or aigcheck. For the
# realizable examples with at most MAX_INPUTS inputs, the controller of knor -a
# must be correct according to Spot (aig2hoa.py, then autfilt --included-in),
# and on the controllers with one output negated, aigcheck and Spot must agree.
#
# Usage: test/spot_crosscheck.sh path/to/knor path/to/aigcheck path/to/autfilt path/to/python3

dir=$(cd "$(dirname "$0")" && pwd)
knor=$1
aigcheck=$2
autfilt=$3
python=$4
max_inputs=6
max_ands=2000

for tool in "$knor" "$aigcheck" "$autfilt" "$python"; do
    if ! [ -x "$tool" ]; then
        echo "not found: $tool"
        exit 2
    fi
done

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

failures=0
fail() {
    echo "FAIL: $1"
    failures=$((failures+1))
}

# verdict CODE: PASS (0), FAIL (1) or ERROR
verdict() {
    case $1 in
        0) echo PASS ;;
        1) echo FAIL ;;
        *) echo "ERROR ($1)" ;;
    esac
}

spot() {   # spot SPEC CONTROLLER: Spot's verdict on the controller
    if ! "$python" -I "$dir/aig2hoa.py" "$2" "$1" >"$tmp/c.hoa"; then
        echo "ERROR (aig2hoa)"
        return
    fi
    "$autfilt" -q "$tmp/c.hoa" --included-in="$1"
    verdict $?
}

controllers=0
mutants=0
wrong=0
while read -r name expect; do
    [ "$expect" = 10 ] || continue
    f="$dir/../examples/$name"
    "$knor" -a "$f" >"$tmp/c.aag" 2>/dev/null
    read -r _ _ i l o a < "$tmp/c.aag"
    [ "$i" -le $max_inputs ] && [ "$a" -le $max_ands ] || continue
    s=$(spot "$f" "$tmp/c.aag")
    [ "$s" = PASS ] || fail "$name: the controller of knor -a is $s according to Spot"
    controllers=$((controllers+1))
    for ((k = 0; k < o && k < 2; k++)); do
        awk -v line=$((i + l + 2 + k)) 'NR == line { $1 = $1 % 2 == 0 ? $1 + 1 : $1 - 1 } { print }' \
            "$tmp/c.aag" >"$tmp/m.aag"
        "$aigcheck" "$f" "$tmp/m.aag" >/dev/null 2>&1
        c=$(verdict $?)
        s=$(spot "$f" "$tmp/m.aag")
        [ "$c" = "$s" ] || fail "$name with output $k negated: aigcheck $c, Spot $s"
        mutants=$((mutants+1))
        [ "$c" = FAIL ] && wrong=$((wrong+1))
    done
done < "$dir/expected_verdicts.txt"

if [ $failures -ne 0 ]; then
    echo "$failures cross-check(s) failed"
    exit 1
fi
echo "Spot accepts all $controllers controllers and agrees with aigcheck on $mutants mutants ($wrong of them wrong)"
