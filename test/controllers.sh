#!/bin/bash
#
# Regression test for the controller options (--auto-encode, --permissive,
# --fpim, --fpjm, --mealy-sat, --mealy-any, --dc, --compress). For every example
# and every combination of options, knor --real and knor -a must both exit with
# the expected verdict, and a realizable example must give an AIGER file with
# one input per uncontrollable and one output per controllable AP, which
# aigcheck accepts. Invalid combinations must exit with 255, and ABC must leave
# no abc.history behind.
#
# Usage: test/controllers.sh [path/to/knor [path/to/aigcheck]]
#        (default: build/knor and build/aigcheck)

dir=$(cd "$(dirname "$0")" && pwd)
knor=${1:-$dir/../build/knor}
aigcheck=${2:-$dir/../build/aigcheck}

for tool in "$knor" "$aigcheck"; do
    if ! [ -x "$tool" ]; then
        echo "not found: $tool"
        exit 2
    fi
done
knor=$(cd "$(dirname "$knor")" && pwd)/$(basename "$knor")

examples=(
    "$dir/../examples/EscalatorNonReactive.tlsf.ehoa 10"
    "$dir/../examples/lilydemo22.tlsf.ehoa 10"
    "$dir/../examples/full_arbiter.tlsf.ehoa 10"
    "$dir/../examples/full_arbiter_3.tlsf.ehoa 10"
    "$dir/../examples/amba_decomposed_arbiter_3.tlsf.ehoa 10"
    "$dir/../examples/ltl2dpa04.tlsf.ehoa 10"
    "$dir/../examples/arbiter.tlsf.ehoa 20"
    "$dir/controllers/no_outputs.ehoa 10"
)

modes=(
    "--dc"
    "--mealy-sat"
    "--mealy-sat --dc"
    "--permissive --mealy-sat"
    "--permissive --mealy-sat --dc"
    "--fpim --mealy-sat"
    "--fpim --permissive --mealy-sat --dc"
    "--auto-encode"
    "--auto-encode --fpjm"
    "--auto-encode --mealy-any --mealy-work 20000000"
    "--auto-encode --compress"
    "--sym --mealy-sat --dc"
    "--bisim --mealy-sat --dc"
    "--permissive --mealy-sat --mealy-any --mealy-work 20000000 --dc"
    # the ablation options, each in one of these
    "--auto-encode --no-relaxed-moves --no-monitor-codes --no-or-gate"
    "--auto-encode --no-cube-codes --no-adjacency-codes --no-rechoose --no-extension"
    "--auto-encode --no-renumber --no-output-codes --no-latch-ite --no-dead-latches"
    "--auto-encode --no-binary-dc --no-guard-complement --no-isop-candidate"
    "--auto-encode --no-onehot-dc --no-all-candidates --reorder-work 0 --gates-only --no-marginal-choice"
    "--fpim --permissive --mealy-sat --dc --no-layered --no-fewer-first --no-strong-bound"
    "--dc --no-unclarified"
)

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

failures=0
fail() {
    echo "FAIL: $1"
    failures=$((failures+1))
}

for example in "${examples[@]}"; do
    f=${example% *}
    expect=${example##* }
    aps=$(sed -n 's/^AP: \([0-9]*\).*/\1/p' "$f")
    outputs=$(sed -n 's/^controllable-AP://p' "$f" | wc -w | tr -d ' ')
    inputs=$((aps - outputs))
    for mode in "${modes[@]}"; do
        # shellcheck disable=SC2086
        "$knor" --real $mode "$f" >/dev/null 2>&1
        res=$?
        [ "$res" = "$expect" ] || fail "$f --real $mode (exit code $res, expected $expect)"
        # shellcheck disable=SC2086
        "$knor" -a $mode "$f" >"$tmp/out.aag" 2>/dev/null
        res=$?
        [ "$res" = "$expect" ] || fail "$f -a $mode (exit code $res, expected $expect)"
        if [ "$expect" = 10 ]; then
            read -r format _ i _ o _ < "$tmp/out.aag"
            if [ "$format" != aag ] || [ "$i" != "$inputs" ] || [ "$o" != "$outputs" ]; then
                fail "$f -a $mode (header \"$(head -1 "$tmp/out.aag")\", expected $inputs inputs and $outputs outputs)"
            elif ! out=$("$aigcheck" "$f" "$tmp/out.aag" 2>&1); then
                fail "$f -a $mode: $out"
            fi
        fi
    done
done

f="$dir/../examples/full_arbiter.tlsf.ehoa"
for mode in "--mealy-any" "--sym --mealy-sat --mealy-any" "--sym --permissive" "--sym --fpim" "--dc --onehot" "--dc --sop" "--auto-encode --onehot"; do
    # shellcheck disable=SC2086
    "$knor" -a $mode "$f" >/dev/null 2>&1
    res=$?
    [ "$res" = 255 ] || fail "$f $mode (exit code $res, expected 255)"
done

mkdir "$tmp/cwd"
(cd "$tmp/cwd" && "$knor" -a --compress "$f" >/dev/null 2>&1)
[ -e "$tmp/cwd/abc.history" ] && fail "$f -a --compress left abc.history in the working directory"

if [ $failures -ne 0 ]; then
    echo "$failures controller test(s) failed"
    exit 1
fi
echo "all controller tests passed"
