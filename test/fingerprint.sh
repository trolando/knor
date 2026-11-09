#!/bin/bash
#
# Fingerprints of knor's results, to check that a refactoring changes nothing:
# for every example of expected_verdicts.txt (or of a list of names), the exit
# code of knor -a with the options, a hash of the AIGER output, the SAT calls
# of the minimiser's searches (s strategy, e extended strategy, a all winning
# moves; per call the states and + suffice, - do not suffice, ? undecided,
# ! too large), and a hash of the statistics without the times (which also
# cover the work of the calls and the encoding candidates). Compare two runs
# with diff. The results are deterministic per platform, not across platforms
# (the order of BDD nodes in memory can break ties).
#
# Usage: test/fingerprint.sh path/to/knor OUTFILE [OPTIONS...]
#        EXAMPLES=file (names as in expected_verdicts.txt) limits the examples.
# Needs python3.

if [ $# -lt 2 ]; then
    echo "usage: $0 path/to/knor OUTFILE [OPTIONS...]"
    exit 2
fi
dir=$(cd "$(dirname "$0")" && pwd)
knor=$1
out=$2
shift 2
if ! [ -x "$knor" ] || ! command -v python3 >/dev/null; then
    echo "not found: $knor or python3"
    exit 2
fi

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

list=${EXAMPLES:-$dir/expected_verdicts.txt}
: > "$out"
while read -r name _; do
    rm -f "$tmp/stats.json"
    "$knor" -a "$@" --stats "$tmp/stats.json" "$dir/../examples/$name" < /dev/null > "$tmp/out.aag" 2>/dev/null
    rc=$?
    python3 - "$name" "$rc" "$tmp/out.aag" "$tmp/stats.json" >> "$out" <<'EOF'
import hashlib, json, os, sys

name, rc, aag, stats = sys.argv[1:5]

def untimed(value):
    if isinstance(value, dict):
        return {k: untimed(v) for k, v in value.items() if k not in ("times", "seconds")}
    if isinstance(value, list):
        return [untimed(v) for v in value]
    return value

def digest(data):
    return hashlib.sha1(data).hexdigest()[:16]

with open(aag, "rb") as f:
    output = digest(f.read())
trace, record = "-", "-"
if os.path.exists(stats):
    with open(stats) as f:
        statistics = untimed(json.load(f))
    outcomes = {"suffice": "+", "do not suffice": "-", "undecided": "?", "too large": "!"}
    trace = " ".join(search["moves"][0] + ":" + ",".join(str(k) + outcomes[outcome] for k, outcome, _ in search["calls"])
                     for search in statistics["searches"]) or "-"
    record = digest(json.dumps(statistics, sort_keys=True, separators=(",", ":")).encode())
print("\t".join([name, rc, output, trace, record]))
EOF
done < "$list"
