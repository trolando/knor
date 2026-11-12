#!/bin/bash
# SYNTCOMP synthesis: read eHOA from a file argument or stdin and write an
# ASCII AIGER controller. Forward all arguments and preserve Knor's exit codes
# (10 realizable, 20 unrealizable; other codes indicate errors).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
exec "$ROOT/scripts/knor.sh" --paritize --auto-encode --compress -a "$@"
