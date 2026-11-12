#!/bin/bash
# Forward arguments and stdin to Knor, preserving its output and exit code.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="$ROOT/build/knor"

if [ ! -x "$BIN" ]; then
    echo "Error: $BIN not found. Run $ROOT/scripts/compile.sh first." >&2
    exit 3
fi

exec "$BIN" "$@"
