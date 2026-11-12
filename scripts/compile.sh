#!/bin/bash
# Build on the machine that will run Knor, using its native instruction set.
# Set KNOR_NATIVE=OFF for a portable build.
# Re-running this script updates the existing build incrementally.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"

NATIVE_MODE="${KNOR_NATIVE:-ON}"
case "$NATIVE_MODE" in
    ON) CPU_FLAGS="-march=native" ;;
    OFF) CPU_FLAGS="" ;;
    *) echo "Error: KNOR_NATIVE must be ON or OFF." >&2; exit 2 ;;
esac
echo "Native CPU optimizations: $NATIVE_MODE (override with KNOR_NATIVE=ON or OFF)."

cmake -S "$ROOT" -B "$ROOT/build" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_FLAGS="$CPU_FLAGS ${CFLAGS:-}" \
    -DCMAKE_CXX_FLAGS="$CPU_FLAGS ${CXXFLAGS:-}" \
    -DOINK_BUILD_NATIVE="$NATIVE_MODE" \
    -DLACE_NATIVE_OPT="$NATIVE_MODE" \
    -DSYLVAN_NATIVE_OPT="$NATIVE_MODE"
cmake --build "$ROOT/build" --target knor --parallel "${CMAKE_BUILD_PARALLEL_LEVEL:-$(getconf _NPROCESSORS_ONLN)}"

echo "Knor built at $ROOT/build/knor. Run $ROOT/scripts/knor.sh with your Knor options."
