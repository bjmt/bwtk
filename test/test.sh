#!/bin/bash
# Smoke test: bedGraph → bigWig → bedGraph round-trip.

set -uo pipefail

BWTK=$(cd "$(dirname "$0")/.." && pwd)/bwtk
TMPBW=$(mktemp /tmp/bwtk_smoke_XXXXX.bw)

trap 'rm -f "$TMPBW"' EXIT

echo "Running smoke test of bwtk output."

if ! "$BWTK" bg2bw \
  -g golden/basic.chrom.sizes \
  -i golden/basic.bedGraph \
  -o "$TMPBW" 2>&1; then
  echo "Test failed: bg2bw encountered an error."
  exit 1
fi

set -e
diff golden/basic.bedGraph.expected \
  <("$BWTK" adjust -B -i "$TMPBW" -o- 2>/dev/null) \
  > /tmp/bwtk_smoke_diff.txt

if [ -s /tmp/bwtk_smoke_diff.txt ]; then
  echo "Test failed, found the following diff:"
  cat /tmp/bwtk_smoke_diff.txt
  rm -f /tmp/bwtk_smoke_diff.txt
  exit 1
else
  rm -f /tmp/bwtk_smoke_diff.txt
  echo "Test succeeded, no changes in output."
  exit 0
fi
