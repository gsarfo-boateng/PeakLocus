#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
"$ROOT/build.sh"
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

java -jar "$ROOT/peaklocus-1.0.0.jar" \
  -i "$ROOT/examples/input/modkit_like.tsv" \
  -g "$ROOT/examples/reference/mini.gtf" \
  -o "$TMP/modkit.tsv" >"$TMP/stdout" 2>"$TMP/stderr"

grep -q "Strand column: 5" "$TMP/stderr"
grep -q $'Gene_ID\tGene_Name\tTranscript_IDs\tGene_Type\tTranscript_Type\tGene_Strand\tLocation' "$TMP/modkit.tsv"
grep -q "3'UTR" "$TMP/modkit.tsv"

diff -u "$ROOT/examples/expected/modkit_like.peaklocus.tsv" "$TMP/modkit.tsv"
echo "INTEGRATION TESTS PASSED"
