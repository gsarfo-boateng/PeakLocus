#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
exec java -jar "$ROOT/peaklocus-1.0.0.jar" "$@"
