#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
BUILD="$ROOT/build"
JAR="$ROOT/peaklocus-1.0.0.jar"
rm -rf "$BUILD"
mkdir -p "$BUILD/classes" "$BUILD/test-classes"
find "$ROOT/src/main/java" -name '*.java' -print0 | xargs -0 javac --release 8 -encoding UTF-8 -d "$BUILD/classes"
cat > "$BUILD/manifest.mf" <<MANIFEST
Manifest-Version: 1.0
Main-Class: org.peaklocus.PeakLocus
Implementation-Title: PeakLocus
Implementation-Version: 1.0.0
MANIFEST
jar cfm "$JAR" "$BUILD/manifest.mf" -C "$BUILD/classes" .
find "$ROOT/src/test/java" -name '*.java' -print0 | xargs -0 javac --release 8 -encoding UTF-8 -cp "$BUILD/classes" -d "$BUILD/test-classes"
java -ea -cp "$BUILD/classes:$BUILD/test-classes" org.peaklocus.PeakLocusSelfTest
printf '\nBuilt: %s\n' "$JAR"
