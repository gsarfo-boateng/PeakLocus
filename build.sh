#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
BUILD="$ROOT/build"
JAR="$ROOT/peaklocus-1.0.0.jar"

rm -rf "$BUILD"
mkdir -p "$BUILD/classes" "$BUILD/test-classes"

echo "Java compiler:"
javac -version

# Java 9+ supports --release.
# Java 8 uses -source/-target instead.
if javac --help 2>&1 | grep -q -- '--release'; then
    JAVA8_FLAGS=(--release 8)
else
    JAVA8_FLAGS=(-source 8 -target 8)
fi

echo "Compiler compatibility flags: ${JAVA8_FLAGS[*]}"
echo

echo "Compiling PeakLocus..."
find "$ROOT/src/main/java" -name '*.java' -print0 | xargs -0 javac "${JAVA8_FLAGS[@]}" -encoding UTF-8 -d "$BUILD/classes"

cat > "$BUILD/manifest.mf" <<MANIFEST
Manifest-Version: 1.0
Main-Class: org.peaklocus.PeakLocus
Implementation-Title: PeakLocus
Implementation-Version: 1.0.0
MANIFEST

echo "Building JAR..."
jar cfm "$JAR" "$BUILD/manifest.mf" -C "$BUILD/classes" .

echo "Compiling tests..."
find "$ROOT/src/test/java" -name '*.java' -print0 | xargs -0 javac  "${JAVA8_FLAGS[@]}" -encoding UTF-8 -cp "$BUILD/classes" -d "$BUILD/test-classes"

echo "Running tests..."
java -ea -cp "$BUILD/classes:$BUILD/test-classes" org.peaklocus.PeakLocusSelfTest
java -ea -cp "$BUILD/classes:$BUILD/test-classes" org.peaklocus.RepresentativeTranscriptRegressionTest

printf '\nBuilt successfully: %s\n' "$JAR"
