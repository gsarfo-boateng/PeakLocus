#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

GRAALVM_HOME="/Library/Java/JavaVirtualMachines/graalvm-21.jdk/Contents/Home"

if [[ ! -d "$GRAALVM_HOME" ]]; then
    echo "ERROR: GraalVM 21 not found:"
    echo "  $GRAALVM_HOME"
    exit 1
fi

export JAVA_HOME="$GRAALVM_HOME"
export PATH="$JAVA_HOME/bin:$PATH"

echo
echo " PeakLocus native build"
echo

echo "JAVA_HOME:"
echo "$JAVA_HOME"

echo
echo "Java:"
java --version

echo
echo "javac:"
javac --version

echo
echo "Native Image:"
native-image --version

echo
echo "Maven:"
mvn --version

echo
echo " Running complete PeakLocus build and regression tests"
echo

bash "$ROOT/build.sh"

echo
echo " Building GraalVM native executable"
echo

mvn clean package -Pnative

NATIVE="$ROOT/target/peaklocus"

if [[ ! -f "$NATIVE" ]]; then
    echo
    echo "ERROR: native executable was not produced:"
    echo "  $NATIVE"
    exit 1
fi

chmod +x "$NATIVE"

echo
echo " Native executable information"
echo

ls -lh "$NATIVE"
file "$NATIVE"

echo
echo " Native executable smoke tests"
echo

"$NATIVE" --version

echo
"$NATIVE" -h >/dev/null

echo
echo " PeakLocus native build PASSED"
echo
echo "Executable:"
echo "  $NATIVE"
echo
