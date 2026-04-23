#!/bin/bash
#
# Compile script for Ptolemy DWBA code (gfortran)
#
# Usage:  ./compile.sh          (compile and link)
#         ./compile.sh clean    (remove .o and executable)
#
# Requirements: gfortran, gcc
#

set -e

EXEC="ptolemy.x"
FFLAGS="-c -g -O0 -w -std=legacy -fallow-argument-mismatch -fno-range-check -fallow-invalid-boz"
CFLAGS="-c"
FC="gfortran"
CC="gcc"

if [ "$1" = "clean" ]; then
    echo "Cleaning..."
    rm -f *.o "$EXEC"
    echo "Done."
    exit 0
fi

echo "=== Compiling Ptolemy ==="

# Fortran source files (order doesn't matter for compilation)
FORTRAN_SOURCES="
    av18.f
    fortlib.f
    gfortran_stuff.f
    keep.f
    keepsub.f
    keptsub.f
    linkule.f
    linkulesfitters.f
    masstable.f
    numbered_store.f
    phiffer.f
    ptolemy-main.f
    source.f
"

# C source files
C_SOURCES="
    dtime.c
    srread.c
"

# Compile Fortran sources
for src in $FORTRAN_SOURCES; do
    echo "  $FC $src"
    $FC $FFLAGS "$src"
    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to compile $src"
        exit 1
    fi
done

# Compile C sources
for src in $C_SOURCES; do
    echo "  $CC $src"
    $CC $CFLAGS "$src"
    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to compile $src"
        exit 1
    fi
done

# Link everything
echo "  Linking -> $EXEC"
$FC -g -o "$EXEC" *.o
if [ $? -ne 0 ]; then
    echo "ERROR: Linking failed"
    exit 1
fi

echo ""
echo "=== Build successful: $EXEC ==="
echo ""
echo "Run with:  ./$EXEC < tests/test01.in"
