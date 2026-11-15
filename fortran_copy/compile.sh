#!/bin/bash
#
# Compile script for Ptolemy DWBA code (gfortran)
#
# Usage:  ./compile.sh          (compile and link)
#         ./compile.sh clean    (remove .o and executable)
#
# Requirements: gfortran, gcc, and 32-bit multilib support
#               (Debian/Ubuntu: gcc-multilib gfortran-multilib)
#
# Compiler flags matter -- the source relies on 32-bit LOC()/pointer
# arithmetic, EQUIVALENCE type punning and default REAL being 8 bytes:
#   -m32                 REQUIRED. 64-bit builds give garbage (1e+37..1e+39)
#                        for many transfer/inelastic/CC inputs.
#   -fdefault-real-8     REQUIRED. The original build (compile-gfortran) used
#                        it; without it 13 of 62 test inputs blow up.
#   -fdefault-double-8   keeps DOUBLE PRECISION at 8 bytes (needed together
#                        with -fdefault-real-8 on modern gfortran).
#   -O0                  REQUIRED. -O3 breaks the type punning (spins read
#                        as 0/2 and the run aborts).
# With these flags ptolemy.x reproduces the Cleopatra 32-bit reference binary
# to full printed precision on all 62 inputs in Cpp/examples + transfer_tests.
#

set -e

EXEC="ptolemy.x"
FFLAGS="-c -m32 -g -O0 -fdefault-real-8 -fdefault-double-8 -w -std=legacy -fallow-argument-mismatch -fno-range-check -fallow-invalid-boz"
# srread.c is K&R-era C (no <unistd.h>); modern gcc rejects it without these
CFLAGS="-c -m32 -std=gnu89 -Wno-implicit-function-declaration"
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
# Link an explicit object list (a stray rcwfn.o duplicates RCWFN in fortlib)
OBJS=""
for src in $FORTRAN_SOURCES $C_SOURCES; do OBJS="$OBJS ${src%.*}.o"; done
$FC -m32 -g -o "$EXEC" $OBJS
if [ $? -ne 0 ]; then
    echo "ERROR: Linking failed"
    exit 1
fi

echo ""
echo "=== Build successful: $EXEC ==="
echo ""
echo "Run with:  ./$EXEC < tests/test01.in"
