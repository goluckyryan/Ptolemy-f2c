# Ptolemy DWBA Code — Working Directory

This repository contains the Fortran source of the Ptolemy DWBA
(Distorted Wave Born Approximation) code and a C++ transliteration.

## Directory Structure

```
fortran_copy/       Clean Fortran sources, verified against the original
                    Cleopatra ptolemy binary. Produces bit-identical
                    cross sections.
Cpp/                C++ transliteration of the full Ptolemy source (~61k
                    lines). Single-channel DWBA verified at 0.0000%
                    deviation across all angles.
reference/          Example input/output and InFileCreator reference code.
```

## Building the Fortran version

```bash
cd fortran_copy
make          # produces ptolemy.x (32-bit, gfortran -m32 -O0)
make clean
```

Requires: `gfortran`, `gcc` with 32-bit multilib support (`gcc-multilib`).

## Building the C++ version

```bash
cd Cpp
make          # produces ptolemy (64-bit, g++ -std=c++17 -O0)
make clean
```

Requires: `g++` with C++17 support.

## Running

```bash
# Fortran
cd fortran_copy
rm -f fort.*              # IMPORTANT: leftover fort.15 causes infinite loop
./ptolemy.x < tests/test01.in
./ptolemy.x < ../test_16O_dp.in

# C++
cd Cpp
rm -f fort.*
./ptolemy < ../test_16O_dp.in
```

## Verification

```bash
bash test_fortran_copy.sh
```

Compares `fortran_copy/ptolemy.x` output against the Cleopatra reference
binary on all test cases. Cross sections are bit-identical.

The C++ version produces 0.0000% deviation vs the Fortran reference across
all 181 angles for 11C(d,p)12C, 16O(d,p)17O, and 48Ca(d,p)49Ca reactions.
