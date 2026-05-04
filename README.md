# Ptolemy DWBA Code — C++ Translation

This repository contains the Fortran source of the Ptolemy DWBA
(Distorted Wave Born Approximation) nuclear reaction code and a faithful
C++ transliteration, verified against the original Cleopatra 32-bit binary.

## Why C++?

The original Fortran source relies on non-standard features (EQUIVALENCE type
punning, COMMON block overlays, HOLLERITH constants, BOZ initializations) that
are compiler-specific. Recompiling with modern gfortran — even with `-m32` and
all compatibility flags — produces **incorrect results** for inelastic
scattering. The C++ translation explicitly maps every COMMON block and
allocator trick to correct C++ structs, making it portable and reliable.

## Directory Structure

```
Cpp/                C++ transliteration of the full Ptolemy source (~61k lines)
fortran_copy/       Clean Fortran sources with Makefile (gfortran -m32)
docs/               Usage guide, keywords, and input file examples
reference/          Example input/output and InFileCreator reference code
inel_tests/         Inelastic scattering test inputs (15 cases)
transfer_tests/     Transfer reaction test inputs (12 cases)
```

## Documentation

See [docs/Usage.md](docs/Usage.md) for a complete reference of all Ptolemy
keywords, control switches, parameter sets, and input file examples for every
supported reaction type.

## Benchmark Results

All benchmarks compare the C++ translation (64-bit, g++) against the
Cleopatra reference binary (32-bit, ifort). Max % error in differential
cross sections across all angles.

### Transfer Reactions (DWBA)

| Reaction | Type | Max Error |
|---|---|---|
| 16O(d,p)17O | 1n stripping | 0.0000% |
| 16O(p,d)15O | 1n pickup | 0.0000% |
| 16O(d,t)15O | 1n pickup | 0.0046% |
| 16O(d,3He)15N | 1p pickup | 0.0012% |
| 12C(d,n)13N | 1p stripping | 0.0014% |
| 18O(p,t)16O | 2n pickup | 0.0000% |
| 16O(t,p)18O | 2n stripping | 0.0000% |
| 12C(3He,n)14O | 2p stripping | 0.0000% |
| 16O(12C,13C)15O | heavy-ion 1n | 0.0000% |
| 12C(8Li,a)16N | heavy-ion 4n | 0.0000% |

### Elastic Scattering

| Reaction | Max Error |
|---|---|
| 48Ca + 16O at 56 MeV | 0.0000% |

### Inelastic Scattering (DWBA, collective model)

| Target | (p,p') 20 MeV | (d,d') 40 MeV | (a,a') 80 MeV |
|---|---|---|---|
| 16O | 0.0000% | 0.0000% | 0.0000% |
| 40Ca | 0.0000% | 0.0000% | 0.0021% |
| 48Ca | 0.0000% | 0.0000% | 0.0000% |
| 130Sn | 0.0000% | 0.0000% | 0.0089% |
| 208Pb | 0.0000% | 0.0000% | 0.0140% |

All 2+ excited states at 4 MeV, beta = 0.10, PARAMETERSET INELOCA3.

### Comparison: C++ vs gfortran recompilation

| Binary | Transfer | Elastic | Inel (p,p') | Inel (a,a') |
|---|---|---|---|---|
| **C++ (g++ 64-bit)** | 0.0000% | 0.0000% | 0.0000% | 0.014% |
| gfortran -m32 | 0.0000% | OK | OK | **BROKEN** |
| gfortran 64-bit | 0.0000% | **BROKEN** | OK | **BROKEN** |

## Building

### C++ version (recommended)

```bash
cd Cpp
make          # produces ptolemy (64-bit, g++ -std=c++17 -O0)
```

Requires: `g++` with C++17 support.

### Fortran version

```bash
cd fortran_copy
make          # produces ptolemy.x (32-bit, gfortran -m32 -O0)
```

Requires: `gfortran`, `gcc` with 32-bit multilib support
(`gcc-multilib`, `gfortran-multilib`).

**Note:** The gfortran build does not reproduce all reaction types correctly.
Use the Cleopatra binary or the C++ translation as reference.

## Running

```bash
cd Cpp
rm -f fort.*              # IMPORTANT: leftover fort.15 causes infinite loop
./ptolemy < ../test_16O_dp.in
./ptolemy < ../inel_tests/inel_208PB_aa.in
```

## Verification

```bash
# Inelastic test suite (15 cases, all < 0.1%)
bash test_inelastic.sh

# Fortran copy vs Cleopatra (transfer + bound state tests)
bash test_fortran_copy.sh
```
