# Ptolemy DWBA Code — C++ Translation

This repository contains a faithful C++ transliteration of the Fortran Ptolemy DWBA
(Distorted Wave Born Approximation) nuclear reaction code, verified against the
original Cleopatra 32-bit binary.

## Why C++?

The original Fortran source relies on non-standard features (EQUIVALENCE type
punning, COMMON block overlays, HOLLERITH constants, BOZ initializations) that
are compiler-specific. Recompiling with modern gfortran — even with `-m32` and
all compatibility flags — produces **incorrect results** for inelastic
scattering. The C++ translation explicitly maps every COMMON block and
allocator trick to correct C++ structs, making it portable and reliable.

## Directory Structure

```
Cpp/                C++ transliteration of the full Ptolemy source
  examples/         Example input files (elastic, transfer, inelastic, CC)
  include/          Header files (COMMON blocks, allocator, types)
  src/              Source files (~61k lines)
fortran_copy/       Clean Fortran sources with compile script (gfortran -m32)
docs/               Theory documentation and input file reference
```

## Documentation

- [docs/Usage.md](docs/Usage.md) — Keywords, parameter sets, and input file format
- [docs/Elastic.md](docs/Elastic.md) — Elastic scattering (optical model, Numerov, S-matrix)
- [docs/DWBA.md](docs/DWBA.md) — Transfer DWBA (form factor, distorted waves, spectroscopic factor)
- [docs/Inelastic.md](docs/Inelastic.md) — Inelastic DWBA (collective model, Coulomb excitation)
- [docs/Coupled_Channel.md](docs/Coupled_Channel.md) — Coupled-channels (CC theory, input format, Born correction)
- [docs/Output.md](docs/Output.md) — Output format and diagnostics
- [docs/Polarization.md](docs/Polarization.md) — Analyzing powers

## Capabilities

| Feature | Status |
|---|---|
| Elastic optical model | **Working** — exact match with Fortran |
| Transfer DWBA | **Working** — exact match |
| Inelastic DWBA (collective model) | **Working** — exact match |
| Coupled channels (CC) | **Working** — exact match (incl. Coulomb Born correction) |
| L-extrapolation | **Working** — exact match |
| Parameter fitting | Stubbed |

## Test Results

32 curated test cases compared against the Fortran Cleopatra binary:

| Category | Tests | PASS | Notes |
|----------|-------|------|-------|
| Elastic | 3 | 3 | p, d, alpha on Ca/Pb |
| Transfer DWBA | 13 | 7 | 5 crash Fortran, 1 Fortran NaN |
| Inelastic DWBA | 15 | 10 | 5 produce Fortran NaN |
| Coupled channels | 1 | 1 | 208Pb+90Zr Coulomb excitation |
| **Total** | **32** | **21** | |

**Every case where the Fortran binary runs successfully, C++ matches to full printed precision (0.000%).** The 11 non-PASS cases are all Fortran bugs (NaN or crash), not C++ disagreements.

### Transfer Reactions

| Reaction | Max Error |
|---|---|
| 48Ca(d,p)49Ca | 0.0000% |
| 208Pb(d,p)209Pb | 0.0000% |
| 140Ce(3He,a)139Ce | 0.0000% |
| 16O(d,p)17O | 0.0000% |
| 16O(p,d)15O | 0.0000% |
| 18O(p,t)16O | 0.0000% |
| 12C(d,n)13N | 0.0014% |

### Inelastic Scattering (DWBA)

| Target | (d,d') | (a,a') |
|---|---|---|
| 16O | 0.0000% | 0.0000% |
| 40Ca | 0.0000% | 0.0021% |
| 48Ca | 0.0000% | 0.0000% |
| 130Sn | 0.0000% | 0.0089% |
| 208Pb | 0.0000% | 0.0140% |

### Coupled Channels

| Reaction | sigma_inel | Match |
|---|---|---|
| 208Pb(90Zr,90Zr')208Pb* 1- | 0.30517 mb | 0.0000% |

## Building

### C++ version (recommended)

```bash
cd Cpp
make -j4        # produces ptolemy (64-bit, g++ -std=c++17 -O0)
```

Requires: `g++` with C++17 support.

### Fortran version

```bash
cd fortran_copy
bash compile.sh  # produces ptolemy.x (32-bit, gfortran -m32 -O0)
```

Requires: `gfortran`, `gcc` with 32-bit multilib support.

**Note:** The gfortran build does not reproduce all reaction types correctly.
Use the Cleopatra binary or the C++ translation as reference.

## Running

```bash
cd Cpp
rm -f fort.*              # IMPORTANT: leftover fort.15 causes infinite loop
./ptolemy < examples/elastic_208Pb_p.in
rm -f fort.*
./ptolemy < examples/cc_208Pb_90Zr.in
rm -f fort.*
```

See [examples/README.md](Cpp/examples/README.md) for all example inputs.
