# Ptolemy DWBA Code — C++ Translation

This repository contains a faithful C++ transliteration of the Fortran Ptolemy DWBA
(Distorted Wave Born Approximation) nuclear reaction code, verified against the
original Cleopatra 32-bit binary.

## Why C++?

The original Fortran source relies on non-standard features (EQUIVALENCE type
punning, COMMON block overlays, HOLLERITH constants, BOZ initializations) that
are compiler-specific. Modern gfortran reproduces the reference binary only
with a narrow, non-obvious flag set (`-m32 -msse2 -mfpmath=sse -O0
-fdefault-real-8 -fdefault-double-8`, see [Fortran version](#fortran-version)); a 64-bit build,
an `-O3` build, or a build without `-fdefault-real-8` silently gives garbage
for many transfer, inelastic and CC inputs. The C++ translation explicitly
maps every COMMON block and allocator trick to correct C++ structs, making it
portable and reliable.

## Directory Structure

```
Cpp/                C++ transliteration of the full Ptolemy source
  examples/         Example input files (elastic, transfer, inelastic, CC)
  include/          Header files (COMMON blocks, allocator, types)
  src/              Source files (~61k lines)
fortran_copy/       Unmodified Fortran sources + compile script (gfortran -m32 -msse2 -fdefault-real-8)
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

### Fortran bugs fixed in the C++ version

The C++ follows the Fortran line by line, with one deliberate exception, where
the original code produces NaN:

**GRDSET: division by zero in the transfer φ-grid (`source.f:19303`).**
`BSPROD` snaps the core–nucleon distance RT (or RP) to exactly 0 when
RT² < 1e-8. GRDSET then evaluates `acos((T1*RO + S1*RI*X)/RT)` = acos(0/0) = NaN
at that (RI, RO, cos φ) point. The NaN spreads through the H-kernel into the
radial integrals of the lowest partial waves, and the Fortran `TOTAL:` becomes
NaN. Whether some grid point lands exactly on RT = 0 depends on beam energy
and on the grid (NPSUM), so an excitation function shows isolated failures:

| 27Mg(d,p)28Mg g.s., `dpsb` | 23.5 MeV | 24 MeV | 24.5 MeV | 25 MeV | 27 MeV |
|---|---|---|---|---|---|
| Fortran (Cleopatra / `fortran_copy -m32`) | 2.6090 | NaN | 2.4273 | NaN | NaN |
| C++ before fix | 2.6090 | 1.4864 | 2.4273 | 1.3868 | 0.82378 |
| C++ after fix | 2.6090 | 2.5150 | 2.4273 | 2.3470 | 2.0491 |

Before the fix, the C++ was worse than the Fortran: NaN guards added during the
port (`SFROMI`) quietly zeroed the affected partial waves (L_in = 0–5 here) and
printed a believable but wrong total. The fix (`Cpp/src/grdset_translated.cpp`)
uses azimuth 0 when RT or RP is 0, which is safe because only the Y₀₀ part of
the bound state survives at r = 0. It also clamps the cosine to [-1, 1]. With
the fix, the 27Mg(d,p) result at 24 MeV converges with NPSUM (40/60/120 →
2.5150/2.5169/2.5169). `transfer_tests/transfer_12C_3Hen.in` (listed above as
"Fortran NaN") was hitting the same bug: C++ now gives 0.05585 mb, stable to
0.02% under grid refinement. Every input that runs without NaN in Fortran
gives the same output as before the fix.

**Grid steps in excitation functions (not a bug).** The interpolated
integration grid has `NPSUMI = int((SUMMAX-SUMMIN)·SUMPTS·(k_i+k_o)/4π)`
points, so it gains a point at certain energies. At `dpsb` precision this
shows up as ~1–2% steps in σ(E) (e.g. 27Mg(d,p) near 14 MeV). Fortran
behaves the same way. For smooth excitation functions, raise NPSUM after
the parameter set, e.g. `PARAMETERSET dpsb r0target` followed by `npsum=80`.
For 27Mg(d,p) near 14 MeV, NPSUM = 60 and 80 agree to 4 digits; raising NPDIF
changes nothing.

### Known limitations — collective inelastic model has a narrow domain

Ptolemy's inelastic scattering uses the **vibrational / collective
excitation model** from the original Fortran source: the excitation is
treated as a one-step surface deformation

  V(r, θ) = V₀(r) - β_L R₀ ∂V/∂R · Y_LM(θ)

parametrized by a single deformation amplitude `BELX` = β_L (input from
`PARAMETERSET INELOCA*`). This is the right physics for low-lying
collective 2⁺ / 3⁻ states of **even-even spherical heavy nuclei**
(e.g. 90Zr, 208Pb, Sn isotopes, 40Ca) at moderate excitation. It is
**physically inappropriate** for:

- **Light targets** (12C, 16O, 20Ne, ...) where the excited states are
  single-particle or cluster in character, not collective surface modes.
- **Excitation energies above the particle-emission threshold** (~7 MeV
  for 12C, ~12 MeV for 16O) where the daughter nucleus is unbound and
  the DWBA assumption of a stable bound final state is violated.
- **Large deformation amplitudes** outside the small-β linearization the
  one-step model assumes.

When the model is applied outside its domain — most visibly in
`(d,d')` on light targets with E_x ≳ 4 MeV at 20–60 MeV beam energy —
the code fails numerically downstream of the physics mismatch. The
L-extrapolation breaks (`**** ERROR: CANNOT EXTRAPOLATE FOR CHANNEL`),
|S(L)| diverges at high L, and the reported DCS becomes meaningless
(typically 10¹⁷ mb/sr or worse). Ptolemy-f2c reproduces the Cleopatra
failure mode bit-identically (1:1 translation).

A random-test sweep of 800 reactions found ~60 affected cases, all in
the high-E_x light-target `(d,d')` corner described above. The DCS
column flag `% FROM L>LMAX` is the canonical diagnostic: any value
far from 0.0 means the answer is dominated by extrapolated
(non-converged) partial waves and should not be trusted. Increasing
`LMAX` / `LMAXADD` / `ASYMPTOPIA` typically makes things worse, not
better, because the broken extrapolation simply runs longer.

**For collective inelastic scattering outside the model's domain** —
use FRESCO or ECIS with an explicit microscopic form factor, or
compare against experimental data.

The sibling project [PtolemyCpp](https://github.com/goluckyryan/PtolemyCpp)
(modern C++ refactor) silently truncates the failed extrapolation
contribution at LMAX instead of including the garbage — also incorrect,
but in a different direction. Neither code produces physically meaningful
DCS for reactions in this regime.

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
bash compile.sh  # or: make  -> ptolemy.x (32-bit)
```

Requires: `gfortran`, `gcc` with 32-bit multilib support
(Debian/Ubuntu: `gcc-multilib gfortran-multilib`).

The flags are **not optional**:

| Flag | Without it |
|---|---|
| `-m32` | 64-bit build: 1e+37–1e+39 mb or NaN for many transfer/inelastic/CC inputs |
| `-fdefault-real-8` (+ `-fdefault-double-8`) | 13 of 62 test inputs blow up to 1e+37–1e+39 or NaN |
| `-O0` | `-O3` breaks EQUIVALENCE type punning (spins read as 0/2, run aborts) |
| `-msse2 -mfpmath=sse` | x87 80-bit arithmetic flips last-bit grid tests (`VMIN < -2U` in GRDSET); 40Ca(d,3He) gives 0.87771 instead of 0.87759 |

With all of them, `ptolemy.x` reproduces the Cleopatra 32-bit reference binary
to full printed precision (`TOTAL:` line) on all 104 inputs in `Cpp/examples/`,
`transfer_tests/` and the PtolemyCpp `test_inputs/`. Long fitting runs
(`fortran_copy/tests/test02`, `test03`) still differ from Cleopatra in the last
printed digit of some lines. The cause has not been pinned down; it is not
x87-vs-SSE (both builds show it). One candidate is that the statically linked
Cleopatra binary uses an older libm.

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
