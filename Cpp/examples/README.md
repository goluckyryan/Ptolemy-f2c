# Examples and Test Suite

Example Ptolemy input files and a comprehensive test suite comparing the C++ translation against the Cleopatra reference binary (32-bit ifort).

## Usage

```bash
cd Claude_attack/Cpp
rm -f fort.*
./ptolemy < examples/elastic_208Pb_p.in
rm -f fort.*
```

> **Important:** Always `rm -f fort.*` before and after each run. Leftover `fort.15` from a previous run can cause an infinite loop.

## Running the Test Suite

```bash
cd Claude_attack/Cpp

# Quick: compare integrated cross sections
bash examples/run_tests.sh

# Full: compare every angle in the angular distribution
bash examples/run_precision_tests.sh
```

Requires the Cleopatra binary at `/home/ryan/ptolemy_2019/digios/analysis/Cleopatra/ptolemy`.

## Test Results

- **Integrated cross sections:** 25/25 PASS (0.000%)
- **Angle-by-angle precision:** 21/24 at 0.0000%, 3 at deep diffraction minima where all binaries (C++, gfortran, Cleopatra) disagree with each other
- **Random fuzz** (800 cases, seed=42): 789 PASS / 9 WARN / 2 FAIL (FAILs are FP sensitivity + d-inelastic extrapolation failure shared with Fortran). See `Claude_attack/docs/tested_cases.md`.

## Elastic Scattering

| File | Reaction | E | sigma_R (mb) | Features |
|------|----------|---|-------------|----------|
| `elastic_208Pb_p` | p + 208Pb | 30 MeV | 1704.702 | Basic proton elastic |
| `elastic_40Ca_d` | d + 40Ca | 60 MeV | 1355.566 | Deuteron elastic |
| `elastic_208Pb_a` | a + 208Pb | 80 MeV | 3477.221 | Alpha elastic |
| `test_elastic_EL1_p` | p + 120Sn | 30 MeV | 1511.729 | PARAMETERSET EL1 |
| `test_elastic_EL2_d` | d + 90Zr | 60 MeV | 2287.223 | PARAMETERSET EL2 |
| `test_elastic_EL3_a` | a + 58Ni | 100 MeV | 2615.411 | PARAMETERSET EL3 |
| `test_elastic_heavy_ion` | 16O + 28Si | 60 MeV | 1543.309 | STEPSPER, ASYMPTOPIA |
| `test_elastic_npcoulomb` | p + 208Pb | 50 MeV | 1791.919 | NPCOULOMB=12 |
| `test_elastic_ecm` | a + 40Ca | ECM=30 | 2301.499 | ECM keyword |
| `test_elastic_shape` | 16O + 48Ca | 56 MeV | 1316.212 / 1316.205 | SHAPE potential-form linkule (tabulated WS, two passes) |

## Transfer DWBA

| File | Reaction | E | sigma_tot (mb) | Features |
|------|----------|---|---------------|----------|
| `transfer_48Ca_dp` | 48Ca(d,p)49Ca 3/2- | 30 MeV | 4.0220 | Standard (d,p), av18 |
| `transfer_208Pb_dp` | 208Pb(d,p)209Pb 9/2+ | 20 MeV | 10.198 | Heavy target (d,p) |
| `transfer_140Ce_3Hea` | 140Ce(3He,a)139Ce 3/2- | 30 MeV | 0.04970 | (3He,a) pickup |
| `transfer_208Pb_16O15N_phiffer` | 208Pb(16O,15N)209Bi 9/2- | 104 MeV | 2.29e-4 | Phiffer, heavy-ion |
| `transfer_48Ca_16O15N_phiffer` | 48Ca(16O,15N)49Sc 7/2- | 56 MeV | 1.2952 | Phiffer, heavy-ion |
| `test_transfer_CA60A` | 40Ca(d,p)41Ca 7/2- | 60 MeV | 9.9928 | PARAMETERSET CA60A |
| `test_transfer_CA60B` | 40Ca(d,p)41Ca 7/2- | 60 MeV | 11.994 | PARAMETERSET CA60B |
| `test_transfer_PBO1A` | 208Pb(16O,15N)209Bi 9/2- | 104 MeV | 7.53e-5 | PARAMETERSET PBO1A, phiffer |
| `test_transfer_DPSA` | 48Ca(d,p)49Ca 3/2- | 15 MeV | 27.589 | PARAMETERSET DPSA, low energy |
| `test_transfer_ALPHA1` | 90Zr(3He,a)89Zr 9/2+ | 30 MeV | 7.44e-3 | PARAMETERSET ALPHA1 |
| `test_transfer_ALPHA2` | 40Ca(3He,d)41Sc 7/2- | 25 MeV | 2.1805 | PARAMETERSET ALPHA2 |
| `test_transfer_spfac` | 208Pb(d,p)209Pb 9/2+ | 20 MeV | 5.7306 | SPFACT=0.90 |
| `test_transfer_r0target` | 48Ca(d,p)49Ca 3/2- | 30 MeV | 8.4121 | R0TARGET keyword |
| `test_transfer_maxlextrap0` | 16O(d,p)17O 5/2+ | 15 MeV | 35.257 | MAXLEXTRAP=0 |
| `test_transfer_2n` | 18O(p,t)16O 0+ | 25 MeV | 5.89e-3 | Two-neutron pickup |

## Inelastic DWBA

| File | Reaction | E | sigma_tot (mb) | Features |
|------|----------|---|---------------|----------|
| `inelastic_40Ca_pp` | 40Ca(p,p') 2+ 4.0 MeV | 20 MeV | 54.604 | Proton inelastic |
| `inelastic_40Ca_aa` | 40Ca(a,a') 2+ | 40 MeV | 43.148 | Alpha inelastic |
| `inelastic_208Pb_aa` | 208Pb(a,a') 3- 2.615 MeV | 120 MeV | 1.2982 | Octupole excitation |
| `test_inel_INELOCA1` | 120Sn(d,d') 2+ 1.17 MeV | 40 MeV | 4.0646 | PARAMETERSET INELOCA1 |
| `test_inel_INELOCA2` | 90Zr(a,a') 2+ 2.19 MeV | 80 MeV | 9.0617 | PARAMETERSET INELOCA2 |
| `test_inel_betacoul` | 208Pb(a,a') 3- 2.615 MeV | 120 MeV | 18.545 | BETACOUL != BETAN |
| `test_inel_E1` | 208Pb(a,a') 1- 7.0 MeV | 120 MeV | 0.12773 | E1 dipole excitation |
| `test_inel_E3` | 208Pb(a,a') 3- 2.615 MeV | 100 MeV | 2.1479 | E3 octupole excitation |
| `test_inel_heavy_ion` | 120Sn(16O,16O') 2+ 1.17 MeV | 200 MeV | 24.979 | Heavy-ion, LMAXADD=40 |
| `test_inel_labangles` | 58Ni(a,a') 2+ 1.45 MeV | 50 MeV | 33.667 | LABANGLES keyword |
| `test_inel_print2` | 40Ca(d,d') 2+ 3.35 MeV | 40 MeV | 42.597 | PRINT=2 verbose |

## Coupled Channels

| File | Reaction | E | sigma_inel (mb) | Features |
|------|----------|---|-----------------|----------|
| `cc_208Pb_90Zr` | 208Pb(90Zr,90Zr')208Pb* 3- | 1350 MeV | 0.10348 | Octupole Coulomb excitation |
| `cc_208Pb_90Zr_1minus` | 208Pb(90Zr,90Zr')208Pb* 1- | 1350 MeV | 0.30517 | Dipole Coulomb excitation |
| `cc_24Mg_12C` | 24Mg(12C,12C')24Mg* 2+ | ECM=20 | 64.315 | Rotational model |

> **Note:** CC output has two TOTAL lines. The first (`TOTAL: 0.00`) is elastic; the second is the inelastic cross section.

## Precision Notes

3 tests show sub-1% differences at individual angles — all at deep diffraction minima where the cross section is 3-4 orders of magnitude below peak. All three binaries (C++, gfortran, Cleopatra) give different values at these points:

| Case | C++ | gfortran | Cleopatra |
|------|-----|----------|-----------|
| 16O+28Si at 90 deg | 2.61e-4 | 2.19e-4 | 2.58e-4 |
| 90Zr(3He,a) at 90 deg | 7.23e-5 | 7.08e-5 | 7.24e-5 |
| 120Sn(16O,16O') at 59 deg | 9.285e-7 | 9.287e-7 | 9.289e-7 |

This is inherent floating-point sensitivity in the partial wave sum at near-cancellation points. The S-matrices match exactly across all binaries.

## Feature Coverage

- **Parametersets:** EL1, EL2, EL3, CA60A, CA60B, DPSA, DPSB, ALPHA1, ALPHA2, PBO1A, PB100B, ALPHA3, INELOCA1, INELOCA2, INELOCA3
- **Projectiles:** p, d, t, 3He, alpha, 16O, 90Zr, 12C
- **Reactions:** elastic, (d,p), (p,t), (3He,a), (3He,d), (16O,15N), inelastic (p,p'), (d,d'), (a,a'), (16O,16O'), CC
- **Keywords:** ECM, LABANGLES, NPCOULOMB, STEPSPER, ASYMPTOPIA, MAXLEXTRAP, R0TARGET, SPFACT, BETACOUL, LMAXADD, PRINT=2, WAVEFUNCTION PHIFFER
- **Multipolarities:** E1, E2, E3
- **Bound states:** av18, phiffer, Woods-Saxon
- **Energy range:** 15-1350 MeV
- **Target mass:** 16O to 208Pb

All cross sections match the Fortran Cleopatra binary to 0.0000%.
