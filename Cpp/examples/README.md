# Example Input Files

Example Ptolemy input files for elastic scattering, transfer DWBA, and inelastic DWBA calculations.

## Usage

```bash
cd /path/to/Cpp
rm -f fort.*
./ptolemy < examples/elastic_208Pb_p.in
rm -f fort.*
```

> **Important:** Always `rm -f fort.*` before and after each run. Leftover `fort.15` from a previous run can cause an infinite loop.

## Elastic Scattering

| File | Reaction | E_lab | σ_R (mb) | Time |
|------|----------|-------|----------|------|
| `elastic_208Pb_p.in` | p + ²⁰⁸Pb | 30 MeV | 1704.702 | 4 ms |
| `elastic_40Ca_d.in` | d + ⁴⁰Ca | 60 MeV | 1355.566 | 7 ms |
| `elastic_208Pb_a.in` | α + ²⁰⁸Pb | 80 MeV | 3477.221 | 4 ms |

## Transfer DWBA

| File | Reaction | E_lab | σ_tot (mb) | Time |
|------|----------|-------|------------|------|
| `transfer_48Ca_dp.in` | ⁴⁸Ca(d,p)⁴⁹Ca g.s. 3/2⁻ | 30 MeV | 4.0220 | 0.36 s |
| `transfer_208Pb_dp.in` | ²⁰⁸Pb(d,p)²⁰⁹Pb g.s. 9/2⁺ | 20 MeV | 10.198 | 1.10 s |
| `transfer_140Ce_3Hea.in` | ¹⁴⁰Ce(³He,α)¹³⁹Ce g.s. 3/2⁻ | 30 MeV | 0.04970 | 0.40 s |

## Inelastic DWBA

| File | Reaction | E_lab | σ_tot (mb) | Time |
|------|----------|-------|------------|------|
| `inelastic_40Ca_pp.in` | ⁴⁰Ca(p,p') 2⁺ 4.0 MeV, β=0.1 | 20 MeV | Coulomb-dominated | 15 ms |
| `inelastic_40Ca_aa.in` | ⁴⁰Ca(α,α') 2⁺, β=0.1 | 40 MeV | 43.148 | 16 ms |
| `inelastic_208Pb_aa.in` | ²⁰⁸Pb(α,α') 3⁻ 2.615 MeV, β=0.11 | 120 MeV | 1.2982 | 34 ms |

> **Note:** The proton inelastic `40Ca(p,p')` total cross section diverges at forward angles due to Coulomb excitation — the angular distribution at specific angles is valid.

## Optical Potential Parameterizations

- **Proton:** Becchetti-Greenlees type (V, R0, A, VSO, etc.)
- **Deuteron:** Lohr-Haeberli type
- **³He:** global parametrization
- **Alpha:** McFadden-Satchler type

All cross sections match the Fortran (Cleopatra) binary to 0.0000%.
