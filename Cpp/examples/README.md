# Example Input Files

Example Ptolemy input files for elastic scattering, transfer DWBA, inelastic DWBA, and coupled-channels calculations.

## Usage

```bash
cd /path/to/Cpp
rm -f fort.*
./ptolemy < examples/elastic_208Pb_p.in
rm -f fort.*
```

> **Important:** Always `rm -f fort.*` before and after each run. Leftover `fort.15` from a previous run can cause an infinite loop.

## Elastic Scattering

| File | Reaction | E_lab | sigma_R (mb) |
|------|----------|-------|--------------|
| `elastic_208Pb_p.in` | p + 208Pb | 30 MeV | 1704.702 |
| `elastic_40Ca_d.in` | d + 40Ca | 60 MeV | 1355.566 |
| `elastic_208Pb_a.in` | alpha + 208Pb | 80 MeV | 3477.221 |

## Transfer DWBA

| File | Reaction | E_lab | sigma_tot (mb) |
|------|----------|-------|----------------|
| `transfer_48Ca_dp.in` | 48Ca(d,p)49Ca g.s. 3/2- | 30 MeV | 4.0220 |
| `transfer_208Pb_dp.in` | 208Pb(d,p)209Pb g.s. 9/2+ | 20 MeV | 10.198 |
| `transfer_140Ce_3Hea.in` | 140Ce(3He,a)139Ce g.s. 3/2- | 30 MeV | 0.04970 |

## Inelastic DWBA

| File | Reaction | E_lab | sigma_tot (mb) |
|------|----------|-------|----------------|
| `inelastic_40Ca_pp.in` | 40Ca(p,p') 2+ 4.0 MeV | 20 MeV | 54.604 |
| `inelastic_40Ca_aa.in` | 40Ca(a,a') 2+, beta=0.1 | 40 MeV | 43.148 |
| `inelastic_208Pb_aa.in` | 208Pb(a,a') 3- 2.615 MeV | 120 MeV | 1.2982 |

## Coupled Channels

| File | Reaction | E | sigma_inel (mb) | Notes |
|------|----------|---|-----------------|-------|
| `cc_208Pb_90Zr.in` | 208Pb(90Zr,90Zr')208Pb* 3- | 1350 MeV lab | 0.10348 | Octupole Coulomb excitation |
| `cc_208Pb_90Zr_1minus.in` | 208Pb(90Zr,90Zr')208Pb* 1- | 1350 MeV lab | 0.30517 | Dipole Coulomb excitation |
| `cc_24Mg_12C.in` | 24Mg(12C,12C')24Mg* 2+ | 20 MeV c.m. | 64.315 | Rotational model coupling |

> **Note:** The first `TOTAL: 0.00` in CC output is the elastic channel (no integrated elastic cross section is printed). The second `TOTAL:` is the inelastic channel cross section.

All cross sections match the Fortran binary to 0.0000%.
