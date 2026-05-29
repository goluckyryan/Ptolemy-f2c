# Ptolemy C++ Translation

Verbatim transliteration of the Fortran DWBA nuclear reaction code Ptolemy (April 2007 version) to C++17.
The translation preserves the original Fortran control flow, variable names, and COMMON block structure.

## Status

| Capability | State |
|---|---|
| Elastic optical model (S-matrix) | **Working** — exact match with Fortran |
| DWBA transfer cross sections | **Working** — exact match |
| Inelastic DWBA (collective model) | **Working** — exact match |
| Coupled-channels (CC) mode | **Working** — exact match (Coulomb Born correction included) |
| L-extrapolation (LXTRP1/2/M, LINLSQ) | **Working** — exact match |
| Parameter fitting engine | Stubbed |
| Analyzing powers (ANAPOW) | Stubbed |

## Test Results

Compared against the Fortran Cleopatra binary across 32 curated test cases:

| Category | Tests | PASS | Notes |
|----------|-------|------|-------|
| Elastic | 3 | 3 | p, d, alpha on Ca/Pb |
| Transfer DWBA | 13 | 7 | 5 crash Fortran (rc=206), 1 Fortran NaN |
| Inelastic DWBA | 15 | 10 | 5 produce Fortran NaN (proton Coulomb cases) |
| Coupled channels | 1 | 1 | 208Pb+90Zr Coulomb excitation |
| **Total** | **32** | **21** | **11 are Fortran bugs, not C++ disagreements** |

**Every case where the Fortran binary runs successfully, C++ matches to full printed precision (0.000%).**

The 11 non-PASS cases are all Fortran bugs:
- 5 inelastic proton cases produce `NaN` in Fortran (known Cleopatra bug)
- 5 transfer cases crash Fortran with rc=206 (array bounds violation)
- 1 transfer case produces `NaN` in Fortran

## Build

```bash
cd Claude_attack/Cpp
make -j4
rm -f fort.*
./ptolemy < examples/elastic_208Pb_p.in
rm -f fort.*
```

Requires g++ with C++17, compiled with `-O0 -g`.

> **Important:** Always `rm -f fort.*` before and after each run. Leftover `fort.15` from a previous run causes an infinite loop.

## Examples

See [examples/README.md](examples/README.md) for a full list of example input files covering elastic, transfer DWBA, inelastic DWBA, and coupled-channels calculations.

## Documentation

- [Usage.md](../docs/Usage.md) — general usage and input format
- [Elastic.md](../docs/Elastic.md) — elastic scattering theory
- [DWBA.md](../docs/DWBA.md) — transfer DWBA theory
- [Inelastic.md](../docs/Inelastic.md) — inelastic DWBA theory
- [Coupled_Channel.md](../docs/Coupled_Channel.md) — coupled-channels theory and input format
- [Output.md](../docs/Output.md) — output format description
- [Polarization.md](../docs/Polarization.md) — analyzing powers

## Calculation Flow

### DWBA

```
CONTRL (main control loop)
  |
  +--> DEFALT (initialize defaults)
  +--> PARAM (load parameter set: DPSB etc.)
  +--> REACTN (parse A(d,p)B reaction)
  |
  +--> SETCHN + SETPOT (set up each channel)
  |     |
  |     +--> PROJECTILE: BOUND -> av18 linkule (deuteron WF) or WOODSX+Numerov
  |     +--> TARGET:     BOUND -> WOODSX + Numerov integration
  |     +--> INCOMING:   WAVSET -> SETSPT (Coulomb) + MAKPOT (optical potential)
  |     +--> OUTGOING:   WAVSET -> SETSPT (Coulomb) + MAKPOT (optical potential)
  |
  +--> PRBPRT (print reaction summary)
  +--> GETSCT (elastic scattering waves: WAVELJ per L, TMATCH, JPTOLX)
  +--> GRDSET (3D integration grids: CUBMAP + BSPROD)
  +--> ANGSET (angular momentum coupling tables)
  +--> INELDC (DWBA radial integrals: A12 + WFGET + SFROMI)
  +--> LINTRP (L-interpolation/extrapolation: INTRCF)
  +--> XSECTN (cross sections: ELDCS + BETCAL + AMPCAL)
```

### Coupled Channels

```
CONTRL
  |
  +--> SETCHN (channels + basis states + coupling definitions)
  +--> WAVSET + MAKPOT (optical potentials for all channels)
  +--> GETSCT (distorted waves)
  +--> COULST (Coulomb integrals + B-factors via SETBFC/COULIN)
  +--> GRDSET (integration grids)
  +--> INRDIN (CC radial integrals, iterative S-matrix solve)
  +--> FFACST (Coulomb Born S-matrix correction)
  |     +--> SETBFC (FF2 Coulomb form factor integrals)
  |     +--> SETBRN (Born subtraction from CC S-matrices)
  +--> LINTRP (L-extrapolation with Born add-back at L960)
  +--> XSECTN (cross sections)
```

## File Organization

```
src/
  ptolemy_main.cpp              Main program dispatch loop
  contrl_translated.cpp         CONTRL — keyword parsing and main dispatch
  defalt_translated.cpp         DEFALT — defaults initialization
  param_reactn_gsinfo_translated.cpp   PARAM, REACTN, GSINFO
  setchn_translated.cpp         SETCHN — channel setup
  chanel_translated.cpp         CHANEL
  clrchn_translated.cpp         CLRCHN — clear channel
  setfg_translated.cpp          SETFG
  bound_translated.cpp          BOUND — bound state Numerov integration
  bound_helpers_translated.cpp  BSPROD, BSSET, CLINTS, JPTOLX
  wavset_translated.cpp         WAVSET, WAVPOT — scattering wave setup
  source_potentials.cpp         SETSPT, MAKPOT, WOODSX, DEFWOO, FFACST, SFROMI, BEBETA
  rcwfn.cpp                     RCWFN — Coulomb wave functions (F, G)
  coulst_translated.cpp         COULST — Coulomb integrals, SETBFC, B-factors
  scattering_partA_translated.cpp   COULIN, COULNG, GETSCT area
  scattering_part2_translated.cpp   TMATCH, LCRITL, SWKB
  scattering_part3_translated.cpp   WAVEF, WAVINH, WFGET
  wavelj_translated.cpp         WAVELJ — Numerov + S-matrix matching
  wavetc_translated.cpp         WAVETC — tensor-coupled channel waves
  prbprt_translated.cpp         PRBPRT — reaction summary print
  grdset_translated.cpp         GRDSET — 3D integration grid setup
  angset_translated.cpp         ANGSET — angular momentum coupling tables
  ineldc_translated.cpp         INELDC — core DWBA transfer integrals
  ineld2_translated.cpp         INELD2 — CC radial integrals
  ingrst_translated.cpp         INGRST — form factor grid setup
  inrdin_translated.cpp         INRDIN — CC radial integral driver
  a12_translated.cpp            A12 — angular momentum transformation
  lintrp_translated.cpp         LINTRP — L-interpolation/extrapolation + Born add-back
  xsectn_translated.cpp         XSECTN — cross section computation
  ampcal_translated.cpp         AMPCAL — transition amplitudes per angle
  angular_part2_translated.cpp  BETCAL, CALANG, ELDCS
  source_angular.cpp            ANAPOW (stub), CALFUN, CALGRA (stubs)
  source_channels.cpp           BASCPL, CCHAN, CCOUP, CTRAN
  source_coupled.cpp            CCDUMP, CCMTCH, CCONV, COUPLN, DRIVE
  source_fitting.cpp            SETBRN, GENBNX, GETBFC, FITEL (stubs), LXTRP1/2, etc.
  source_scattering.cpp         Scattering helpers
  source_io.cpp                 NXINT, NXWORD, NEWCD, NXVAL, MSCAN, DATAIN, etc.
  source_misc.cpp               INTRCF, CUBMAP, AITLAG, EPSLON, MUELCO, CUPSPN, etc.
  source_bound.cpp              Bound state helpers
  numbered_store.cpp            NALLOC, IALLOC, NAMLOC — Fortran-style allocator
  fortlib_part1.cpp             CLEBSH, THREEJ, RACAH, SIXJ, WIG9J, DSGMAL, etc.
  fortlib_part2.cpp             SPLNCB, INTRPC, LIN, MATINV, GAUSSL, etc.
  masstable.cpp                 EXCESS, MASEXX, AZCODE — nuclear mass table
  av18.cpp                      Argonne v18 deuteron wave function linkule
  linkule.cpp                   LINKUL — linkule dispatcher
  linkulesfitters_fitters.cpp   Optimizer linkules (DLSMIN, LMCHOL, etc.) — stubs
  linkulesfitters_potentials.cpp  Potential linkules
  linkulesfitters_wavefn.cpp    Wave function linkules
  phiffer.cpp                   PHIFFER linkule
  phsprt_translated.cpp         PHSPRT — phase shift printing
  rcasym_translated.cpp         RCASYM — asymptotic Coulomb expansion
  baslbl_translated.cpp         BASLBL — basis state labels
  keep.cpp, keepsub.cpp, keptsub.cpp   KEEP/KEPT — Speakeasy object I/O
  srread.cpp                    SRREAD — string reader
  dtime.cpp                     CPU timing
  gfortran_stuff.cpp            gfortran runtime compatibility shims

include/
  ptolemy_types.h               Type aliases (real8, integer, char8, etc.)
  ptolemy_commons.h             All COMMON block structs + 1-based accessors
  ptolemy_alloc.h               ALLOC/ILLOC/ALLOC4 accessors and ALLOC_base helpers
  ptolemy_forward.h             Forward declarations of all subroutines
  ptolemy_intrinsics.h          Fortran intrinsic wrappers (DABS, DEXP, DLOG, etc.)
  ptolemy_io.h                  I/O helpers (NXWORD interface, char8 utilities)
```

## Key Translation Notes

### COMMON Block Layout
All Fortran COMMON blocks are mapped to C++ structs in `ptolemy_commons.h`.
1-based Fortran array indexing is preserved via inline accessor functions:
- `FLOAT_arr(I)` — padding-aware access into FloatCommon flat array
- `ALLOC(I)` / `ALLOC_base(ptr)` — numbered allocator access

### Numerov Integration (WAVELJ)
Uses the Raynal complex form: `z(I) = u(I-1) - 10*z(I-1) - z(I-2)` with `u(I) = 12*z(I)/W(I)`.
Integration proceeds in sub-loops of N1ADD steps with overflow rescaling at `|u| > 1e30`.

### Allocator (numbered_store)
The Fortran "numbered storage" (Z array) is a named-block dynamic allocator over a 50 MB pool.
`NALLOC(name, size)` returns a handle; `Z[handle]` is the base index; `ALLOC(Z[h]+i)` accesses element i.
`ALLOC_base(Z[h])` returns a pointer where `ptr[1]` = first element (1-based C array pointer).

### Remaining Stubs
- **Fitting engine**: FITEL, FITEL2, FITINP, SETFIT, LMCFUN, MAKDER, SDERIV, SECDER
- **CC advanced**: CCDUMP, CCMTCH, CCONV, COUPLN, DRIVE (multi-step CC iteration)
- **Other**: CALFUN/CALGRA, PRTCHI
