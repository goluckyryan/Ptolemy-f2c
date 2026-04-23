# Ptolemy C++ Translation

Verbatim transliteration of the Fortran DWBA nuclear reaction code Ptolemy (April 2007 version) to C++17.
The translation preserves the original Fortran control flow, variable names, and COMMON block structure.

## Status

| Capability | State |
|---|---|
| Elastic optical model (S-matrix) | **Working** — matches Fortran to full printed precision for all L |
| DWBA transfer cross sections | **Working** — pipeline complete, cross sections produced |
| Coupled-channels (CC) mode | Stubbed |
| Parameter fitting engine | Stubbed |
| Analyzing powers (ANAPOW) | Stubbed |
| Phase shift printing (PHSPRT) | Stubbed |

## Build

```bash
cd Claude_attack/Cpp
make
./ptolemy < ../elastic_incoming.in      # elastic test
./ptolemy < ../reference/DWBA.in        # DWBA test
```

Requires g++ with C++17, compiled with `-O0 -g`.

## DWBA Calculation Flow

```
CONTRL (main control loop)
  |
  +--> DEFALT (initialize defaults)
  +--> PARAM (load parameter set: DPSB etc.)
  +--> REACTN (parse A(d,p)B reaction)
  |
  +--> SETCHN + SETPOT (set up each channel)
  |     |
  |     +--> PROJECTILE: BOUND → av18 linkule (deuteron WF) or WOODSX+Numerov
  |     +--> TARGET:     BOUND → WOODSX + Numerov integration
  |     +--> INCOMING:   WAVSET → SETSPT (Coulomb) + MAKPOT (optical potential)
  |     +--> OUTGOING:   WAVSET → SETSPT (Coulomb) + MAKPOT (optical potential)
  |
  +--> PRBPRT (print reaction summary)
  +--> GETSCT (elastic scattering waves: WAVELJ per L, TMATCH, JPTOLX)
  +--> GRDSET (3D integration grids: CUBMAP + BSPROD)
  +--> ANGSET (angular momentum coupling tables)
  +--> INELDC (DWBA radial integrals: A12 + WFGET + SFROMI)
  +--> LINTRP (L-interpolation/extrapolation: INTRCF)
  +--> XSECTN (cross sections: ELDCS + BETCAL + AMPCAL)
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
  scattering_partA_translated.cpp   GETSCT, COULIN, COULNG, LCRITL area
  scattering_part2_translated.cpp   TMATCH, LCRITL, SWKB
  scattering_part3_translated.cpp   WAVEF, WAVINH, WFGET
  wavelj_translated.cpp         WAVELJ — Numerov + S-matrix matching
  wavetc_translated.cpp         WAVETC — tensor-coupled channel waves
  prbprt_translated.cpp         PRBPRT — reaction summary print
  grdset_translated.cpp         GRDSET — 3D integration grid setup
  angset_translated.cpp         ANGSET — angular momentum coupling tables
  ineldc_translated.cpp         INELDC — core DWBA transfer integrals
  a12_translated.cpp            A12 — angular momentum transformation
  lintrp_translated.cpp         LINTRP — L-interpolation and extrapolation
  xsectn_translated.cpp         XSECTN — cross section computation
  ampcal_translated.cpp         AMPCAL — transition amplitudes per angle
  angular_part2_translated.cpp  BETCAL, CALANG, ELDCS
  source_angular.cpp            ANAPOW (stub), CALFUN, CALGRA (stubs)
  source_channels.cpp           BASCPL, CCHAN, CCOUP, CTRAN (CC stubs)
  source_coupled.cpp            CCDUMP, CCMTCH, CCONV, COUPLN, DRIVE, INELD2, INGRST, INRDIN (CC stubs)
  source_fitting.cpp            FITEL, FITEL2, FITINP, SETFIT, LMCFUN, LXTRP1/2, MAKDER, SDERIV, SECDER, DERCHK, GENBNX, GETBFC (fitting stubs)
  source_scattering.cpp         COULST (stub)
  source_io.cpp                 NXINT, NXWORD, NEWCD, NXVAL, MSCAN, DATAIN, LSTKEY, GETNUM, etc.
  source_misc.cpp               INTRCF, CUBMAP, AITLAG, EPSLON, MUELCO, VCSQ12, etc.
  source_bound.cpp              SETBRN, SETBFC, GETBFC helpers
  numbered_store.cpp            NALLOC, IALLOC, NAMLOC — Fortran-style allocator
  fortlib_part1.cpp             CLEBSH, THREEJ, RACAH, SIXJ, WIG9J, DSGMAL, CCNFRC, CCONTF, etc.
  fortlib_part2.cpp             SPLNCB, INTRPC, LIN, MATINV, GAUSSL, etc.
  masstable.cpp                 EXCESS, MASEXX, AZCODE — nuclear mass table
  av18.cpp                      Argonne v18 deuteron wave function linkule
  linkule.cpp                   LINKUL — linkule dispatcher
  linkulesfitters_fitters.cpp   Optimizer linkules (DLSMIN, LMCHOL, etc.) — stubs
  linkulesfitters_potentials.cpp  Potential linkules
  linkulesfitters_wavefn.cpp    Wave function linkules
  phiffer.cpp                   PHIFFER linkule
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
- `VTEN()[I-1]`, `RTEN()[I-1]`, `ATEN()[I-1]` — 0-based C++ pointers for 1-based Fortran arrays

### Numerov Integration (WAVELJ)
Uses the Raynal complex form: `z(I) = u(I-1) - 10*z(I-1) - z(I-2)` with `u(I) = 12*z(I)/W(I)`.
The effective potential W is stored in `WAVR[I+4]` during the potential loop, then overwritten
with `z` values during the XSI initialization and Numerov sweep.
Integration proceeds in sub-loops of N1ADD steps with overflow rescaling at `|u| > 1e30`.

### Allocator (numbered_store)
The Fortran "numbered storage" (Z array) is a named-block dynamic allocator over a 50 MB pool.
`NALLOC(name, size)` → returns handle; `Z[handle]` → base index; `ALLOC(Z[h]+i)` → element i.
`ALLOC_base(Z[h])` returns a pointer where `ptr[1]` = first element (1-based C array pointer).

### Fortran Linkules
Linkules are Fortran-callable shared-object plugins (av18, phiffer, etc.).
The `LINKUL` dispatcher resolves linkule names to registered C++ function objects.

## Known Issues and Fixes Applied

### Applied Fixes
1. **Tensor potential off-by-one** (`wavset_translated.cpp`, `source_potentials.cpp`):
   `VTEN()`, `RTEN()`, `ATEN()` return 0-based pointers; code used 1-based index `I`.
   Changed `VTEN()[I]` → `VTEN()[I-1]` everywhere. Previously caused NaN in WAVELJ via
   IEEE 754 `Inf × 0 = NaN` when `A=0` tensor parameter was passed to WOODSX.
   **Result:** Elastic S-matrix now matches Fortran reference exactly for all L.

2. **INELDC ALLOC pointer off-by-one**: `&ALLOC(LWAVR)` → `ALLOC_base(LWAVR)` for double* args to WFGET.

3. **TEMPVS not populated in WAVSET**: Added explicit `TEMPVS.TVR = TVR_l` etc. assignments.

4. **VCSQ12 Coulomb formula**: Replaced dimensionless `ZZ = 2*eta*h` with physical
   `VC0 = Z1*Z2*(HBARC/AFINE)` (MeV·fm).

5. **LINTRP INTRCF pointer**: `&ILLOC(LLIS+1)` → `ILLOC_base(LLIS+1)`.

6. **XFACS index transposition in WAVELJ**: `XFACS[1][NWP]` → `XFACS[NWP][1]`.

7. **SFROMI NaN guard**: Added `isnan` check on radial integrals for forbidden transitions.

### Remaining Stubs
- **Coupled-channels mode**: `BASCPL`, `CCHAN`, `CCOUP`, `CTRAN`, `CCDUMP`, `CCMTCH`, `CCONV`, `COUPLN`, `DRIVE`, `INELD2`, `INGRST`, `INRDIN`
- **Fitting engine**: `FITEL`, `FITEL2`, `FITINP`, `SETFIT`, `LMCFUN`, `LXTRP1/2`, `MAKDER`, `SDERIV`, `SECDER`, `DERCHK`, `GENBNX`, `GETBFC`
- **COULST**: Coulomb scattering amplitude (source.f L11210–11873)
- **ANAPOW**, **CALFUN**, **CALGRA**: Analyzing power routines
- **PHSPRT**: Phase shift table printing
- **INELDC scratch file I/O**: Binary Fortran WRITE/READ to unit 1/2 not yet implemented

### Struct Padding
`FloatCommon` has padding elements (`EXS[6]`, `AMXCS[6]`, etc. for 1-based Fortran arrays)
that shift flat-array offsets. Use `FLOAT_arr(I)` for Fortran-indexed access, not raw pointer arithmetic.
