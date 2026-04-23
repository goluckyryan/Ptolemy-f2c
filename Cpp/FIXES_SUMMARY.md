# C++ Ptolemy Fixes Summary

## Verified Fixes (elastic S-matrix matches Fortran to 0.0000%)

### 1. INELDC &ALLOC() off-by-one for double arrays
**File:** `src/ineldc_translated.cpp`
Changed `&ALLOC(LWAVR)` to `ALLOC_base(LWAVR)` for WFGET double* arguments.
A12 call arguments kept as `&ALLOC()` because A12 does internal `-1` adjustment.
SPLNCB/INTRPC calls use `ALLOC_base()` (direct 1-based indexing).

### 2. INELDC &ALLOC4() off-by-one for float arrays
**File:** `src/ineldc_translated.cpp`
Changed `&ALLOC4(LRO+1)` to `ALLOC4_base(LRO+1)` for WFGET float* arguments.

### 3. TEMPVS not populated in WAVSET
**File:** `src/wavset_translated.cpp`
Added: `TEMPVS.TVR = TVR_l;` etc. and `TEMPVS.TVSOR = TVSOR_l;` etc.

### 4. VCSQ12 wrong Coulomb potential formula
**File:** `src/source_misc.cpp`
Replaced dimensionless `ZZ = 2*eta*h` with physical `VC0 = Z1*Z2*(HBARC/AFINE)` (MeV fm).

### 5. LINTRP INTRCF LIS argument
**File:** `src/lintrp_translated.cpp`
Changed `&ILLOC(LLIS + 1)` to `ILLOC_base(LLIS + 1)`.

### 6. XFACS index transposition in WAVELJ Coulomb extension
**Files:** `src/wavelj_translated.cpp`, `src/wavetc_translated.cpp`
Changed `XFACS[1][NWP]` to `XFACS[NWP][1]` etc. The array was accessed with
transposed indices, causing the outgoing channel Coulomb extension to diverge
(max_WF = 3.2×10^260 → 0.81 after fix, matching Fortran exactly).

### 7. SFROMI NaN guard for forbidden transitions
**File:** `src/source_potentials.cpp`
Added NaN check on XIREAL/XIIMAG before S-matrix accumulation.

### 8. DW product NaN/inf guard in INELDC
**File:** `src/ineldc_translated.cpp`
Added NaN/inf check on DW products after distorted wave computation.

### 9. ELDCS BETCAL call indexing
**File:** `src/angular_part2_translated.cpp`
Reverted `ALLOC_base`/`ALLOC4_base` to `&ALLOC`/`&ALLOC4` for BETCAL calls,
because BETCAL uses F1/F2/F3 macros with `(i)-1` (internal 0-based adjustment).

## Remaining Issues

1. **GRDSET zero weights:** The RIROWTS array (bound state form factor × integration
   weights) is all zero, preventing transfer cross sections. The second pass of GRDSET
   produces zero ALLOC4(LWIO) values. Root cause TBD — likely in the phi-integrated
   bound state product (WOW) or the CUBMAP grid weights.

2. **PHSPRT stub:** Phase shift printing is not implemented.

3. **ELDCS AMPCAL not called:** TODO comment in ELDCS prevents elastic cross section computation.
