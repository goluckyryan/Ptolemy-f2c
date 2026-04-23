// linkulesfitters_wavefn.cpp -- translated from linkulesfitters.f
//
// Wavefunction linkule routines: CK, CKSET
// Stubs for: OHTA, RAWITSch, REID (complex routines, translated later)
//
// Translated from Fortran to C++ preserving all logic.

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include "ptolemy_forward.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_alloc.h"
#include "ptolemy_io.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <complex>

// ============================================================================
// FUNCTION CK(R)
// Computes the effective complex K for Rawitscher boundary conditions.
// ============================================================================

// CKBLOK.CPARAM is double[5] in the struct, but Fortran uses COMPLEX*16 CPARAM(4).
// Each COMPLEX*16 = 2 doubles. We reinterpret CPARAM as complex16* for access.
static inline complex16* CPARAM_cx() {
    return reinterpret_cast<complex16*>(&CKBLOK.CPARAM[1]);
    // CPARAM_cx()[0] = CPARAM(1), etc.
}

complex16 CK_func(double R)
{
    complex16* CP = CPARAM_cx();
    double X = R - CKBLOK.R1;
    complex16 CSUM = CP[0] + X * (CP[1] + X * (CP[2] + X * CP[3]));
    return std::sqrt(CKBLOK.FACTOR * CSUM);
}

// ============================================================================
// SUBROUTINE CKSET
// Set up common block for CK function.
// ============================================================================

void CKSET(double RR1, double STEPSZ, double A1, double B1, double DC1, double D1,
           double A2, double B2, double DC2, double D2)
{
    complex16* CP = CPARAM_cx();
    CKBLOK.R1 = RR1;
    CP[0] = complex16(A1 - 1.0, A2);
    CP[1] = complex16(B1, B2);
    CP[2] = complex16(DC1, DC2);
    CP[3] = complex16(D1, D2);
    CKBLOK.FACTOR = 12.0 / (STEPSZ * STEPSZ);
}

// ============================================================================
// Stub routines for wavefunction linkules
// These are complex integration routines that will be fully translated
// when needed for wavefunction calculations.
// ============================================================================

void OHTA(char8 ALIAS, int* MYINTS, int IPOTYP, int IREQUE, int& IRETUR,
          int IUNIT, int& NUMOUT,
          int L, double J, double RSTART, double STEPSZ, int NUMPTS,
          double* WAVR, double* WAVI,
          char* ID, int* IPARAM)
{
    printf("\n**** OHTA LINKULE: wavefunction linkule (Phase 7 stub).\n");
    IRETUR = -1;
    NUMOUT = 1;
}

void RAWITSch(char8 ALIAS, int* MYINTS, int IPOTYP, int IREQUE, int& IRETUR,
              int IUNIT, int& NUMOUT,
              int L, double J, double RSTART, double STEPSZ, int NUMPTS,
              double* WAVR, double* WAVI,
              char* ID, int* IPARAM)
{
    printf("\n**** RAWITSch LINKULE: wavefunction linkule (Phase 7 stub).\n");
    IRETUR = -1;
    NUMOUT = 1;
}

void REID(char8 ALIAS, real4* SAVINT, int IPOTYP, int IREQUE, int& IRETUR,
          int IUNIT, int& NUMOUT, int L, double& JP, double RSTART,
          double STEPSZ, int NUMPTS, double* ARRAY1, double* ARRAY2,
          double* FLT, double* TEMPVS, double* PARAM, int* INTGER_arr,
          double* JBLOCK_arr, int* ISWTCH, double* INTRNL_arr, double* FINTRN,
          double* CNSTNT_arr, double* WAVBL, double* FKANDM,
          char* ID, double* ALLOC_arr, int* ILLOC_arr, int FACFR4,
          integer4* LOC_arr, int* LENGTH_arr,
          int (*NALLOC_fn)(const char*, int), int (*NAMLOC_fn)(const char*),
          int* IPARAM)
{
    printf("\n**** REID LINKULE: deuteron wavefunction linkule (Phase 7 stub).\n");
    IRETUR = -1;
    NUMOUT = 1;
}
