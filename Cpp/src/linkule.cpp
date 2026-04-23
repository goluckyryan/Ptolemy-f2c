// linkule.cpp -- translated from linkule.f
//
// GETLNK - processes LINKULE keywords and loads linkule
// LINKUL - linkule caller for Ptolemy (dispatches to specific linkule routines)
//
// 11/24/77 - first version
// 3/15/03 - add AV18
// 2/15/07 - add PHIFFER
//
// Translated from Fortran to C++ preserving all logic and labels.

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include "ptolemy_forward.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_alloc.h"
#include "ptolemy_io.h"
#include <cstdio>
#include <cstring>

// ============================================================================
// SUBROUTINE GETLNK ( AKEY, IRET )
//
// Processes linkule keywords and loads linkule.
// Determines what (real potential, imag potential, wave functions) is to be
// set by the linkule and also the linkule name.
// Then loads the linkule and saves its address.
//
// IRET: 0 = error; 1 = O.K.
// ============================================================================

void GETLNK(char8 AKEY, int& IRET)
{
    // implicit real*8 (a-h, o-z)

    // COMMON /LNKBLK/
    auto& NUMLNK = LNKBLK.NUMLNK;
    auto& IUNIQU = LNKBLK.IUNIQU;
    auto* LNKADR = LNKBLK.LNKADR;    // 1-based: LNKADR[i][j]

    // character*8 LNKNAM equivalenced to LNKADR
    // LNKNAM(1,I) <=> LNKADR(1,I) as character*8

    char8 WORD;

    // Key list for linkule usage types
    static const char8 AKEYS[14] = {
        char8(),              // [0] unused
        char8("REALPOTE"),    // [1]
        char8("IMAGPOTE"),    // [2]
        char8("REALSOPO"),    // [3]
        char8("IMAGSOPO"),    // [4]
        char8("COULOMBP"),    // [5]
        char8("WAVEFUNC"),    // [6]
        char8("REALTRPO"),    // [7]
        char8("IMAGTRPO"),    // [8]
        char8("REALTLPO"),    // [9]
        char8("IMAGTLPO"),    // [10]
        char8("REALTPRP"),    // [11]
        char8("IMAGTPRP"),    // [12]
        char8("SIPOTENT")     // [13]
    };

    // List of built-in linkule names
    constexpr int NUMNAM = 17;
    static const char8 NAMES[18] = {
        char8(),              // [0] unused
        char8("BKGPTELP"),    // [1]
        char8("FIXEDWOO"),    // [2]
        char8("GAUSSIAN"),    // [3]
        char8("LAGRANGE"),    // [4]
        char8("LTSTELP "),    // [5]
        char8("RAWITSCH"),    // [6]
        char8("REID"),        // [7]
        char8("SHAPE"),       // [8]
        char8("SPLINE"),      // [9]
        char8("TWOSHAPE"),    // [10]
        char8("DEFORMED"),    // [11]
        char8("JDEPEN"),      // [12]
        char8("JDEPENWS"),    // [13]
        char8("OHTA"),        // [14]
        char8("PARITWOO"),    // [15]
        char8("AV18"),        // [16]
        char8("PHIFFER")      // [17]
    };

    static const char8 INTERN("INTERNAL");

    int IKEY, ILOC;

    // Find the key
    for (IKEY = 1; IKEY <= NUMLNK; IKEY++) {
        if (AKEY == AKEYS[IKEY]) goto L50;
    }

    // Get the linkule name
L50:
    if (NXWORD(WORD.data) != 0) goto L900;

    // Store linkule name: LNKNAM(1,IKEY) = WORD
    // LNKNAM is equivalenced to LNKADR, so LNKNAM(1,IKEY) = first 8 bytes of LNKADR(1,IKEY)
    std::memcpy(&LNKADR[IKEY][1], &WORD, 8);

    ILOC = 0;
    if (WORD == INTERN) goto L200;

    // In the non-IBM version, linkules are built-in
    for (int I = 1; I <= NUMNAM; I++) {
        if (WORD == NAMES[I]) goto L140;
    }
    printf("\n*** %.8s IS NOT IN THIS VERSION OF PTOLEMY.\n", WORD.data);
    goto L950;

L140:
    // 'I' is the loop variable from the for loop above — but we exited via goto
    // We need to capture the matching I. Let's redo this properly:
    {
        int I;
        for (I = 1; I <= NUMNAM; I++) {
            if (WORD == NAMES[I]) break;
        }
        ILOC = I;
    }

L200:
    LNKADR[IKEY][3] = ILOC;
    printf(" LINKULE %.8s IS LINKULE NUMBER/location%8d\n", WORD.data, ILOC);
    IRET = 1;
    return;

L900:
    printf("\n**** A LINKULE NAME MUST FOLLOW THE %.8s KEYWORD.\n", AKEY.data);
L950:
    IRET = 0;
    return;
}


// ============================================================================
// SUBROUTINE LINKUL ( LOC, ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,
//                     L, J, RSTART, STEPSZ, NUMPTS, ARRAY1, ARRAY2, ID )
//
// Linkule caller for Ptolemy.
// Adds the fixed part of a linkule argument list and calls the linkule.
// Non-IBM version: linkules are part of the processor and we just call them.
// ============================================================================

void LINKUL(int LOC, char8 ALIAS, int* MYINTS, int IPOTYP, int IREQUE,
            int& IRETUR, int L, double J, double RSTART, double STEPSZ,
            int NUMPTS, double* ARRAY1, double* ARRAY2, char* ID)
{
    // implicit real*8 (a-h, o-z)

    // COMMON blocks
    auto* Z    = LOCPTRS.Z;        // /LOCptrs/
    auto& FACFR4 = ALLOCS.FACFR4;  // /ALLOCS/
    auto* LENG = LENGTH.LENG;      // /LENGTH/

    // Create FLAT copies of padded COMMON blocks for linkule access.
    // The C++ structs have padding (EXS[6], JS[6], etc.) that shifts
    // flat array offsets. Linkules use raw array indexing matching Fortran
    // layout (no padding). We copy to flat arrays, call linkule, copy back.
    double F_flat[NUMFLOAT + 1]; // 1-based: F_flat[1]..F_flat[153]
    for (int i = 1; i <= NUMFLOAT; i++) F_flat[i] = FLOAT_arr(i);
    auto* F = F_flat;  // Fortran F(1) = F_flat[1], accessed as F[N] for Fortran F(N)

    // INTGER: padding at pos 19 (IZS[0]), 38 (PARITS[0]), 43 (PARIPT[0])
    int I_flat[NUMINTEGER + 1]; // 1-based
    for (int i = 1; i <= NUMINTEGER; i++) I_flat[i] = INTGER_arr_f(i);
    auto* I_arr = I_flat;

    auto& IPRINT = INTGER.IPRINT;

    // JBLOCK: padding at pos 2 (JS[0])
    double JB_flat[NUMJWORD + 1]; // 1-based
    for (int i = 1; i <= NUMJWORD; i++) JB_flat[i] = JBLOCK_arr_f(i);
    auto* JB = JB_flat;

    auto* IS   = SWITCH_arr();     // /SWITCH/  IS(1) — no padding
    auto* C    = reinterpret_cast<double*>(&CNSTNT);   // /CNSTNT/ C(1) — no padding
    auto* T    = TEMPVS_arr();     // /TEMPVS/  T(1)
    auto* WAV  = reinterpret_cast<double*>(&WAVCOM);   // /WAVCOM/
    auto* KM   = reinterpret_cast<double*>(&KANDM);    // /KANDM/
    auto* IN   = reinterpret_cast<double*>(&INTRNL);   // /INTRNL/

    // For CNSTNT we need 1-based, so we shift pointer by -1
    double* C1 = C - 1;
    double* T1 = T - 1;
    double* WAV1 = WAV - 1;
    double* KM1  = KM - 1;
    double* IN1  = IN - 1;

    int IUNIT = 6;
    int NUMOUT;

    bool DBUGSW = (IPRINT % 10) == 9;

    if (DBUGSW) {
        printf(" CALLING LINKULE AT%8d %.8s%8d%8d%4d%4d %.4s\n",
               LOC, ALIAS.data, MYINTS[0], MYINTS[1], IPOTYP, IREQUE, ID);
    }

    // NOTE: F(112) is beginning of PARAM's
    //       I(46) is beginning of IPARAM's

    // The linkule calling convention passes many COMMON block arrays.
    // All linkule routines take the same long argument list.
    // We dispatch based on LOC (1-17) to the appropriate built-in linkule.

    // Arguments common to all linkule calls:
    // (ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,
    //  IUNIT, NUMOUT,
    //  L, J, RSTART, STEPSZ, NUMPTS, ARRAY1, ARRAY2,
    //  F, T, F(112), I, JB, IS, IN, IN, C, WAV, KM,
    //  ID, ALLOC, ILLOC, FACFR4, Z, LENG, NALLOC, NAMLOC,
    //  I(46))

    // Note: MYINTS is real*4[2] in Fortran equivalenced to savint.
    // In the linkule routines (av18, phiffer, etc.) savint is real*4.
    // For the dispatch we pass the raw int* pointer.

    switch (LOC) {
        case 1:  // bkgptelp
            BKGPTElp(ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,
                     IUNIT, NUMOUT, L, J, RSTART, STEPSZ, NUMPTS,
                     ARRAY1, ARRAY2, ID, &INTGER.IPARM1);
            break;
        case 2:  // fixedwoo
            FIXEDWOo(ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,
                     IUNIT, NUMOUT, L, J, RSTART, STEPSZ, NUMPTS,
                     ARRAY1, ARRAY2, ID);
            break;
        case 3:  // gaussian
            GAUSSIAn(ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,
                     IUNIT, NUMOUT, L, J, RSTART, STEPSZ, NUMPTS,
                     ARRAY1, ARRAY2, ID);
            break;
        case 4:  // lagrange
            LAGRANGE(ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,
                     IUNIT, NUMOUT, L, J, RSTART, STEPSZ, NUMPTS,
                     ARRAY1, ARRAY2, ID);
            break;
        case 5:  // ltstelp
            LTSTELp(ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,
                    IUNIT, NUMOUT, L, J, RSTART, STEPSZ, NUMPTS,
                    ARRAY1, ARRAY2, ID, &INTGER.IPARM1);
            break;
        case 6:  // rawitsch
            RAWITSch(ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,
                     IUNIT, NUMOUT, L, J, RSTART, STEPSZ, NUMPTS,
                     ARRAY1, ARRAY2, ID, &INTGER.IPARM1);
            break;
        case 7:  // reid
            REID(ALIAS, reinterpret_cast<real4*>(MYINTS), IPOTYP, IREQUE, IRETUR,
                 IUNIT, NUMOUT, L, J, RSTART, STEPSZ, NUMPTS, ARRAY1, ARRAY2,
                 F + 1, T1, &F[112], I_arr + 1, JB + 1, IS + 1, IN1, IN1,
                 C1, WAV1, KM1,
                 ID, ALLOC_base(1), ILLOC_base(1), FACFR4, Z + 1, LENG + 1,
                 NALLOC, NAMLOC, &I_arr[46]);
            break;
        case 8:  // shape
            SHAPE(ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,
                  IUNIT, NUMOUT, L, J, RSTART, STEPSZ, NUMPTS,
                  ARRAY1, ARRAY2, ID);
            break;
        case 9:  // spline
            SPLINE_linkule(ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,
                           IUNIT, NUMOUT, L, J, RSTART, STEPSZ, NUMPTS,
                           ARRAY1, ARRAY2, ID);
            break;
        case 10: // twoshape
            TWOSHApe(ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,
                     IUNIT, NUMOUT, L, J, RSTART, STEPSZ, NUMPTS,
                     ARRAY1, ARRAY2, ID);
            break;
        case 11: // deformed
            DEFORMEd(ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,
                     IUNIT, NUMOUT, L, J, RSTART, STEPSZ, NUMPTS,
                     ARRAY1, ARRAY2, ID, &INTGER.IPARM1);
            break;
        case 12: // jdepen
            JDEPEN(ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,
                   IUNIT, NUMOUT, L, J, RSTART, STEPSZ, NUMPTS,
                   ARRAY1, ARRAY2, ID);
            break;
        case 13: // jdepenws
            JDEPENWS(ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,
                     IUNIT, NUMOUT, L, J, RSTART, STEPSZ, NUMPTS,
                     ARRAY1, ARRAY2, ID);
            break;
        case 14: // ohta
            OHTA(ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,
                 IUNIT, NUMOUT, L, J, RSTART, STEPSZ, NUMPTS,
                 ARRAY1, ARRAY2, ID, &INTGER.IPARM1);
            break;
        case 15: // paritwoo
            PARITWOO(ALIAS, MYINTS, IPOTYP, IREQUE, IRETUR,
                     IUNIT, NUMOUT, L, J, RSTART, STEPSZ, NUMPTS,
                     ARRAY1, ARRAY2, ID);
            break;

        case 16: // av18
            // Fortran: CALL av18(..., F, T, F(112), I, JB, IS, IN, IN, ...)
            // F, I, JB, IS start at position 1 of their COMMON blocks (1-based Fortran)
            // C++ raw pointers already start at position 1, so pass directly
            // av18 uses Fortran 1-based indexing: flt(N) → flt[N] in C++ with base offset
            av18(ALIAS, reinterpret_cast<real4*>(MYINTS), IPOTYP, IREQUE, IRETUR,
                 IUNIT, NUMOUT, L, J, RSTART, STEPSZ, NUMPTS, ARRAY1, ARRAY2,
                 F, T1, &F[112], I_arr, JB, IS, IN1, IN1,
                 C1, WAV1, KM1,
                 ID, ALLOC_base(1), ILLOC_base(1), FACFR4, Z + 1, LENG + 1,
                 NALLOC, NAMLOC, &I_arr[46]);
            break;

        case 17: // phiffer
            phiffer(ALIAS, reinterpret_cast<real4*>(MYINTS), IPOTYP, IREQUE, IRETUR,
                    IUNIT, NUMOUT, L, J, RSTART, STEPSZ, NUMPTS, ARRAY1, ARRAY2,
                    F, T1, &F[112], I_arr, JB, IS, IN1, IN1,
                    C1, WAV1, KM1,
                    ID, ALLOC_base(1), ILLOC_base(1), FACFR4, Z + 1, LENG + 1,
                    NALLOC, NAMLOC, &I_arr[46]);
            break;

        default:
            goto L900;
    }

    // Copy flat arrays back to padded structs
    for (int i = 1; i <= NUMFLOAT; i++) FLOAT_arr(i) = F_flat[i];
    for (int i = 1; i <= NUMINTEGER; i++) INTGER_arr_f(i) = I_flat[i];
    for (int i = 1; i <= NUMJWORD; i++) JBLOCK_arr_f(i) = JB_flat[i];

    // L190:
    if (DBUGSW) {
        printf(" LINKULE RETURNS%10d%10d%10d%10d%10d\n",
               IRETUR, MYINTS[0], MYINTS[1], IRETUR, NUMOUT);
    }
    return;

L900:
    printf("\n**** LINKULE%10d NOT AVAILABLE.\n", LOC);
    FSTOP(1122);
}
