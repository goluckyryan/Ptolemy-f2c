// linkulesfitters_potentials.cpp -- translated from linkulesfitters.f
//
// Potential linkule routines for Ptolemy.
// Contains: AITKEN, BKGPTElp, DEFORMEd, FIXEDWOo, GAUSSIAn,
//           JDEPEN, JDEPENWS, LAGRANGE, LTSTELp, PARITWOO,
//           SHAPE, SPLINE, TWOSHApe
//
// NOTE: These routines access COMMON blocks via global structs
//       (FLOAT_common, INTGER, CNSTNT, etc.) rather than through
//       the passed array parameters, to avoid struct padding issues.
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
#include <cmath>

// Forward declarations for source.f routines (Phase 8)
extern void WOODSX(int NUMPTS, double RSTART, double STEPSZ, double* ARRAY,
                   int& N1, int& N2, int ITYPE, double V, double R, double A, double POWER);
extern void DEFWOO(int NUMPTS, double RSTART, double STEPSZ, double* ARRAY,
                   int& N1, int& N2, int ITYPE, double V, double R, double A, double POWER,
                   int IPORT, int LXMIN, int* LDEFMX,
                   double* RBETAS, int* LXS, int NCOSIN);


// ============================================================================
// SUBROUTINE AITKEN
//
// Aitken-Lagrange interpolation on a nonuniform grid.
// Nearest surrounding points are used, biased toward center.
//
// May 15, 1976 - first version - S. Pieper
// Nov 5, 1976 - always use surrounding points
// ============================================================================

void AITKEN(int NDEGRE, double DELTA, double EPSLON, int NIN,
            double* XIN, double* FIN,
            int NOUT, double* XOUT, double* FOUT,
            int& NFAIL, double& WORST, double* WORK)
{
    // WORK is dimensioned (2, NDEGRE+1) minimum — treated as WORK[2*(NDEGRE+1)]
    // WORK(1,j) = WORK[2*(j-1)], WORK(2,j) = WORK[2*(j-1)+1]
    // Using 1-based: WORK1(j) = WORK[2*(j-1)], WORK2(j) = WORK[2*(j-1)+1]
    #define WORK1(j) WORK[2*((j)-1)]
    #define WORK2(j) WORK[2*((j)-1)+1]

    bool DELTSW = DELTA > 0.0;
    bool UPSW;

    int NI = 1;
    NFAIL = 0;
    WORST = 0.0;
    double XMID = XIN[(NIN + 1) / 2];
    int NDEGR = MIN0(NDEGRE, NIN - 1);

    for (int NO = 1; NO <= NOUT; NO++) {
        double X = XOUT[NO];
        double F, ERROR;

        // Find nearest point
        UPSW = false;
        double D = DABS(X - XIN[NI]);

    L120:
        if (NI == NIN) goto L150;
        {
            double D2 = DABS(X - XIN[NI + 1]);
            if (D2 >= D) goto L150;
            UPSW = true;
            NI = NI + 1;
            D = D2;
        }
        goto L120;

    L150:
        if (UPSW) goto L200;
    L160:
        if (NI == 1) goto L200;
        {
            double D2 = DABS(X - XIN[NI - 1]);
            if (D2 >= D) goto L200;
            NI = NI - 1;
            D = D2;
        }
        goto L160;

        // Have found nearest XIN
    L200:
        {
            int NL = NI - 1;
            int NU = NI + 1;
            F = FIN[NI];
            WORK1(1) = X - XIN[NI];
            UPSW = WORK1(1) < 0.0;
            WORK2(1) = F;
            ERROR = 0.0;
            if (NDEGR <= 0) goto L700;

            // Do up to NDEGR iterations
            for (int ND = 1; ND <= NDEGR; ND++) {
                double FLAST = F;
                int I;

                // Bias choice of 3rd point towards center
                if (ND == 2) UPSW = X < XMID;
                if (UPSW) goto L350;

                // Last point was lower - now use upper
            L320:
                I = NU;
                NU = NU + 1;
                UPSW = true;
                if (I <= NIN) goto L400;

                // Last point was upper - now use lower
            L350:
                I = NL;
                NL = NL - 1;
                UPSW = false;
                if (I < 1) goto L320;

                // Now use next point
            L400:
                {
                    double DEL = X - XIN[I];
                    F = FIN[I];
                    for (int JJ = 1; JJ <= ND; JJ++) {
                        F = (WORK1(JJ) * F - DEL * WORK2(JJ)) / (WORK1(JJ) - DEL);
                    }

                    // Have we converged
                    if (!DELTSW) goto L550;
                    ERROR = DABS(F);
                    if (ERROR <= EPSLON) goto L700;
                    ERROR = DABS((F - FLAST) / F);
                    if (ERROR <= DELTA) goto L700;
                L550:
                    WORK1(ND + 1) = DEL;
                    WORK2(ND + 1) = F;
                }
            } // end ND loop

            // Stopped by the degree
            if (DELTSW) NFAIL = NFAIL + 1;
        }

        // Converged
    L700:
        FOUT[NO] = F;
        WORST = DMAX1(WORST, ERROR);
    } // end NO loop

    #undef WORK1
    #undef WORK2
    return;
}


// ============================================================================
// SUBROUTINE GAUSSIAn
//
// Generalized Gaussian potential:
//   POT(X) = -V * EXP( (R^(2P) - X^(2P)) / A^(2P) )
//
// 5/30/78 - made from FIXEDWOOD - S.P.
// ============================================================================

void GAUSSIAn(char8 ALIAS, int* MYINTS, int IPOTYP, int IREQUE, int& IRETUR,
              int IUNIT, int& NUMOUT,
              int L, double J, double RSTART, double STEPSZ, int NUMPTS,
              double* ARRAY1, double* ARRAY2,
              char* ID)
{
    double UNDEF = INTRNL.UNDEF;

    NUMOUT = 0;
    IRETUR = 0;

    // Get parameters based on potential type
    double V = TEMPVS_arr()[IPOTYP - 1];          // TEMPVS(IPOTYP)
    double A = TEMPVS_arr()[3 + IPOTYP];           // TEMPVS(4+IPOTYP)
    double R = FLOAT_common.R;                     // FLOAT(43)
    if (IPOTYP == 2) R = FLOAT_common.RI;          // FLOAT(47)
    double POWER;
    if (IPOTYP == 1) POWER = FLOAT_common.POWRL;   // FLOAT(102)
    else POWER = FLOAT_common.POWIM;               // FLOAT(103)

    if (R == UNDEF || A == UNDEF || V == UNDEF || POWER == UNDEF) goto L950;
    if (R <= 0.0 || A <= 0.0 || POWER <= 0.0) goto L960;

    if (IREQUE == 1) goto L100;
    if (IREQUE == 2) goto L500;
    if (IREQUE == 3) goto L600;
    return;

    // Generate radius array for INTRPC
L100:
    MYINTS[0] = 0;
    MYINTS[1] = 0;
    if (NUMOUT != 0) goto L900;
    return;

    // Print it out
L500:
    printf(" THE POTENTIAL IS A GENERALIZED GAUSSIAN:\n");
    printf(" V =%13.5g AT R =%6.2f FM;   A =%8.3f FM,   POWER =%8.4f\n",
           V, R, A, POWER);
    goto L900;

    // Do it
L600:
    {
        double TWOP = POWER + POWER;
        // In the following 46 is LN(1E+20) and 23 is LN(1E+10)
        int NMAX = (int)(std::pow(std::pow(R, TWOP) + std::pow(A, TWOP) *
                   (46.0 + DLOG(DABS(V))), 1.0 / TWOP) / STEPSZ);
        NMAX = MIN0(NUMPTS, NMAX);
        int NMIN = 1;
        double TEMP = std::pow(R, TWOP) - std::pow(A, TWOP) * (23.0 - DLOG(DABS(V)));
        if (TEMP > 0.0) NMIN = (int)(std::pow(TEMP, 1.0 / TWOP) / STEPSZ);
        NMIN = MAX0(NMIN, 1);

        if (NMIN != 1) {
            for (int I = 1; I <= NMIN; I++) {
                ARRAY1[I] = -1.0e+10;
            }
        }

        double RVAL = (NMIN - 1) * STEPSZ / A;
        for (int I = NMIN; I <= NMAX; I++) {
            ARRAY1[I] = -V * DEXP(std::pow(R / A, TWOP) - std::pow(RVAL, TWOP));
            RVAL = RVAL + (STEPSZ / A);
        }

        if (NMAX >= NUMPTS) return;
        NMAX = NMAX + 1;
        for (int I = NMAX; I <= NUMPTS; I++) {
            ARRAY1[I] = 0.0;
        }
    }
    return;

L900:
    NUMOUT = 1;
    return;

L950:
    printf("\n**** R, V, AND A MUST BE DEFINED.\n");
    IRETUR = -1;
    goto L900;

L960:
    IRETUR = -1;
    return;
}


// ============================================================================
// SUBROUTINE FIXEDWOo
//
// Woods-Saxon potential fixed at a point.
// V(R=PARAM1) = VFIX, with shape from R, A, POWER.
//
// 4/6/78 - based on SHAPE - S.P.
// ============================================================================

void FIXEDWOo(char8 ALIAS, int* MYINTS, int IPOTYP, int IREQUE, int& IRETUR,
              int IUNIT, int& NUMOUT,
              int L, double J, double RSTART, double STEPSZ, int NUMPTS,
              double* ARRAY1, double* ARRAY2,
              char* ID)
{
    double UNDEF = INTRNL.UNDEF;

    NUMOUT = 0;
    IRETUR = 0;

    double RFIX;
    if (IPOTYP == 1) RFIX = FLOAT_common.PARAM1;
    else RFIX = FLOAT_common.PARAM2;
    double VFIX = TEMPVS_arr()[IPOTYP - 1];
    double A = TEMPVS_arr()[3 + IPOTYP];
    double R = FLOAT_common.R;
    if (IPOTYP == 2) R = FLOAT_common.RI;
    double POWER;
    if (IPOTYP == 1) POWER = FLOAT_common.POWRL;
    else POWER = FLOAT_common.POWIM;

    if (R == UNDEF || A == UNDEF || VFIX == UNDEF || RFIX == UNDEF) {
        printf("\n**** R, V, A, AND PARAM1 MUST BE DEFINED.\n");
        IRETUR = -1;
        NUMOUT = 1;
        return;
    }
    if (R <= 0.0 || A <= 0.0 || RFIX <= 0.0) {
        IRETUR = -1;
        return;
    }

    double V = VFIX * std::pow(1.0 + DEXP((RFIX - R) / A), POWER);

    if (IREQUE == 1) {
        MYINTS[0] = 0;
        MYINTS[1] = 0;
        if (NUMOUT != 0) { NUMOUT = 1; return; }
        return;
    }

    if (IREQUE == 2) {
        printf(" THE POTENTIAL IS A WOODS-SAXON:%15.5g%15.5g%15.5g     POWER =%14.5g\n",
               V, R, A, POWER);
        printf(" WITH V =%15.5g     AT R =%15.5g\n", VFIX, RFIX);
        NUMOUT = 1;
        return;
    }

    if (IREQUE == 3) {
        int N1, N2;
        WOODSX(NUMPTS, RSTART, STEPSZ, ARRAY1, N1, N2, 1, V, R, A, POWER);
        return;
    }
}


// ============================================================================
// SUBROUTINE JDEPEN
//
// J-dependent depth Woods-Saxon potential:
//   V = (1 + PARAM1*J + PARAM2*J^2) * WS(V, R, A)
//
// 11/16/81 - JDEPEN based on PARITWOOD - S.P.
// ============================================================================

void JDEPEN(char8 ALIAS, int* MYINTS, int IPOTYP, int IREQUE, int& IRETUR,
            int IUNIT, int& NUMOUT,
            int L, double J, double RSTART, double STEPSZ, int NUMPTS,
            double* ARRAY1, double* ARRAY2,
            char* ID)
{
    double UNDEF = INTRNL.UNDEF;
    double BIGLOG = CNSTNT.BIGLOG;

    NUMOUT = 0;
    IRETUR = 0;

    double PARAM1_val = FLOAT_common.PARAM1;
    double PARAM2_val = FLOAT_common.PARAM2;
    if (PARAM2_val == UNDEF) PARAM2_val = 0.0;
    double V = TEMPVS_arr()[IPOTYP - 1];
    double A = TEMPVS_arr()[3 + IPOTYP];
    double R = FLOAT_common.R;
    if (IPOTYP == 2) R = FLOAT_common.RI;
    double POWER;
    if (IPOTYP == 1) POWER = FLOAT_common.POWRL;
    else POWER = FLOAT_common.POWIM;

    if (R == UNDEF || A == UNDEF || V == UNDEF || PARAM1_val == UNDEF) goto L950;
    if (R <= 0.0 || A <= 0.0) goto L960;

    if (IREQUE == 1) goto L100;
    if (IREQUE == 2) goto L500;
    if (IREQUE == 3) goto L600;
    if (IREQUE == 4) goto L700;
    return;

    // Generate work array to save the W.S.
L100:
    {
        char worknm[9];
        std::memcpy(worknm, ID, 4);
        worknm[4] = '\0';
        int I = NALLOC(worknm, NUMPTS);
        MYINTS[0] = I;
        MYINTS[1] = 0;
    }
    // This is an L-dependent potential
    IRETUR = +1;
    if (NUMOUT != 0) goto L900;
    return;

L500:
    printf(" THE POTENTIAL IS A J-DEPENDANT WOODS-SAXON:%15.5g%15.5g%15.5g     POWER =%14.5g\n",
           V, R, A, POWER);
    printf(" WITH THE J-DEPENDANT DEPTH FACTOR  1 + J *%13.5g + J**2 *%13.5g\n",
           PARAM1_val, PARAM2_val);
    goto L900;

    // Set up the two pieces
L600:
    {
        int LL = LOCPTRS.Z[MYINTS[0]];
        int N1, N2;
        WOODSX(NUMPTS, RSTART, STEPSZ, ALLOC_base(LL), N1, N2, 1, V, R, A, POWER);
        for (int I = 1; I <= NUMPTS; I++) {
            ARRAY1[I] = ALLOC(LL - 1 + I);
        }
    }
    return;

    // Add in the J-dependent part
L700:
    {
        int LL = LOCPTRS.Z[MYINTS[0]];
        double DJ = 0.5 * JBLOCK.J;
        double FAC = ARRAY2[1] * DJ * (PARAM1_val + DJ * PARAM2_val);
        int N1 = (int)(RSTART / STEPSZ + 0.5);
        LL = LL - 1 + N1;
        for (int I = 1; I <= NUMPTS; I++) {
            ARRAY1[I] = ARRAY1[I] + FAC * ALLOC(LL + I);
        }
    }
    return;

L900:
    NUMOUT = 1;
    return;

L950:
    printf("\n**** R, V, A, AND PARAM1 MUST BE DEFINED.\n");
    IRETUR = -1;
    goto L900;

L960:
    IRETUR = -1;
    return;
}


// ============================================================================
// SUBROUTINE JDEPENWS
//
// J-dependent W.S. depth potential:
//   V = F(J) * WS(V, R, A)
//   F(J) = 1 / (1 + EXP((J - PARM100)/PARM101))
//
// 5/11/82 - JDEPENWS made from JDEPEN - S.P.
// ============================================================================

void JDEPENWS(char8 ALIAS, int* MYINTS, int IPOTYP, int IREQUE, int& IRETUR,
              int IUNIT, int& NUMOUT,
              int L, double J, double RSTART, double STEPSZ, int NUMPTS,
              double* ARRAY1, double* ARRAY2,
              char* ID)
{
    double UNDEF = INTRNL.UNDEF;
    double BIGLOG = CNSTNT.BIGLOG;

    NUMOUT = 0;
    IRETUR = 0;

    // Fortran PARAM(10) and PARAM(11) via padding-aware accessor.
    double PARM10 = PARAMS_at(10);
    double PARM11 = PARAMS_at(11);
    double V = TEMPVS_arr()[IPOTYP - 1];
    double A = TEMPVS_arr()[3 + IPOTYP];
    double R = FLOAT_common.R;
    if (IPOTYP == 2) R = FLOAT_common.RI;
    double POWER;
    if (IPOTYP == 1) POWER = FLOAT_common.POWRL;
    else POWER = FLOAT_common.POWIM;

    if (R == UNDEF || A == UNDEF || V == UNDEF ||
        PARM10 == UNDEF || PARM11 == UNDEF) goto L950;
    if (R <= 0.0 || A <= 0.0) goto L960;

    if (IREQUE == 1) goto L100;
    if (IREQUE == 2) goto L500;
    if (IREQUE == 3) goto L600;
    if (IREQUE == 4) goto L700;
    return;

L100:
    {
        char worknm[9];
        std::memcpy(worknm, ID, 4);
        worknm[4] = '\0';
        int I = NALLOC(worknm, NUMPTS);
        MYINTS[0] = I;
        MYINTS[1] = 0;
    }
    IRETUR = +1;
    if (NUMOUT != 0) goto L900;
    return;

L500:
    printf(" THE POTENTIAL IS A J-DEPENDANT WOODS-SAXON:%15.5g%15.5g%15.5g     POWER =%14.5g\n",
           V, R, A, POWER);
    printf(" WITH THE J-DEPENDANT DEPTH FACTOR  1 / ( 1 + EXP( ( J - %12.5g) / %12.5g) )\n",
           PARM10, PARM11);
    goto L900;

L600:
    {
        int LL = LOCPTRS.Z[MYINTS[0]];
        int N1, N2;
        WOODSX(NUMPTS, RSTART, STEPSZ, ALLOC_base(LL), N1, N2, 1, V, R, A, POWER);
        for (int I = 1; I <= NUMPTS; I++) {
            ARRAY1[I] = ALLOC(LL - 1 + I);
        }
    }
    return;

    // Add in the J-dependent part
L700:
    {
        int LL = LOCPTRS.Z[MYINTS[0]];
        double DJ = 0.5 * JBLOCK.J;
        double FAC = 1.0;
        double X = (DJ - PARM10) / PARM11;
        if (X < -30.0) goto L710;
        FAC = 0.0;
        if (X > 30.0) goto L710;
        FAC = 1.0 / (1.0 + DEXP(X));
    L710:
        FAC = ARRAY2[1] * (FAC - 1.0);
        int N1 = (int)(RSTART / STEPSZ + 0.5);
        LL = LL - 1 + N1;
        for (int I = 1; I <= NUMPTS; I++) {
            ARRAY1[I] = ARRAY1[I] + FAC * ALLOC(LL + I);
        }
    }
    return;

L900:
    NUMOUT = 1;
    return;

L950:
    printf("\n**** R, V, A, PARAM10 & PARAM11 MUST BE DEFINED.\n");
    IRETUR = -1;
    goto L900;

L960:
    IRETUR = -1;
    return;
}


// ============================================================================
// SUBROUTINE PARITWOO
//
// Parity-dependent Woods-Saxon:
//   V = (1 + PARAM1*(-1)^L) * WS(V, R, A)
//
// 6/13/78 - based on FIXEDWOOD - S.P.
// ============================================================================

void PARITWOO(char8 ALIAS, int* MYINTS, int IPOTYP, int IREQUE, int& IRETUR,
              int IUNIT, int& NUMOUT,
              int L, double J, double RSTART, double STEPSZ, int NUMPTS,
              double* ARRAY1, double* ARRAY2,
              char* ID)
{
    double UNDEF = INTRNL.UNDEF;

    NUMOUT = 0;
    IRETUR = 0;

    double PCOEF;
    if (IPOTYP == 1) PCOEF = FLOAT_common.PARAM1;
    else PCOEF = FLOAT_common.PARAM2;
    double V = TEMPVS_arr()[IPOTYP - 1];
    double A = TEMPVS_arr()[3 + IPOTYP];
    double R = FLOAT_common.R;
    if (IPOTYP == 2) R = FLOAT_common.RI;
    double POWER;
    if (IPOTYP == 1) POWER = FLOAT_common.POWRL;
    else POWER = FLOAT_common.POWIM;

    if (R == UNDEF || A == UNDEF || V == UNDEF || PCOEF == UNDEF) goto L950;
    if (R <= 0.0 || A <= 0.0) goto L960;

    if (IREQUE == 1) goto L100;
    if (IREQUE == 2) goto L500;
    if (IREQUE == 3) goto L600;
    if (IREQUE == 4) goto L700;
    return;

L100:
    {
        char worknm[9];
        std::memcpy(worknm, ID, 4);
        worknm[4] = '\0';
        int I = NALLOC(worknm, NUMPTS);
        MYINTS[0] = I;
        MYINTS[1] = 0;
    }
    IRETUR = +1;
    if (NUMOUT != 0) goto L900;
    return;

L500:
    printf(" THE POTENTIAL IS A PARITY WOODS-SAXON:%15.5g%15.5g%15.5g     POWER =%14.5g\n",
           V, R, A, POWER);
    printf(" WITH THE PARITY-DEPENDANT DEPTH FACTOR  1 + (-PARAM1)**L:  PARAM1 =%15.5g\n", PCOEF);
    goto L900;

L600:
    {
        int LL = LOCPTRS.Z[MYINTS[0]];
        int N1, N2;
        WOODSX(NUMPTS, RSTART, STEPSZ, ALLOC_base(LL), N1, N2, 1, V, R, A, POWER);
        for (int I = 1; I <= NUMPTS; I++) {
            ARRAY1[I] = ALLOC(LL - 1 + I);
        }
    }
    return;

    // Add in the L-dependent part
L700:
    {
        int LL = LOCPTRS.Z[MYINTS[0]];
        double FAC = ARRAY2[1] * PCOEF;
        if (BTEST(L, 0)) FAC = -FAC;
        for (int I = 1; I <= NUMPTS; I++) {
            ARRAY1[I] = ARRAY1[I] + FAC * ALLOC(LL - 1 + I);
        }
    }
    return;

L900:
    NUMOUT = 1;
    return;

L950:
    printf("\n**** R, V, A, AND PARAM1 MUST BE DEFINED.\n");
    IRETUR = -1;
    goto L900;

L960:
    IRETUR = -1;
    return;
}


// ============================================================================
// SUBROUTINE SHAPE
//
// Simple linkule: fixed shape times depth parameter.
// Shape is read from named objects (REALSHAP, IMAGSHAP, etc.)
//
// 12/8/77 - first version - S.P.
// ============================================================================

void SHAPE(char8 ALIAS, int* MYINTS, int IPOTYP, int IREQUE, int& IRETUR,
           int IUNIT, int& NUMOUT,
           int L, double J, double RSTART, double STEPSZ, int NUMPTS,
           double* ARRAY1, double* ARRAY2,
           char* ID)
{
    static const char NAMES[6][9] = {
        "", "REALSHAP", "IMAGSHAP", "REALSOSH", "IMAGSOSH", "SHAPE   "
    };
    static const char SCALNM[6][9] = {
        "", "REALSCAL", "IMAGSCAL", "REALSOSC", "IMAGSOSC", "SHAPESCA"
    };

    double SMLNUM = CNSTNT.SMLNUM;
    NUMOUT = 0;
    IRETUR = 0;

    double CONS = TEMPVS_arr()[IPOTYP - 1];
    if (CONS != INTRNL.UNDEF) goto L45;
    printf("\n**** THE WELL DEPTH PARAMETER (V, VI, VSO, VSOI) MUST BE DEFINED FOR THE SHAPE LINKULE.\n");
    IRETUR = -1;
    goto L900;
L45:
    CONS = -CONS;

    if (IREQUE == 1) goto L100;
    if (IREQUE == 2) goto L500;
    if (IREQUE == 3) goto L600;
    return;

    // Locate the shape
L100:
    {
        int II = IPOTYP;
        int I = NAMLOC(NAMES[II]);
        if (I != 0) goto L150;
        II = 5;
        I = NAMLOC(NAMES[II]);
        if (I != 0) goto L150;
        printf("\n**** AN OBJECT WITH THE NAME %.8s OR SHAPE MUST BE DEFINED TO USE THE SHAPE LINKULE.\n",
               NAMES[IPOTYP]);
        IRETUR = -1;
        goto L900;

    L150:
        {
            int ISCAL = NAMLOC(SCALNM[II]);
            double RMAX = RSTART + (NUMPTS - 1) * STEPSZ;
            int NIN = LENGTH.LENG[I];

            if (ISCAL != 0) goto L200;

            // No scale array — check number of points
            if (NUMPTS == NIN) goto L180;
            printf("\n**** ASYMPTOPIA AND STEPSIZE =%15.5g%15.5g AND REQUIRE%6d POINTS.\n",
                   RMAX, STEPSZ, NUMPTS);
            printf("      HOWEVER, %.8s HAS%6d POINTS AND %.8s IS NOT DEFINED.\n",
                   NAMES[II], NIN, SCALNM[II]);
            IRETUR = -1;
            goto L900;

        L180:
            printf(" %.8s IS BEING ASSUMED TO BE DEFINED FOR%15.5g =< R =<%15.5g WITH STEPSIZE =%15.5g\n",
                   NAMES[II], RSTART, RMAX, STEPSZ);
            NUMOUT = 1;
            goto L400;

        L200:
            {
                int LISCAL = LENGTH.LENG[ISCAL];
                bool RSW;

                if (LISCAL == 2) goto L210;
                if (LISCAL == NIN) goto L220;
                printf("\n**** SCALLING ARRAY %.8s SHOULD HAVE 2 OR%4d ELEMENTS, BUT IT HAS%6d ELEMENTS.\n",
                       SCALNM[II], NIN, LISCAL);
                IRETUR = -1;
                goto L900;

            L210:
                {
                    double R1 = ALLOC(LOCPTRS.Z[ISCAL]);
                    double R2 = ALLOC(LOCPTRS.Z[ISCAL] + 1);
                    RSW = false;
                    if (R1 == RSTART && NIN == NUMPTS &&
                        DABS(R2 - RMAX) / STEPSZ < 1.0e-5)
                        goto L400;
                    goto L240;
                }

            L220:
                RSW = true;
                goto L240;

            L240:
                {
                    // Must rescale — use AITKEN interpolation
                    double R1, R2;
                    if (RSW) {
                        R1 = ALLOC(LOCPTRS.Z[ISCAL]);
                        R2 = ALLOC(LOCPTRS.Z[ISCAL] + NIN - 1);
                    } else {
                        R1 = ALLOC(LOCPTRS.Z[ISCAL]);
                        R2 = ALLOC(LOCPTRS.Z[ISCAL] + 1);
                    }

                    int IOLD = I;
                    int ISIZE = NIN + NUMPTS;
                    char worknm[9] = "    SHAP";
                    std::memcpy(worknm, ID, 4);
                    I = NALLOC(worknm, NUMPTS);
                    char worknm2[9] = "    WORK";
                    std::memcpy(worknm2, ID, 4);
                    int IWORK = NALLOC(worknm2, ISIZE);

                    int LXIN = LOCPTRS.Z[IWORK] - 1;
                    int LXOUT = LXIN + NIN;
                    int LYIN = LOCPTRS.Z[IOLD] - 1;
                    int LYOUT = LOCPTRS.Z[I] - 1;

                    if (RSW) LXIN = LOCPTRS.Z[ISCAL] - 1;
                    if (!RSW) {
                        // Setup input X's
                        double RR = R1;
                        double STEPIN = (R2 - R1) / (NIN - 1);
                        for (int N = 1; N <= NIN; N++) {
                            ALLOC(LXIN + N) = RR;
                            RR = RR + STEPIN;
                        }
                    }

                    // Get extrapolating function
                    double VSTEP = 0.0;
                    double AVAL = 0.0;
                    double YNM1 = ALLOC(LYIN + NIN - 1);
                    double YN = ALLOC(LYIN + NIN);
                    if (YN * YNM1 > 0.0) {
                        AVAL = YN;
                        double BVAL = DLOG(YNM1 / AVAL) / (R2 - ALLOC(LXIN + NIN - 1));
                        VSTEP = DEXP(-BVAL * STEPSZ);
                    }

                    // Setup output X's
                    double RR = RSTART;
                    int NSTRT = 1;
                    int NLAST = 1;
                    double VAL0 = ALLOC(LYIN + 1);
                    for (int N = 1; N <= NUMPTS; N++) {
                        if (RR >= R1) {
                            if (RR <= R2) {
                                ALLOC(LXOUT + N) = RR;
                                NLAST = N;
                            } else {
                                AVAL = AVAL * VSTEP;
                                if (DABS(AVAL) < SMLNUM) AVAL = 0.0;
                                ALLOC(LYOUT + N) = AVAL;
                            }
                        } else {
                            ALLOC(LYOUT + N) = VAL0;
                            NSTRT = N + 1;
                        }
                        RR = RR + STEPSZ;
                    }
                    int NOUT = NLAST + 1 - NSTRT;

                    // Do the interpolation
                    double WORK[20];
                    int NFAIL;
                    double WORST;
                    AITKEN(5, 0.0, 0.0, NIN, ALLOC_base(LXIN + 1), ALLOC_base(LYIN + 1),
                           NOUT, ALLOC_base(LXOUT + NSTRT), ALLOC_base(LYOUT + NSTRT),
                           NFAIL, WORST, WORK);

                    // Free work area
                    LOCPTRS.Z[IWORK] = -LOCPTRS.Z[IWORK];
                }
            }
        }

    L400:
        MYINTS[0] = I;
        MYINTS[1] = II;
        if (NUMOUT != 0) goto L900;
        return;
    }

    // Print it out
L500:
    printf(" THE POTENTIAL SHAPE STORED IN OBJECT %.8s IS BEING MULTIPLIED BY%15.5g\n",
           NAMES[MYINTS[1]], CONS);
    goto L900;

    // Do it
L600:
    {
        int LL = LOCPTRS.Z[MYINTS[0]] - 1;
        for (int I = 1; I <= NUMPTS; I++) {
            ARRAY1[I] = CONS * ALLOC(LL + I);
        }
    }
    return;

L900:
    NUMOUT = 1;
    return;
}


// ============================================================================
// SUBROUTINE SPLINE
//
// Construct potential by cubic spline interpolation of log|V|.
// Uses PARAM pairs (R1,V1)...(R10,V10).
//
// 4/6/78 - based on SHAPE - S.P.
// ============================================================================

void SPLINE_linkule(char8 ALIAS, int* MYINTS, int IPOTYP, int IREQUE, int& IRETUR,
                    int IUNIT, int& NUMOUT,
                    int L, double J, double RSTART, double STEPSZ, int NUMPTS,
                    double* ARRAY1, double* ARRAY2,
                    char* ID)
{
    double UNDEF = INTRNL.UNDEF;

    // Local arrays for sorting points
    double R_arr[21], V_arr[21], B[21], C_arr[21], D_arr[21];

    NUMOUT = 0;
    IRETUR = 0;

    // Get the defined points and values and order them.
    // Fortran PARAM is a 2x10 array equivalenced to PARAM1..PARAM20:
    //   PARAM(1,IV) = R-coord (Fortran PARAMS(2*IV-1))
    //   PARAM(2,IV) = V-coord (Fortran PARAMS(2*IV))
    // Use PARAMS_at to handle PAR620[16] 1-based padding correctly.

    int NV = 0;
    int IV = 0;

    // Find first valid pair
    while (IV < 10) {
        double r_val = PARAMS_at(2 * IV + 1);
        double v_val = PARAMS_at(2 * IV + 2);
        if (r_val != UNDEF && v_val != UNDEF) {
            NV = 1;
            R_arr[1] = r_val;
            V_arr[1] = v_val;
            IV++;
            break;
        }
        IV++;
    }
    if (NV == 0) goto L950;

    // Find remaining valid pairs and insert in order
    while (IV < 10) {
        double r_val = PARAMS_at(2 * IV + 1);
        double v_val = PARAMS_at(2 * IV + 2);
        IV++;
        if (r_val == UNDEF || v_val == UNDEF) continue;

        // Find insertion point
        int ins = NV + 1;
        for (int I = 1; I <= NV; I++) {
            if (r_val < R_arr[I]) {
                ins = I;
                break;
            }
        }
        // Shift right
        for (int II = NV; II >= ins; II--) {
            R_arr[II + 1] = R_arr[II];
            V_arr[II + 1] = V_arr[II];
        }
        NV++;
        R_arr[ins] = r_val;
        V_arr[ins] = v_val;
    }

    if (NV <= 2) goto L950;

    if (IREQUE == 1) goto L100;
    if (IREQUE == 2) goto L500;
    if (IREQUE == 3) goto L600;
    return;

L100:
    {
        char worknm[9] = "    SHAP";
        std::memcpy(worknm, ID, 4);
        int I = NALLOC(worknm, NUMPTS);
        double RR = RSTART;
        int LOC_I = LOCPTRS.Z[I];
        for (int II = 1; II <= NUMPTS; II++) {
            ALLOC(LOC_I - 1 + II) = RR;
            RR = RR + STEPSZ;
        }
        MYINTS[0] = I;
        MYINTS[1] = 0;
    }
    if (NUMOUT != 0) goto L900;
    return;

L500:
    printf(" THE POTENTIAL IS DEFINED BY THE POINTS:\n");
    for (int I = 1; I <= NV; I++) {
        printf("%10.3f%13.5g", R_arr[I], V_arr[I]);
        if (I % 5 == 0 || I == NV) printf("\n");
    }
    goto L900;

L600:
    {
        int LL = LOCPTRS.Z[MYINTS[0]] - 1;

        // Setup output X's
        double RR = RSTART;
        int NSTRT = 1;
        int NLAST = 1;
        for (int N = 1; N <= NUMPTS; N++) {
            if (RR < R_arr[1]) {
                ARRAY1[N] = -V_arr[1];
                NSTRT = N + 1;
            } else if (RR > R_arr[NV]) {
                ARRAY1[N] = 0.0;
            } else {
                NLAST = N;
            }
            RR = RR + STEPSZ;
        }
        int NOUT = NLAST + 1 - NSTRT;

        // Interpolate the logs
        for (int I = 1; I <= NV; I++) {
            V_arr[I] = DLOG(V_arr[I]);
        }

        SPLNCB(NV, R_arr, V_arr, B, C_arr, D_arr);
        INTRPC(NV, R_arr, V_arr, B, C_arr, D_arr,
               NOUT, ALLOC_base(LL + NSTRT), &ARRAY1[NSTRT]);

        for (int I = NSTRT; I <= NLAST; I++) {
            ARRAY1[I] = -DEXP(ARRAY1[I]);
        }

        // Fill in the end with exponential decay
        double TERM = ARRAY1[NLAST];
        NLAST = NLAST + 1;
        if (NLAST > NUMPTS) return;
        double STEP = DEXP(B[NV] * STEPSZ);
        TERM = TERM * STEP;
        for (int I = NLAST; I <= NUMPTS; I++) {
            ARRAY1[I] = TERM;
            if (TERM > -1.0e-30) return;
            TERM = TERM * STEP;
        }
    }
    return;

L900:
    NUMOUT = 1;
    return;

L950:
    printf("\n**** POTENTIAL MUST BE DEFINED ON AT LEAST TWO POINTS FOR SPLINE.\n");
    IRETUR = -1;
    goto L900;
}


// ============================================================================
// SUBROUTINE LAGRANGE
//
// Construct potential by Lagrange interpolation.
// Uses PARAM pairs and AITKEN interpolation.
//
// 4/10/78 - based on SPLINE - S.P.
// ============================================================================

void LAGRANGE(char8 ALIAS, int* MYINTS, int IPOTYP, int IREQUE, int& IRETUR,
              int IUNIT, int& NUMOUT,
              int L, double J, double RSTART, double STEPSZ, int NUMPTS,
              double* ARRAY1, double* ARRAY2,
              char* ID)
{
    double UNDEF = INTRNL.UNDEF;
    double R_arr[21], V_arr[21], A_work[31];

    NUMOUT = 0;
    IRETUR = 0;

    // Interpolation order from INTGER — INTGER(38) = MAXFUN
    int IORDER = INTGER.MAXFUN;

    // Get the defined points and values and order them.
    // Fortran PARAM(1,IV)/PARAM(2,IV) ↔ PARAMS_at(2*IV-1)/PARAMS_at(2*IV).

    int NV = 0;
    int IV = 0;

    while (IV < 10) {
        double r_val = PARAMS_at(2 * IV + 1);
        double v_val = PARAMS_at(2 * IV + 2);
        if (r_val != UNDEF && v_val != UNDEF) {
            NV = 1;
            R_arr[1] = r_val;
            V_arr[1] = v_val;
            IV++;
            break;
        }
        IV++;
    }
    if (NV == 0) goto L950;

    while (IV < 10) {
        double r_val = PARAMS_at(2 * IV + 1);
        double v_val = PARAMS_at(2 * IV + 2);
        IV++;
        if (r_val == UNDEF || v_val == UNDEF) continue;

        int ins = NV + 1;
        for (int I = 1; I <= NV; I++) {
            if (r_val < R_arr[I]) { ins = I; break; }
        }
        for (int II = NV; II >= ins; II--) {
            R_arr[II + 1] = R_arr[II];
            V_arr[II + 1] = V_arr[II];
        }
        NV++;
        R_arr[ins] = r_val;
        V_arr[ins] = v_val;
    }

    if (NV <= 2) goto L950;

    if (IREQUE == 1) goto L100;
    if (IREQUE == 2) goto L500;
    if (IREQUE == 3) goto L600;
    return;

L100:
    {
        char worknm[9] = "    SHAP";
        std::memcpy(worknm, ID, 4);
        int I = NALLOC(worknm, NUMPTS);
        double RR = RSTART;
        int LOC_I = LOCPTRS.Z[I];
        for (int II = 1; II <= NUMPTS; II++) {
            ALLOC(LOC_I - 1 + II) = RR;
            RR = RR + STEPSZ;
        }
        MYINTS[0] = I;
        MYINTS[1] = 0;
    }
    if (NUMOUT != 0) goto L900;
    return;

L500:
    printf(" THE POTENTIAL IS DEFINED BY LAGRANGE INTERPOLATION OF ORDER%2d BASED ON THE POINTS:\n",
           IORDER);
    for (int I = 1; I <= NV; I++) {
        printf("%10.3f%13.5g", R_arr[I], V_arr[I]);
        if (I % 5 == 0 || I == NV) printf("\n");
    }
    goto L900;

L600:
    {
        int LL = LOCPTRS.Z[MYINTS[0]] - 1;

        double RR = RSTART;
        int NSTRT = 1;
        int NLAST = 1;
        for (int N = 1; N <= NUMPTS; N++) {
            if (RR < R_arr[1]) {
                ARRAY1[N] = -V_arr[1];
                NSTRT = N + 1;
            } else if (RR > R_arr[NV]) {
                ARRAY1[N] = 0.0;
            } else {
                NLAST = N;
            }
            RR = RR + STEPSZ;
        }
        int NOUT = NLAST + 1 - NSTRT;

        // Does the potential change sign
        int LNORDR = IORDER;
        double VSIGN = -DSIGN(1.0, V_arr[NV]);
        int LNSTRT = 1;

        int change_i = 0;
        for (int II = 2; II <= NV; II++) {
            int I = NV + 2 - II;
            if (VSIGN * V_arr[I - 1] < 0.0) {
                change_i = I;
                break;
            }
        }

        if (change_i > 0) {
            // Yes, use straight interpolation up to the last change
            LNSTRT = change_i;
            int NSTRT2 = (int)((R_arr[LNSTRT] - RSTART) / STEPSZ + 2);
            int NFAIL;
            double WORST;
            AITKEN(MIN0(IORDER, LNSTRT - 1), 0.0, 0.0, LNSTRT, R_arr, V_arr,
                   NSTRT2 - NSTRT, ALLOC_base(LL + NSTRT), &ARRAY1[NSTRT],
                   NFAIL, WORST, A_work);
            for (int I = NSTRT; I <= NSTRT2; I++) {
                ARRAY1[I] = -ARRAY1[I];
            }

            LNORDR = MIN0(IORDER, NV - LNSTRT);
            if (LNORDR <= 0) return;
            NSTRT = NSTRT2;
            NOUT = NLAST + 1 - NSTRT2;
        }

        // Interpolate the logs
        for (int I = LNSTRT; I <= NV; I++) {
            V_arr[I] = DLOG(DMAX1(DABS(V_arr[I]), 1.0e-15));
        }

        {
            int NFAIL;
            double WORST;
            AITKEN(LNORDR, 0.0, 0.0, NV - LNSTRT + 1,
                   &R_arr[LNSTRT], &V_arr[LNSTRT],
                   NOUT, ALLOC_base(LL + NSTRT), &ARRAY1[NSTRT],
                   NFAIL, WORST, A_work);
        }

        for (int I = NSTRT; I <= NLAST; I++) {
            double X = ARRAY1[I];
            X = DMAX1(X, -69.0);
            X = DMIN1(X, 69.0);
            ARRAY1[I] = VSIGN * DEXP(X);
        }

        // Fill in the end with exponential decay
        double TERM = ARRAY1[NLAST];
        NLAST = NLAST + 1;
        if (NLAST > NUMPTS) return;
        double STEP = TERM / ARRAY1[NLAST - 2];
        TERM = TERM * STEP;
        if (1.0 - STEP < 1.0e-5) return;
        for (int I = NLAST; I <= NUMPTS; I++) {
            ARRAY1[I] = TERM;
            if (TERM > -1.0e-30) return;
            TERM = TERM * STEP;
        }
    }
    return;

L900:
    NUMOUT = 1;
    return;

L950:
    printf("\n**** POTENTIAL MUST BE DEFINED ON AT LEAST TWO POINTS FOR LAGRANGE.\n");
    IRETUR = -1;
    goto L900;
}


// ============================================================================
// Stub routines for linkules that need WOODSX from Phase 8
// These are complete translations but depend on WOODSX/DEFWOO.
// ============================================================================

// BKGPTElp, DEFORMEd, LTSTELp, and TWOSHApe are more complex and will be
// added when their dependencies (WOODSX, DEFWOO) are available from Phase 8.
// For now, they print a message and return an error.

void BKGPTElp(char8 ALIAS, int* MYINTS, int IPOTYP, int IREQUE, int& IRETUR,
              int IUNIT, int& NUMOUT,
              int L, double J, double RSTART, double STEPSZ, int NUMPTS,
              double* ARRAY1, double* ARRAY2,
              char* ID, int* IPARAM)
{
    printf("\n**** BKGPTElp LINKULE NOT YET FULLY OPERATIONAL (needs WOODSX from Phase 8).\n");
    IRETUR = -1;
    NUMOUT = 1;
}

void DEFORMEd(char8 ALIAS, int* MYINTS, int IPOTYP, int IREQUE, int& IRETUR,
              int IUNIT, int& NUMOUT,
              int L, double J, double RSTART, double STEPSZ, int NUMPTS,
              double* ARRAY1, double* ARRAY2,
              char* ID, int* IPARAM)
{
    printf("\n**** DEFORMEd LINKULE NOT YET FULLY OPERATIONAL (needs DEFWOO from Phase 8).\n");
    IRETUR = -1;
    NUMOUT = 1;
}

void LTSTELp(char8 ALIAS, int* MYINTS, int IPOTYP, int IREQUE, int& IRETUR,
             int IUNIT, int& NUMOUT,
             int L, double J, double RSTART, double STEPSZ, int NUMPTS,
             double* ARRAY1, double* ARRAY2,
             char* ID, int* IPARAM)
{
    printf("\n**** LTSTELp LINKULE NOT YET FULLY OPERATIONAL (needs WOODSX from Phase 8).\n");
    IRETUR = -1;
    NUMOUT = 1;
}

void TWOSHApe(char8 ALIAS, int* MYINTS, int IPOTYP, int IREQUE, int& IRETUR,
              int IUNIT, int& NUMOUT,
              int L, double J, double RSTART, double STEPSZ, int NUMPTS,
              double* ARRAY1, double* ARRAY2,
              char* ID)
{
    printf("\n**** TWOSHApe LINKULE NOT YET FULLY OPERATIONAL.\n");
    IRETUR = -1;
    NUMOUT = 1;
}
