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
        // INTRPC expects 1-based ptrs (YS[1] = first element). ARRAY1 is already 1-based
        // (caller passes ALLOC_base form), so &ARRAY1[NSTRT] - 1 yields YS[1] = ARRAY1[NSTRT].
        INTRPC(NV, R_arr, V_arr, B, C_arr, D_arr,
               NOUT, ALLOC_base(LL + NSTRT), &ARRAY1[NSTRT] - 1);

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

    // Fortran: IORDER = INTGER(38) = PARITS(1). The LAGRANGE linkule co-opts the
    // PARITS(1) slot to carry the interpolation order (not INTGER.MAXFUN as the
    // original comment suggested). With no explicit order set, it defaults to 0
    // (constant interpolation).
    int IORDER = INTGER.PARITS[1];

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

        // Fortran: IF (VSIGN*V(I-1) .GE. 0) GO TO 660 — exits when V(I-1) has
        // OPPOSITE sign to V(NV) (i.e., a sign change has been found).
        int change_i = 0;
        for (int II = 2; II <= NV; II++) {
            int I = NV + 2 - II;
            if (VSIGN * V_arr[I - 1] >= 0.0) {
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
            // AITKEN expects 1-based XOUT/FOUT (FOUT[1] = first element).
            // ARRAY1 is already 1-based, so &ARRAY1[NSTRT] - 1 makes FOUT[1] = ARRAY1[NSTRT].
            AITKEN(MIN0(IORDER, LNSTRT - 1), 0.0, 0.0, LNSTRT, R_arr, V_arr,
                   NSTRT2 - NSTRT, ALLOC_base(LL + NSTRT), &ARRAY1[NSTRT] - 1,
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
            // Same 1-based shift: ARRAY1, R_arr, V_arr are all 1-based; pass pointer - 1
            // so AITKEN sees XIN[1]=R_arr[LNSTRT], FIN[1]=V_arr[LNSTRT], FOUT[1]=ARRAY1[NSTRT].
            AITKEN(LNORDR, 0.0, 0.0, NV - LNSTRT + 1,
                   &R_arr[LNSTRT] - 1, &V_arr[LNSTRT] - 1,
                   NOUT, ALLOC_base(LL + NSTRT), &ARRAY1[NSTRT] - 1,
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

// ============================================================================
// SUBROUTINE BKGPTElp — Baltz-Kauffmann-Glendenning-Pruess 2+ TELP linkule
// Translated from linkulesfitters.f L190-L501.
// Computes the imaginary potential as WS + surface + 2+ TELP (Baltz et al.).
// Note: Fortran FKANDM(20) is a REAL*8 view of /KANDM/; FKANDM(ICHANW) = AKS(ICHANW),
// FKANDM(ICHANW+16) = ETAS(ICHANW). Fortran INTRNL(5) as INTEGER view = ICHANB.
// L-dependent linkule — IRETUR=+1 on IREQUE=1.
// NOT yet exercised by any regression test (untested translation).
// ============================================================================

void BKGPTElp(char8 ALIAS, int* MYINTS, int IPOTYP, int IREQUE, int& IRETUR,
              int IUNIT, int& NUMOUT,
              int L, double J, double RSTART, double STEPSZ, int NUMPTS,
              double* ARRAY1, double* ARRAY2,
              char* ID, int* IPARAM)
{
    static const char WORDS[3][9] = { "        ", "PROJCTIL", "  TARGET" };

    NUMOUT = 0;
    IRETUR = 0;

    auto ATERM_fn = [](double DL, double ETA, double ETASQ, double AK) {
        return ETA * AK*AK / (DL*DL) *
               (ETA * (3.0*DL*DL + ETASQ) / std::pow(DL*DL + ETASQ, 2.0)
                - DATAN(DL / ETA) / DL);
    };
    auto BTERM_fn = [](double DL, double ETA, double ETASQ, double AK) {
        return 4.0 * ETA * AK * DL*DL / std::pow(DL*DL + ETASQ, 2.0);
    };
    auto CTERM_fn = [](double DL, double ETASQ) {
        return 2.0 * std::pow(DL, 4.0) / std::pow(DL*DL + ETASQ, 2.0);
    };

    int    LWRK = 0, IWRK = 0;
    double WP = 0.0, WPP = 0.0, AK = 0.0, ETA = 0.0;
    double AC = 0.0, BC = 0.0, CC = 0.0;
    double RC = 0.0, BETA = 0.0, RD = 0.0, FATRC = 0.0, DLC = 0.0;
    int    IPORT = 0;
    double R0MASS = 0.0;
    double VI = 0.0, AI = 0.0, RI = 0.0, POWER = 0.0;
    double VSI = 0.0, RSI = 0.0, ASI = 0.0;
    double G2_val = 1.0;

    if (IREQUE == 4) goto L700;

    {
        double UNDEF  = INTRNL.UNDEF;
        int    NOTDEF = NOTDEF_INT;

        VI    = TEMPVS_arr()[1];          // TEMPVS(2)
        AI    = TEMPVS_arr()[5];          // TEMPVS(6)
        RI    = FLOAT_common.RI;          // FLOAT(47)
        POWER = FLOAT_common.POWIM;       // FLOAT(103)
        VSI   = FLOAT_common.VSI;         // FLOAT(71)
        RSI   = FLOAT_common.RSI;         // FLOAT(69)
        ASI   = FLOAT_common.ASI;         // FLOAT(68)
        G2_val = PARAMS_at(20);
        if (G2_val == UNDEF) G2_val = 1.0;
        int ICHANW_local = INTRNL.ICHANB; // Fortran INTRNL(5) as INTEGER view = ICHANB
        AK   = AKS()[ICHANW_local];       // FKANDM(ICHANW)
        RC   = FLOAT_common.RC;           // FLOAT(45)
        int ZP = INTGER.IZP;
        int ZT = INTGER.IZT;
        ETA  = KANDM.ETAS[ICHANW_local];  // FKANDM(ICHANW+16) = ETAS(ICHANW)
        double ALPHA = 1.0 / CNSTNT.AFINE;
        double PI    = CNSTNT.PI;
        double HBARC = CNSTNT.HBARC;
        double E     = FLOAT_common.E;
        double AM    = FLOAT_common.AM;
        R0MASS = INTRNL.R0MASS;

        if (IREQUE == 1) {
            IPORT = IPARAM[4];            // Fortran IPARAM(5)
            if (IPORT == NOTDEF) {
                IPORT = 0;
                if ((int)JBLOCK.JS[2] != (int)JBLOCK.JS[3]) IPORT = 1;
                if ((int)JBLOCK.JS[4] != (int)JBLOCK.JS[5]) IPORT = 2;
            }
            if (IPORT != 1 && IPORT != 2) goto L930;
            double R2FAC = std::pow(FLOAT_arr(39 + IPORT), 1.0/3.0) / R0MASS;
            int IB = NAMLOC("BETACOUL");
            if (IB == 0) goto L940;
            if (VI != 0.0 && (RI == UNDEF || AI == UNDEF)) goto L950;
            if (VSI != 0.0 && (RSI == UNDEF || ASI == UNDEF)) goto L945;

            char wname[9] = "    SHAP";
            std::memcpy(wname, ID, 4);
            IWRK = NALLOC(wname, NUMPTS + 10);
            LWRK = LOCPTRS.Z[IWRK];
            ALLOC(LWRK)     = R2FAC;
            ALLOC(LWRK + 1) = ALLOC(LOCPTRS.Z[IB]);
            goto L300;
        }

        // IREQUE = 2 or 3 — use cached MYINTS
        IWRK = MYINTS[0];
        LWRK = LOCPTRS.Z[IWRK];
        if (VI  != 0.0 && (AI  <= 0.0 || RI  <= 0.0)) goto L960;
        if (VSI != 0.0 && (ASI <= 0.0 || RSI <= 0.0)) goto L960;
        double RC2  = RC * ALLOC(LWRK);
        BETA = ALLOC(LWRK + 1);
        RD   = ZP * ZT * ALPHA * HBARC / E;
        WP   = 3.0 * std::pow(ZP * ZT * ALPHA * BETA, 2.0) * AM * G2_val
             * std::pow(RC2, 4.0) / (50.0 * PI * AK);
        WPP  = 3.0 * WP / 8.0;
        ALLOC(LWRK + 2) = ETA;
        ALLOC(LWRK + 3) = WPP;
        ALLOC(LWRK + 8) = AK;
        double ETASQ = ETA * ETA;

        double DLCSQ = AK * RC * (AK * RC - 2.0 * ETA);
        if (DLCSQ < 0.25) DLCSQ = 0.25;
        DLC = DSQRT(DLCSQ);
        AC = ATERM_fn(DLC, ETA, ETASQ, AK);
        BC = BTERM_fn(DLC, ETA, ETASQ, AK);
        CC = CTERM_fn(DLC, ETASQ);
        ALLOC(LWRK + 4) = AC;
        ALLOC(LWRK + 5) = BC;
        ALLOC(LWRK + 6) = CC;
        FATRC = ((CC / RC + BC) / RC + AC) / std::pow(RC, 3.0);
    }

L300:
    if (IREQUE == 1) {
        MYINTS[0] = IWRK;
        MYINTS[1] = 0;
        IRETUR = +1;
        if (NUMOUT != 0) goto L900;
        return;
    }
    if (IREQUE == 2) goto L500;
    if (IREQUE == 3) goto L600;
    return;

L500:
    if (VI != 0.0) {
        double R0 = RI / R0MASS;
        std::printf(" W.S. WELL:        %14.3g     %7.4f     %7.4f       %7.4f      %7.3f\n",
                    VI, RI, AI, R0, POWER);
    }
    if (VSI != 0.0) {
        double R0 = RSI / R0MASS;
        std::printf(" SURFACE ABSORPTION:%14.3g     %7.4f     %7.4f       %7.4f\n",
                    VSI, RSI, ASI, R0);
    }
    {
        double WPRC5 = WP / std::pow(RC, 5.0);
        std::printf(" %.8s 2+ BKGP TELP WITH BETACOULOMB =%7.4f     G2 =%12.5g"
                    "     RC =%7.4f     RD =%7.4f     WP/RC**5 =%12.5g MEV\n",
                    WORDS[IPORT], BETA, G2_val, RC, RD, WPRC5);
        double FFC = FATRC * std::pow(RC, 5.0);
        double AAC = AC * RC * RC;
        double BBC = BC * RC;
        std::printf(" SIMPLE L-HAT CRITICAL =%8.1f     F(LCRIT, RC)*RC**5 =%13.5g"
                    "     TERMS =%14.5g%14.5g%14.5g\n",
                    DLC, FFC, AAC, BBC, CC);
    }
    goto L900;

L600:
    {
        int LL = LWRK + 10;
        int N1 = 0, N2 = 0;
        WOODSX(NUMPTS, RSTART, STEPSZ, ALLOC_base(LL), N1, N2, 1, VI, RI, AI, POWER);
        for (int I = 1; I <= NUMPTS; I++) ARRAY1[I] = ALLOC(LL - 1 + I);

        if (VSI != 0.0) {
            WOODSX(NUMPTS, RSTART, STEPSZ, ALLOC_base(LL), N1, N2, 3, VSI, RSI, ASI, 1.0);
            for (int I = N1; I <= N2; I++) ARRAY1[I] = ARRAY1[I] + ALLOC(LL - 1 + I);
        }

        N1 = (int)((RC - RSTART) / STEPSZ + 1.5);
        N1 = MIN0(N1, NUMPTS);
        ALLOC(LWRK + 7) = (double)N1;
        ALLOC(LWRK + 9) = RC;
        for (int I = 1; I <= N1; I++) ARRAY1[I] = ARRAY1[I] - WPP * FATRC;

        N1 = N1 + 1;
        if (N1 > NUMPTS) return;
        double R = RSTART + (N1 - 1) * STEPSZ;
        for (int I = N1; I <= NUMPTS; I++) {
            double RINV = 1.0 / R;
            ARRAY1[I] = ARRAY1[I]
                - WPP * ((CC * RINV + BC) * RINV + AC) * RINV * RINV * RINV;
            R += STEPSZ;
        }
    }
    return;

L700:
    {
        LWRK = LOCPTRS.Z[MYINTS[0]];
        ETA = ALLOC(LWRK + 2);
        WPP = ALLOC(LWRK + 3) * ARRAY2[1];   // ARRAY2(1) Fortran = ARRAY2[1] (caller pre-shifts)
        AC  = ALLOC(LWRK + 4);
        BC  = ALLOC(LWRK + 5);
        CC  = ALLOC(LWRK + 6);
        int N1 = (int)ALLOC(LWRK + 7);
        AK  = ALLOC(LWRK + 8);
        RC  = ALLOC(LWRK + 9);
        double ETASQ = ETA * ETA;
        double DL = L + 0.5;
        double A_v = ATERM_fn(DL, ETA, ETASQ, AK);
        double B_v = BTERM_fn(DL, ETA, ETASQ, AK);
        double C_v = CTERM_fn(DL, ETASQ);
        FATRC = (((C_v - CC) / RC + (B_v - BC)) / RC + (A_v - AC)) / std::pow(RC, 3.0);
        for (int I = 1; I <= N1; I++) ARRAY1[I] = ARRAY1[I] - WPP * FATRC;
        double R = RSTART + N1 * STEPSZ;
        N1 = N1 + 1;
        if (N1 > NUMPTS) return;
        for (int I = N1; I <= NUMPTS; I++) {
            double RINV = 1.0 / R;
            ARRAY1[I] = ARRAY1[I]
                - WPP * (((C_v - CC) * RINV + (B_v - BC)) * RINV + (A_v - AC)) * RINV * RINV * RINV;
            R += STEPSZ;
        }
    }
    return;

L900:
    NUMOUT = 1;
    return;

L930:
    std::printf("0**** IPARAM5 =%3d MUST BE 1 OR 2 TO INDICATE PROJECTILE OR TARGET EXCITATION FOR BKGPTELP.\n", IPORT);
    IRETUR = -1; goto L900;
L940:
    std::printf("0**** BETACOULOMB MUST BE DEFINED FOR BKGPTELP.\n");
    IRETUR = -1; goto L900;
L945:
    std::printf("0**** RSI (OR RSI0) AND ASI MUST BE DEFINED.\n");
    IRETUR = -1; goto L900;
L950:
    std::printf("0**** RI (OR RI0) AND AI MUST BE DEFINED.\n");
    IRETUR = -1; goto L900;
L960:
    IRETUR = -1;
    return;
}

// ============================================================================
// SUBROUTINE DEFORMEd — Deformed Woods-Saxon optical potential
// Translated from linkulesfitters.f L504-L729.
// ============================================================================

void DEFORMEd(char8 ALIAS, int* MYINTS, int IPOTYP, int IREQUE, int& IRETUR,
              int IUNIT, int& NUMOUT,
              int L, double J, double RSTART, double STEPSZ, int NUMPTS,
              double* ARRAY1, double* ARRAY2,
              char* ID, int* IPARAM)
{
    static const char DEFNAM[3][9] = { "PROJECTI", "TARGET  ", "SUMMED  " };

    double UNDEF  = INTRNL.UNDEF;
    int    NOTDEF = NOTDEF_INT;
    double PI     = CNSTNT.PI;
    double BIGLOG = CNSTNT.BIGLOG;

    NUMOUT = 0;
    IRETUR = 0;

    // Local arrays: BETAS(2,20), RBETAS(2,20), IBETA(2), LDEFMX(2), RDEFS(2), LXS(2)
    // (Fortran 1-based, dim2 indexed by I=1..2, dim1 by LX=LXMIN..LXMAX up to 20)
    double BETAS[3][21];   // BETAS[I][LX], I=1..2, LX=1..20
    double RBETAS[3][21];
    int    IBETA[3]   = {0,0,0};
    int    LDEFMX[3]  = {0,0,0};
    double RDEFS[3]   = {0,0,0};
    int    LXS[3]     = {0,0,0};

    double V, A, R, POWER;
    if (IPOTYP == 13) {
        // Deformed surface potential (derivative of WS)
        V = FLOAT_arr(71);
        R = FLOAT_arr(69);
        A = FLOAT_arr(68);
        POWER = 1.0;
    } else {
        V = TEMPVS_arr()[IPOTYP - 1];          // TEMPVS(IPOTYP)
        A = TEMPVS_arr()[3 + IPOTYP];          // TEMPVS(4+IPOTYP)
        R = FLOAT_common.R;                    // FLOAT(43)
        if (IPOTYP == 2) R = FLOAT_common.RI;  // FLOAT(47)
        if (IPOTYP == 1) POWER = FLOAT_common.POWRL;  // FLOAT(102)
        else             POWER = FLOAT_common.POWIM;  // FLOAT(103)
    }

    double VV = 0.0;
    int    LXOUT = 0;
    int    LXMIN = 2;
    int    NCOSIN = 10;
    int    IPORT = 0;

    if (V != 0.0) {
        if (R == UNDEF || A == UNDEF || V == UNDEF) goto L950;
        if (R <= 0.0 || A <= 0.0) goto L960;

        LXOUT = INTGER.LX;                     // INTGER(5)
        if (LXOUT == -1 || LXOUT == NOTDEF) LXOUT = 0;
        LXS[1] = 0;
        LXS[2] = 0;
        LXMIN = IPARAM[3];                     // Fortran IPARAM(4) = IPARM4 (caller passes &IPARM1)
        if (LXMIN == NOTDEF) LXMIN = 2;
        IBETA[1] = MYINTS[0];                  // MYINTS(1) — C++ 0-based MYINTS[0]
        IBETA[2] = MYINTS[1];                  // MYINTS(2)

        if (IREQUE == 1) {
            // Find BETAP, BETAT
            IBETA[1] = NAMLOC("BETAP   ");
            IBETA[2] = NAMLOC("BETAT   ");
            if (IBETA[1] + IBETA[2] == 0) goto L970;
        }

        // Get the betas and compute R*beta
        double R0MASS = INTRNL.R0MASS;
        for (int I = 1; I <= 2; I++) {
            LDEFMX[I] = 0;
            if (IBETA[I] == 0) continue;
            LDEFMX[I] = LXMIN + LENGTH.LENG[IBETA[I]] - 1;
            int LXMAX_loc = LDEFMX[I];
            double RDEF = R * std::pow(FLOAT_arr(39 + I), 1.0/3.0) / R0MASS;
            RDEFS[I] = RDEF;
            int LBETA = LOCPTRS.Z[IBETA[I]];
            for (int LX = LXMIN; LX <= LXMAX_loc; LX++) {
                double BETA = ALLOC(LBETA);
                BETAS[I][LX]  = BETA;
                RBETAS[I][LX] = RDEF * BETA;
                LBETA++;
            }
        }

        NCOSIN = INTGER.NCOSIN;                // INTGER(11)
        if (NCOSIN == NOTDEF) NCOSIN = 10;
        if (NCOSIN > 100)     NCOSIN = 10;
        VV = V;
        IPORT = 0;
        if (INTGER.LX == -1 || INTGER.LX == NOTDEF) goto L100;

        // Determin what projection is needed for excitation
        if (IREQUE != 3) goto L100;
        IPORT = (INTGER_arr_f(58) % 100) - 2;  // JBLOCK(1) — actually JB(1); let me use JBLOCK.J directly
        // (JBLOCK(1) is J — see ptolemy_commons.h JblockCommon: COMMON /JBLOCK/ J, JS(5), ...)
        // Re-do using direct JBLOCK access:
        IPORT = ((int)JBLOCK.J % 100) - 2;
        if (IPORT >= 3) goto L80;
        if (IPORT < 1 || IBETA[IPORT] == 0) goto L965;
        LXS[IPORT] = LXOUT;
        VV = -V / RDEFS[IPORT];
        goto L100;

    L80:
        LXS[1] = ((int)JBLOCK.J / 100)   % 100;
        LXS[2] = ((int)JBLOCK.J / 10000) % 100;
        VV = -V / (RDEFS[1] * RDEFS[2]);
    }

L100:
    if (IREQUE == 1) {
        MYINTS[0] = IBETA[1];
        MYINTS[1] = IBETA[2];
        return;
    }
    if (IREQUE == 2) goto L500;
    if (IREQUE == 3) goto L600;
    return;

L500:
    if (V == 0.0) {
        std::printf(" THE POTENTIAL IS IDENTICALLY ZERO.\n");
        goto L900;
    }
    if (IPOTYP != 13)
        std::printf(" THE POTENTIAL IS A DEFORMED WOODS-SAXON:\n");
    else
        std::printf(" THE POTENTIAL IS A DEFORMED DERIVATIVE OF A WOODS-SAXON:\n");
    std::printf(" V =%7.2f   R =%7.4f   A =%7.4f     POWER =%8.4f\n", V, R, A, POWER);
    for (int I = 1; I <= 2; I++) {
        if (IBETA[I] == 0) continue;
        int LXMAX_loc = LDEFMX[I];
        std::printf("0THE %.8s RADIUS OF%8.4f FM. IS BEING DEFORMED.\n"
                    " LX   BETA   DEFORMATION LENGTH\n", DEFNAM[I-1], RDEFS[I]);
        for (int LX = LXMIN; LX <= LXMAX_loc; LX++) {
            std::printf("%3d%8.4f%15.4f\n", LX, BETAS[I][LX], RBETAS[I][LX]);
        }
    }
    std::printf("%5d-POINT GAUSS INTEGRATION IS BEING USED.\n", NCOSIN);
    goto L900;

L600:
    {
        int ITYPE = 1;
        if (IPOTYP == 13) ITYPE = 3;
        int N1, N2;
        // DEFWOO's internal macro `RBETAS(i,lx) RBETAS_flat[((lx)-1)*2 + (i)]`
        // uses 1-based access — RBETAS_flat[1] is RBETAS(1,1). Same for LDEFMX/LXS[I=1..2].
        // expects 1-based access: RBETAS_flat[1] = RBETAS(1,1), [2] = (2,1), [3] = (1,2), etc.
        // Allocate one extra slot at index 0 for the padding so the math lines up.
        double RBETAS_flat[2 * 20 + 1] = {0};
        for (int LX = 1; LX <= 20; LX++) {
            RBETAS_flat[(LX-1)*2 + 1] = RBETAS[1][LX];   // RBETAS(1, LX) — projectile
            RBETAS_flat[(LX-1)*2 + 2] = RBETAS[2][LX];   // RBETAS(2, LX) — target
        }
        // LDEFMX_arr and LXS_arr — DEFWOO indexes LDEFMX[I] for I=1..2 so use 1-based with padding
        int LDEFMX_arr[3] = { 0, LDEFMX[1], LDEFMX[2] };
        int LXS_arr[3]    = { 0, LXS[1],    LXS[2]    };
        DEFWOO(NUMPTS, RSTART, STEPSZ, ARRAY1, N1, N2,
               ITYPE, VV, R, A, POWER, IPORT, LXMIN, LDEFMX_arr,
               RBETAS_flat, LXS_arr, NCOSIN);

        if (VV == 0.0) goto L900;
        VV = VV / V;
        if (IPORT == 0) {
            std::printf(" DEFORMED OPTICAL POTENTIAL WEIGHTED BY%15.5g"
                        " WAS COMPUTED WITH LXP AND LXT MAXIMUM =%3d%3d\n",
                        VV, LDEFMX[1], LDEFMX[2]);
        } else if (IPORT < 3) {
            std::printf(" THE LAMBDA =%2d PROJECTION OF THE DEFORMED POTENTIAL WEIGHTED BY"
                        "%15.5g WAS COMPUTED WITH A DEFORMED %.8s RADIUS.\n",
                        LXOUT, VV, DEFNAM[IPORT-1]);
        } else {
            std::printf(" THE MUTUAL DEFORMED POTENTIAL WEIGHTED BY%15.5g"
                        " WAS COMPUTED FOR LXP, LXT =%3d%3d\n",
                        VV, LXS[1], LXS[2]);
        }
        // Indicate that no derivative is needed in the form factor
        if (IPORT > 0) SWITCH.KFRMTP = 1;
    }
    // fallthrough to L900

L900:
    NUMOUT = 1;
    return;

L950:
    std::printf("0**** R, V, AND A MUST BE DEFINED.\n");
    IRETUR = -1;
    goto L900;

L960:
    std::printf("0*** R OR A IS INVALID:%20.10g%20.10g\n", R, A);
    IRETUR = -1;
    goto L900;

L965:
    std::printf("0*** WRONG NUCLEUS DEFORMED:%15d%8d%8d%8d%8d%8d%8d\n",
                (int)JBLOCK.J, IPORT, IREQUE, LXS[1], LXS[2], IBETA[1], IBETA[2]);
    IRETUR = -1;
    goto L900;

L970:
    std::printf("0*** BETAP OR BETAT MUST BE DEFINED\n");
    IRETUR = -1;
    goto L900;
}

// ============================================================================
// SUBROUTINE LTSTELp — Love-Terasawa-Satchler 2+ TELP linkule
// Translated from linkulesfitters.f L1442-L1696.
// Like BKGPTELP but with the LTS TELP formula; NO L-dependence (IRETUR=0).
// NOT yet exercised by any regression test (untested translation).
// ============================================================================

void LTSTELp(char8 ALIAS, int* MYINTS, int IPOTYP, int IREQUE, int& IRETUR,
             int IUNIT, int& NUMOUT,
             int L, double J, double RSTART, double STEPSZ, int NUMPTS,
             double* ARRAY1, double* ARRAY2,
             char* ID, int* IPARAM)
{
    static const char WORDS[3][9] = { "        ", "PROJCTIL", "  TARGET" };

    NUMOUT = 0;
    IRETUR = 0;

    double UNDEF  = INTRNL.UNDEF;
    int    NOTDEF = NOTDEF_INT;

    double VI    = TEMPVS_arr()[1];
    double AI    = TEMPVS_arr()[5];
    double RI    = FLOAT_common.RI;
    double POWER = FLOAT_common.POWIM;
    double VSI   = FLOAT_common.VSI;
    double RSI   = FLOAT_common.RSI;
    double ASI   = FLOAT_common.ASI;
    double G2_val = PARAMS_at(20);
    if (G2_val == UNDEF) G2_val = 1.0;
    int ICHANW_local = INTRNL.ICHANB;
    double AK   = AKS()[ICHANW_local];
    double RC   = FLOAT_common.RC;
    int    ZP   = INTGER.IZP;
    int    ZT   = INTGER.IZT;
    double ALPHA = 1.0 / CNSTNT.AFINE;
    double PI    = CNSTNT.PI;
    double HBARC = CNSTNT.HBARC;
    double E     = FLOAT_common.E;
    double AM    = FLOAT_common.AM;
    double R0MASS = INTRNL.R0MASS;

    int IWRK = 0, LWRK = 0;
    int IPORT = 0;
    double WP = 0.0, RC2 = 0.0, BETA = 0.0, RD = 0.0;

    if (IREQUE == 1) {
        IPORT = IPARAM[4];
        if (IPORT == NOTDEF) {
            IPORT = 0;
            if ((int)JBLOCK.JS[2] != (int)JBLOCK.JS[3]) IPORT = 1;
            if ((int)JBLOCK.JS[4] != (int)JBLOCK.JS[5]) IPORT = 2;
        }
        if (IPORT != 1 && IPORT != 2) goto LL930;
        double R2FAC = std::pow(FLOAT_arr(39 + IPORT), 1.0/3.0) / R0MASS;
        int IB = NAMLOC("BETACOUL");
        if (IB == 0) goto LL940;
        if (VI != 0.0 && (RI == UNDEF || AI == UNDEF)) goto LL950;
        if (VSI != 0.0 && (RSI == UNDEF || ASI == UNDEF)) goto LL945;
        char wname[9] = "    SHAP";
        std::memcpy(wname, ID, 4);
        IWRK = NALLOC(wname, NUMPTS + 10);
        LWRK = LOCPTRS.Z[IWRK];
        ALLOC(LWRK)     = R2FAC;
        ALLOC(LWRK + 1) = ALLOC(LOCPTRS.Z[IB]);
        goto LL300;
    }

    IWRK = MYINTS[0];
    LWRK = LOCPTRS.Z[IWRK];
    if (VI  != 0.0 && (AI  <= 0.0 || RI  <= 0.0)) goto LL960;
    if (VSI != 0.0 && (ASI <= 0.0 || RSI <= 0.0)) goto LL960;
    RC2  = RC * ALLOC(LWRK);
    BETA = ALLOC(LWRK + 1);
    RD   = ZP * ZT * ALPHA * HBARC / E;
    WP   = 3.0 * std::pow(ZP * ZT * ALPHA * BETA, 2.0) * AM * G2_val
         * std::pow(RC2, 4.0) / (50.0 * PI * AK);

LL300:
    if (IREQUE == 1) {
        MYINTS[0] = IWRK;
        MYINTS[1] = 0;
        IRETUR = 0;
        if (NUMOUT != 0) goto LL900;
        return;
    }
    if (IREQUE == 2) goto LL500;
    if (IREQUE == 3) goto LL600;
    return;

LL500:
    if (VI != 0.0) {
        double R0 = RI / R0MASS;
        std::printf(" W.S. WELL:        %14.3g     %7.4f     %7.4f       %7.4f      %7.3f\n",
                    VI, RI, AI, R0, POWER);
    }
    if (VSI != 0.0) {
        double R0 = RSI / R0MASS;
        std::printf(" SURFACE ABSORPTION:%14.3g     %7.4f     %7.4f       %7.4f\n",
                    VSI, RSI, ASI, R0);
    }
    {
        double WPRC5 = WP / std::pow(RC, 5.0);
        std::printf(" %.8s 2+ LTS TELP WITH BETACOULOMB =%7.4f     G2 =%12.5g"
                    "     RC =%7.4f     RD =%7.4f     WP/RC**5 =%12.5g MEV\n",
                    WORDS[IPORT], BETA, G2_val, RC, RD, WPRC5);
    }
    goto LL900;

LL600:
    {
        int LL = LWRK + 10;
        int N1 = 0, N2 = 0;
        WOODSX(NUMPTS, RSTART, STEPSZ, ALLOC_base(LL), N1, N2, 1, VI, RI, AI, POWER);
        for (int I = 1; I <= NUMPTS; I++) ARRAY1[I] = ALLOC(LL - 1 + I);

        if (VSI != 0.0) {
            WOODSX(NUMPTS, RSTART, STEPSZ, ALLOC_base(LL), N1, N2, 3, VSI, RSI, ASI, 1.0);
            for (int I = N1; I <= N2; I++) ARRAY1[I] = ARRAY1[I] + ALLOC(LL - 1 + I);
        }

        // LTS TELP: Coulomb braking term
        double X = 1.0 / DSQRT(0.10);
        N1 = (int)((RD - RSTART) / (0.90 * STEPSZ) + 1.5);
        N1 = MIN0(N1, NUMPTS);
        for (int I = 1; I <= N1; I++) ALLOC(LL - 1 + I) = X;
        double R = N1 * STEPSZ + RSTART;
        for (int I = N1; I <= NUMPTS; I++) {
            ALLOC(LL - 1 + I) = 1.0 / DSQRT(1.0 - RD / R);
            R += STEPSZ;
        }

        // LTS TELP term (two regions)
        N1 = (int)((RC - RSTART) / STEPSZ + 1.5);
        N1 = MIN0(N1, NUMPTS);
        R = RSTART;
        for (int I = 1; I <= N1; I++) {
            ARRAY1[I] = ARRAY1[I]
                - (2.0 * WP / (3.0 * std::pow(RC, 9.0))) * ALLOC(LL - 1 + I) * std::pow(R, 4.0);
            R += STEPSZ;
        }

        N1 = N1 + 1;
        if (N1 > NUMPTS) return;
        for (int I = N1; I <= NUMPTS; I++) {
            double RCoR = RC / R;
            ARRAY1[I] = ARRAY1[I]
                - WP * ALLOC(LL - 1 + I) * (1.0 - RCoR*RCoR * (2.0/7.0 + RCoR*RCoR * (1.0/21.0)))
                  / std::pow(R, 5.0);
            R += STEPSZ;
        }
    }
    return;

LL900:
    NUMOUT = 1;
    return;

LL930:
    std::printf("0**** IPARAM5 =%3d MUST BE 1 OR 2 TO INDICATE PROJECTILE OR TARGET EXCITATION FOR LTSTELP.\n", IPORT);
    IRETUR = -1; goto LL900;
LL940:
    std::printf("0**** BETACOULOMB MUST BE DEFINED FOR LTSTELP.\n");
    IRETUR = -1; goto LL900;
LL945:
    std::printf("0**** RSI (OR RSI0) AND ASI MUST BE DEFINED.\n");
    IRETUR = -1; goto LL900;
LL950:
    std::printf("0**** RI (OR RI0) AND AI MUST BE DEFINED.\n");
    IRETUR = -1; goto LL900;
LL960:
    IRETUR = -1;
    return;
}

// ============================================================================
// SUBROUTINE TWOSHApe — Use two fixed shapes scaled by depth parameters
// Translated from linkulesfitters.f L3438-L3691.
// Potential = CONS * SHAPE1 + CONS2 * SHAPE2 where SHAPE1/2 are stored arrays.
// ============================================================================

void TWOSHApe(char8 ALIAS, int* MYINTS, int IPOTYP, int IREQUE, int& IRETUR,
              int IUNIT, int& NUMOUT,
              int L, double J, double RSTART, double STEPSZ, int NUMPTS,
              double* ARRAY1, double* ARRAY2,
              char* ID)
{
    static const char NAMES[6][9]  = { "", "REALSHAP", "IMAGSHAP", "REALSOSH", "IMAGSOSH", "SHAPE   " };
    static const char NAMES2[6][9] = { "", "REALADD ", "IMAGADD ", "REALSOAD", "IMAGSOAD", "SHAPEADD" };
    static const char SCALNM[6][9] = { "", "REALSCAL", "IMAGSCAL", "REALSOSC", "IMAGSOSC", "SHAPESCA" };

    double UNDEF = INTRNL.UNDEF;
    NUMOUT = 0;
    IRETUR = 0;

    double CONS  = TEMPVS_arr()[IPOTYP - 1];   // Fortran TEMPVS(IPOTYP)
    double CONS2 = PARAMS_at(IPOTYP);          // Fortran PARAM(IPOTYP)
    if (CONS == UNDEF || CONS2 == UNDEF) {
        std::printf("0**** THE WELL DEPTH PARAMETER (V, VI, VSO, VSOI) AND CORRESPONDINGLY"
                    " PARAM1 - PARAM4 MUST BE DEFINED FOR THE TWOSHAPE LINKULE.\n");
        IRETUR = -1;
        goto L900;
    }
    CONS  = -CONS;
    CONS2 = -CONS2;

    if (IREQUE == 1) goto L100;
    if (IREQUE == 2) goto L500;
    if (IREQUE == 3) goto L600;
    return;

L100:
    {
        int II = IPOTYP;
        int I  = NAMLOC(NAMES[II]);
        if (I == 0) {
            II = 5;
            I = NAMLOC(NAMES[II]);
            if (I == 0) {
                std::printf("0**** AN OBJECT WITH THE NAME %.8s OR TWOSHAPE MUST BE DEFINED"
                            " TO USE THE SHAPE LINKULE.\n", NAMES[IPOTYP]);
                IRETUR = -1;
                goto L900;
            }
        }
        int I2 = NAMLOC(NAMES2[II]);
        if (I2 == 0) {
            std::printf("0**** AN OBJECT WITH THE NAME %.8s MUST BE DEFINED"
                        " FOR THE TWOSHAPE LINKULE.\n", NAMES2[II]);
            IRETUR = -1;
            goto L900;
        }
        int NIN = LENGTH.LENG[I];
        if (NIN != LENGTH.LENG[I2]) {
            std::printf("0**** %.8s AND %.8s MUST BE DEFINED ON THE SAME GRID.\n",
                        NAMES[II], NAMES2[II]);
            IRETUR = -1;
            goto L900;
        }

        int ISCAL = NAMLOC(SCALNM[II]);
        double RMAX = RSTART + (NUMPTS - 1) * STEPSZ;
        int LISCAL = (ISCAL == 0) ? 0 : LENGTH.LENG[ISCAL];

        bool need_interp = false;
        bool skip_to_400 = false;
        bool RSW = false;
        double R1 = 0, R2 = 0;

        if (ISCAL == 0) {
            // No scale array — just check point count
            if (NUMPTS != NIN) {
                std::printf("0**** ASYMPTOPIA AND STEPSIZE =%15.5g%15.5g AND REQUIRE%6d POINTS.\n"
                            "      HOWEVER, %.8s HAS%6d POINTS AND %.8s IS NOT DEFINED.\n",
                            RMAX, STEPSZ, NUMPTS, NAMES[II], NIN, SCALNM[II]);
                IRETUR = -1;
                goto L900;
            }
            std::printf("0%.8s AND %.8s ARE BEING ASSUMED TO BE DEFINED FOR%15.5g =< R =<%15.5g WITH STEPSIZE =%15.5g\n",
                        NAMES[II], NAMES2[II], RSTART, RMAX, STEPSZ);
            NUMOUT = 1;
            skip_to_400 = true;
        } else if (LISCAL == 2) {
            // Scaling array = (start, stop)
            R1 = ALLOC(LOCPTRS.Z[ISCAL]);
            R2 = ALLOC(LOCPTRS.Z[ISCAL] + 1);
            RSW = false;
            if (R1 == RSTART && NIN == NUMPTS && DABS(R2 - RMAX) / STEPSZ < 1.0e-5)
                skip_to_400 = true;
            else
                need_interp = true;
        } else if (LISCAL == NIN) {
            R1 = ALLOC(LOCPTRS.Z[ISCAL]);
            R2 = ALLOC(LOCPTRS.Z[ISCAL] + NIN - 1);
            RSW = true;
            need_interp = true;
        } else {
            std::printf("0**** SCALLING ARRAY %.8s SHOULD HAVE 2 OR%4d ELEMENTS, BUT IT HAS%6d ELEMENTS.\n",
                        SCALNM[II], NIN, LISCAL);
            IRETUR = -1;
            goto L900;
        }

        if (need_interp && !skip_to_400) {
            // Must rescale using AITKEN
            int IOLD  = I;
            int IOLD2 = I2;
            int ISIZE = NIN + NUMPTS;
            char wname[9]  = "    SHAP";
            char wname2[9] = "    SHP2";
            char wname3[9] = "    WORK";
            std::memcpy(wname,  ID, 4);
            std::memcpy(wname2, ID, 4);
            std::memcpy(wname3, ID, 4);
            I  = NALLOC(wname,  NUMPTS);
            I2 = NALLOC(wname2, NUMPTS);
            int IWORK = NALLOC(wname3, ISIZE);

            int LXIN   = LOCPTRS.Z[IWORK] - 1;
            int LXOUT  = LXIN + NIN;
            int LYIN   = LOCPTRS.Z[IOLD]  - 1;
            int LYOUT  = LOCPTRS.Z[I]     - 1;
            int LYIN2  = LOCPTRS.Z[IOLD2] - 1;
            int LYOUT2 = LOCPTRS.Z[I2]    - 1;

            if (RSW) {
                LXIN = LOCPTRS.Z[ISCAL] - 1;
            } else {
                // Setup input X's evenly between R1 and R2
                double R = R1;
                double STEPIN = (R2 - R1) / (NIN - 1);
                for (int N = 1; N <= NIN; N++) {
                    ALLOC(LXIN + N) = R;
                    R += STEPIN;
                }
            }

            // Setup output X's, no extrapolation (use F(1) for R<R1, 0 for R>R2)
            double R = RSTART;
            int NSTRT = 1;
            int NLAST = 1;
            double VAL0  = ALLOC(LYIN  + 1);
            double VAL02 = ALLOC(LYIN2 + 1);
            for (int N = 1; N <= NUMPTS; N++) {
                if (R < R1) {
                    ALLOC(LYOUT  + N) = VAL0;
                    ALLOC(LYOUT2 + N) = VAL02;
                    NSTRT = N + 1;
                } else if (R > R2) {
                    ALLOC(LYOUT  + N) = 0.0;
                    ALLOC(LYOUT2 + N) = 0.0;
                } else {
                    ALLOC(LXOUT + N) = R;
                    NLAST = N;
                }
                R += STEPSZ;
            }
            int NOUT = NLAST + 1 - NSTRT;

            double WORK[20];
            int NFAIL; double WORST;
            AITKEN(5, 0.0, 0.0, NIN,
                   ALLOC_base(LXIN + 1), ALLOC_base(LYIN + 1),
                   NOUT, ALLOC_base(LXOUT + NSTRT), ALLOC_base(LYOUT + NSTRT),
                   NFAIL, WORST, WORK);
            AITKEN(5, 0.0, 0.0, NIN,
                   ALLOC_base(LXIN + 1), ALLOC_base(LYIN2 + 1),
                   NOUT, ALLOC_base(LXOUT + NSTRT), ALLOC_base(LYOUT2 + NSTRT),
                   NFAIL, WORST, WORK);

            // Free work area
            LOCPTRS.Z[IWORK] = -LOCPTRS.Z[IWORK];
        }

        MYINTS[0] = 256 * I + I2;
        MYINTS[1] = II;
        if (NUMOUT != 0) goto L900;
        return;
    }

L500:
    std::printf(" THE POTENTIAL IS%15.5g TIMES THE SHAPE IN %.8s\n"
                "               PLUS%15.5g TIMES THE SHAPE IN %.8s\n",
                CONS, NAMES[MYINTS[1]], CONS2, NAMES2[MYINTS[1]]);
    goto L900;

L600:
    {
        int I  = MYINTS[0] / 256;
        int I2 = MYINTS[0] - 256 * I;
        int LL  = LOCPTRS.Z[I]  - 1;
        int LL2 = LOCPTRS.Z[I2] - 1;
        for (int i = 1; i <= NUMPTS; i++) {
            ARRAY1[i] = CONS * ALLOC(LL + i) + CONS2 * ALLOC(LL2 + i);
        }
    }
    return;

L900:
    NUMOUT = 1;
    return;
}
