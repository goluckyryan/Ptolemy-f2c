// ampcal_translated.cpp — Verbatim C++ translation of AMPCAL
// Translated from source.f lines 222-580
//
// CALCULATES THE TRANSITION AMPLITUDES F AT ONE ANGLE.

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_forward.h"
#include <cstdio>
#include <cmath>

// Helper: Fortran 1-based access for raw pointer arrays.
// F1(ptr, i) accesses Fortran index i (1-based) from a 0-based C pointer.
#define F1(ptr, i) (ptr)[(i) - 1]

// Helper: Fortran column-major 2D access (1-based) for raw pointer.
// F2(ptr, dim1, r, c) accesses element (r,c) in array with first dim = dim1.
#define F2(ptr, dim1, r, c) (ptr)[((c) - 1) * (dim1) + (r) - 1]

// Helper: Fortran column-major 3D access (1-based) for raw pointer.
// F3(ptr, d1, d2, r, c, k) accesses element (r,c,k) with dims (d1,d2,*).
#define F3(ptr, d1, d2, r, c, k) (ptr)[((k) - 1) * (d1) * (d2) + ((c) - 1) * (d1) + (r) - 1]


// ============================================================================
// SUBROUTINE AMPCAL — lines 222-580
//
// CALCULATES THE TRANSITION AMPLITUDES F AT ONE ANGLE.
//
// INPUT SCALARS:
//   ANGLE = C.M. ANGLE (DEGREES).
//   ELSW = .TRUE. FOR ELASTIC SCATTERING.
//   NSPL = NUMBER OF AMPLITUDES F.
//   LMN, LMX, LSKP = MINIMUM, MAXIMUM, INCREMENT IN INCOMING L.
//   LXMX = MAXIMUM LX.
//   IDNPAR = 0 IF THERE ARE NO IDENTICAL PARTICLES
//          = 1 IF ONLY ONE PAIR (IN OR OUT) IS IDENTICAL
//          = 2 IF BOTH PAIRS ARE IDENTICAL
//   ETA = COULOMB PARAMETER (ELASTIC ONLY).
//   UK = INCIDENT WAVENUMBER (ELASTIC ONLY).
//   SIGZRO = L=0 COULOMB PHASE (ELASTIC ONLY).
//   LOHISW = .TRUE. IF FLO AND FHI ARE TO BE RETURNED.
//   LHI = MAXIMUM L TO BE INCLUDED IN FHI.
//   LEBACK = NUMBER OF L VALUES TO USE IN THE EPSILON ALGORITHM; 0 = NOT USED.
//
// OUTPUT SCALAR:
//   NFLOP = APPROX NUMBER OF FLOATING OPS IN THIS CALL
//
// INPUT ARRAYS:
//   JTOCS(4,NSPL)
//   BETAS(2,NSPL,LMX-LBASE+1) — REAL*4
//   ALOWFC(2,2)
//
// OUTPUT ARRAYS:
//   F(2,NSPL), FLO(2,NSPL), FHI(2,NSPL), FERROR(NSPL), FCOUL(2,3)
//   PLM(LMX+1+((2*LMX+1-LXMX)*LXMX)/2)
//   CONTR(LMX+1), CONTI(LMX+1), FEPSLO(2,LEBACK)
// ============================================================================
void AMPCAL(double ANGLE, int ELSW, int NSPL, int LMN, int LMX, int LSKP,
            int LXMX, int IDNPAR, double ETA, double UK, double SIGZRO,
            int LOHISW, int LHI,
            int* JTOCS, float* BETAS, double* ALOWFC,
            double* F, double* FLO, double* FHI, double* FERROR,
            double* FCOUL, double* PLM, double* CONTR, double* CONTI,
            int LEBACK, double* FEPSLO, int& NFLOP)
{
    // COMMON /CNSTNT/
    double& PI      = CNSTNT.PI;
    double& DEGREE  = CNSTNT.DEGREE;
    double& BIGNUM  = CNSTNT.BIGNUM;
    double& SMLNUM  = CNSTNT.SMLNUM;

    // Local variables
    int IDENSW;   // LOGICAL
    int IDELSW;   // LOGICAL
    double FT[3];    // 1-based: FT[1], FT[2]
    double FTS[3];   // 1-based: FTS[1], FTS[2]
    double FTC[3];   // 1-based: FTC[1], FTC[2]
    double FTSC[3];  // 1-based: FTSC[1], FTSC[2]
    double FOUT[3];  // 1-based: FOUT[1], FOUT[2]

    int NNFLOP, IODD, LOMN, LOMX, LOHI, LBASE;
    int N, K, KOFFS, LX, MX, LOMNMX, LPLM, LTEMP, LO, LPARIT;
    int IRI, II, NN, IER;
    double AN, ACOS_val, TEMP, PHASE, PLMFAC, FACMBL;
    double CONR, CONI, DUMMY;

    //
    // INITIALIZATION
    //
    NFLOP = 0;
    NNFLOP = 0;
    IODD = MOD(IABS(F2(JTOCS, 4, 1, 1)), 2);
    LOMN = LMN;
    if (LSKP == 2) LOMN = LOMN + IODD;
    LOMX = LMX;
    if (LSKP == 2) LOMX = LOMX - MOD(LOMX + IODD, 2);
    LOHI = LHI;
    if (LSKP == 2) LOHI = LOHI - MOD(LOHI + IODD, 2);
    LBASE = LMN;
    IDENSW = (IDNPAR != 0) ? TRUE_F : FALSE_F;
    FACMBL = 1;
    if (IDENSW && !ELSW) FACMBL = DSQRT(2.0);
    if (IDNPAR == 2 && !ELSW) FACMBL = 2;
    IDELSW = (IDENSW && ELSW) ? TRUE_F : FALSE_F;
    N = 2 * NSPL;
    if (IDELSW) N = 6;
    for (K = 1; K <= N; K++) {
        F1(F, K) = 0;
        if (!LOHISW) goto L39;
        F1(FLO, K) = 0;
        F1(FHI, K) = 0;
    L39:;
    }

    //
    // CALCULATE THE LEGENDRE FUNCTIONS.
    //
    AN = DEGREE * ANGLE;
    ACOS_val = DCOS(AN);
    PLMSUB(LMX, LXMX, ACOS_val, PLM-1);
    NFLOP = 16 + 7 * (LMX + 1) * (LXMX + 1);

    //
    // ELASTIC: CALCULATE THE COULOMB AMPLITUDE.
    // IDENTICAL PARTICLES NEED IT FOR BOTH ANGLE AND PI-ANGLE.
    //
    if (!ELSW) goto L100;
    LOMN = 0;
    LBASE = 0;
    TEMP = DSIN(0.5 * AN);
    TEMP = DMAX1(TEMP, 1.0e-10);
    PHASE = 2 * (SIGZRO - ETA * DLOG(TEMP));
    TEMP = -0.5 * ETA / (UK * TEMP * TEMP);
    if (DABS(TEMP) > BIGNUM) TEMP = DSIGN(BIGNUM, TEMP);
    FTC[1] = TEMP * DCOS(PHASE);
    FTC[2] = TEMP * DSIN(PHASE);
    FTSC[1] = 0;
    FTSC[2] = 0;
    NFLOP = NFLOP + 72;
    if (!IDENSW) goto L80;
    AN = PI - AN;
    TEMP = DSIN(0.5 * AN);
    TEMP = DMAX1(TEMP, 1.0e-10);
    PHASE = 2 * (SIGZRO - ETA * DLOG(TEMP));
    TEMP = -0.5 * ETA / (UK * TEMP * TEMP);
    if (DABS(TEMP) > BIGNUM) TEMP = DSIGN(BIGNUM, TEMP);
    FTSC[1] = TEMP * DCOS(PHASE);
    FTSC[2] = TEMP * DSIN(PHASE);
    NFLOP = NFLOP + 72;
    FTS[1] = 0;
    FTS[2] = 0;

    //
    // STILL ELASTIC: STORE SYMMETRIC, ANTISYMMETRIC, AND
    // UNSYMMETRIZED COULOMB AMPLITUDES.
    //
L80:
    for (IRI = 1; IRI <= 2; IRI++) {
        F2(FCOUL, 2, IRI, 1) = FTC[IRI] + FTSC[IRI];
        F2(FCOUL, 2, IRI, 2) = FTC[IRI] - FTSC[IRI];
        F2(FCOUL, 2, IRI, 3) = FTC[IRI];
    }
    NFLOP = NFLOP + 4;

    //
    // LOOP THROUGH KOFFS = JT, JP, LX, MX.
    //
L100:
    for (KOFFS = 1; KOFFS <= NSPL; KOFFS++) {
        if (F2(JTOCS, 4, 4, KOFFS) < 0) goto L299;
        LX = F2(JTOCS, 4, 2, KOFFS);
        MX = (F2(JTOCS, 4, 1, KOFFS) + LX + 1) / 2;
        LOMNMX = MAX0(LOMN, MX);
        if (LSKP == 2) LOMNMX = LOMNMX + MOD(LOMNMX + IODD, 2);
        FT[1] = 0;
        FT[2] = 0;
        LPLM = MX * (2 * LMX + 1 - MX) / 2 + 1;

        //
        // SUM OVER LO.
        //
        // FLO = AMPLITUDE FROM LO = LOMN, LOMN+1.
        // FOR ELASTIC AND LX=0, IT DOES NOT INCLUDE THE 1 IN 1-S(L)
        //
        if (LX != 0 || !ELSW) goto L110;
        if (!LOHISW) goto L110;
        LTEMP = LOMN + 1;
        for (LO = LOMN; LO <= LTEMP; LO += LSKP) {
            PLMFAC = FACMBL * F1(PLM, LO + 1);
            LPARIT = MOD(LO, 2) + 1;
            F2(FLO, 2, 1, 1) = F2(FLO, 2, 1, 1) + PLMFAC * F2(ALOWFC, 2, 1, LPARIT);
            F2(FLO, 2, 2, 1) = F2(FLO, 2, 2, 1) + PLMFAC * F2(ALOWFC, 2, 2, LPARIT);
        }
        NNFLOP = 5 * (2 + LSKP);

        //
        // COMPUTE ALL CONTRIBUTIONS FROM NON-ZERO S-MATRICES
        //
    L110:
        for (LO = LOMNMX; LO <= LOMX; LO += LSKP) {
            PLMFAC = FACMBL * F1(PLM, LO + LPLM);
            F1(CONTR, LO + 1) = PLMFAC * F3(BETAS, 2, NSPL, 1, KOFFS, LO - LBASE + 1);
            F1(CONTI, LO + 1) = PLMFAC * F3(BETAS, 2, NSPL, 2, KOFFS, LO - LBASE + 1);
        }
        NNFLOP = NNFLOP + 5 * (LOMX - LOMNMX + LSKP);

        //
        // FIRST SUM AMPLITUDE UP TO LMIN+1 TO GET FLO
        //
        LTEMP = MIN0(LOMN + 1, LOMX);
        if (LTEMP < LOMNMX) goto L125;
        for (LO = LOMNMX; LO <= LTEMP; LO += LSKP) {
            FT[1] = FT[1] + F1(CONTR, LO + 1);
            FT[2] = FT[2] + F1(CONTI, LO + 1);
        }

        //
    L125:
        if (!LOHISW) goto L130;
        if (ELSW && LX == 0) goto L130;
        F2(FLO, 2, 1, KOFFS) = FT[1];
        F2(FLO, 2, 2, KOFFS) = FT[2];

        //
        // SUM CONTRIBUTION FOR LMIN+2 TO LMAX (LMAX-2 FOR ELASTIC)
        //
    L130:
        LTEMP = MAX0(LOMN + 2, LOMNMX);
        if (LTEMP > LOHI) goto L140;
        for (LO = LTEMP; LO <= LOHI; LO += LSKP) {
            FT[1] = FT[1] + F1(CONTR, LO + 1);
            FT[2] = FT[2] + F1(CONTI, LO + 1);
        }

        //
        // FHI = AMPLITUDE FROM LO <= LOHI.
        //
    L140:
        F2(FHI, 2, 1, KOFFS) = FT[1];
        F2(FHI, 2, 2, KOFFS) = FT[2];

        //
        // NOW THE LAST TERMS > LMAX OR LMAX-1
        //
        LTEMP = LOHI + LSKP;
        if (LTEMP > LOMX) goto L150;
        for (LO = LTEMP; LO <= LOMX; LO += LSKP) {
            FT[1] = FT[1] + F1(CONTR, LO + 1);
            FT[2] = FT[2] + F1(CONTI, LO + 1);
        }

        //
        // IF NECESSARY APPLY EPSILON ALGORITHM
        //
    L150:
        F1(FERROR, KOFFS) = (DABS(F1(CONTR, LOMX + 1)) + DABS(F1(CONTI, LOMX + 1)))
            / (DABS(FT[1]) + DABS(FT[2]) + SMLNUM);
        if (LEBACK <= 0) goto L200;
        N = MIN0((LOMX - LOMNMX) / LSKP, LEBACK - 1);
        if (N <= 5) goto L200;
        NN = N + 1;
        F2(FEPSLO, 2, 1, NN) = FT[1];
        F2(FEPSLO, 2, 2, NN) = FT[2];
        for (int I = 1; I <= N; I++) {
            II = NN - I;
            F2(FEPSLO, 2, 1, II) = F2(FEPSLO, 2, 1, II + 1) - F1(CONTR, LOMX + 1 - LSKP * (I - 1));
            F2(FEPSLO, 2, 2, II) = F2(FEPSLO, 2, 2, II + 1) - F1(CONTI, LOMX + 1 - LSKP * (I - 1));
        }
        {
            double DEPS_val = 1.0e-5;
            EPSLON(FEPSLO, NN, &FT[1], DEPS_val, F1(FERROR, KOFFS), IER);
        }

        N = MIN0((LOHI - LOMNMX) / LSKP, LEBACK - 1);
        if (N <= 5) goto L200;
        NN = N + 1;
        F2(FEPSLO, 2, 1, NN) = F2(FHI, 2, 1, KOFFS);
        F2(FEPSLO, 2, 2, NN) = F2(FHI, 2, 2, KOFFS);
        for (int I = 1; I <= N; I++) {
            II = NN - I;
            F2(FEPSLO, 2, 1, II) = F2(FEPSLO, 2, 1, II + 1) - F1(CONTR, LOHI + 1 - LSKP * (I - 1));
            F2(FEPSLO, 2, 2, II) = F2(FEPSLO, 2, 2, II + 1) - F1(CONTI, LOHI + 1 - LSKP * (I - 1));
        }
        {
            double DEPS_val = 1.0e-5;
            EPSLON(FEPSLO, NN, &FOUT[1], DEPS_val, DUMMY, IER);
            F2(FHI, 2, 1, KOFFS) = FOUT[1];
            F2(FHI, 2, 2, KOFFS) = FOUT[2];
        }

        //
        // FOR ELASTIC SCATTERING OF IDENTICAL PARTICLES,
        // DO IT AGAIN FOR PI - ANGLE.
        //
    L200:
        if (!IDELSW) goto L270;
        for (LO = LOMNMX; LO <= LOMX; LO += LSKP) {
            CONR = F1(CONTR, LO + 1);
            CONI = F1(CONTI, LO + 1);
            if (MOD(LO, 2) == 0) goto L220;
            CONR = -CONR;
            CONI = -CONI;
        L220:
            FTS[1] = FTS[1] + CONR;
            FTS[2] = FTS[2] + CONI;
            if (!LOHISW) goto L259;
            if (LO > LOMN + 1) goto L250;
            LPARIT = MOD(LO, 2) + 1;
            PLMFAC = FACMBL * F1(PLM, LO + LPLM);
            if (LPARIT == 2) PLMFAC = -PLMFAC;
            F2(FLO, 2, 1, 2) = F2(FLO, 2, 1, 2) + PLMFAC * F2(ALOWFC, 2, 1, LPARIT);
            F2(FLO, 2, 2, 2) = F2(FLO, 2, 2, 2) + PLMFAC * F2(ALOWFC, 2, 2, LPARIT);
        L250:
            if (LO != LOHI) goto L259;
            F2(FHI, 2, 1, 2) = FTS[1];
            F2(FHI, 2, 2, 2) = FTS[2];

        L259:;
        }
        NNFLOP = NNFLOP + 5 * (LOMX - LOMNMX + LSKP);

        //
        // END OF LO SUM. TRANSFER TO OUTPUT ARRAYS.
        // FLO IS ALL EXCEPT THE FIRST 2 L'S.
        //
    L270:
        F2(F, 2, 1, KOFFS) = FT[1];
        F2(F, 2, 2, KOFFS) = FT[2];

        if (!LOHISW) goto L299;
        F2(FLO, 2, 1, KOFFS) = FT[1] - F2(FLO, 2, 1, KOFFS);
        F2(FLO, 2, 2, KOFFS) = FT[2] - F2(FLO, 2, 2, KOFFS);
    L299:;
    }

    //
    // END OF KOFFS LOOP.
    //
    NFLOP = NFLOP + NNFLOP / LSKP;

    //
    // FOR ELASTIC SCATTERING, ADD THE COULOMB AMPLITUDES.
    //
    if (!ELSW) return;
    for (IRI = 1; IRI <= 2; IRI++) {
        F2(F, 2, IRI, 1) = F2(F, 2, IRI, 1) + FTC[IRI];
        if (!LOHISW) goto L330;
        F2(FHI, 2, IRI, 1) = F2(FHI, 2, IRI, 1) + FTC[IRI];
        F2(FLO, 2, IRI, 1) = F2(FLO, 2, IRI, 1) + FTC[IRI];

        //
        // FOR ELASTIC SCATTERING OF IDENTICAL PARTICLES,
        // FORM SYMMETRIC, ANTISYMMETRIC, AND UNSYMMETRIZED COMBINATIONS.
        //
    L330:
        if (!IDELSW) goto L359;

        FTS[IRI] = FTS[IRI] + FTSC[IRI];
        F2(F, 2, IRI, 3) = F2(F, 2, IRI, 1);
        F2(F, 2, IRI, 1) = F2(F, 2, IRI, 3) + FTS[IRI];
        F2(F, 2, IRI, 2) = F2(F, 2, IRI, 3) - FTS[IRI];
        if (!LOHISW) goto L359;
        F2(FLO, 2, IRI, 2) = FTS[IRI] - F2(FLO, 2, IRI, 2);
        F2(FLO, 2, IRI, 3) = F2(FLO, 2, IRI, 1);
        F2(FLO, 2, IRI, 1) = F2(FLO, 2, IRI, 3) + F2(FLO, 2, IRI, 2);
        F2(FLO, 2, IRI, 2) = F2(FLO, 2, IRI, 3) - F2(FLO, 2, IRI, 2);
        F2(FHI, 2, IRI, 2) = F2(FHI, 2, IRI, 2) + FTSC[IRI];
        F2(FHI, 2, IRI, 3) = F2(FHI, 2, IRI, 1);
        F2(FHI, 2, IRI, 1) = F2(FHI, 2, IRI, 3) + F2(FHI, 2, IRI, 2);
        F2(FHI, 2, IRI, 2) = F2(FHI, 2, IRI, 3) - F2(FHI, 2, IRI, 2);
    L359:;
    }
    NFLOP = NFLOP + 6;
    return;
}
