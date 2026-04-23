// angular_part2_translated.cpp — Verbatim C++ translation of BETCAL, CALANG, ELDCS
// Translated from source.f lines 3796-4080, 5590-5804, 14845-15232

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include "ptolemy_alloc.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_forward.h"
#include "ptolemy_io.h"
#include <cstdio>
#include <cstring>
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
// SUBROUTINE BETCAL — lines 3796-4080
//
// CALCULATES THE ANGLE-INDEPENDENT PARTS (BETA) OF THE AMPLITUDES F.
// ============================================================================
void BETCAL(int ELSW, double UK, int JSP, int NSPL, int LMN, int LMX, int LSKP,
            int LXMN, int LXMX, int LDELMX, int NTEMPS, int ISTAT, int PRNTSW,
            double& SIGTOT, double& SIGREA,
            int* JTOCS, double* S, float* SMAG, float* SPHASE, double* SIGIN, double* SIGOT,
            float* BETAS, double* TOTLX, double* TOTMX, double* ALOWFC, double* TEMPS)
{
    // JTOCS(4,NSPL), S(2,NSPL,LMX+1), SMAG(NSPL,*), SPHASE(NSPL,*)
    // SIGIN(LMX+1), SIGOT(LMX+1)
    // BETAS(2,NSPL,LMX-LBASE+1) — REAL*4
    // TOTLX(LXMX+1), TOTMX(NSPL), ALOWFC(2,2), TEMPS(NTEMPS)

    // COMMON /CNSTNT/
    double& PI = CNSTNT.PI;

    double FOPT[3];   // 1-based: FOPT[1], FOPT[2]
    double REACTN_arr[3]; // 1-based: REACTN_arr[1], REACTN_arr[2]
    int IZERO = 0;

    //
    // INITIALIZATION
    //
    int IODD = MOD(LDELMX, 2);
    double FACTOR = 0.5 / UK;
    if (IODD != 0) FACTOR = -FACTOR;
    int NMX = LXMX + 1;
    int NMLX = LXMX - LXMN + 1;
    int LBASE, LOMN;
    if (ELSW) goto L30;
    // REACTION ONLY
    for (int KOFFS = 1; KOFFS <= NSPL; KOFFS++) {
        F1(TOTMX, KOFFS) = 0;
    }
    for (int LX = IZERO; LX <= LXMX; LX++) {
        F1(TOTLX, LX + 1) = 0;
    }
    SIGTOT = 0;
    LBASE = LMN;
    LOMN = LMN;
    if (LSKP == 2) LOMN = LOMN + IODD;
    goto L40;
    // ELASTIC SCATTERING ONLY
L30:
    LBASE = 0;
    LOMN = 0;
    FOPT[1] = 0;
    FOPT[2] = 0;
    REACTN_arr[1] = 0;
    REACTN_arr[2] = 0;
    for (int I = 1; I <= 4; I++) {
        // ALOWFC(I,1) = 0  — linear fill of (2,2) array, elements 1..4
        F1(ALOWFC, I) = 0;
    }
    // BOTH ELASTIC AND REACTION
L40:
    {
        int N = 2 * NSPL * (LMX - LBASE + 1);
        for (int I = 1; I <= N; I++) {
            // BETAS(I,1,1) = 0 — linear fill
            F1(BETAS, I) = 0;
        }
    }

    //
    // THE FOLLOWING FORCES THE COMMON BLOCKS USED BY CLEBSH TO BE
    // LOADED INTO CORE.
    //
    DUMMY1();

    //
    // THE OUTER LOOP IS OVER LO.
    //
    for (int LO = LOMN; LO <= LMX; LO += LSKP) {
        double DLO = LO;
        double D2L1 = 2 * DLO + 1.;
        int ILO = LO - LBASE + 1;
        int LPARIT = MOD(LO, 2) + 1;
        int LIMN = LO - LDELMX;
        int LIMX = LO + LDELMX;
        int MXMX = MIN0(LXMX, LO);

        //
        // INSERT CLEBSCH-GORDANS ETC. INTO THE TEMP ARRAY.
        //
        for (int I = 1; I <= NTEMPS; I++) {
            F1(TEMPS, I) = 0;
        }
        for (int LI = LIMN; LI <= LIMX; LI += 2) {
            if (LI < LBASE) goto L189;
            {
                int LX1 = MAX0(IABS(LI - LO), LXMN);
                for (int LX = LX1; LX <= LXMX; LX++) {
                    int MXZ = MOD(LX + LI - LO, 2);
                    int I = 1 + NMX * (LX - LXMN + NMLX * (LI - LIMN) / 2);
                    for (int MX = MXZ; MX <= LX; MX++) {
                        F1(TEMPS, I + MX) = FACTOR * (2 * LI + 1)
                            * CLEBSH(2 * LI, 2 * LX, 0, 2 * MX, 2 * LO, 2 * MX);
                    }
                }
            }
L189:       ;
        }

        //
        // START THE MAIN LOOP OVER JT, JP, LX, LI-LO.
        //
        int LIPREV = -100000;
        int MXZ = 0, KOFFZ = 0;
        for (int KOFFS = 1; KOFFS <= NSPL; KOFFS++) {
            if (F2(JTOCS, 4, 4, KOFFS) < 0) goto L299;
            { int LX = F2(JTOCS, 4, 2, KOFFS);
            int LI = LO - F2(JTOCS, 4, 1, KOFFS);

            //
            // A INCREASE OR NO CHANGE IN LI MEANS A NEW (JT,JP,LX) SET, AND
            // HENCE A NEW SET OF MX VALUES TO INCLUDE IN THE SUM.
            //
            if (LI < LIPREV) goto L210;
            MXZ = LX + F2(JTOCS, 4, 1, KOFFS);
            KOFFZ = KOFFS - MXZ;
L210:       LIPREV = LI;
            if (LI < LBASE || LI > LMX) goto L299;

            //
            // INSERT THE COULOMB PHASES INTO THE S-MATRIX ELEMENT, AND
            // DIVIDE BY I.
            //
            // FOR REACTIONS, S COMES FROM SMAG AND SPHASE.
            //
            {
                int I = LI - LBASE + 1;
                double SMATR, SMATI, PHASE, AMAG, TR, TI, SR, SI;
                if (ELSW) goto L230;
                AMAG = F2(SMAG, NSPL, KOFFS, I);
                if (AMAG == 0) goto L299;
                PHASE = F2(SPHASE, NSPL, KOFFS, I) + F1(SIGIN, LI + 1) + F1(SIGOT, LO + 1);
                SMATR = AMAG * DSIN(PHASE);
                SMATI = -AMAG * DCOS(PHASE);
                goto L250;

                //
                // FOR ELASTIC SCATTERING, S COMES FROM S.
                // ALSO UPDATE REACTION CROSS SECTION.
                //
L230:           PHASE = F1(SIGIN, LI + 1) + F1(SIGIN, LO + 1);
                TR = DCOS(PHASE);
                TI = DSIN(PHASE);
                if (LO >= LMN) goto L235;
                SMATR = 0;
                SMATI = 0;
                if (LX < 0) goto L250;
                if (LX == 0) goto L240;
                goto L250;
L235:           SR = F3(S, 2, NSPL, 1, KOFFS, I);
                SI = F3(S, 2, NSPL, 2, KOFFS, I);
                SMATR = SR * TI + SI * TR;
                SMATI = SI * TI - SR * TR;
                REACTN_arr[LPARIT] = REACTN_arr[LPARIT] - D2L1 * (SR * SR + SI * SI);
                if (LX != 0) goto L250;
                if (LO > LMN + 1) goto L240;
                // ALOWFC(1,LPARIT) and ALOWFC(2,LPARIT)
                F2(ALOWFC, 2, 1, LPARIT) = FACTOR * D2L1 * SMATR;
                F2(ALOWFC, 2, 2, LPARIT) = FACTOR * D2L1 * SMATI;
L240:           SMATR = SMATR - TI;
                SMATI = SMATI + TR;
                REACTN_arr[LPARIT] = REACTN_arr[LPARIT] + D2L1;

                //
                // MULTIPLY S BY THE CLEBSCH-GORDAN ETC, AND SUM INTO BETAS
                // FOR EACH RELEVANT MX.
                //
L250:           I = 1 + NMX * (LX - LXMN + NMLX * (LI - LIMN) / 2);
                for (int MX = MXZ; MX <= LX; MX++) {
                    // BETAS(1,KOFFZ+MX,ILO) and BETAS(2,KOFFZ+MX,ILO)
                    F3(BETAS, 2, NSPL, 1, KOFFZ + MX, ILO) =
                        F3(BETAS, 2, NSPL, 1, KOFFZ + MX, ILO)
                        + (float)(F1(TEMPS, I + MX) * SMATR);
                    F3(BETAS, 2, NSPL, 2, KOFFZ + MX, ILO) =
                        F3(BETAS, 2, NSPL, 2, KOFFZ + MX, ILO)
                        + (float)(F1(TEMPS, I + MX) * SMATI);
                }
            }
            } // end scope for LX, LI

L299:       ;
        }

        //
        // NOW CALCULATE THE SQRT(FACTORIALS) INTO TEMPS.
        //
        F1(TEMPS, 1) = 1;
        if (MXMX == 0) goto L320;
        for (int MX = 1; MX <= MXMX; MX++) {
            F1(TEMPS, MX + 1) = F1(TEMPS, MX) / DSQRT((DLO + MX) * (DLO - MX + 1));
        }

        //
        // CROSS SECTIONS MUST BE DOUBLED FOR IDENTICAL PARTICLES.
        //
L320:   {
            double TEMP = 40 * PI / (2 * DLO + 1);
            if (ISTAT != 3) TEMP = 2 * TEMP;

            //
            // SECOND PASS THROUGH KOFFS LOOP:
            // INCREMENT THE TOTAL AND PARTIAL CROSS SECTIONS,
            // AND PUT THE SQRT(FACTORIALS) INTO BETA.
            //
            for (int KOFFS = 1; KOFFS <= NSPL; KOFFS++) {
                if (F2(JTOCS, 4, 4, KOFFS) < 0) goto L399;
                { int LX = F2(JTOCS, 4, 2, KOFFS);
                int MX = (F2(JTOCS, 4, 1, KOFFS) + LX + 1) / 2;
                float BETAR = F3(BETAS, 2, NSPL, 1, KOFFS, ILO);
                float BETAI = F3(BETAS, 2, NSPL, 2, KOFFS, ILO);
                if (ELSW) goto L350;
                // REACTION ONLY.
                {
                    double TERM = TEMP * ((double)BETAR * BETAR + (double)BETAI * BETAI);
                    F1(TOTMX, KOFFS) = F1(TOTMX, KOFFS) + TERM;
                    if (MX > 0) TERM = 2 * TERM;
                    SIGTOT = SIGTOT + TERM;
                    F1(TOTLX, LX + 1) = F1(TOTLX, LX + 1) + TERM;
                }
                goto L360;
                // ELASTICS ONLY
L350:           if (MX == 0) FOPT[LPARIT] = FOPT[LPARIT] + (double)BETAI;

L360:           BETAR = (float)(F1(TEMPS, MX + 1) * BETAR);
                BETAI = (float)(F1(TEMPS, MX + 1) * BETAI);
                F3(BETAS, 2, NSPL, 1, KOFFS, ILO) = BETAR;
                F3(BETAS, 2, NSPL, 2, KOFFS, ILO) = BETAI;

                if (PRNTSW) std::printf(" JP, JT, LX, MX, LO, B =%3d/2%3d/2%3d%3d%4d%17.5G +%13.5GI\n",
                    F2(JTOCS, 4, 3, KOFFS), F2(JTOCS, 4, 4, KOFFS),
                    LX, MX, LO, (double)BETAR, (double)BETAI);
                } // end scope for LX, MX, BETAR, BETAI

L399:           ;
            }
        }
    }

    //
    // FOR ELASTICS, FINISH CALCULATING THE REACTION AND
    // NUCLEAR CROSS SECTIONS.
    //
    if (!ELSW) return;
    SIGREA = REACTN_arr[1] + REACTN_arr[2];
    SIGTOT = FOPT[1] + FOPT[2];
    if (ISTAT == 3) goto L520;
    {
        double TEMP = JSP;
        double TERM = TEMP + 1;
        SIGREA = (TEMP * SIGREA + 2 * REACTN_arr[ISTAT]) / TERM;
        SIGTOT = (TEMP * SIGTOT + 2 * FOPT[ISTAT]) / TERM;
    }
L520:
    SIGREA = 10 * PI * SIGREA / (UK * UK);
    SIGTOT = 40 * PI * SIGTOT / UK;
    return;
}


// ============================================================================
// SUBROUTINE CALANG — lines 5590-5804
//
// EVALUATES CROSS SECTIONS AND GRADIENTS FOR ONE ANGLE
// ============================================================================
void CALANG(int GRADSW, int LLMIN, int LLMAX, int LSKIP, int MPARA,
            double ANGLE,
            double* PL, double* COULPH, double* DFDX, double* SINCOS,
            double* SMAT, double* DSDP,
            double& FR, double& FI, double& RUTH, double& HIGHR, double& HIGHI)
{
    // PL(LLMAX), COULPH(LLMAX), DFDX(2,MPARA),
    // SINCOS(2,LLMAX), SMAT(2,LLMAX), DSDP(2,MPARA,LLMAX)

    // COMMON blocks
    double& PI     = CNSTNT.PI;

    // /KANDM/
    // AKS(1)=AKI, AKS(2)=AKO; ETAS(1), ETAS(2)

    int& ITEMP_loc = LOCFIT.full.ITEMP;

    //
    // GET LEGENDRE POLYNOMIALS
    //
    PLSUB(LLMAX, DCOS(ANGLE), PL);

    //
    // RUTHERFORD CROSS SECTION FIRST
    //
    double RATJAC = KANDM.ETAS[1] / (2 * KANDM.AKI);
    double TEMP = ANGLE;
    FR = 0;
    FI = 0;

    //
    // FOLLOWING WORKS ONLY FOR:
    //  1) NON-IDENTICAL PARTICLES (LSKIP=1)
    //  2) IDENTICAL SPIN-0 PARTICLES (LSKIP=2)
    //
    for (int I = 1; I <= LSKIP; I++) {
        double S2 = DSIN(TEMP * .5);
        FR = FR - RATJAC / (S2 * S2) * DCOS(2 * (F1(COULPH, 1)
            - KANDM.ETAS[1] * DLOG(S2)));
        FI = FI - RATJAC / (S2 * S2) * DSIN(2 * (F1(COULPH, 1)
            - KANDM.ETAS[1] * DLOG(S2)));
        TEMP = PI - TEMP;
    }
    // CROSS SECTIONS IN MILLI BARNS
    RUTH = 10.0e0 * (FR * FR + FI * FI);
    if (RUTH == 0) RUTH = 1;

    //
    // NOW NUCLEAR AMPS
    // FOLLOWING FACTORS CONTAIN FACTOR OF 2 FOR IDENTICAL PARTICLES
    //
    double CINCR = LSKIP / KANDM.AKI;
    double C = .5 * CINCR;
    CINCR = LSKIP * CINCR;

    if (!GRADSW) goto L310;
    for (int I = 1; I <= MPARA; I++) {
        F2(DFDX, 2, 1, I) = 0;
        F2(DFDX, 2, 2, I) = 0;
    }

L310:
    if (LLMIN == 0) goto L330;

    //
    // SUM FOR L < LMIN; S = 0
    //
    {
        int LMINM1 = LLMIN - LSKIP;
        int IZERO = 0;
        for (int LL = IZERO; LL <= LMINM1; LL += LSKIP) {
            double AFAC = C * F1(PL, LL + 1);
            C = C + CINCR;
            FR = FR - AFAC * F2(SINCOS, 2, 1, LL + 1);
            FI = FI + AFAC * F2(SINCOS, 2, 2, LL + 1);
        }
    }

    //
    // NOW SUM FOR LLMIN =< L =< LLMAX
    //
L330:
    {
        // ASSIGN GOTO replacement: integer labels
        // IGOTO/JGOTO can be 355, 360, or 369
        int IGOTO = 360;
        int JGOTO = 360;
        if (!GRADSW) {
            IGOTO = 369;
            JGOTO = 369;
        }
        if (ITEMP_loc != 2) IGOTO = 355;

        for (int LL = LLMIN; LL <= LLMAX; LL += LSKIP) {
            double AFAC = C * F1(PL, LL + 1);
            double SINSIG = AFAC * F2(SINCOS, 2, 1, LL + 1);
            double COSSIG = AFAC * F2(SINCOS, 2, 2, LL + 1);
            C = C + CINCR;
            FR = FR + SINSIG * (F2(SMAT, 2, 1, LL + 1) - 1)
                + COSSIG * F2(SMAT, 2, 2, LL + 1);
            FI = FI + SINSIG * F2(SMAT, 2, 2, LL + 1)
                - COSSIG * (F2(SMAT, 2, 1, LL + 1) - 1);

            // GO TO IGOTO, ( 355, 360, 369 )
            switch (IGOTO) {
            case 355: goto L355;
            case 360: goto L360;
            case 369: goto L369;
            }
L355:       if (LL != (LLMAX - 2)) goto L359;
            HIGHR = FR;
            HIGHI = FI;
L359:       // GO TO JGOTO, ( 360, 369 )
            switch (JGOTO) {
            case 360: goto L360;
            case 369: goto L369;
            }

            //
            // COMPUTE THE GRADIENT OF F IF NEEDED
            //
L360:       for (int I = 1; I <= MPARA; I++) {
                F2(DFDX, 2, 1, I) = F2(DFDX, 2, 1, I)
                    + SINSIG * F3(DSDP, 2, MPARA, 1, I, LL + 1)
                    + COSSIG * F3(DSDP, 2, MPARA, 2, I, LL + 1);
                F2(DFDX, 2, 2, I) = F2(DFDX, 2, 2, I)
                    + SINSIG * F3(DSDP, 2, MPARA, 2, I, LL + 1)
                    - COSSIG * F3(DSDP, 2, MPARA, 1, I, LL + 1);
            }

L369:       ;
        }
    }

    return;
}


// ============================================================================
// SUBROUTINE ELDCS — lines 14845-15232
//
// CALCULATES DIFFERENTIAL CROSS SECTIONS FOR ELASTIC SCATTERING
// ============================================================================
void ELDCS(double UK, double ETA, double ANGMIN_arg, double ANGMAX_arg, double ANGSTP_arg,
           int LMN, int LMX, int LSKP, int ISTAT, int JSP,
           int ISMAT, int ITOC, int NSPL, int ISIG, double ELAB, double TAU,
           int PRNTSW, int& ITORUT, const char* RUTNAM,
           int FKEPSW, int& IFEL, const char* FNAME, double& SIGREA,
           int& IRUTH, const char* RUTHNM)
{
    // COMMON blocks
    auto& PI     = CNSTNT.PI;
    auto& RADIAN = CNSTNT.RADIAN;
    auto& DEGREE = CNSTNT.DEGREE;
    auto& IPRINT = INTGER.IPRINT;
    auto& LABANG = SWITCH.LABANG;
    auto& LEBACK = INTGER.LEBACK;
    auto* Z = LOCPTRS.Z;
    auto& FACFR4_v = ALLOCS.FACFR4;
    auto& ABAR  = KANDM.ABAR;

    int PBETSW;  // LOGICAL
    int IDNTSW;  // LOGICAL

    double FCOUL[7];     // flat 1-based: Fortran FCOUL(2,3), access via F2(FCOUL,2,r,c)
    double ALOWFC_arr[5]; // flat 1-based: Fortran ALOWFC(2,2), 4 elements
    float DUM4[2];       // DUM4(1,1) scratch
    double DUM8[2];      // DUM8(1) scratch

    double ANGMIN = ANGMIN_arg;
    double ANGMAX = ANGMAX_arg;
    double ANGSTP = ANGSTP_arg;

    PBETSW = (MOD(IPRINT, 10) >= 4) ? TRUE_F : FALSE_F;

    //
    // PRINT HEADER
    //
    if (PRNTSW) {
        std::printf("1");
        for (int i = 0; i < 47; i++) std::printf(" ");
        std::printf("P T O L E M Y\n");
        std::printf(" ");
        for (int i = 1; i <= 45; i++) std::printf("%c", HEDCOM.REACT[i]);
        std::printf("%8.2f MEV     ", ELAB);
        for (int i = 1; i <= 65; i++) std::printf("%c", HEDCOM.HEADER[i]);
        std::printf("\n");
        std::printf("0    ANGLE");
        std::printf("%*sSIGMA/", 11, "");
        std::printf("%*sSIGMA", 19, "");
        std::printf("%*sRUTHERFORD", 19, "");
        std::printf("%*s%% PER", 19, "");
        std::printf("%*s%% PER\n", 9, "");
        std::printf("  C.M.    LAB     RUTHERFORD");
        std::printf("%*sC.M.", 6, "");
        std::printf("%*sLAB", 13, "");
        std::printf("%*sC.M.", 15, "");
        std::printf("%*sLOW L", 16, "");
        std::printf("%*sHIGH L\n", 8, "");
        std::printf("\n");
    }

    double ANMIN = ANGMIN;
    double ANMAX = ANGMAX;
    if (PRNTSW) ABAR = 180;

    //
    // FOR LAB ANGLES WE HAVE FUNNY STUFF IF TAU > 1
    // FOR XSECTN CALLS, IT HAS ALL BEEN WORKED OUT ALREADY
    //
    if (LABANG == 0) goto L50;
    if (!PRNTSW) goto L50;
    if (TAU <= 1) goto L50;

    //
    // OUR CONVENTION IS THAT NEGATIVE VALUES FOR ANGLEMIN OR ANGLEMAX
    // CORRESPOND TO THE 2ND BRANCH.
    //
    {
        double ANGBAR = RADIAN * std::asin(1.0 / TAU);
        ANGSTP = DABS(ANGSTP);
        ABAR = ANGSTP * (std::trunc((ANGBAR - ANMIN) / ANGSTP + 1.e-10)
            + ANMIN);
    }

L50:
    {
        double TABAR = 2 * ABAR;
        ANMAX = DMIN1(ANMAX, TABAR);
        ANGSTP = DSIGN(ANGSTP, ANMAX - ANMIN);
        int NUMANG = (int)((ANMAX - ANMIN) / ANGSTP + 1.5);

        //
        // MAKE A SPACE EVERY ONCE AND A WHILE.
        //
        double ANGBLK = 52;
        if (NUMANG <= 50) goto L80;
        ANGBLK = 10;
        if (ANGSTP >= 1) goto L80;
        ANGBLK = 5;
        if (ANGSTP >= .5) goto L80;
        ANGBLK = 1;
L80:    int LINECT = 1;
        int NUMBLK = (int)(ANGBLK / ANGSTP + .5);
        int LINEMX = 53;

        if (NSPL == 1 || ISTAT == 3) goto L100;
        std::printf("0**** CANNOT EVALUATE ELASTIC SCATTERING FOR "
                    "IDENTICAL PARTICLES WITH S.O. FORCE\n");
        return;

        //
        // ALLOCATE SPACE FOR THETA-INDEPENDENT ARRAYS
        //
L100:   int LXMX = 0;
        if (NSPL > 1) LXMX = JSP;
        int LDELMX = LXMX - MOD(LXMX, 2);
        int NTEMPS = (LXMX + 1) * (LXMX + 1) * (LDELMX + 1);
        int ITEMPS = NALLOC("TEMPS   ", NTEMPS);
        int IBETA = NALLOC("BETAELAS", 2 * NSPL * (LMX + 1));

        SETLOG(2 * (LMX + LXMX));

        //
        // CALCULATE "BETAELAS", THE ANGLE-INDEPENDENT PARTS OF F.
        // ALSO CALCULATE TOTAL REACTION AND NUCLEAR CROSS SECTIONS.
        //
        int LTOC = FACFR4_v * Z[ITOC] - FACFR4_v + 1;
        int LBETA = FACFR4_v * Z[IBETA] - FACFR4_v + 1;
        double SIGNUC;
        BETCAL(TRUE_F, UK, JSP, NSPL, LMN, LMX, LSKP,
            0, LXMX, LDELMX, NTEMPS, ISTAT, PBETSW,
            SIGNUC, SIGREA,
            &ILLOC(LTOC), &ALLOC(Z[ISMAT]), DUM4, DUM4, &ALLOC(Z[ISIG]),
            DUM8, &ALLOC4(LBETA), DUM8, DUM8, ALOWFC_arr, &ALLOC(Z[ITEMPS]));

        //
        // RELEASE TEMPORARY ARRAY.
        //
        Z[ITEMPS] = -Z[ITEMPS];

        //
        // ALLOCATE SPACE FOR LEGENDRE FUNCTIONS, AMPLITUDES, CROSS SECTIONS.
        //
        int ISPACE = LMX + 1 + ((2 * LMX + 1 - LXMX) * LXMX) / 2;
        int IPLM = NALLOC("PLM     ", ISPACE);
        int ITMPCN = NALLOC("TMPCONTR", 2 * (LMX + 1));
        IDNTSW = (ISTAT != 3) ? TRUE_F : FALSE_F;
        int IDNPAR = 0;
        if (IDNTSW) IDNPAR = 2;
        int NF = 2 * NSPL;
        if (IDNTSW) NF = 6;
        int IFTEMP = NALLOC("FTEMP   ", NF);
        int IFLO = NALLOC("FLO     ", NF);
        int IFHI = NALLOC("FHI     ", NF);
        int IFERR = NALLOC("FERROR  ", NF);
        int IFEPSL = NALLOC("FEPSILON", 2 * LEBACK + 1);
        IFEL = 0;
        if (FKEPSW) IFEL = NALLOC(FNAME, NF * NUMANG);
        int ICROSS = NALLOC("CROSSSEC", NUMANG);
        ITORUT = NALLOC(RUTNAM, NUMANG);
        int ILAB = 0;
        if (PRNTSW) ILAB = NALLOC("LABCROSS", NUMANG);
        IRUTH = NALLOC(RUTHNM, NUMANG);

        //
        // ALLOCATIONS ARE DONE, GET THE BASE ADDRESSES.
        //
        int LFTEMP = Z[IFTEMP] - 2;
        int LFLO = Z[IFLO] - 2;
        int LFHI = Z[IFHI] - 2;
        int LFERR = Z[IFERR] - 1;
        LTOC = Z[ITOC] * FACFR4_v - FACFR4_v;
        int LF = Z[IFEL] - 2 * NUMANG - 2;
        LBETA = FACFR4_v * Z[IBETA] - FACFR4_v + 1;
        int LPLM = Z[IPLM];
        int LTMPCN = Z[ITMPCN];

        double STATFC = 0;
        if (ISTAT < 3) STATFC = (ISTAT - 1.50e0) * JSP / (JSP + 1.e0);
        int LHI = LMX - 2 - MOD(LMX, LSKP);
        double SIGZRO = ALLOC(Z[ISIG]);

        //
        // LOOP OVER ALL ANGLES
        //
        for (int NANG = 1; NANG <= NUMANG; NANG++) {

            double ANGLE = ANMIN + (NANG - 1) * ANGSTP;

            if (LABANG == 0) goto L120;
            {
                double AN = DEGREE * ANGLE;
                if (ANGLE > ABAR) AN = DEGREE * (TABAR - ANGLE);
                double TERM2 = 1 - (TAU * DSIN(AN)) * (TAU * DSIN(AN));
                if (TERM2 < 0) TERM2 = 0;
                TERM2 = DCOS(AN) * DSQRT(TERM2);
                if (ANGLE > ABAR) TERM2 = -TERM2;
                ANGLE = RADIAN * std::acos(-TAU * (DSIN(AN)) * (DSIN(AN)) + TERM2);
            }

L120:       {
                double COSAN = DCOS(ANGLE * DEGREE);
                double ANGLAB = RADIAN * DATAN2(DSIN(ANGLE * DEGREE), COSAN + TAU);

                double TERM2 = 1 + TAU * (TAU + 2 * COSAN);
                double TEMP = 1 + TAU * COSAN;
                double AJACOB = 0;
                if (TERM2 < 0) TERM2 = 0;
                if (TEMP != 0) AJACOB = DABS(TERM2 * DSQRT(TERM2) / TEMP);

                //
                // CALCULATE THE AMPLITUDES FOR THIS ANGLE.
                //
                int NFLPAM = 0;
                AMPCAL(ANGLE, TRUE_F, NSPL, LMN, LMX, LSKP,
                       LXMX, IDNPAR, ETA, UK, SIGZRO, PRNTSW, LHI,
                       &ILLOC(LTOC + 1), &ALLOC4(LBETA), ALOWFC_arr,
                       &ALLOC(LFTEMP + 2), &ALLOC(LFLO + 2), &ALLOC(LFHI + 2),
                       &ALLOC(LFERR + 1), FCOUL,
                       &ALLOC(LPLM), &ALLOC(LTMPCN), &ALLOC(LTMPCN + LMX + 1),
                       LEBACK, &ALLOC(Z[IFEPSL]), NFLPAM);

                //
                // CALCULATE THE CROSS SECTIONS.
                //
                double SIGMA = 0;
                double SIGLOW = 0;
                double SIGHI = 0;
                double RUTH;
                if (IDNTSW) goto L250;
                for (int KOFFS = 1; KOFFS <= NSPL; KOFFS++) {
                    if (ILLOC(LTOC + 4 * KOFFS) < 0) goto L229;
                    { int MX = (ILLOC(LTOC + 4 * KOFFS - 3) + ILLOC(LTOC + 4 * KOFFS - 2) + 1) / 2;
                    TEMP = 10.;
                    if (MX != 0) TEMP = 20.;
                    SIGMA = SIGMA + TEMP *
                        (ALLOC(LFTEMP + 2 * KOFFS) * ALLOC(LFTEMP + 2 * KOFFS)
                       + ALLOC(LFTEMP + 2 * KOFFS + 1) * ALLOC(LFTEMP + 2 * KOFFS + 1));
                    if (!PRNTSW) goto L229;
                    SIGLOW = SIGLOW + TEMP *
                        (ALLOC(LFLO + 2 * KOFFS) * ALLOC(LFLO + 2 * KOFFS)
                       + ALLOC(LFLO + 2 * KOFFS + 1) * ALLOC(LFLO + 2 * KOFFS + 1));
                    SIGHI = SIGHI + TEMP *
                        (ALLOC(LFHI + 2 * KOFFS) * ALLOC(LFHI + 2 * KOFFS)
                       + ALLOC(LFHI + 2 * KOFFS + 1) * ALLOC(LFHI + 2 * KOFFS + 1));
                    } // end scope for MX
L229:               ;
                }
                RUTH = (F2(FCOUL, 2, 1, 3) * F2(FCOUL, 2, 1, 3)
                      + F2(FCOUL, 2, 2, 3) * F2(FCOUL, 2, 2, 3)) * 10.;
                goto L260;

                //
                // IDENTICAL PARTICLES (NO SPIN-DEPENDENT FORCES).
                //
                // IDENTICAL PARTICLES: FCOUL(r,c) via F2(FCOUL,2,r,c)
L250:           RUTH = (F2(FCOUL, 2, 1, ISTAT) * F2(FCOUL, 2, 1, ISTAT)
                    + F2(FCOUL, 2, 2, ISTAT) * F2(FCOUL, 2, 2, ISTAT)
                    + STATFC * (F2(FCOUL, 2, 1, 1) * F2(FCOUL, 2, 1, 1)
                    + F2(FCOUL, 2, 2, 1) * F2(FCOUL, 2, 2, 1)
                    - F2(FCOUL, 2, 1, 2) * F2(FCOUL, 2, 1, 2)
                    - F2(FCOUL, 2, 2, 2) * F2(FCOUL, 2, 2, 2))) * 10.;
                SIGMA = (ALLOC(LFTEMP + 2 * ISTAT) * ALLOC(LFTEMP + 2 * ISTAT)
                    + ALLOC(LFTEMP + 2 * ISTAT + 1) * ALLOC(LFTEMP + 2 * ISTAT + 1)
                    + STATFC * (ALLOC(LFTEMP + 2) * ALLOC(LFTEMP + 2) + ALLOC(LFTEMP + 3) * ALLOC(LFTEMP + 3)
                    - ALLOC(LFTEMP + 4) * ALLOC(LFTEMP + 4) - ALLOC(LFTEMP + 5) * ALLOC(LFTEMP + 5))) * 10.;
                if (!PRNTSW) goto L260;
                SIGLOW = (ALLOC(LFLO + 2 * ISTAT) * ALLOC(LFLO + 2 * ISTAT)
                    + ALLOC(LFLO + 2 * ISTAT + 1) * ALLOC(LFLO + 2 * ISTAT + 1)
                    + STATFC * (ALLOC(LFLO + 2) * ALLOC(LFLO + 2) + ALLOC(LFLO + 3) * ALLOC(LFLO + 3)
                    - ALLOC(LFLO + 4) * ALLOC(LFLO + 4) - ALLOC(LFLO + 5) * ALLOC(LFLO + 5))) * 10.;
                SIGHI = (ALLOC(LFHI + 2 * ISTAT) * ALLOC(LFHI + 2 * ISTAT)
                    + ALLOC(LFHI + 2 * ISTAT + 1) * ALLOC(LFHI + 2 * ISTAT + 1)
                    + STATFC * (ALLOC(LFHI + 2) * ALLOC(LFHI + 2) + ALLOC(LFHI + 3) * ALLOC(LFHI + 3)
                    - ALLOC(LFHI + 4) * ALLOC(LFHI + 4) - ALLOC(LFHI + 5) * ALLOC(LFHI + 5))) * 10.;

                //
                // SAVE F IF NECESSARY.
                //
L260:           if (!FKEPSW) goto L280;
                {
                    //
                    // NF/2 = NSPL EXCEPT =3 FOR IDENTICAL PARTICLES
                    //
                    int NFBY2 = NF / 2;
                    for (int K = 1; K <= NFBY2; K++) {
                        ALLOC(LF + 2 * NANG + 2 * NUMANG * K) = ALLOC(LFTEMP + 2 * K);
                        ALLOC(LF + 2 * NANG + 2 * NUMANG * K + 1) = ALLOC(LFTEMP + 2 * K + 1);
                    }
                }

                //
                // COMPUTE RELATIVE ERRORS FROM LOWER AND UPPER L-TRUNCATIONS
                //
L280:           {
                    double EL = 50 * (1 - SIGLOW / SIGMA);
                    double EH = 50 * (1 - SIGHI / SIGMA);
                    double SIGMAR = 0;
                    if (RUTH != 0) SIGMAR = SIGMA / RUTH;

                    if (!PRNTSW) goto L400;

                    { double SIGLAB = AJACOB * SIGMA;

                    //
                    // PRINT EVERYTHING IN MILLIBARNS
                    //
                    LINECT = LINECT + 1;
                    if (LINECT <= LINEMX || NANG == NUMANG) goto L300;
                    // Reprint header
                    std::printf("1");
                    for (int i = 0; i < 47; i++) std::printf(" ");
                    std::printf("P T O L E M Y\n");
                    std::printf(" ");
                    for (int i = 1; i <= 45; i++) std::printf("%c", HEDCOM.REACT[i]);
                    std::printf("%8.2f MEV     ", ELAB);
                    for (int i = 1; i <= 65; i++) std::printf("%c", HEDCOM.HEADER[i]);
                    std::printf("\n");
                    std::printf("0    ANGLE");
                    std::printf("%*sSIGMA/", 11, "");
                    std::printf("%*sSIGMA", 19, "");
                    std::printf("%*sRUTHERFORD", 19, "");
                    std::printf("%*s%% PER", 19, "");
                    std::printf("%*s%% PER\n", 9, "");
                    std::printf("  C.M.    LAB     RUTHERFORD");
                    std::printf("%*sC.M.", 6, "");
                    std::printf("%*sLAB", 13, "");
                    std::printf("%*sC.M.", 15, "");
                    std::printf("%*sLOW L", 16, "");
                    std::printf("%*sHIGH L\n", 8, "");
                    std::printf("\n");
                    LINECT = 1;

L300:               std::printf("%7.2f%7.2f%14.4G%14.4G%12.4G%15.4G%13.2f%8.2f\n",
                        ANGLE, ANGLAB, SIGMAR, SIGMA, SIGLAB, RUTH, EL, EH);

                    if (DMOD(ANGLE + ANGSTP + 1.e-7, ANGBLK) > 1.e-5) goto L380;
                    std::printf(" \n");
                    LINECT = LINECT + 1;
                    if (LINECT + NUMBLK > LINEMX) LINECT = LINEMX + 1;
L380:               ALLOC(Z[ILAB] - 1 + NANG) = SIGLAB;
                    } // end scope for SIGLAB
L400:               ALLOC(Z[ITORUT] - 1 + NANG) = SIGMAR;
                    ALLOC(Z[IRUTH] - 1 + NANG) = RUTH;
                    if (LABANG != 0 && !PRNTSW)
                        SIGMA = AJACOB * SIGMA;
                    ALLOC(Z[ICROSS] - 1 + NANG) = SIGMA;
                }
            }
        }

        if (PRNTSW) {
            std::printf("\n0TOTAL REACTION CROSS SECTION =%13.3f MB\n", SIGREA);
            std::printf("  NUCLEAR TOTAL CROSS SECTION =%13.3f MB\n", SIGNUC);
        }

        //
        // FREE ALLOCATED SPACE
        //
        Z[IPLM] = -Z[IPLM];
        Z[IBETA] = -Z[IBETA];
        Z[IFTEMP] = -Z[IFTEMP];
        Z[IFLO] = -Z[IFLO];
        Z[IFHI] = -Z[IFHI];
    }
    return;
}
