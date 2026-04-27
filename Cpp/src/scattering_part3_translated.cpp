// scattering_part3_translated.cpp — WAVEF, WAVINH, WFGET
// Translated from source.f lines 34863-35430, 36902-36989, 38022-38103
// Verbatim transliteration from Fortran to C++

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include "ptolemy_alloc.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_forward.h"
#include "ptolemy_io.h"
#include <cstdio>
#include <cstring>
#include <cmath>

// ============================================================================
// SUBROUTINE WAVEF ( IRET )
//
//  WAVEF CALCULATES RADIAL WAVEFUNCTION FOR ALL PARTIAL WAVES
//  USING WOOD-SAXON POTENTIAL
//
//     THIS ROUTINE IS USED ONLY FOR STAND-ALONE CALCULATIONS OF
//     ELASTIC SCATTERING.  IT IS AN INTERFACE TO WAVSET, WAVPOT
//     AND WAVELJ WITH A COUPLE OF BELLS AND WHISTLES ADDED.
//
//     IF L IS DEFINED THEN THE CALC IS DONE FOR ONLY ONE VALUE OF
//     L; ELSE THE RANGE (LMIN, LMAX) IS USED.
//     IF THERE IS A SPIN-ORBIT FORCE THEN JP IS CHECKED.  IF IT IS
//     DEFINED THEN ONLY THAT JP (AND THE ASSOCIATED L'S) ARE DONE.
//     IF IT IS NOT DEFINED THEN ALL JP'S FOR EACH L ARE DONE.
//
// Translated from source.f L34863-35430
// ============================================================================
void WAVEF(int& IRET)
{
    // implicit real*8 ( a-h, o-z )
    // Local variables (I-N are integer, A-H,O-Z are double)

    // References to COMMON blocks
    auto& L      = INTGER.L;
    auto& LMAX   = INTGER.LMAX;
    auto& LMIN   = INTGER.LMIN;
    auto& LX     = INTGER.LX;
    auto& IPRINT = INTGER.IPRINT;
    auto& LMAXAD = INTGER.LMAXAD;
    auto& LMINSB = INTGER.LMINSB;

    auto& JP     = JBLOCK.JP;
    auto& JSP    = JBLOCK.JSP;
    auto& JST    = JBLOCK.JST;

    auto& IASYMP = SWITCH.IASYMP;
    auto& ICHECK = SWITCH.ICHECK;
    auto& IELAST = SWITCH.IELAST;
    auto& IRLWAV = SWITCH.IRLWAV;
    auto& ISAVB  = SWITCH.ISAVB;

    auto& ICHANW = INTRNL.ICHANW;

    auto& ANGMIN = FLOAT_common.ANGMIN;
    auto& ANGMAX = FLOAT_common.ANGMAX;
    auto& ANGSTP = FLOAT_common.ANGSTP;
    auto& ELAB   = FLOAT_common.ELAB;
    auto& E      = FLOAT_common.E;
    auto& STEPSZ = FLOAT_common.STEPSZ;
    auto& WRITES = FLOAT_common.WRITES;
    auto& ALMNMT = FLOAT_common.ALMNMT;
    auto& ALMXMT = FLOAT_common.ALMXMT;

    auto& RADIAN = CNSTNT.RADIAN;
    auto& DEGREE = CNSTNT.DEGREE;

    auto* Z = LOCPTRS.Z;   // 1-based

    int NOTDEF = NOTDEF_INT;
    double UNDEF = INTRNL.UNDEF;

    // CHARACTER*8 locals
    char8 RUTNAM("TORUTHER");
    char8 RUTHNM("RUTHERFO");

    // LOGICAL locals
    int ONEJSW, SOSW, TCSW, THISTC, FKEPSW;

    // Other locals
    int IPRNT, LLMAX, LLMIN, LSAVE;
    double JPSAVE;  // must be double: preserves NOTDEF bit pattern (int cast loses it)
    int LSKIP, ISTAT, NSTEP, NSTEPS, NSPL;
    int LOCWRK, NFJ, NGJ;
    int LWAVR, LWAVI, LSMATS, LTOC;
    int LFJ, LGJ, LWORK;
    int LINES, ISPST, ISPND, ISP;
    int I, II, KOFFS, LL, NWRITE;
    int IFLAG, ISTRT, LREALV, LIMAGV, LCENTR;
    int ITORUT, IFEL, IRUTH, IAPOW;
    double SPNAV, ETA, H, TL, SR, SI, TEMP, PHASE;
    double ROTR, ROTI, FR, FI, RHOO, RVAL, REALF, AMAGF, RATR, RR;
    double SIGREA;

    //
    IPRNT = MOD(IPRINT, 10);
    if (IPRNT > 0) {
        std::printf("1%47s%s\n", "", "P T O L E M Y");
        std::printf("0");
        for (int k = 1; k <= 45; k++) std::printf("%c", HEDCOM.REACT[k]);
        std::printf("%8.2f MEV     ", ELAB);
        for (int k = 1; k <= 65; k++) std::printf("%c", HEDCOM.HEADER[k]);
        std::printf("\n");
    }

    //
    //     TEST THE INPUT AND RETURN  0  COMPLETION CODE IF BAD
    //
    IRET = 0;
    LLMAX = (int)LMAX;
    LLMIN = (int)LMIN;
    LSAVE = L;
    JPSAVE = JP;  // save as double to preserve NOTDEF bit pattern

    //
    if (L == NOTDEF) goto L60;

    //
    //     DOING ONLY ONE L
    //
    LLMAX = L;
    LLMIN = L;

    //
L60:
    if (ICHANW > 2) ICHANW = 1;

    //
    //     SET UP ARRAYS AND PRINT PARAMETERS.
    //
    WAVSET(IRET, TRUE_F, FALSE_F);
    if (IRET == 0) goto L9900;
    IRET = 0;
    LSKIP = WAVCOM.LSKIPS[ICHANW];
    ISTAT = WAVCOM.ISTATS[ICHANW];

    //
    if (LLMIN == NOTDEF || LLMAX == NOTDEF) goto L80;

    //
    //     LSKIP = 2 ONLY FOR IDENTICAL SPIN 0 SCATTERING
    //
    if (LSKIP == 2) LLMIN = IABS(LLMIN - MOD(LLMIN, 2));
    if (LSKIP == 2) LLMAX = LLMAX + MOD(LLMAX, 2);
    goto L100;

    //
    //     GET ESTIMATE OF L CRITICAL
    //
L80:
    KANDM.LCRIT = KANDM.LCRITS[ICHANW];

    //
    if (LLMIN != NOTDEF) goto L85;
    LLMIN = (int)(KANDM.LCRIT * ALMNMT);
    LLMIN = MIN0(LLMIN, KANDM.LCRIT - LMINSB);
    LLMIN = MAX0(LLMIN, 0);
    if (LSKIP == 2) LLMIN = IABS(LLMIN - MOD(LLMIN, 2));
L85:
    if (LLMAX == NOTDEF)
        LLMAX = MAX0(KANDM.LCRIT + LMAXAD, INT(ALMXMT * KANDM.LCRIT));
    if (LSKIP == 2) LLMAX = LLMAX + MOD(LLMAX, 2);

    //
L100:
    std::printf("0%10d =< L =<%4d\n", LLMIN, LLMAX);

    //     ARE WE DOING ONLY ONE JP
    //
    SOSW = WAVCOM.SOSWS[ICHANW];
    TCSW = WAVCOM.TCSWS[ICHANW];
    ONEJSW = (JP != INTRNL.NOTDEF) ? TRUE_F : FALSE_F;
    ISPST = 0;
    ISPND = 0;
    SPNAV = 1;
    if (!ftobool(SOSW)) goto L200;
    SPNAV = 1.0 / (JSP + 1);
    if (ftobool(ONEJSW)) goto L150;
    //     DOING ALL JP FOR EACH L
    ISPST = -(int)JSP;
    ISPND = (int)JSP;
    goto L200;
L150:
    if (MOD((int)JP + (int)JSP, 2) == 0) goto L160;
    std::printf("0**** JP AND SP ARE AN INVALID HALF-INTEGER AND "
                "INTEGER COMBINATION:%6d/2%6d/2\n", (int)JP, (int)JSP);
    goto L9900;
    //     ONLY DOING ONE JP, REDUCE L-RANGE
L160:
    LLMIN = MAX0(LLMIN, ((int)JP - (int)JSP) / 2);
    LLMAX = MIN0(LLMAX, ((int)JP + (int)JSP) / 2);
    if (LLMIN <= LLMAX) goto L200;
    std::printf("0**** THERE ARE NO VALID COMBINATIONS OF JP, SP AND L:"
                "   JP, SP =%5d/2%5d/2\n LMIN, LMAX =%5d%5d  OR  L =%5d\n",
                (int)JP, (int)JSP, LMIN, LMAX, L);
    goto L9900;

    //
L200:
    KANDM.LOMOST = LLMAX;
    if (ftobool(TCSW) && L != NOTDEF) KANDM.LOMOST = LLMAX + 2;

    //
    //     SET UP S-MATRIX AND COULOMB WF ARRAYS.
    //
    WAVPOT(IRET);
    if (IRET == 0) goto L9900;
    IRET = 0;
    NSTEP = WAVCOM.NSTPSS[ICHANW];
    NSTEPS = NSTEP + 1;
    ETA = KANDM.ETAS[ICHANW];
    H = WAVCOM.HS[ICHANW];
    NSPL = WAVCOM.NUMJS[ICHANW];

    //
    //     IF NECESSARY, GET WORK AREAS FOR THE ASYMPTOTIC CHECK
    //
    if (ICHECK == 0) goto L210;
    LOCWRK = NALLOC("ASYWRK  ", 2 * (LLMAX + 1));
    NFJ = NALLOC("ASYF    ", LLMAX + 1);
    NGJ = NALLOC("ASYG    ", LLMAX + 1);

    //
L210:
    ;

    //
    //     ALL ALLOCATIONS DONE, GET LOCATIONS
    //
    LWAVR = Z[WAVCOM.IWAVR];
    LWAVI = Z[WAVCOM.IWAVI];
    LSMATS = Z[KANDM.ISMATS[ICHANW]] - 1;
    LTOC = Z[WAVCOM.ITOCE[ICHANW]] * ALLOCS.FACFR4 - ALLOCS.FACFR4;
    if (ICHECK == 0) goto L220;
    LFJ = Z[NFJ];
    LGJ = Z[NGJ];
    LWORK = Z[LOCWRK];

    //
    //  SOME DERIVED DATA
    //
    //
    //     INITIALIZE S-MATRIX ELEMENTS TO ZERO.
    //
L220:
    INIT8(KANDM.ISMATS[ICHANW], 0.0);

    //
    if (LLMIN == LLMAX) IPRNT = 1;

    //
    //     PRINT HEADINGS IF NECESSARY.
    //
    LINES = 35;
    if (IPRNT > 0) {
        std::printf("\n0 L    L'   LX%19sS%24s|S|%10sPHASE SHIFT%11sFRACTION\n",
                    "", "", "", "");
        std::printf("%70sDEGREES%11sABSORBED\n", "", "");
        std::printf(" \n");
    }

    //
    // L1260:
    ;

    //
    //     THIS BRINGS IN THE TABLES NEEDED BY RACAH.
    //
    DUMMY1();
    LREALV = Z[WAVCOM.IRLVS[ICHANW]];
    LIMAGV = Z[WAVCOM.IIMVS[ICHANW]];
    LCENTR = Z[WAVCOM.ICENTR[ICHANW]];

    //
    //     LOOP THROUGH INCIDENT L.
    //
    THISTC = FALSE_F;
    for (L = LLMIN; L <= LLMAX; L += LSKIP) {

        //
        //        LOOP THROUGH JP = 2*J(PROJECTILE).  SPIN OF TARGET IS IGNORED.
        //        S-MATRIX ELEMENTS FOR ALL LX ARE ACCUMULATED AS SUMS OVER JP.
        //
        for (ISP = ISPST; ISP <= ISPND; ISP += 2) {
            if (!ftobool(ONEJSW)) JP = 2 * L + ISP;
            if (ICHANW == 2 && L <= 2)
            THISTC = (ftobool(TCSW) && (int)JP != 2 * L && (int)JP - L >= 0)
                     ? TRUE_F : FALSE_F;

            //
            //           SUBROUTINE WFGET CALLS EITHER WAVELJ (NO TENSOR) OR WAVETC
            //           (TENSOR COUPLING) TO DO THE ACTUAL CALCULATION.
            //
            WFGET(L, L, (int)JP, ICHANW, 0, &ALLOC4(1), &ALLOC4(1), &ALLOC4(1),
                  ALLOC_base(LWAVR), ALLOC_base(LWAVI), ALLOC_base(LREALV), ALLOC_base(LIMAGV),
                  ALLOC_base(LCENTR));

            //
        } // end ISP loop (259)

        //
        //        END OF CALCULATION LOOP THROUGH JP.
        //        IF NECESSARY, RETRIEVE AND PRINT THE S-MATRIX ELEMENTS.
        //
        if (IPRNT == 0) goto L399;

        //
        //        PRINT HEADINGS IF NECESSARY.
        //
        if (LINES + NSPL < 59) goto L270;
        std::printf("1%47s%s\n", "", "P T O L E M Y");
        std::printf("0");
        for (int k = 1; k <= 45; k++) std::printf("%c", HEDCOM.REACT[k]);
        std::printf("%8.2f MEV     ", ELAB);
        for (int k = 1; k <= 65; k++) std::printf("%c", HEDCOM.HEADER[k]);
        std::printf("\n");
        std::printf("\n0 L    L'   LX%19sS%24s|S|%10sPHASE SHIFT%11sFRACTION\n",
                    "", "", "", "");
        std::printf("%70sDEGREES%11sABSORBED\n", "", "");
        std::printf(" \n");
        LINES = 7;

        //
        //        CALCULATE THE FRACTION ABSORBED FROM THIS L.
        //
L270:
        TL = 1.0;
        I = LSMATS + 2 * NSPL * L - 1;
        for (KOFFS = 1; KOFFS <= NSPL; KOFFS++) {
            I = I + 2;
            TL = TL - ALLOC(I) * ALLOC(I) - ALLOC(I + 1) * ALLOC(I + 1);
        } // end KOFFS loop (279)

        //
        //        LOOP THROUGH KOFFS = LX, LAS.  QUANTUM NUMBERS ARE
        //        IN THE TABLE-OF CONTENTS ARRAY.
        //
        for (KOFFS = 1; KOFFS <= NSPL; KOFFS++) {
            LX = ILLOC(LTOC + 4 * KOFFS - 2);
            LL = L + ILLOC(LTOC + 4 * KOFFS - 3);
            if (LX > L + LL) goto L299;

            //
            //           RETRIEVE THE S-MATRIX ELEMENT, AND CALCULATE THE
            //           MAGNITUDE, PHASE, AND FRACTION ABSORBED.
            //
            I = LSMATS + 2 * (NSPL * L + KOFFS) - 1;
            SR = ALLOC(I);
            SI = ALLOC(I + 1);
            TEMP = DSQRT(SR * SR + SI * SI);
            if (TEMP == 0) goto L299;
            PHASE = 0.5 * RADIAN * DATAN2(SI, SR);

            //
            //           PRINT IT.
            //
            if (KOFFS == 1) goto L290;
            std::printf(" %3d%5d%5d%17.5G +%12.5G I%17.5G%12.2f\n",
                        L, LL, LX, SR, SI, TEMP, PHASE);
            goto L295;
L290:
            std::printf(" %3d%5d%5d%17.5G +%12.5G I%17.5G%12.2f%16.5G\n",
                        L, LL, LX, SR, SI, TEMP, PHASE, TL);
L295:
            LINES = LINES + 1;
L299:
            ;
        } // end KOFFS loop (299)

        //
        //        END OF KOFFS LOOP FOR PRINTING.  INSERT A BLANK LINE IF NEEDED.
        //
        if (NSPL == 1) goto L399;
        std::printf(" \n");
        LINES = LINES + 1;

        //
L399:
        ;
    } // end L loop (399)

    //
    //     END OF LOOP OVER ALL L'S AND J'S
    //
    L = LLMAX;

    //
    //     FOR SINGLE WAVEFUNCTION, REDEFINE WAVER AND WAVEI TO BE
    //     THE CORRECT SIZE WITH NO SLOP AT THE END
    //
    if (LLMAX == LLMIN && !ftobool(THISTC))
        IREDEF(NSTEP + 1, WAVCOM.IWAVR);
    if (LLMAX == LLMIN && !ftobool(THISTC))
        IREDEF(NSTEP + 1, WAVCOM.IWAVI);

    //
    LWAVR = Z[WAVCOM.IWAVR];
    LWAVI = Z[WAVCOM.IWAVI];
    ISTRT = WAVCOM.LASTNF[ICHANW] - 2;

    //
    //     IF REQUESTED, ATTEMPT TO MAKE SCATTERING WAVE REAL
    //
    ROTR = 1;
    ROTI = 0;
    if (IRLWAV == 0) goto L500;
    ROTR = DCOS(DEGREE * PHASE);
    ROTI = -DSIN(DEGREE * PHASE);
    for (II = ISTRT; II <= NSTEP; II++) {
        FR = ROTR * ALLOC(LWAVR + II) - ROTI * ALLOC(LWAVI + II);
        ALLOC(LWAVI + II) = ROTR * ALLOC(LWAVI + II) + ROTI * ALLOC(LWAVR + II);
        ALLOC(LWAVR + II) = FR;
    } // end II loop (459)
    if (!ftobool(THISTC)) goto L500;
    for (II = ISTRT; II <= NSTEP; II++) {
        FR = ROTR * ALLOC(LWAVR + II + WAVCOM.NWAVEF)
             - ROTI * ALLOC(LWAVI + II + WAVCOM.NWAVEF);
        ALLOC(LWAVI + II + WAVCOM.NWAVEF) = ROTR * ALLOC(LWAVI + II + WAVCOM.NWAVEF)
             + ROTI * ALLOC(LWAVR + II + WAVCOM.NWAVEF);
        ALLOC(LWAVR + II + WAVCOM.NWAVEF) = FR;
    } // end II loop (469)

    //
    //     WF TABULATION INTERVAL IS 'WRITESTEP' KEYWORD.
    //
L500:
    NWRITE = 1;
    if (WRITES != UNDEF) NWRITE = (int)DMAX1(WRITES / STEPSZ + 0.5, 1.0);

    //
    //     IF REQUESTED ('CHECKASYM'), DO ASYMPTOTIC CHECK.
    //
    if (ICHECK == 0) goto L600;
    std::printf("         ** ASYMTOTIC CHECK ** \n"
                "  SOLUTION SHOULD GO TO ASYMPTOTIC FORM\n\n"
                "    R%16sF(R)%15sASYMPTOTIC FORM%7s"
                "REL. ERROR\n\n", "", "", "");

    //
    for (II = NWRITE; II <= NSTEP; II += NWRITE) {
        if (II < ISTRT) goto L559;
        RHOO = H * II;
        RVAL = STEPSZ * II;
        RCWFN(RHOO, ETA, L, L, &ALLOC(LFJ), &ALLOC(LWORK),
              &ALLOC(LGJ), &ALLOC(LWORK + LLMAX + 1), 1.0e-8, IFLAG);
        if (IFLAG != 0) goto L9920;
        REALF = (ALLOC(LFJ + L) * (1.0 + SR) + SI * ALLOC(LGJ + L)) / 2.0;
        AMAGF = (ALLOC(LGJ + L) * (1.0 - SR) + SI * ALLOC(LFJ + L)) / 2.0;
        FR = ROTR * REALF - ROTI * AMAGF;
        AMAGF = ROTR * AMAGF + ROTI * REALF;
        REALF = FR;
        FR = ALLOC(LWAVR + II);
        FI = ALLOC(LWAVI + II);
        RATR = (DABS(FR - REALF) + DABS(FI - AMAGF))
             / (DABS(FR) + DABS(FI));
        std::printf("%7.2f  %11.4G%11.4G    %11.4G%11.4G%14.2G\n",
                    RVAL, FR, FI, REALF, AMAGF, RATR);
L559:
        ;
    } // end II loop (559)

    //
    //     FREE ALLOCATED SPACE
    //
    Z[NFJ] = -Z[NFJ];
    Z[NGJ] = -Z[NGJ];
    Z[LOCWRK] = -Z[LOCWRK];
    goto L700;

    //
    //
    //     ARE WE TO PRINT THE WAVEFUNCTION
    //
L600:
    if (WRITES == UNDEF) goto L700;

    //
    if (!ftobool(SOSW)) std::printf("0OPTICAL MODEL WAVEFUNCTION FOR L =%4d\n\n", L);
    if (ftobool(SOSW))
        std::printf("0OPTICAL MODEL WAVE FUNCTION FOR L =%4d  J =%4d/2\n\n", L, (int)JP);
    if (ftobool(THISTC)) goto L670;
    std::printf("     R       REAL F       IMAG F\n");
    for (II = 1; II <= NSTEPS; II += NWRITE) {
        if (II < ISTRT) goto L659;
        RR = (II - 1) * STEPSZ;
        std::printf("%8.3f  %13.5G%13.5G\n", RR, ALLOC(LWAVR + II - 1), ALLOC(LWAVI + II - 1));
L659:
        ;
    } // end II loop (659)
    goto L700;

    //
    //     IF TENSOR COUPLING EXISTS, PRINT BOTH PARTS OF COUPLED WF.
    //
L670:
    {
        int LAS = (int)JP - L;
        std::printf("%19sLAS =%3d%3s%19sLAS =%3d%3s\n", "", L, "", "", LAS, "");
        std::printf("     R%7sREAL F%7sIMAG F%4s%7sREAL F%7sIMAG F\n",
                    "", "", "", "", "");
        for (II = 1; II <= NSTEPS; II += NWRITE) {
            if (II < ISTRT) goto L679;
            RR = (II - 1) * STEPSZ;
            std::printf("%8.3f  %13.5G%13.5G     %13.5G%13.5G\n",
                        RR, ALLOC(LWAVR + II - 1), ALLOC(LWAVI + II - 1),
                        ALLOC(LWAVR + II - 1 + WAVCOM.NWAVEF),
                        ALLOC(LWAVI + II - 1 + WAVCOM.NWAVEF));
L679:
            ;
        } // end II loop (679)
    }

    //
    //     IF DESIRED, COMPUTE THE ELASTIC DIFFERENTIAL CROSS SECTIONS
    //
L700:
    if (IELAST == 0) goto L850;
    FKEPSW = (ISAVB == 1 || ftobool(SOSW)) ? TRUE_F : FALSE_F;
    ELDCS(H/STEPSZ, ETA, ANGMIN, ANGMAX, ANGSTP,
          LLMIN, LLMAX, LSKIP, ISTAT, (int)JSP,
          KANDM.ISMATS[ICHANW], WAVCOM.ITOCE[ICHANW],
          NSPL, KANDM.ISIGS[ICHANW], ELAB, ELAB/E - 1.0, TRUE_F,
          ITORUT, RUTNAM.data, FKEPSW,
          IFEL, "F       ", SIGREA,
          IRUTH, RUTHNM.data);

    //
    //     NOW COMPUTE THE ANALYZING POWERS (IF NONZERO).
    //
    if (!ftobool(SOSW)) goto L850;
    WAVCOM.PWBGSW = (MOD(IPRINT, 10) >= 3) ? TRUE_F : FALSE_F;
    ANAPOW(ANGMIN, ANGMAX, ANGSTP, (int)JSP, (int)JSP, (int)JST, (int)JST,
           NSPL, 1, TRUE_F, WAVCOM.PWBGSW,
           "ANPOW   ", ELAB, "ELASTIC ",
           IFEL, WAVCOM.ITOCE[ICHANW], IAPOW);

    //
    //     A L L   D O N E
    //
    //
L850:
    ;
    IRET = 1;
    goto L9900;

    //
    //
    //     ERROR RETURNS
    //
    //
    //     ERRORS IN RCWFN
L9920:
    std::printf("-***** ERROR RETURN FROM RCWFN: %18.8G%18.8G%18.8G%18.8G\n",
                (double)IFLAG, ETA, RHOO, (double)L);
    Z[LOCWRK] = -Z[LOCWRK];
    IRET = 0;

    //
    //     RESTORE COMMON
    //
L9900:
    L = LSAVE;
    JP = JPSAVE;
    return;

    //
}

// ============================================================================
// SUBROUTINE WAVINH ( NFIRST, NSTEP, WAVR, WAVI, RHSR, RHSI, FLP, NSTPDM )
//
//     NUMEROV LOOP FOR INHOMOGENEOUS WAVEFUNCTIONS
//
// Translated from source.f L36902-36989
// ============================================================================
void WAVINH(int NFIRST, int NSTEP, double* WAVR, double* WAVI,
            double* RHSR, double* RHSI, double& FLP, int NSTPDM)
{
    // implicit real*8 ( a-h, o-z )
    // WAVR, WAVI, RHSR, RHSI are 1-based arrays (dimension 1)

    //
    //      NUMEROV METHOD (A LA RAYNAL) IS USED IN THE FOLLOWING LOOP
    //
    //     WE DEFINE
    //       W(I) = 1 - (H**2/12)*(V(I)/E - 1)
    //       XSI(I) = W(I)*U(I)/12
    //     WHERE
    //       V(I) = V(I*STEPSZ) = V(I*H/K)
    //     CONTAINS THE COMPLETE NUCLEAR, COULOMB AND CENTRIFIGAL POTENTIALS.
    //     W(I) IS COMPUTED BY COUPLN BEFORE CALLING US.
    //
    //     THE DESIRED WAVEFUNCTION IS
    //       U(I) = U(I*STEPSZ)
    //     EXCEPT THAT HERE WE REALLY COMPUTE 12*U WHICH IS WHAT
    //     IS DESIGNATED AS U BELOW  ( U = 12*PSI*R ).
    //
    //     AT STEP I IN THE LOOP, THE WAVR AND WAVI ARRAYS CONTAIN:
    //
    //     WAV_(...)   CONTENTS
    //          I-1    U(I-2)     NOT USED
    //          I      U(I-1)
    //          I+1    XSI(I-2),   CHANGED TO U(I) DURING STEP
    //          I+2    XSI(I-1)
    //          I+3    W(I-1), NOT USED, CHANGED TO XSI(I) DURING STEP
    //          I+4    W(I)
    //
    //     THE SMALLEST POSSIBLE NFIRST IS 2 WHICH RESULTS IN ISTRT=0.
    //     THE U AND XSI FOR NFIRST-2, NFIRST-1 HAVE ALREADY BEEN STORED.
    //

    double D = 12.0 / (WAVR[NFIRST + 4] * WAVR[NFIRST + 4]
                      + WAVI[NFIRST + 4] * WAVI[NFIRST + 4]);

    //
    for (int I = NFIRST; I <= NSTEP; I++) {

        //
        WAVR[I + 3] = 10.0 * (RHSR[I] - WAVR[I + 2]) + (WAVR[I]
                    - WAVR[I + 1] + RHSR[I + 1] + RHSR[I - 1]);
        WAVI[I + 3] = 10.0 * (RHSI[I] - WAVI[I + 2]) + (WAVI[I]
                    - WAVI[I + 1] + RHSI[I + 1] + RHSI[I - 1]);
        WAVR[I + 1] = (WAVR[I + 4] * WAVR[I + 3]
                    + WAVI[I + 4] * WAVI[I + 3]) * D;
        WAVI[I + 1] = (WAVR[I + 4] * WAVI[I + 3]
                    - WAVI[I + 4] * WAVR[I + 3]) * D;
        D = (12.0 / (WAVR[I + 5] * WAVR[I + 5] + WAVI[I + 5] * WAVI[I + 5]));

        //
    } // end I loop (210)

    //
    FLP = 24 * (NSTEP - NFIRST + 1);

    //
    return;
}

// ============================================================================
// SUBROUTINE WFGET ( L, LAS, JP, NWP, NUMPTS, RGRID, WAVER, WAVEI,
//    WAVR, WAVI, VREAL, VIMAG, VCENT )
//
//     INTERFACE TO CALL WAVELJ OR WAVETC AS APPROPRIATE.
//
// Translated from source.f L38022-38103
// ============================================================================
void WFGET(int L, int LAS, int JP, int NWP, int NUMPTS,
           float* RGRID, float* WAVER, float* WAVEI,
           double* WAVR, double* WAVI, double* VREAL, double* VIMAG,
           double* VCENT)
{
    // implicit real*8 ( a-h, o-z )

    //
    if (WAVCOM.PWBGSW)
        std::printf(" WFGET:  L,LAS,JP,NWP,NUMPTS=%8d%8d%8d%8d%8d\n",
                    L, LAS, JP, NWP, NUMPTS);

    //
    //     CHECK FOR LEGALITY AND COUPLING.
    //
    if (L < 0) goto L200;
    if (!ftobool(WAVCOM.SOSWS[NWP])) goto L100;
    if (JP < IABS(2 * L - WAVCOM.JSPS[NWP])) goto L200;
    if (!ftobool(WAVCOM.TCSWS[NWP])) goto L100;
    if (LAS < 0 || LAS > KANDM.LOMOST) goto L200;
    if (JP == 2 * L || JP == 0 || JP >= 2 * KANDM.LOMOST)
        goto L100;

    //
    //     TENSOR COUPLED
    //
    WAVETC(L, LAS, JP, NWP, WAVCOM.NWAVEF, NUMPTS, RGRID,
           WAVER, WAVEI, WAVR, WAVI, VREAL, VIMAG, VCENT);
    return;

    //
    //     NOT TENSOR COUPLED
    //
L100:
    WAVELJ(L, JP, NWP, NUMPTS, RGRID, WAVER, WAVEI,
           WAVR, WAVI, VREAL, VIMAG, VCENT);
    return;

    //
    //     ILLEGAL - SET TO ZERO.
    //
L200:
    if (NUMPTS == 0) return;
    for (int I = 1; I <= NUMPTS; I++) {
        WAVER[I] = 0;
        WAVEI[I] = 0;
    } // end I loop (239)
    return;
}
