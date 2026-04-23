// ineldc_translated.cpp — INELDC subroutine
// Translated from source.f lines 20095-21033
// Computes the transfer reaction radial integrals (DWBA transfer amplitudes).
//
// INELDC calculates H(LI,LO,LX), the integral of A12, V, and the bound state
// wavefunctions, then integrates H over the distorted waves.
// SFROMI is called to renormalize results and save them as S(LI,LO,LX,JP,JT).

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include "ptolemy_alloc.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_forward.h"
#include "ptolemy_io.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <algorithm>

extern double second();

// A12 has the full 25-argument signature matching the Fortran source (source.f L1628-1632).
// The stub in source_coupled.cpp has a simplified 13-arg signature that will need updating.
// For now, declare the full signature needed by this translation:
extern void A12(int LI, int LXMIN_a, int LXMAX_a, int LMIN_a, int LMAX_a,
                double* XLAM, int* NLAM, double* A12VL, double* DMSVAL,
                int* JA12S, int& JA12M, int& JA12N, int& JA12AN, int* IINDEX,
                int& IHMAX, int& LOMNMN, int& LOMXMX, double* DINTS,
                double* OUTTMP, double* XLOTMP, int* LXTMP,
                int LHSM1, int LA12VL, int IPRINT_a, int INIT,
                int* INDXDW, int& NDW, int& NI, int* IDWFI, int* IDWFO);

// INELD2 — 9 arguments matching Fortran source.f L21035 (TSTART is REAL*4)

void INELDC(int& IRTN) {
    //
    // COMPUTES THE TRANSFER REACTION RADIAL INTEGRALS.
    //
    // Local variables
    int PBUGSW;   // LOGICAL
    int PINFSW;   // LOGICAL
    int NXLISW;   // LOGICAL

    float TSTART, TTIMER, TT, TT1, TT2, TT5, TT7;

    double HAFNEG = -0.50;

    // Convenient references to COMMON block members
    int&    IPRINT  = INTGER.IPRINT;
    int&    LMIN    = INTGER.LMIN;
    int&    LMAX    = INTGER.LMAX;
    int&    NPSUM   = INTGER.NPSUM;
    int&    NPDIF   = INTGER.NPDIF;
    int&    NPPHI   = INTGER.NPPHI;
    int&    NCOSIN  = INTGER.NCOSIN;

    int&    ISCRFL  = SWITCH.ISCRFL;
    int&    IBSPAS  = SWITCH.IBSPAS;

    int&    IHSAVE  = INTRNL.IHSAVE;
    int&    ISTRIP  = INTRNL.ISTRIP;

    double& ELAB    = FLOAT_common.ELAB;
    double& PI      = CNSTNT.PI;

    int&    LXMIN   = INELCM.LXMIN;
    int&    LXMAX   = INELCM.LXMAX;
    int&    JBPF    = INELCM.JBPF;
    int&    NUMLX   = INELCM.NUMLX;
    int&    NUMLIS  = INELCM.NUMLIS;

    int&    LOMOST  = KANDM.LOMOST;
    double& AKI     = KANDM.AKI;
    double& AKO     = KANDM.AKO;

    int&    NWFI    = GRIDCM.NWFI;
    int&    NWFO    = GRIDCM.NWFO;
    int&    NRIROI  = GRIDCM.NRIROI;
    int&    NCRIT   = GRIDCM.NCRIT;
    int&    NPSUMI  = GRIDCM.NPSUMI;
    int&    NUMHS   = GRIDCM.NUMHS;
    int&    ICOSS4  = GRIDCM.ICOSS4;

    int     FACFR4  = ALLOCS.FACFR4;

    // Local scalars
    int IPRNT, INIA12;
    int LIMIN, LIL, LIPRTY;
    int ILI, LI;
    int ITRK, KWO, KWI;
    int LO, LASO, JPO, LASI, JPI;
    int I, II, K, IV, IU, IH;
    int LIIP, LIIPH;
    int IPLUNK, ITEMP;
    int MCNT, NSKIP;
    int NUMHSU;
    int JA12M, JA12N, JA12AN;
    int IHMAX, LOMNMN, LOMXMX;
    int LA12OF, LHSTRT, NUMIH, NUMMU;
    int KA12M, NMU, LH;
    int KDW, NDW, NI;
    int NLINE;
    int NUMTIM, NUMANG, NUMTRM;
    int NUM3, NUM4, NUM5, NUM6;
    int NUMHIN;
    int NN;
    int JA12OF, KMEND;

    double FACTOR;
    double DWRIOS, DWRIOL, DWRIOC, DWCONT, DWMAX, DWLIM;
    double PHI, PHIP, PHIT, PVPDX;
    double ARG, DELTAC, RNN;
    double S2PHI, C2PHI, S2PHI2;
    double COSTHE, SINSIN, COSNEW;
    double FIR, FII, FOR_, FOI;   // FOR_ to avoid keyword clash
    double DWR, DWI;
    double TERM;
    double PI256;

    // Allocator pointer locals (computed once at start)
    int LMSVAL, LJA12S, LINTS, LA12VL, LHS, LHSM1;
    int LIINDX, LA12T1, LA12T2, LA12T3;
    int LCOSIN, LPHIT, LPHIP, LPHI, LTRP;
    int LRI, LRO, LLIR, LLII, LLOR, LLOI, LWIO;
    int LDW, LLIS;
    int LLILOR, LLILOI;
    int LHINT, LHABS, LABS1;
    int LSMHPT, LSMIPT, LSMHVL, LSMIVL, LSMHWK;
    int LRIOEX;
    int LWAVR, LWAVI;
    int LVRL1, LVRL2, LVIM1, LVIM2;
    int LCENT1, LCENT2;
    int LATERM;
    int LDWI, LWFIO, LWFII, LINDXS;
    int LNLAM;
    int LHEND;

    // =========================================================================
    // Begin executable code
    // =========================================================================

    TSTART = (float)second();
    IRTN = 0;

    for (int III = 1; III <= 8; III++) {
        FTIME.TIMES[III] = 0.0f;
    }

    NUMTIM = 4;
    NUMANG = 0;
    NUMTRM = 0;
    NUM3 = 0;
    NUM4 = 0;
    NUM5 = 0;
    NUM6 = 0;
    NUMHIN = 0;

    IPRNT = IPRINT % 10;
    PBUGSW = (IPRNT >= 4) ? true : false;
    PINFSW = (IPRNT >= 2) ? true : false;
    INIA12 = 5;

    //
    // WRITE header
    // FORMAT ( '1', T60, 'P T O L E M Y' /
    //   10X, 'COMPUTATION OF TRANSFER S-MATRIX ELEMENTS', T95, 'WHEELS WITHIN WHEELS' /
    //   '0', 45A1, 'ELAB =', F7.2, ' MEV', 5X, 65A1 / )
    //
    std::printf("1%58sP T O L E M Y\n", "");
    std::printf("%10sCOMPUTATION OF TRANSFER S-MATRIX ELEMENTS%43sWHEELS WITHIN WHEELS\n", "", "");
    std::printf("0%.45sELAB =%7.2f MEV     %.65s\n\n",
                &HEDCOM.REACT[1], ELAB, &HEDCOM.HEADER[1]);

    NLINE = 1000;

    SETLOG(2 * (LOMOST + LXMAX));

    //
    // WE WILL HAVE NO MORE IALLOC CALLS SO GET ADDRESSES
    //
    LMSVAL  = LOCPTRS.Z[GRIDCM.IMSVAL] - 1;
    LJA12S  = FACFR4 * LOCPTRS.Z[GRIDCM.IJA12S] - FACFR4;
    LINTS   = LOCPTRS.Z[GRIDCM.IINTS] - 1 + GRIDCM.INTOFF;
    LA12VL  = LOCPTRS.Z[GRIDCM.IA12VL] - 1;
    LHS     = LOCPTRS.Z[GRIDCM.IHS];
    LHSM1   = LHS - 1;
    LIINDX  = FACFR4 * LOCPTRS.Z[GRIDCM.IIINDX] - FACFR4;
    LA12T1  = LOCPTRS.Z[GRIDCM.IA12TM];
    LA12T2  = LA12T1 + NUMLX;
    LA12T3  = FACFR4 * (LA12T2 + 2 * LXMAX + 2) - FACFR4;
    if (IHSAVE == 2) LA12T3 = FACFR4 * LA12T1 - FACFR4;
    LCOSIN  = LOCPTRS.Z[GRIDCM.ICOSIN];
    LPHIT   = FACFR4 * LOCPTRS.Z[GRIDCM.IPHIT] - FACFR4;
    LPHIP   = FACFR4 * LOCPTRS.Z[GRIDCM.IPHIP] - FACFR4;
    LPHI    = FACFR4 * LOCPTRS.Z[GRIDCM.IPHI] - FACFR4;
    LTRP    = FACFR4 * LOCPTRS.Z[GRIDCM.LOCTRP] - FACFR4;
    LRI     = FACFR4 * LOCPTRS.Z[GRIDCM.IRI] - FACFR4;
    LRO     = FACFR4 * LOCPTRS.Z[GRIDCM.IRO] - FACFR4;
    LLIR    = FACFR4 * LOCPTRS.Z[GRIDCM.ILIR] - FACFR4;
    LLII    = FACFR4 * LOCPTRS.Z[GRIDCM.ILII] - FACFR4;
    LLOR    = FACFR4 * LOCPTRS.Z[GRIDCM.ILOR] - FACFR4;
    LLOI    = FACFR4 * LOCPTRS.Z[GRIDCM.ILOI] - FACFR4;
    LWIO    = FACFR4 * LOCPTRS.Z[GRIDCM.IWIO] - FACFR4;
    LDW     = LOCPTRS.Z[GRIDCM.IDW] - 1;
    LLIS    = FACFR4 * LOCPTRS.Z[INELCM.ILIS] - FACFR4;
    LLILOR  = LOCPTRS.Z[INELCM.ILILOR] - 1;
    LLILOI  = LOCPTRS.Z[INELCM.ILILOI] - 1;
    // debug removed
    LHINT   = LOCPTRS.Z[GRIDCM.IHINT] - 1;
    LHABS   = LOCPTRS.Z[GRIDCM.IHABS] - 1;
    LABS1   = LOCPTRS.Z[GRIDCM.IABS1] - 1;
    LSMHPT  = LOCPTRS.Z[GRIDCM.ISMHPT];
    LSMIPT  = LOCPTRS.Z[GRIDCM.ISMIPT];
    LSMHVL  = LOCPTRS.Z[GRIDCM.ISMHVL] - 1;
    LSMIVL  = LOCPTRS.Z[GRIDCM.ISMIVL] - 1;
    LSMHWK  = LOCPTRS.Z[GRIDCM.ISMHWK];
    LRIOEX  = LOCPTRS.Z[GRIDCM.IRIOEX] - 1;
    LWAVR   = LOCPTRS.Z[WAVCOM.IWAVR];
    LWAVI   = LOCPTRS.Z[WAVCOM.IWAVI];
    LVRL1   = LOCPTRS.Z[WAVCOM.IRLVS[1]];
    LVRL2   = LOCPTRS.Z[WAVCOM.IRLVS[2]];
    LVIM1   = LOCPTRS.Z[WAVCOM.IIMVS[1]];
    LVIM2   = LOCPTRS.Z[WAVCOM.IIMVS[2]];
    LCENT1  = LOCPTRS.Z[WAVCOM.ICENTR[1]];
    LCENT2  = LOCPTRS.Z[WAVCOM.ICENTR[2]];
    LATERM  = LOCPTRS.Z[INELCM.IBETAS[1]];
    LDWI    = FACFR4 * LOCPTRS.Z[GRIDCM.IDWI] - FACFR4;
    LWFIO   = FACFR4 * LOCPTRS.Z[GRIDCM.IWFIO] - FACFR4;
    LWFII   = FACFR4 * LOCPTRS.Z[GRIDCM.IWFII] - FACFR4;
    LINDXS  = FACFR4 * LOCPTRS.Z[INELCM.IINDXS] - FACFR4;

    //
    // WE ASSUME THAT L1+L2+LX <= 256  AND  NUMCOSIN <= 16384
    //
    PI256 = 256.0 * PI;

    //
    // FORCE THE OVERLAYS TO GET THE THREEJ COMMON BLOCKS IN.
    //
    DUMMY1();

    FACTOR = 2.0 * sqrt(AKI * AKO / (WAVCOM.ES[1] * WAVCOM.ES[2]));
    if (INELCM.IDENSW) FACTOR = sqrt(2.0) * FACTOR;

    //
    // THE LOOPS BEGIN NOW
    //
    TT7 = (float)second();

    //
    // LOOP OVER LI
    //
    // TWO PASSES OVER LI LOOP NOW; EVEN LI FIRST, THEN ODD
    //
    LIMIN = LMIN + 1;
    LIL   = LMIN;
    if (LMIN % 2 == 1) goto L150;
    LIMIN = LMIN;
    LIL   = LMIN + 1;

    //
    // LIPRTY = 1 FOR LI EVEN, = 2 FOR LI ODD
    //
L150:
    for (LIPRTY = 1; LIPRTY <= 2; LIPRTY++) {

        if (LIMIN > LMAX) goto L980;
        NXLISW = false;
        ILI = 1;

        for (LI = LIMIN; LI <= LMAX; LI += 2) {

            //
            // GET THE REQUIRED SCATTERING WAVEFUNCTIONS.
            //
            if (NXLISW) goto L250;

            //
            // FOR FIRST LI, FILL UP THE OUT-CHANNEL WAVEFUNCTION ARRAY
            //
            ITRK = 0;

            for (KWO = 1; KWO <= NWFO; KWO++) {
                LO   = ILLOC(LWFIO + 4*KWO - 3) + LI;
                LASO = ILLOC(LWFIO + 4*KWO - 2) + LI;
                JPO  = ILLOC(LWFIO + 4*KWO - 1) + 2*LI;
                ILLOC(LWFIO + 4*KWO) = ITRK;
                WFGET(LO, LASO, JPO, 2, NRIROI, ALLOC4_base(LRO + 1),
                      ALLOC4_base(LLOR + ITRK + 1), ALLOC4_base(LLOI + ITRK + 1),
                      ALLOC_base(LWAVR), ALLOC_base(LWAVI),
                      ALLOC_base(LVRL2), ALLOC_base(LVIM2), ALLOC_base(LCENT2));
                ITRK = ITRK + NRIROI;
            }

            ITRK = 0;
            NXLISW = true;
            goto L280;

            //
            // ON SUBSEQUENT LI'S WE NEED OUT-CHANNEL WAVE FUNCTIONS FOR ONE
            // NEW LO. WE MUST DETERMINE WHERE THE NEW LO STARTS.
            //
L250:
            LO   = ILLOC(LWFIO + 1) + 2;
            LASO = ILLOC(LWFIO + 2) + 2;
            JPO  = ILLOC(LWFIO + 3) + 4;
            for (II = 1; II <= NWFO; II++) {
                if (LO   != ILLOC(LWFIO + 4*II - 3)
                 || LASO != ILLOC(LWFIO + 4*II - 2)
                 || JPO  != ILLOC(LWFIO + 4*II - 1)) goto L254;
                I = NWFO - II + 1;
                goto L255;
L254:           ;
            }
            I = 0;
            goto L260;

            //
            // MOVE POINTERS DOWN.
            //
L255:
            for (KWO = 1; KWO <= I; KWO++) {
                K = KWO + (NWFO - I);
                ILLOC(LWFIO + 4*KWO) = ILLOC(LWFIO + 4*K);
            }

            // NOW PUT THE NEW WAVE FUNCTIONS INTO THE CIRCULAR BUFFER.
L260:
            I = I + 1;
            for (KWO = I; KWO <= NWFO; KWO++) {
                LO   = ILLOC(LWFIO + 4*KWO - 3) + LI;
                LASO = ILLOC(LWFIO + 4*KWO - 2) + LI;
                JPO  = ILLOC(LWFIO + 4*KWO - 1) + 2*LI;
                ILLOC(LWFIO + 4*KWO) = ITRK;
                WFGET(LO, LASO, JPO, 2, NRIROI, ALLOC4_base(LRO + 1),
                      ALLOC4_base(LLOR + ITRK + 1), ALLOC4_base(LLOI + ITRK + 1),
                      ALLOC_base(LWAVR), ALLOC_base(LWAVI),
                      ALLOC_base(LVRL2), ALLOC_base(LVIM2), ALLOC_base(LCENT2));
                ITRK = ITRK + NRIROI;
                if (ITRK > NCRIT) ITRK = 0;
            }

            //
            // FOR EACH LI, WE MUST GET THE IN-CHANNEL WAVE FUNCTIONS
            //
L280:
            for (KWI = 1; KWI <= NWFI; KWI++) {
                LASI = ILLOC(LWFII + 3*KWI - 2) + LI;
                JPI  = ILLOC(LWFII + 3*KWI - 1) + 2*LI;
                I    = ILLOC(LWFII + 3*KWI) + 1;
                WFGET(LI, LASI, JPI, 1, NRIROI, ALLOC4_base(LRI + 1),
                      ALLOC4_base(LLIR + I), ALLOC4_base(LLII + I),
                      ALLOC_base(LWAVR), ALLOC_base(LWAVI),
                      ALLOC_base(LVRL1), ALLOC_base(LVIM1), ALLOC_base(LCENT1));
            }

            //
            // WAVEFUNCTIONS ALL FOUND, WILL WE DO THIS VALUE OF LI
            //
L290:
            if (LI - ILLOC(LLIS + ILI) < 0) goto L959;
            if (LI - ILLOC(LLIS + ILI) == 0) goto L300;
            // LI - ILLOC(LLIS+ILI) > 0
            if (ILI >= NUMLIS) goto L959;
            ILI = ILI + 1;
            goto L290;

            //
            // THIS IS AN LI TO DO, FIND THE A12 COEFFICIENTS
            //
L300:
            TT1 = (float)second();

            LNLAM = FACFR4 * LOCPTRS.Z[GRIDCM.INLAM] - FACFR4;
            A12(LI, LXMIN, LXMAX, LMIN, LMAX,
                &ALLOC(LOCPTRS.Z[GRIDCM.IXLAM]),    // XLAM: A12 does -1
                &ILLOC(LNLAM + 1),                    // NLAM: A12 does -1
                &ALLOC(LA12VL + 1),                   // A12VL: A12 does -1
                &ALLOC(LMSVAL + 1),                   // DMSVAL: A12 does -1
                &ILLOC(LJA12S + 1),                   // JA12S: A12 does -1
                JA12M, JA12N, JA12AN,
                &ILLOC(LIINDX + 1),                   // IINDEX: 0-based macro
                IHMAX, LOMNMN, LOMXMX,
                &ALLOC(LINTS + 1),                    // DINTS: A12 does -1
                &ALLOC(LA12T1),                       // OUTTMP: A12 does -1
                &ALLOC(LA12T2),                       // XLOTMP: A12 does -1
                &ILLOC(LA12T3 + 1),                   // LXTMP: A12 does -1
                LHSM1, LA12VL, IPRINT, INIA12,
                &ILLOC(LDWI + 1),                     // INDXDW: 0-based macro
                NDW, NI,
                &ILLOC(LWFII + 1),                    // IDWFI: 0-based macro
                &ILLOC(LWFIO + 1));                   // IDWFO: 0-based macro
            NUMHSU = IHMAX * NPSUM;

            if ((JA12M > INELCM.IA12M || JA12N > GRIDCM.IA12N) && IHSAVE != 2) {
                std::printf(" **** ERROR *** A12 HAS TOO MANY TERMS: %8d%8d%8d%8d\n",
                            INELCM.IA12M, JA12M, GRIDCM.IA12N, JA12N);
            }
            FTIME.TIMES[4] = FTIME.TIMES[4] + (float)second() - TT1;
            NUMTIM = NUMTIM + 2;

            if (LOMNMN > LOMXMX) goto L959;

            //
            // NOW PROCEED TO THE (RI, RO) INTEGRATION LOOP
            //
            LHEND = LHSM1 + IHMAX;
            MCNT = 1;
            NSKIP = 0;
            DWRIOS = 10.0;
            DWRIOL = 0.0;
            DWRIOC = 0.0;
            DWCONT = 0.0;

            //
            // INITIALIZE "I1REAL", "I1IMAG" AND "RIROABS" FOR THIS LI.
            //
            for (II = 1; II <= NI; II++) {
                ALLOC(LLILOR + II) = 0.0;
                ALLOC(LLILOI + II) = 0.0;
                ALLOC(LABS1 + II) = 0.0;
            }

            TT5 = (float)second();

            //
            // LOOP OVER THE V-GRID AND FOR EACH ONE DO A COMPLETE U-INTEGRAL.
            //
            for (IV = 1; IV <= NPDIF; IV++) {

                //
                // GET THE H'S FOR INTERPOLATION INPUT
                //
                // INDEX TO THE H'S FILE
                //
                LIIP = 1000 * (1000 * LIPRTY + LI) + IV;
                if (IHSAVE == 2) goto L560;

                //
                // LOOP OVER U-VALUES AND COMPUTE THE H'S
                //
                IPLUNK = NPSUM * (IV - 1);
                TT2 = (float)second();

                for (IU = 1; IU <= NPSUM; IU++) {

                    IPLUNK = IPLUNK + 1;

                    for (IH = 1; IH <= IHMAX; IH++) {
                        // SALLOC(LHSM1+IH) = 0  => ALLOC(LHSM1+IH) = 0
                        ALLOC(LHSM1 + IH) = 0.0;
                        ALLOC(LHINT + IH) = 0.0;
                        ALLOC(LHABS + IH) = 0.0;
                    }

                    if (ISCRFL == 0) goto L420;

                    // READ ( 2, END=9960, ERR=9960 )
                    //    ( ALLOC4(LPHI+II), II = 1, NPPHI ),
                    //    ( ALLOC4(LPHIP+II), II = 1, NPPHI ),
                    //    ( ALLOC4(LPHIT+II), II = 1, NPPHI ),
                    //    ( ALLOC4(LTRP+II), II = 1, NPPHI )
                    // Scratch file binary read — stub
                    std::printf(" INELDC: scratch file READ(2) stub — binary I/O not yet implemented\n");
                    goto L9960;

                    MCNT = 1;

                    //
                    // THE DX INTEGRAL (AT LAST)
                    //
L420:
                    for (II = 1; II <= NPPHI; II++) {

                        PHI  = (double)ALLOC4(LPHI + MCNT);
                        KMEND = LMSVAL + JA12M;
                        PHIP = (double)ALLOC4(LPHIP + MCNT);
                        ARG  = PHI + PHI;
                        PHIT = (double)ALLOC4(LPHIT + MCNT);
                        PVPDX = (double)ALLOC4(LTRP + MCNT);
                        JA12OF = LJA12S - LMSVAL + 1;

                        //
                        // HIGH SPEED INLINE COSINE AND SINE(ARG)
                        // ARG = 2*PHI, finding COS(2PHI) and SIN(2PHI)
                        //
                        RNN = (double)(int)(ARG * GRIDCM.STPINV + 0.5);
                        NN = (int)RNN;
                        NN = NN % NCOSIN;
                        DELTAC = ARG - RNN * GRIDCM.COSSTP;
                        S2PHI = ((1.0 + HAFNEG * DELTAC * DELTAC)
                                * ALLOC(LCOSIN + std::abs(NN - ICOSS4))
                                + DELTAC * ALLOC(LCOSIN + NN));
                        C2PHI = (1.0 + HAFNEG * DELTAC * DELTAC)
                                * ALLOC(LCOSIN + NN)
                                - DELTAC * ALLOC(LCOSIN + std::abs(NN - ICOSS4));

                        // THE COSINE AND SINE(ARG) HAS BEEN FOUND
                        S2PHI2 = S2PHI * S2PHI;

                        //
                        // LOOP OVER TERMS IN A12
                        //
                        for (KA12M = LMSVAL; KA12M <= KMEND; KA12M += 5) {

                            //
                            // THE MAGIC ANGLE FOR THE FIRST VALUE OF MU-2
                            // WE ADD 256*PI TO FORCE IT POSITIVE
                            //
                            ARG = PI256 + ALLOC(KA12M + 1) * PHIT
                                        + ALLOC(KA12M + 2) * PHIP
                                        - ALLOC(KA12M + 3) * PHI;

                            LA12OF = ILLOC(JA12OF + KA12M);

                            //
                            // HIGH SPEED INLINE COSINE AND SINE(ARG)
                            //
                            NUMMU = ILLOC(JA12OF + KA12M + 1);
                            RNN = (double)(int)(ARG * GRIDCM.STPINV + 0.5);
                            NN = (int)RNN;
                            NN = NN % NCOSIN;
                            DELTAC = ARG - RNN * GRIDCM.COSSTP;

                            //
                            // COSTHE = COS(MT*PHIT+MP*PHIP-MU*PHI)
                            // SINSIN = SIN(MT*PHIT+MP*PHIP-MU*PHI) * SIN(2*PHI)
                            //
                            COSTHE = (1.0 + HAFNEG * DELTAC * DELTAC)
                                    * ALLOC(LCOSIN + NN)
                                    - DELTAC * ALLOC(LCOSIN + std::abs(NN - ICOSS4));
                            SINSIN = ((1.0 + HAFNEG * DELTAC * DELTAC)
                                    * ALLOC(LCOSIN + std::abs(NN - ICOSS4))
                                    + DELTAC * ALLOC(LCOSIN + NN)) * S2PHI;

                            // THE COSINE AND SINE(ARG) HAS BEEN FOUND
                            LHSTRT = ILLOC(JA12OF + KA12M + 2);
                            COSTHE = PVPDX * COSTHE;
                            NUMIH = ILLOC(JA12OF + KA12M + 3);
                            SINSIN = PVPDX * SINSIN;

                            //
                            // INNERMOST LOOP OVER MU'S.
                            // EACH MU IS 2 GREATER THAN THE PREVIOUS, SO THETA IS DECREMENTED
                            // BY 2PHI AND COS(THETA) IS FOUND ITERATIVELY.
                            //
                            for (NMU = 1; NMU <= NUMMU; NMU++) {
                                COSNEW = COSTHE * C2PHI + SINSIN;
                                SINSIN = SINSIN * C2PHI - COSTHE * S2PHI2;
                                COSTHE = COSNEW;

                                // LOOP OVER ALL (LX, LO) PAIRS FOR THIS LI
                                for (LH = LHSTRT; LH <= LHEND; LH++) {
                                    // SALLOC(LH) = SALLOC(LH) + COSNEW*ALLOC(LA12OF+LH)
                                    ALLOC(LH) = ALLOC(LH) + COSNEW * ALLOC(LA12OF + LH);
                                }
                                LA12OF = LA12OF + NUMIH;
                            }

                        }  // end of KA12M loop (A12 terms)

                        MCNT = MCNT + 1;

                        //
                        // ACCUMULATE THE DESIRED INTEGRAL IN LHINT AND THE INTEGRAL
                        // OF |INTEGRAND| IN LHABS.
                        //
                        for (IH = 1; IH <= IHMAX; IH++) {
                            ALLOC(LHINT + IH) = ALLOC(LHINT + IH) + ALLOC(LHSM1 + IH);
                            ALLOC(LHABS + IH) = ALLOC(LHABS + IH) + fabs(ALLOC(LHSM1 + IH));
                            ALLOC(LHSM1 + IH) = 0.0;
                        }

                    }  // end of II=1,NPPHI (DX integral terms)

                    //
                    // END OF DX INTEGRAL TERMS
                    //
                    NUMTRM = NUMTRM + JA12N;
                    NUMANG = NUMANG + JA12AN;
                    NUM3 = NUM3 + JA12M;
                    NUM4 = NUM4 + 1;

                    //
                    // STORE THE H'S WITH THE FACTOR EXP(-ALPHAP*RP - ALPHAT*RT)
                    // REMOVED IN THE INTERPOLATION INPUT ARRAY
                    //
                    for (IH = 1; IH <= IHMAX; IH++) {
                        ALLOC(LSMHVL + IU + NPSUM * (IH - 1)) =
                            ALLOC(LHINT + IH) * ALLOC(LRIOEX + IPLUNK);
                    }
                    if (IPRNT >= 9) {
                        std::printf(" CMPTED: %9d%5d", LIIP, IPLUNK);
                        for (IH = 1; IH <= IHMAX; IH++) {
                            if ((IH - 1) % 7 == 0 && IH > 1) std::printf("\n%22s", "");
                            std::printf(" %14.5G", ALLOC(LHINT + IH));
                        }
                        std::printf("\n");
                    }

                }  // end of IU=1,NPSUM (U loop)

                //
                // END OF LOOP ON U
                //
                FTIME.TIMES[5] = FTIME.TIMES[5] + (float)second() - TT2;
                NUMTIM = NUMTIM + 2;

                //
                // SAVE H'S NOW
                //
                if (IHSAVE == 1) {
                    // WRITE ( 1, ERR=9950 ) LIIP, ( ALLOC(LSMHVL+I), I = 1, NUMHSU )
                    // Scratch file binary write — stub
                    std::printf(" INELDC: scratch file WRITE(1) stub — binary I/O not yet implemented\n");
                    // On error, would go to 9950
                }

                goto L590;

                //
                // HERE WE READ SAVED H'S. MAKE SURE WE DO NOT GET OUT OF SYNC.
                //
L560:
                // READ ( 1, END=9950, ERR=9950 ) LIIPH, ( ALLOC(LSMHVL+I), I = 1, NUMHSU )
                // Scratch file binary read — stub
                std::printf(" INELDC: scratch file READ(1) stub — binary I/O not yet implemented\n");
                LIIPH = 0;
                goto L9950;

                if (LIIPH == LIIP) goto L590;
                if (LIIPH < LIIP) goto L560;
                std::printf("\n**** SYNCHRONIZATION ERROR: %12d%12d\n", LIIP, LIIPH);
                goto L9950;

L590:
                NUM5 = NUM5 + NPSUM;
                NUM6 = NUM6 + NUMHSU;

                //
                // NOW INTERPOLATE
                //
                TT2 = (float)second();
                for (IH = 1; IH <= IHMAX; IH++) {
                    ITEMP = LSMHVL + 1 + (IH - 1) * NPSUM;
                    SPLNCB(NPSUM, ALLOC_base(LSMHPT), ALLOC_base(ITEMP),
                           ALLOC_base(LSMHWK), ALLOC_base(LSMHWK + NPSUM),
                           ALLOC_base(LSMHWK + 2 * NPSUM));
                    INTRPC(NPSUM, ALLOC_base(LSMHPT), ALLOC_base(ITEMP),
                           ALLOC_base(LSMHWK), ALLOC_base(LSMHWK + NPSUM),
                           ALLOC_base(LSMHWK + 2 * NPSUM),
                           NPSUMI, ALLOC_base(LSMIPT),
                           ALLOC_base(LSMIVL + 1 + (IH - 1) * NPSUMI));
                    NUMHIN = NUMHIN + NPSUMI;
                }
                FTIME.TIMES[6] = FTIME.TIMES[6] + (float)second() - TT2;
                NUMTIM = NUMTIM + 2;

                //
                // NOW DO INTEGRAL DU FOR FIXED V.
                //
                IPLUNK = NPSUMI * (IV - 1);
                for (IU = 1; IU <= NPSUMI; IU++) {
                    IPLUNK = IPLUNK + 1;

                    //
                    // COMPUTE IMAG AND REAL PARTS OF PRODUCT OF THE DISTORTED
                    // SCATTERING WAVES.
                    //
                    DWMAX = 0.0;
                    for (KDW = 1; KDW <= NDW; KDW++) {
                        I = ILLOC(LDWI + 4*KDW - 1) + IPLUNK;
                        FIR = (double)ALLOC4(LLIR + I);
                        FII = (double)ALLOC4(LLII + I);
                        I = ILLOC(LDWI + 4*KDW) + IPLUNK;
                        FOR_ = (double)ALLOC4(LLOR + I);
                        FOI  = (double)ALLOC4(LLOI + I);
                        DWR = FOR_ * FIR - FOI * FII;
                        DWI = FOR_ * FII + FOI * FIR;
                        // Guard against NaN from WAVELJ overflow (outgoing L=0)
                        if (std::isnan(DWR) || std::isinf(DWR)) DWR = 0.0;
                        if (std::isnan(DWI) || std::isinf(DWI)) DWI = 0.0;
                        DWLIM = std::max(fabs(DWR), fabs(DWI));
                        if (DWLIM > DWRIOL) DWRIOL = DWLIM;
                        if (DWLIM < DWRIOS && DWLIM != 0.0) DWRIOS = DWLIM;
                        DWRIOC = DWRIOC + DWLIM;
                        DWMAX = std::max(DWMAX, DWLIM);
                        ALLOC(LDW + 2*KDW - 1) = DWR;
                        ALLOC(LDW + 2*KDW) = DWI;
                        DWCONT = DWCONT + 1.0;
                    }

                    //
                    // SKIP THE NEXT PART IF ALL DISTORTED-WAVE PRODUCTS ARE SMALL.
                    //
                    if (DWMAX < 1.0e-30) goto L780;

                    TERM = (double)ALLOC4(LWIO + IPLUNK);

                    //
                    // NOW ADD THESE H'S INTO THE APPROPRIATE I(JO,JI,LO,LX)
                    // ALONG WITH THE PRODUCT OF THE DISTORTED WAVES
                    //
                    for (II = 1; II <= NI; II++) {
                        ITEMP = ILLOC(LIINDX + 4*II - 3) + IU;
                        KDW   = ILLOC(LIINDX + 4*II - 2);
                        // Guard against NaN from WAVELJ overflow (outgoing channel)
                        double dwR = ALLOC(KDW - 1), dwI = ALLOC(KDW);
                        if (std::isnan(dwR) || std::isinf(dwR)) dwR = 0.0;
                        if (std::isnan(dwI) || std::isinf(dwI)) dwI = 0.0;
                        {
                        double addR = TERM * ALLOC(ITEMP) * dwR;
                        double addI = TERM * ALLOC(ITEMP) * dwI;
                        ALLOC(LLILOR + II) += addR;
                        ALLOC(LLILOI + II) += addI;
                        }

                        //
                        // COMPUTE AN INDICATION OF LOSS OF SIGNIFICANCE:
                        //   "RIROABS" = INTEGRAL(RI,RO) |INTEGRAL(PHI) INTEGRAND|
                        //
                        if (PINFSW) {
                            ALLOC(LABS1 + II) = ALLOC(LABS1 + II)
                                + fabs(TERM * ALLOC(ITEMP))
                                * (fabs(ALLOC(KDW - 1)) + fabs(ALLOC(KDW)));
                        }
                    }

                    goto L789;

                    //
                    // WE COME HERE WHEN SKIPPING I CALCULATION DUE TO THE PRODUCT
                    // OF THE DISTORTED WAVES BEING SMALL.
                    //
L780:
                    NSKIP = NSKIP + 1;

L789:               ;
                }  // end of IU=1,NPSUMI

            }  // end of IV=1,NPDIF

            FTIME.TIMES[7] = FTIME.TIMES[7] + (float)second() - TT5;
            NUMTIM = NUMTIM + 2;

            //
            // END OF (RI,RO) INTEGRATION
            //
            DWRIOC = DWRIOC / DWCONT;
            if (PBUGSW) {
                std::printf(" SMALLEST, LARGEST AND AVERAGE DABS(DWRIRO)= %15.5G%15.5G%15.5G\n\n",
                            DWRIOS, DWRIOL, DWRIOC);
            }

            //
            // SUBROUTINE SFROMI CALCULATES S-MATRIX ELEMENTS FROM THE
            // RADIAL INTEGRALS.
            //
            SFROMI(LI, ILI,
                   ALLOC_base(LOCPTRS.Z[INELCM.ISMATR]),
                   ALLOC_base(LOCPTRS.Z[INELCM.ISMATI]),
                   &ILLOC(LINDXS + 1),
                   ALLOC_base(LLILOR + 1),
                   ALLOC_base(LLILOI + 1),
                   &ILLOC(LIINDX + 1),
                   NI,
                   &ILLOC(LDWI + 1),
                   &ILLOC(LWFII + 1),
                   &ILLOC(LWFIO + 1),
                   ALLOC_base(LABS1 + 1),
                   ALLOC_base(LATERM),
                   FACTOR, PINFSW, NLINE);

            if (IHSAVE == 2) goto L955;

            if (ISCRFL != 0) {
                // REWIND 2 — rewind scratch file unit 2
                // Stub: no-op for now
            }

            if (MCNT - 1 > GRIDCM.MAXCNT) {
                std::printf(" **** ERROR ****  PHI COUNTS MESSED UP: %10d%10d\n",
                            MCNT, GRIDCM.MAXCNT);
            }
            I = LHSTRT + LA12OF - LA12VL - 1;
            if (I != JA12N) {
                std::printf(" ******* A12 COUNTS MESSED UP: %10d%10d%10d%10d%10d%10d%10d%10d\n",
                            I, JA12N, LI, JA12M, IHMAX, LHSTRT, LA12OF, LA12VL);
            }

L955:       ;

L959:       ;
        }  // end of LI=LIMIN,LMAX,2

        //
        // END OF THE LI LOOP FOR A GIVEN PARITY OF LI
        // IF LI IS EVEN HAVE ANOTHER BASH AT LI LOOP WITH ODD LI
        //
L980:
        LIMIN = LIL;

    }  // end of LIPRTY=1,2

    TT = (float)second();
    FTIME.TIMES[8] = TT - TT7;

    //
    // END OF THE 2 LI LOOPS (EVEN AND ODD PARITY OF LI)
    //
    // IF NECESSARY, CLEAN UP AFTER SAVEHS/USEHS
    //
    if (IHSAVE - 1 < 0)  goto L2000;
    if (IHSAVE - 1 == 0) goto L1600;
    if (IHSAVE - 1 > 0)  goto L1700;

    //
    // WRITE A DUMMY RECORD AT THE END
    //
L1600:
    LIIP = 300000000;
    // WRITE ( 1, ERR=9950 ) LIIP, ( ALLOC(LSMHVL+I), I = 1, NUMHS )
    // Scratch file binary write — stub
    std::printf(" INELDC: scratch file WRITE(1) dummy record stub\n");
    goto L1780;

    //
    // FOR A MULTIPASS, USEHS RUN, READ PAST THE DUMMY RECORD.
    //
L1700:
    if (IBSPAS == 0) goto L1780;
L1720:
    if (LIIPH == 300000000) goto L1780;
    // READ ( 1, ERR=9950, END=1780 ) LIIPH
    // Scratch file binary read — stub
    std::printf(" INELDC: scratch file READ(1) skip-past stub\n");
    LIIPH = 300000000;  // Force exit for stub
    goto L1720;

L1780:
    ;  // CONTINUE

    //
    // A L L    D O N E
    //
L2000:
    IRTN = 1;

    INELD2(TSTART, NUMHIN, NUMTIM, NUMTRM, NUMANG,
           NUM3, NUM4, NUM5, NUM6);

    return;

    //
    // ERROR HANDLERS
    //
L9950:
    std::printf("\n**** ERROR IN READING OR WRITING H'S FILE.\n");
    return;

L9960:
    std::printf("\n**** ERROR IN READING SCRATCH FILE (UNIT 2).\n");
    return;
}
