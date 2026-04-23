// angset_part1.cpp — Verbatim C++ translation of ANGSET (first half)
// Translated from source.f lines 916-1545
//
// SETS UP ARRAYS FOR THE ANGULAR TRANSFORMS IN THE RADIAL INTS.
//
// THIS ROUTINE IS CALLED AFTER THE WAVPOT CALLS AND BEFORE INELDC.
// IT ALLOCATES ARRAYS USED BY A12 AND OTHER ARRAYS USED BY
// INELDC.  IT MUST BE CALLED AFTER WAVPOT TO ALLOW THE COULOMB
// ROUTINE IN WAVPOT TO USE LOTS OF ALLOCATOR AS TEMPORARY SCRATCH.
//
// IRTN IS SET TO 0 IF AN ERROR OCCURS, 1 OTHERWISE.

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include "ptolemy_alloc.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_forward.h"
#include "ptolemy_io.h"
#include <cstdio>
#include <cstring>
#include <cmath>

void ANGSET(int& IRTN)
{
    // =====================================================================
    // Local variables
    // =====================================================================
    double TSTART;
    int IPRNT;
    int LOC;
    int N, I, K, LI, LIH;
    int ILI, ILIH;
    int NUMLIH, ITEMP;
    int LLIS;
    int LCHNDF, NCHNDF, MCHNDF;
    int LINDXS, LTOCS;
    int INEED;
    int LAS, LASMIN, LASMAX;
    int JPMN, JPMX, JPI, JPO;
    int LWFII, LWFIO;
    int LOMAX, LOMIN, LASR;
    int LO;

    logical PBUGSW, PINFSW, COULSW, MULTSW, COUPSW;

    // CHARACTER*8 DATA initializations
    char8 NAMER[3];  // 1-based
    NAMER[1] = char8("WAVSAVIR");
    NAMER[2] = char8("WAVSAVOR");
    char8 NAMEI[3];  // 1-based
    NAMEI[1] = char8("WAVSAVII");
    NAMEI[2] = char8("WAVSAVOI");

    int NUMFRE = 21;
    char8 FRENAM[22];  // 1-based
    FRENAM[1]  = char8("TORUTHIN");
    FRENAM[2]  = char8("SMAG    ");
    FRENAM[3]  = char8("SPHASE  ");
    FRENAM[4]  = char8("SINMAG  ");
    FRENAM[5]  = char8("SINPHASE");
    FRENAM[6]  = char8("SOUTMAG ");
    FRENAM[7]  = char8("SOUTPHAS");
    FRENAM[8]  = char8("TORUTHOT");
    FRENAM[9]  = char8("CROSSSEC");
    FRENAM[10] = char8("LXCROSSS");
    FRENAM[11] = char8("LXTOTS  ");
    FRENAM[12] = char8("MXTOTS  ");
    FRENAM[13] = char8("BETAS   ");
    FRENAM[14] = char8("FIN     ");
    FRENAM[15] = char8("FOUT    ");
    FRENAM[16] = char8("F       ");
    FRENAM[17] = char8("LABCROSS");
    FRENAM[18] = char8("MXCROSSS");
    FRENAM[19] = char8("ANPOW   ");
    FRENAM[20] = char8("ANPOWA  ");
    FRENAM[21] = char8("ANPOWB  ");

    double DUMMY[2];  // DIMENSION DUMMY(1)
    int IDUMMY[2][2]; // DIMENSION IDUMMY(1,1)

    // =====================================================================
    // COMMON block aliases (dot access)
    // =====================================================================
    auto& FACFR4    = ALLOCS.FACFR4;
    auto* Z         = LOCPTRS.Z;

    auto& IPRINT    = INTGER.IPRINT;
    auto& LMAX      = INTGER.LMAX;
    auto& LMIN      = INTGER.LMIN;
    auto& LSTEP     = INTGER.LSTEP;
    auto& LBACK     = INTGER.LBACK;
    auto& LX        = INTGER.LX;

    auto& PROBLM    = SWITCH.PROBLM;
    auto& IBSPAS    = SWITCH.IBSPAS;

    auto& IHSAVE    = INTRNL.IHSAVE;
    auto& ISTRIP    = INTRNL.ISTRIP;
    auto& NOTDEF    = INTRNL.NOTDEF;

    auto& LBP       = FORMF.LBP;
    auto& LBT       = FORMF.LBT;

    auto& IINDXS    = INELCM.IINDXS;
    auto& ISMATR    = INELCM.ISMATR;
    auto& ISMATI    = INELCM.ISMATI;
    auto& NSMAT     = INELCM.NSMAT;
    auto& NLX       = INELCM.NLX;
    auto& NSPLI     = INELCM.NSPLI;
    auto& ITOCS     = INELCM.ITOCS;
    auto& NASPLI    = INELCM.NASPLI;
    auto& LXMIN     = INELCM.LXMIN;
    auto& LXMAX     = INELCM.LXMAX;
    auto& MSTOP     = INELCM.MSTOP;
    auto& NMLOLX    = INELCM.NMLOLX;
    auto& NUMLX     = INELCM.NUMLX;
    auto& NUMLIS    = INELCM.NUMLIS;
    auto& ILIS      = INELCM.ILIS;
    auto& LSKIP     = INELCM.LSKIP;
    auto& LIFIT     = INELCM.LIFIT;
    auto& ILIFIT    = INELCM.ILIFIT;
    auto& NOLFIT    = INELCM.NOLFIT;
    auto& IELCUP    = INELCM.IELCUP;
    auto& IDELSR    = INELCM.IDELSR;
    auto& IDELSI    = INELCM.IDELSI;
    auto& IA12M     = INELCM.IA12M;
    auto& IRDINT    = INELCM.IRDINT;

    auto& JPMIN     = INELCM.JPMIN;
    auto& JPMAX     = INELCM.JPMAX;
    auto& JPBASE    = INELCM.JPBASE;
    auto& NJP       = INELCM.NJP;
    auto& JTMIN     = INELCM.JTMIN;
    auto& JTMAX     = INELCM.JTMAX;
    auto& JTBASE    = INELCM.JTBASE;
    auto& NJT       = INELCM.NJT;

    auto& ISIZE     = INTGER.ISIZE;

    auto& JACOB     = GRIDCM.JACOB;
    auto& IHS       = GRIDCM.IHS;
    auto& IHINT     = GRIDCM.IHINT;
    auto& IHABS     = GRIDCM.IHABS;
    auto& INLAM     = GRIDCM.INLAM;
    auto& IXLAM     = GRIDCM.IXLAM;
    auto& IA12VL    = GRIDCM.IA12VL;
    auto& IMSVAL    = GRIDCM.IMSVAL;
    auto& IJA12S    = GRIDCM.IJA12S;
    auto& IA12TM    = GRIDCM.IA12TM;
    auto& IA12N     = GRIDCM.IA12N;
    auto& IDW       = GRIDCM.IDW;
    auto& IDWI      = GRIDCM.IDWI;
    auto& NOFLO     = GRIDCM.NOFLO;
    auto& NRIROI    = GRIDCM.NRIROI;
    auto& NCRIT     = GRIDCM.NCRIT;
    auto& ILIR      = GRIDCM.ILIR;
    auto& ILII      = GRIDCM.ILII;
    auto& ILOR      = GRIDCM.ILOR;
    auto& ILOI      = GRIDCM.ILOI;
    auto& IWFII     = GRIDCM.IWFII;
    auto& IWFIO     = GRIDCM.IWFIO;
    auto& NWFI      = GRIDCM.NWFI;
    auto& NWFO      = GRIDCM.NWFO;

    auto& ICHNDF    = CCBLK.ICHNDF;
    auto& IBASDF    = CCBLK.IBASDF;

    auto* TCSWS     = WAVCOM.TCSWS;
    auto* SOSWS     = WAVCOM.SOSWS;
    auto* NUMJS     = WAVCOM.NUMJS;
    auto* JSPS      = WAVCOM.JSPS;
    auto* IWAVRS    = WAVCOM.IWAVRS;
    auto* IWAVIS    = WAVCOM.IWAVIS;

    // =====================================================================
    // BEGIN EXECUTABLE CODE
    // =====================================================================
    IRTN = 0;
//
    IPRNT = MOD(IPRINT, 10);
    PBUGSW = IPRNT >= 3;
    PINFSW = IPRNT >= 2;
//
    COUPSW = PROBLM == 24;
//
//     ON PASS>=2, FREE AREAS NO LONGER NEEDED.
//
    if (IBSPAS <= 1) goto L300;
    for (K = 1; K <= NUMFRE; K++) {
        I = NAMLOC(FRENAM[K].data);
        if (I > 0) Z[I] = -Z[I];
    } // 259
//
//
//     WE DO NOT NECESSARILY COMPUTE EVERY I( LI, LO, LX ).  RATHER
//     WE FIND ONLY THOSE FOR
//        LI = LMIN, LMIN+LSTEP, LMIN+2LSTEP, ..., LMAX
//     AND
//        LI = LMIN, LMIN+1, ..., LXMAX
//     WHERE THE SERIES TERMINATES AT OR BEFORE LMAX.  ALL LO AND LX
//     FOR THE ABOVE LI'S ARE FOUND.  THE ARRAY "LIS" CONTAINS THOSE
//     LI'S FORWHICH THE CALCULATION IS TO BE DONE
//     NOTE THAT IF LMIN < LXMAX, WE EFFECTIVELY START OUT WITH
//     LSTEP=1 TO POVIDE A SMOOTH STARTING POINT FOR INTERPOLATION
//
//**************************************************
//
//     NOTE:  FOR COUPLED CHANNELS CALCULATIONS THESE ARE
//            REALLY THE TOTAL J VALUES.  HOWEVER WE DO
//            NOT YET DOUBLE THEM.
//
//**************************************************
//
//
L300:
    if (LSTEP != 1 && TCSWS[1])
        std::printf("\n***** WARNING:  LSTEP MUST = 1 FOR TR OR TP"
                    " POTENTIALS IN INCIDENT CHANNEL.\n");
    if (TCSWS[1]) LSTEP = 1;
//
//     COMPUTE THE NUMBER OF LI'S INCLUDING THE LI<LXMAX PART
//
    NUMLIS = (LMAX - LMIN + LSTEP) / LSTEP;
    N = 0;
    if (LXMAX > LMIN) N = (LXMAX - LMIN) / LSKIP
        - (LXMAX - LMIN) / LSTEP;
    NUMLIS = NUMLIS + N;
    ILIS = NALLOC("LIS     ", (NUMLIS + 1) / FACFR4);
    for (I = 1; I <= NUMLIS; I++) {
        int L = LMIN + (I - 1) * LSKIP;
        if (L > LXMAX) L = (I - N - 1) * LSTEP + LMIN;
        ILLOC(FACFR4 * Z[ILIS] - FACFR4 + I) = L;
    } // 309
//
    if (IHSAVE == 1) {
        // WRITE (1, ERR=9950) NUMLIS
        // WRITE (1, ERR=9950) ( ILLOC(FACFR4*Z(ILIS)-FACFR4+I), I=1,NUMLIS )
        // TODO: unformatted binary I/O to unit 1 (SAVEHS)
    }
    if (IHSAVE != 2) goto L400;
//
//     RECOVER LI'S FROM THE H'S FILE AND MAKE SURE THE PRESENT SET
//     IS A SUBSET.
//
    // READ (1, ERR=9950, END=9950) NUMLIH
    // TODO: unformatted binary I/O from unit 1 (USEHS)
    NUMLIH = 0; // placeholder
    ITEMP = IALLOC((NUMLIH + 1) / FACFR4);
    // READ (1, ERR=9950, END=9950) ( ILLOC(FACFR4*Z(ITEMP)-FACFR4+I), I=1,NUMLIH )
    // TODO: unformatted binary I/O from unit 1 (USEHS)
    ILIH = 1;
    for (ILI = 1; ILI <= NUMLIS; ILI++) { // 359
        LI = ILLOC(FACFR4 * Z[ILIS] - FACFR4 + ILI);
L330:
        LIH = ILLOC(FACFR4 * Z[ITEMP] - FACFR4 + ILIH);
        if (LI - LIH < 0) goto L380;
        if (LI - LIH == 0) continue; // goto 359
        // LI - LIH > 0:
        // L340:
        ILIH = ILIH + 1;
        if (ILIH <= NUMLIH) goto L330;
        goto L380;
    } // 359
    Z[ITEMP] = -Z[ITEMP];
    goto L400;
L380:
    std::printf("\n***** ERROR: THIS SPECIFICATION OF LMIN, LSTEP,"
                " LMAX (%5d%5d%5d ) REQUIRES LI'S THAT WERE NOT SAVED"
                " IN THE SAVEHS RUN.\n", LMIN, LSTEP, LMAX);
    std::printf("\n***** REQUIRED LI'S -\n               ");
    for (ILI = 1; ILI <= NUMLIS; ILI++)
        std::printf("%5d", ILLOC(FACFR4 * Z[ILIS] - FACFR4 + ILI));
    std::printf("\n");
    std::printf("\n***** AVAILABLE LI'S -\n               ");
    for (ILI = 1; ILI <= NUMLIH; ILI++)
        std::printf("%5d", ILLOC(FACFR4 * Z[ITEMP] - FACFR4 + ILI));
    std::printf("\n");
    IRTN = 0;
    return;
//
//     DETERMINE THE RANGE OF LI'S TO BE USED FOR EXTRAPOLATION.
//
L400:
    LLIS = FACFR4 * Z[ILIS] - FACFR4;
    if (LBACK == NOTDEF_INT)
        LBACK = LMAX - ILLOC(LLIS + NUMLIS - 3);
    LIFIT = LMAX - LBACK;
    for (ILI = 1; ILI <= NUMLIS; ILI++) { // 429
        if (LIFIT < ILLOC(LLIS + ILI)) goto L430;
    } // 429
L430:
    ILIFIT = ILI - 1;
    LIFIT = ILLOC(LLIS + ILIFIT);
    NOLFIT = NUMLIS - ILIFIT + 1;
//
//
//     "INDXS" WILL CONTAIN ALL INFORMATION NECESSARY TO FIND A
//     GIVEN S-MATRIX ELEMENT.  IT HAS 3 VALUES FOR EACH
//     (LX, JP, JT) TRIPLE, INDEXED BY
//        K = 1 + LX + NLX*( JP-JPBASE + NJP*(JT-JTBASE) )/2
//
//        (3*K-2) = KOFFS = OFFSET INTO S-MATRIX (=0 IF NONEXISTENT).
//        (3*K-1) = LDELMN = MINIMUM (LO-LI).
//        (3*K  ) = NDEL = NUMBER OF (LO-LI) VALUES.
//
//
//     "TOCS" WILL CONTAIN A TABLE OF CONTENTS FOR THE S MATRIX.
//     IT HAS 4 VALUES IDENTIFYING EACH S-MATRIX ELEMENT USED, INDEXED
//     BY ITS POSITION (KOFFS + (LO-LI+LDELMN)/2 ) IN THE S MATRIX:
//
//        (4*KOFFS-3) = LO - LI
//        (4*KOFFS-2) = LX
//        (4*KOFFS-1) = JP
//        (4*KOFFS  ) = JT
//
//
//     "SMATR", "SMATI" WILL HOLD THE CALCULATED (NOT INTERPOLATED)
//     REACTION S MATRIX.  FOR A GIVEN (LO, LX, JP, JT, LI), THE
//     S MATRIX IS INDEXED AS
//
//        (LI-LMIN)*NSPLI/LSTEP + KOFFS + ( LO - LI - LDELMN )/2
//
//     WHERE KOFFS AND LDELMN ARE TAKED FROM "INDXS" AS DESCRIBED ABOVE.
//
//
//     FIRST PASS, ALLOCATE "INDXS" AND "TOCS" AND INITIALIZE POINTERS.
//
    if (IBSPAS > 1) goto L1450;
//
//     FIRST FIND INDXS SIZE
//
    ISIZE = 0;
    if (COUPSW) goto L1410;
    SETSPT(ISIZE, LOC, 1, (int)JBLOCK.JS[4], (int)JBLOCK.JS[2], PBUGSW,
           &IDUMMY[0][0], 1, 1);
    goto L1420;
//
L1410:
    LCHNDF = FACFR4 * (Z[ICHNDF] - 1) + 1;
    NCHNDF = ILLOC(LCHNDF);
    MCHNDF = ILLOC(LCHNDF + 1);
    LCHNDF = LCHNDF + ILLOC(LCHNDF + 2);
    for (I = 2; I <= NCHNDF; I++) { // 1419
        LCHNDF = LCHNDF + MCHNDF;
        ILLOC(LCHNDF + 11) = ISIZE;
        SETSPT(ISIZE, LOC, 1,
               IABS(ILLOC(LCHNDF + 1)), IABS(ILLOC(LCHNDF + 3)), PBUGSW,
               &IDUMMY[0][0], 1, 1);
    } // 1419
//
//     DEFINE THE INDXS ARRAY
//
L1420:
    IINDXS = NALLOC("INDXS   ", (ISIZE + 1) / FACFR4);
    INIT4i(IINDXS, 0);
    LINDXS = Z[IINDXS] * FACFR4 - FACFR4;
//     NASPLI IS NSPLI SUMMED OVER ALL CHANNELS
    NASPLI = 0;
//
//     NOW FILL IN "INDXS" FOR ALL S-MATRIX ELEMENTS NEEDED ON THIS
//     AND SUBSEQUENT PASSES.
//
    if (COUPSW) goto L1430;
    SETSPT(NASPLI, LINDXS, 2, (int)JBLOCK.JS[4], (int)JBLOCK.JS[2], PBUGSW,
           &IDUMMY[0][0], 1, 1);
    goto L1440;
//
L1430:
    LCHNDF = FACFR4 * (Z[ICHNDF] - 1) + 1;
    LCHNDF = LCHNDF + ILLOC(LCHNDF + 2);
    {
        int LBASDF = FACFR4 * (Z[IBASDF] - 1) + 1;
        int MBASDF = ILLOC(LBASDF + 1);
        LBASDF = LBASDF + ILLOC(LBASDF + 2);
        for (I = 2; I <= NCHNDF; I++) { // 1439
            LCHNDF = LCHNDF + MCHNDF;
            int IB = ILLOC(LCHNDF + 7);
            int NBASDF = ILLOC(LCHNDF + 8);
            ILLOC(LCHNDF + 9) = NUMLIS * NASPLI;
            ILLOC(LCHNDF + 12) = 4 * NASPLI;
            SETSPT(NASPLI, LINDXS + ILLOC(LCHNDF + 11), 2,
                   IABS(ILLOC(LCHNDF + 1)), IABS(ILLOC(LCHNDF + 3)), PBUGSW,
                   &ILLOC(LBASDF + MBASDF * (IB - 1)), MBASDF, NBASDF);
            ILLOC(LCHNDF + 10) = NSPLI;
        } // 1439
    }
//
//     ALLOCATE AND INITIALIZE THE S-MATRIX AND ITS TABLE OF CONTENTS.
//     DELTASR, DELTASI  ARE FOR S - S(BORN).  WE START OUT WITH
//     THEM DEFINED AS THE FULL S-MATRIX ARRAYS.
//
L1440:
    NSMAT = NASPLI * NUMLIS;
    ISMATR = NALLOC("SMATR   ", NSMAT);
    ISMATI = NALLOC("SMATI   ", NSMAT);
    IDELSR = ISMATR;
    IDELSI = ISMATI;
    ITOCS = NALLOC("TOCS    ", 4 * NASPLI / FACFR4);
    INIT8(ISMATR, 0.0);
    INIT8(ISMATI, 0.0);
    INIT4i(ITOCS, -1);
    if (PROBLM >= 23) IELCUP = NALLOC("ELCOUPS ", 2 * NUMLIS);
//
//
//     ALL PASSES:  FILL IN "TOCS" FOR THOSE S-MATRIX ELEMENTS WHICH
//     ARE USED ON THIS PASS.
//
L1450:
    LINDXS = Z[IINDXS] * FACFR4 - FACFR4;
    LTOCS = Z[ITOCS] * FACFR4 - FACFR4;
    if (COUPSW) goto L1460;
    SETSPT(LTOCS, LINDXS, 3, (int)JBLOCK.JS[4], (int)JBLOCK.JS[2], PBUGSW,
           &IDUMMY[0][0], 1, 1);
//
//     THIS IS THE NUMBER OF (LO', LX') PAIRS PER LI IN THE INNER
//     INTEGRATION LOOPS.  LX'S DUE ONLY TO SPIN-ORBIT DISTORTIONS
//     ARE NOT INCLUDED.
//
    NMLOLX = ((LXMAX + LXMIN + 1) * (LXMAX - LXMIN + 1)
              + MOD(LXMIN + LBP + LBT + 1, 2)) / 2;
    goto L1500;
//
L1460:
    LCHNDF = FACFR4 * (Z[ICHNDF] - 1) + 1;
    LCHNDF = LCHNDF + ILLOC(LCHNDF + 2);
    for (I = 2; I <= NCHNDF; I++) { // 1469
        LCHNDF = LCHNDF + MCHNDF;
        { int SETSPT_arg1 = LTOCS + ILLOC(LCHNDF + 12);
        int SETSPT_arg2 = LINDXS + ILLOC(LCHNDF + 11);
        SETSPT(SETSPT_arg1, SETSPT_arg2, 3,
               IABS(ILLOC(LCHNDF + 1)), IABS(ILLOC(LCHNDF + 3)), PBUGSW,
               &IDUMMY[0][0], 1, 1);
        } } // 1469
//
//     FOR C.C. WE USE NMLOLX = NASPLI
//
    NMLOLX = NASPLI;
//
//
//     STORAGE FOR THE H'S (TRANSFER ONLY)
//
L1500:
    IHS = 0;
    IHINT = 0;
    IHABS = 0;
    INLAM = 0;
    IXLAM = 0;
    IA12VL = 0;
    IMSVAL = 0;
    IJA12S = 0;
    IA12TM = 0;
    if (ISTRIP == 0) goto L1800;
    IA12TM = (3 * NUMLX + 1) / FACFR4;
    if (IHSAVE == 2) goto L1750;
//
    IHS = NALLOC("H       ", NMLOLX);
    IHINT = NALLOC("HINTEGRL", NMLOLX);
    IHABS = NALLOC("HABSINT ", NMLOLX);
//
//     ALLOCATE ARRAYS FOR  A12
//
    MSTOP = MAX0(LXMAX + LBP, LXMAX + LBT, LBP + LBT + 1, 2 * LXMAX);
    MSTOP = 2 * (MSTOP / 2);
    INEED = (MSTOP / 2 + 1) * (KANDM.LOMOST + 1);
    INLAM = NALLOC("NLAMBDA ", (KANDM.LOMOST + 1) / FACFR4 + 1);
    IXLAM = NALLOC("XLAMBDA ", INEED);
//
//     FOLLOWING IS AN UPPER LIMIT ON THE NUMBER OF TERMS IN THE
//     A12 SUMS.
//     FIRST THE MAXIMUM NUMBER OF M1 AND M2
//     COMBINATIONS.  WE USE
//        M1  <= (LBP, LX+LBT)    M2  <= (LBT, LX+LBP)
//        MX  <= (LX,LBP+LBT)   M1+LBP AND M2+LBT  EVEN
//
    IA12M = (MIN0(LXMAX, LBP + LBT) + 1) * (MIN0(LBP, LBT) + 1);
//
//     THEN FOR EACH OF THE ABOVE MU TAKES ON VALUES SUCH THAT
//       0 <= MU <= LO <= LOMOST,    MU+LO EVEN
//     AND WE HAVE  NMLOLX  SETS OF SUCH A12 VALUES.
//
    IA12N = IA12M * (KANDM.LOMOST / 2 + 1) * NMLOLX;
//
//
//     FOR EACH (M1, M2) PAIR WE HAVE UP TO FIVE ITEMS TO STORE
//
    IA12M = 5 * IA12M;
    IA12VL = NALLOC("A12VALS ", IA12N);
    IMSVAL = NALLOC("A12MSVAL", IA12M);
    IJA12S = NALLOC("JA12S   ", (IA12M + 1) / FACFR4);
//
//     FOLLOWING IS 3 ARRAYS; ONLY THE THIRD IS USED FOR USEHS RUNS.
//     REAL*8   NUMLX
//     REAL*8   2*LXMAX+2
//     INTEGER  3*NUMLX  ( WE TREAT AS REAL*8 TO ALLOW FOR CDC )
//
    IA12TM = 4 * NUMLX + 2 * LXMAX + 2;
L1750:
    IA12TM = NALLOC("A12TEMPS", IA12TM);
//
//     "IDWFI" WILL CONTAIN INFORMATION ABOUT THE INCIDENT
//     SCATTERING WAVE FUNCTIONS, FOR A GIVEN LI:
//
//        (3*K-2) = LASI - LI
//        (3*K-1) = JPI - 2*LI  (MAX. VALUE IF NO S. O.)
//        (3*K  ) = INDEX TO BASE OF WF IN ALLOCATOR
//
L1800:
    IWFII = NALLOC("IDWFI   ", (3 * NUMJS[1] + 1) / FACFR4);
    LWFII = Z[IWFII] * FACFR4 - FACFR4;
    I = 0;
    if (TCSWS[1]) I = JSPS[1];
    LASMIN = -I;
    LASMAX = I;
    NWFI = 0;
//     THE OUTER LOOP IS OVER LAS = ASYMPTOTIC L (TENSOR ONLY).
    for (LAS = LASMIN; LAS <= LASMAX; LAS += 2) { // 1819
//        THE INNER LOOP IS OVER JP = TOTAL PROJECTILE J
        JPMX = MIN0(2 * LAS, 0) + JSPS[1];
        JPMN = JPMX;
        if (SOSWS[1]) JPMN = MAX0(2 * LAS, 0) - JSPS[1];
        for (JPI = JPMN; JPI <= JPMX; JPI += 2) { // 1809
            NWFI = NWFI + 1;
            ILLOC(LWFII + 3 * NWFI - 2) = LAS;
            ILLOC(LWFII + 3 * NWFI - 1) = JPI;
            ILLOC(LWFII + 3 * NWFI) = (NWFI - 1) * NRIROI;
        } // 1809
    } // 1819
    I = 3 * NWFI;
    if (PBUGSW) {
        for (K = 1; K <= I; K++)
            std::printf(" IDWFI%10d\n", ILLOC(LWFII + K));
    }
//
    if (COUPSW) goto L2100;
//
//     "IDWFO" WILL CONTAIN INFORMATION ABOUT THE OUTGOING
//     SCATTERING WAVE FUNCTIONS, FOR A GIVEN LI:
//
//        (4*K-3) = LO - LI
//        (4*K-2) = LASO - LI
//        (4*K-1) = JPO - 2*LI
//        (4*K  ) = INDEX TO BASE OF WF ARRAY (SET IN INELDC).
//
    I = 0;
    if (ISTRIP != 0) I = MOD(LBP + LBT + LXMAX, 2);
    NOFLO = LXMAX + 1 - I;
    IWFIO = NALLOC("IDWFO   ", 4 * NOFLO * NUMJS[2] / FACFR4);
    LWFIO = Z[IWFIO] * FACFR4 - FACFR4;
    LOMAX = NOFLO - 1;
    LOMIN = -LOMAX;
    LASR = 0;
    if (TCSWS[2]) LASR = JSPS[2];
    NWFO = 0;
    for (LO = LOMIN; LO <= LOMAX; LO += 2) { // 1869
        LASMIN = LO - LASR;
        LASMAX = LO + LASR;
        for (LAS = LASMIN; LAS <= LASMAX; LAS += 2) { // 1859
            JPMX = 2 * MIN0(LO, LAS) + JSPS[2];
            JPMN = JPMX;
            if (SOSWS[2]) JPMN = 2 * MAX0(LO, LAS) - JSPS[2];
            for (JPO = JPMN; JPO <= JPMX; JPO += 2) { // 1849
//     ELIMINATE INITIAL JPO'S TOO LOW TO COUPLE TO LI.
                if (SOSWS[2] && NWFO == 0 && JPO + JSPS[1] + FORMF.JBT < 0) continue; // goto 1849
                NWFO = NWFO + 1;
                ILLOC(LWFIO + 4 * NWFO - 3) = LO;
                ILLOC(LWFIO + 4 * NWFO - 2) = LAS;
                ILLOC(LWFIO + 4 * NWFO - 1) = JPO;
            } // 1849
        } // 1859
    } // 1869
//     ELIMINATE FINAL JP'S TOO LARGE TO COUPLE TO LI.
    if (!SOSWS[2]) goto L1900;
    I = NWFO;
    for (K = 1; K <= I; K++) { // 1879
        if (ILLOC(LWFIO + 4 * NWFO - 1) <= JSPS[1] + FORMF.JBT) goto L1890;
        NWFO = NWFO - 1;
    } // 1879
L1890:
    if (NWFO < NOFLO * NUMJS[2])
        IREDEF(4 * NWFO / FACFR4, IWFIO);
L1900:
    I = 4 * NWFO;
    if (PBUGSW) {
        for (K = 1; K <= I; K++)
            std::printf(" IDWFO%10d\n", ILLOC(LWFIO + K));
    }
//
//     NOW ALLOCATE SPACE FOR THE WAVE FUNCTIONS (AT GAUSS POINTS)
//
    I = (NWFI * NRIROI + 1) / FACFR4;
    ILIR = NALLOC("WAVEAR  ", I);
    ILII = NALLOC("WAVEAI  ", I);
    I = (NWFO * NRIROI + 1) / FACFR4;
    ILOR = NALLOC("WAVEBR  ", I);
    ILOI = NALLOC("WAVEBI  ", I);
//     FOR TENSOR COUPLING, WAVETC NEEDS SOME SCRATCH STORAGE SPACE.
    I = 4 * NRIROI / FACFR4;
    for (K = 1; K <= 2; K++) { // 1909
        if (!TCSWS[K]) continue; // goto 1909
        IWAVRS[K] = NALLOC(NAMER[K].data, I);
        IWAVIS[K] = NALLOC(NAMEI[K].data, I);
    } // 1909
//     THIS POINTER CONTROLS WF STORAGE IN THE CIRCULAR BUFFER.
    NCRIT = NWFO * NRIROI - NRIROI;
//
//     "DW" WILL CONTAIN PRODUCTS OF SCATTERING WAVE FUNCTIONS.
//     "INDXDW" WILL CONTAIN POINTERS USED TO CALCULATE AND USE "DW".
//
    IDW = NALLOC("DW      ", NWFI * NWFO * 2);
    IDWI = NALLOC("INDXDW  ", NWFI * NWFO * 4 / FACFR4);
    goto L2200;
//
//     FOR C.C. THERE ARE NO OUTGOING WAVE ARRAYS BUT WE
//     STILL STORE INCOMING S SEPARATELY.
//
L2100:
    LCHNDF = FACFR4 * (Z[ICHNDF] - 1) + 1;
    LCHNDF = LCHNDF + ILLOC(LCHNDF + 2);
    ILLOC(LCHNDF + 9) = -1;
    ILLOC(LCHNDF + 10) = NWFI;
    ILLOC(LCHNDF + 11) = 0;
    ILLOC(LCHNDF + 12) = 0;

    // --- ANGSET lines 1545-1627 (L2200 through END) ---

    // Compute max number of H's
L2200:
    GRIDCM.NUMHS = INELCM.NMLOLX;

    if (COUPSW) goto L3000;
    I = GRIDCM.NUMHS * WAVCOM.NUMJS[1] * WAVCOM.NUMJS[2];

    // "RIROABS" for round-off error checking
    GRIDCM.IABS1 = NALLOC("RIROABS ", I);

    if (ISTRIP == 0) goto L3000;

    // "SUMHVALS" and "SUMIVALS" — H's before/after interpolation
    GRIDCM.ISMHVL = NALLOC("SUMHVALS", INTGER.NPSUM * GRIDCM.NUMHS);
    GRIDCM.ISMIVL = NALLOC("SUMIVALS", GRIDCM.NPSUMI * GRIDCM.NUMHS);

    // "I1REAL", "I1IMAG" — final 3D integrals
    // "IINDEX" — pointers into H and DW
    INELCM.ILILOR = NALLOC("I1REAL  ", I);
    INELCM.ILILOI = NALLOC("I1IMAG  ", I);
    GRIDCM.IIINDX = NALLOC("IINDEX  ", 4 * I / FACFR4);

    if (IHSAVE == 2) goto L3000;

    // Setup for high speed cosine
    I = (int)(DLOG(DFLOAT(INTGER.NCOSIN)) / DLOG(2.0) + 0.5);
    INTGER.NCOSIN = 1 << I;
    {
        double DICOSS = INTGER.NCOSIN;
        GRIDCM.COSSTP = 2.0 * CNSTNT.PI / DICOSS;
        GRIDCM.ICOSS4 = INTGER.NCOSIN / 4;
        GRIDCM.STPINV = DICOSS / (2.0 * CNSTNT.PI);
        GRIDCM.MSKCOS = INTGER.NCOSIN - 1;
        int ICOSSP = INTGER.NCOSIN + 1;
        GRIDCM.ICOSIN = NALLOC("COSINES ", ICOSSP);
        double X = -GRIDCM.COSSTP;
        for (int II = 1; II <= ICOSSP; II++) {
            X = X + GRIDCM.COSSTP;
            ALLOC(LOCPTRS.Z[GRIDCM.ICOSIN] + II - 1) = DCOS(X);
        }
    }

    // Table of integers in double precision form
    {
        I = INTGER.LMAX + INELCM.LXMAX + FORMF.LBP + FORMF.LBT;
        int IMIN = MAX0(I, 2 * (FORMF.LBP + FORMF.LBT) + 20);
        GRIDCM.INTOFF = IMIN + 1;
        int IMAX_val = MAX0(IMIN, 2 * I + 2);
        GRIDCM.IINTS = NALLOC("INTEGERS", IMAX_val + IMIN + 1);
        int LINTS = LOCPTRS.Z[GRIDCM.IINTS] - 1 + GRIDCM.INTOFF;
        IMIN = -IMIN;
        double X = IMIN;
        for (I = IMIN; I <= IMAX_val; I++) {
            ALLOC(LINTS + I) = X;
            X = X + 1.0;
        }
    }

L3000:
    IRTN = 1;
    return;

L9950:
    std::printf("\n**** ERROR IN READING OR WRITING H'S FILE.\n");
    return;
}
