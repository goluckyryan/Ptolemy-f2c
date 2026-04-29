// coulst_translated.cpp — COULST subroutine
// Translated from source.f lines 11210-11873
//
// COULST computes Coulomb integrals (FF, FG, GF, GG) for inelastic DWBA and
// coupled-channels calculations.  It:
//   1. Determines LIMOST (the maximum L needed for Coulomb extrapolation)
//      using asymptotic formulae (Alder et al., Rev.Mod.Phys. 28, 432, 1956).
//   2. For non-CC calculations, allocates ICL1FF/FG/GF/GG and integrates
//      INTEGRAL(SUMMAX -> inf) FF, FG, GF, GG by calling COULIN for each LX.
//   3. Allocates ICL2FF and computes pure-Coulomb (FF-only) integrals out to
//      LOMOST for the Coulomb extrapolation correction.
//   4. For coupled-channels (PROBLM==24), delegates to GENBNX, SETBFC, SETFG
//      and allocates CC-specific arrays (HOMO, INHR, PADE, SMATITER …).
//   Returns IRTN=1 on success, 0 on error.

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

// ---------------------------------------------------------------------------
// Local extern declarations for routines whose forward.h signatures differ
// from what COULST actually calls (the forward.h stubs were written for
// other callers).
// ---------------------------------------------------------------------------

// 4-argument RTXLNX: finds real solution of  A*X + B*LN(X) + C = 0
// (forward.h only declares the 1-arg stub used elsewhere)
extern double RTXLNX(double A, double B, double C, double ACC);

// GIVEAL / CHOPIT — memory allocator helpers
extern int GIVEAL(int IWRDS);
extern int CHOPIT(int IWRDS);

// COULIN is declared in ptolemy_forward.h (real signature, source.f L10647)
// and implemented in coulin_translated.cpp.

// Full GENBNX signature (source.f L17285-17287):
extern void GENBNX_full(
    int IPASS, int& NBINDX, int& MBINDX, int* BINDEX,
    int NBASCP, int MBASCP, int* BASCUP,
    int MBASDF, int* BASDEF,
    int MCHNDF, float* CHNDEF,
    int& NMBFAC, int& NMFFAC, int& MXLXGS, int& MXLXBF,
    int COULSW, int PBUGSW);

// Full SETBFC signature (source.f L31512-31516):
extern void SETBFC_full(
    int NBINDX, int MBINDX, int* BINDEX,
    int MCHNDF, int* CHNDEF_i, float* RCHNDF,
    int NUMLIS, int* LIS, double* R2S_arg, int LMXMX,
    double* SIG1, double* SIG2,
    double RLOWER, int ALLSW, int IRORC,
    int NMBFAC, float* BFAC,
    double* FFWORK, int IDIM1, int LDLDIM,
    double* WRKWK, double* FIWORK, int IDIM2,
    double* WRKST, int& IRET, double& CLTIME);

// Full SETFG signature (source.f L32481-32483):
extern void SETFG_full(
    int NCHNDF, int MCHNDF, int* CHNDEF, float* RCHNDF,
    int NBASDF, int MBASDF, int* BASDEF,
    int MCHNVL, double* CHNVAL,
    int NUMLIS, double* LIS,
    double* FG,
    int L1, int L2P1,
    double* FGWORK, double R1, double R2, int& IRET);

// ---------------------------------------------------------------------------
// Stub implementations for routines not yet fully translated.
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------

void COULST(int& IRTN) {
    //
    // COULOMB INTEGRALS FOR INELASTIC DWBA & COUPLED CHANNELS
    //

    // Local logicals (Fortran LOGICAL = C++ int)
    int PBUGSW;  // LOGICAL
    int PINFSW;  // LOGICAL
    int COULSW;  // LOGICAL
    int COUPSW;  // LOGICAL

    // Local scalars
    int    IPRNT;
    int    NCHNDF;
    int    LBETAR;
    int    NCHN, LX, LI, LO, I, K, IDEL, ILI;
    int    LBINDX, MBINDX_loc, NBINDX_loc;
    int    LBASCP, NBASCP, MBASCP;
    int    LBASDF, NBASDF, MBASDF;
    int    LCHNDF, MCHNDF;
    int    IBINDX_loc;
    int    MAXDL1, LDLDIM, LMNMN1, LMXMX, LMXMX1;
    int    IWRKST;
    int    IWRKFF, IWRKFI, IWRKFG, IWRKWK;
    int    IWRKFO, IWRKGI, IWRKGO;
    int    IWRKFG_CC;
    int    LWRKFF, LWRKFI, LWRKFG, LWRKGG;
    int    LWRKFO, LWRKGI, LWRKGO, LWRKWK, LWRKST;
    int    LSIG1, LSIG2;
    int    LLIS;
    int    LINDXS;
    int    KBASE, LDELMN, NDEL, MAXDEL;
    int    IDELMN, LDEL, ID, IL, I1;
    int    IDIM1, IDIM2;
    int    LMNMN;
    int    LMP1;
    int    MLX;
    int    L;
    int    IFG;
    int    IWRKFG2;
    int    L2P1, NWRDS;
    int    IRET;
    int    LCL1FF;
    int    LCHNDF2;
    int    LBASDF2;

    // double temporaries
    double TT, T2, T3, T4, TALL;
    double CLTIME;
    double DLCRIT;
    double ATERM, BTERM, CTERM;
    double CI, CO;

    // timing
    double T_TT4;

    // -------------------------------------------------------------------------
    // Convenient references to COMMON block members
    // -------------------------------------------------------------------------
    int&    IPRINT  = INTGER.IPRINT;
    int&    LMIN    = INTGER.LMIN;
    int&    LMAX    = INTGER.LMAX;
    int&    MAXLEX  = INTGER.MAXLEX;
    int&    NPCOUL  = INTGER.NPCOUL;
    int&    MXCOUL  = INTGER.MXCOUL;
    int&    MAXITR  = INTGER.MAXITR;
    int&    NBACK   = INTGER.NBACK;

    int&    PROBLM  = SWITCH.PROBLM;
    int&    NOBFAC  = SWITCH.NOBFAC;
    int&    ISAVSM  = SWITCH.ISAVSM;

    int&    LIMOST  = INTRNL.LIMOST;
    int&    IEXCIT  = INTRNL.IEXCIT;

    double& SUMMAX  = FLOAT_common.SUMMAX;
    double& ACCINE  = FLOAT_common.ACCINE;
    double& COULML  = FLOAT_common.COULML;
    double& DWCUT   = FLOAT_common.DWCUT;

    // KANDM
    int*    ISIGS   = KANDM.ISIGS;    // 1-based
    int&    LOMOST  = KANDM.LOMOST;
    int&    LCRIT   = KANDM.LCRIT;
    double* ETAS    = KANDM.ETAS;     // 1-based

    // GRIDCM — using the inelastic-code C++ member names that correspond to
    // the COULST Fortran /gridcmf/ view.  Mapping derived from comparing the
    // two COMMON /GRIDCM/ declarations (source.f L1089 vs L11345).
    //
    // COULST Fortran name → C++ GridcmCommon field (same physical slot)
    //   IBINDX → IMSVAL,  MBINDX → IINTS,   NBINDX → INTOFF
    //   NMFFAC  → IA12VL,  MXLXGS → ICOSIN,  NMBFAC  → IPHIT
    //   IINHR   → IPHIP,   IINHI  → IPHI
    //   IRHSR   → LOGIC,   IRHSI  → LOCTRP
    //   NPTMIN  → IWIO,    NUMPT  → NRIROI
    //   ISMX    → IXLAM,   ISMSQ  → INLAM,   IPADE  → IA12N
    //   IHOMO   → MSKCOS,  ISMTSV → MAXCNT
    //   ICL1FF  → IGRDC2[1], ICL1FG → IGRDC2[2]
    //   ICL1GF  → IGRDC2[3], ICL1GG → IGRDC2[4]
    //
    int&    IBINDX  = GRIDCM.IMSVAL;
    int&    MBINDX  = GRIDCM.IINTS;
    int&    NBINDX  = GRIDCM.INTOFF;
    int&    NMFFAC  = GRIDCM.IA12VL;
    int&    MXLXGS  = GRIDCM.ICOSIN;
    int&    NMBFAC  = GRIDCM.IPHIT;
    int&    IINHR   = GRIDCM.IPHIP;
    int&    IINHI   = GRIDCM.IPHI;
    int&    IRHSR   = GRIDCM.LOGIC;
    int&    IRHSI   = GRIDCM.LOCTRP;
    int&    NPTMIN  = GRIDCM.IWIO;    // Fortran NPTMIN = position 17 = C++ IWIO
    int&    NUMPT   = GRIDCM.NRIROI;  // Fortran NUMPT  = position 19 = C++ NRIROI
    int&    ISMX    = GRIDCM.IXLAM;
    int&    ISMSQ   = GRIDCM.INLAM;
    int&    IPADE   = GRIDCM.IA12N;
    int&    IHOMO   = GRIDCM.MSKCOS;
    int&    ISMTSV  = GRIDCM.MAXCNT;
    // ICL1FF/FG/GF/GG live in GRIDCM.IGRDC2[1..4]
    int&    ICL1FF  = GRIDCM.IGRDC2[1];
    int&    ICL1FG  = GRIDCM.IGRDC2[2];
    int&    ICL1GF  = GRIDCM.IGRDC2[3];
    int&    ICL1GG  = GRIDCM.IGRDC2[4];

    // WAVCOM
    double* RSTEPS  = WAVCOM.RSTEPS;  // 1-based
    int*    NFS     = WAVCOM.NFS;     // 1-based
    int*    NUMJS   = WAVCOM.NUMJS;   // 1-based

    // INELCM
    int&    LXMIN   = INELCM.LXMIN;
    int&    LXMAX   = INELCM.LXMAX;
    int&    NUMLIS  = INELCM.NUMLIS;
    int&    ILIS    = INELCM.ILIS;
    int&    ICL2FF  = INELCM.ICL2FF;
    int&    NSMAT   = INELCM.NSMAT;
    int&    NLX     = INELCM.NLX;
    int&    NSPLI   = INELCM.NSPLI;
    int&    IINDXS  = INELCM.IINDXS;
    int*    IBETAS  = INELCM.IBETAS;  // 1-based
    int&    NMLOLX  = INELCM.NMLOLX;

    // CCBLK
    int&    ICHNDF_loc = CCBLK.ICHNDF;
    int&    ICHNCP_loc = CCBLK.ICHNCP;
    int&    IBASDF_loc = CCBLK.IBASDF;
    int&    IBASCP_loc = CCBLK.IBASCP;
    int&    MCHNVL_loc = CCBLK.MCHNVL;
    int&    ICHNVL_loc = CCBLK.ICHNVL;
    int&    IHOMOA_loc = CCBLK.IHOMOA;
    int&    IHOMOB_loc = CCBLK.IHOMOB;
    int&    IINHR8_loc = CCBLK.IINHR8;
    int&    IINHI8_loc = CCBLK.IINHI8;

    // FORMF: R2S(4) EQUIVALENCE (RBNDS(1))
    // R2S(4) = RBNDS(4) = FORMF.RBNDS[4]
    // In C++: FORMF.RBNDS[4] (1-based → index 4)
    // R2S(I) = FORMF.RBNDS[I]

    // FTIME
    float* TIMES    = FTIME.TIMES;    // 1-based [1..8]

    // ALLOCS
    int&    FACFR4  = ALLOCS.FACFR4;

    // =========================================================================
    // Begin executable code
    // =========================================================================

    IRTN = 0;

    IPRNT  = IPRINT % 10;
    PBUGSW = (IPRNT >= 3) ? 1 : 0;
    PINFSW = (IPRNT >= 2) ? 1 : 0;

    //
    // COUPSW = ( PROBLM .EQ. 24 )
    //
    COUPSW = (PROBLM == 24) ? 1 : 0;

    //
    // THESE INTEGRALS WILL BE JUST ZERO IF THE BETA(COULOMB) IS ZERO.
    // IN THIS CASE THE LARGE COULOMB LMAX IS ALSO NOT NEEDED.
    //
    LIMOST        = LMAX;
    TIMES[5]      = 0.0f;
    TIMES[6]      = 0.0f;
    TIMES[7]      = 0.0f;
    TIMES[8]      = 0.0f;
    CLTIME        = 0.0;
    COULSW        = 0;   // .FALSE.
    int MXLXBF    = LXMAX;

    if (COUPSW) goto L100;

    //
    // NON-CC CASE: CHECK WHETHER ANY COULOMB BETA IS NONZERO
    //
    NCHNDF = 2;
    // LBETAR = Z(IBETAS(3)) - LXMIN/2
    // In Fortran: LBETAR is used as base pointer so ALLOC(LBETAR + LX/2) = beta(LX)
    LBETAR = LOCPTRS.Z[IBETAS[2]] - LXMIN / 2;

    for (LX = LXMIN; LX <= LXMAX; LX += 2) {
        if (ALLOC(LBETAR + LX / 2) == 0.0) goto L69;
        COULSW = 1;   // .TRUE.
        goto L200;
    L69:
        ;
    }

    //
    // NO COULOMB TERMS WILL BE NEEDED
    //
    COULSW = 0;
    NOBFAC = 1;
    goto L300;

    // -----------------------------------------------------------------------
    // FOR COUPLED CHANNELS: go through the basis coupling array,
    // find required storage size for B-factors and FF2 integrals,
    // construct B index array and pointers to it.
    // -----------------------------------------------------------------------

L100:
    //
    // IBINDX = GIVEAL(I) — allocate a named slot for the B-index array
    //
    IBINDX_loc = GIVEAL(0);   // 'I' is uninitialized in Fortran — pass 0
    IBINDX     = IBINDX_loc;
    // ANAMES(IBINDX) = BNDXNM  (character assignment — skip in C++)
    // LBINDX = FACFR4*(Z(IBINDX)-1) + 1
    LBINDX = FACFR4 * (LOCPTRS.Z[IBINDX] - 1) + 1;
    MBINDX = 9;

    //
    // Read base array pointers out of the LOCptrs table
    //
    LBASCP = FACFR4 * (LOCPTRS.Z[IBASCP_loc] - 1) + 1;
    NBASCP = ILLOC(LBASCP);
    MBASCP = ILLOC(LBASCP + 1);
    LBASCP = LBASCP + ILLOC(LBASCP + 2);

    LBASDF = FACFR4 * (LOCPTRS.Z[IBASDF_loc] - 1) + 1;
    NBASDF = ILLOC(LBASDF);
    MBASDF = ILLOC(LBASDF + 1);
    LBASDF = LBASDF + ILLOC(LBASDF + 2);

    LCHNDF = FACFR4 * (LOCPTRS.Z[ICHNDF_loc] - 1) + 1;
    NCHNDF = ILLOC(LCHNDF);
    MCHNDF = ILLOC(LCHNDF + 1);
    LCHNDF = LCHNDF + ILLOC(LCHNDF + 2);

    //
    // CALL GENBNX ( 1, NBINDX, MBINDX, ILLOC(LBINDX), NBASCP,
    //   MBASCP, ILLOC(LBASCP), MBASDF, ILLOC(LBASDF), MCHNDF,
    //   ALLOC4(LCHNDF), NMBFAC, NMFFAC, MXLXGS, MXLXBF,
    //   COULSW, PBUGSW )
    //
    GENBNX_full(1, NBINDX, MBINDX, &ILLOC(LBINDX),
                NBASCP, MBASCP, &ILLOC(LBASCP),
                MBASDF, &ILLOC(LBASDF),
                MCHNDF, &ALLOC4(LCHNDF),
                NMBFAC, NMFFAC, MXLXGS, MXLXBF,
                COULSW, PBUGSW);

    {
        int I_loc = MAX0((NBINDX * MBINDX + 1) / FACFR4, 1);
        int ixx = CHOPIT(I_loc);
        (void)ixx;
    }

    if (NBINDX != 0) goto L150;

    //
    // B-FACTORS ARE NOT NEEDED OR ARE INHIBITED.
    //
    NOBFAC = 1;

L150:
    if (MAXLEX == 0 || !COULSW) goto L300;

    // -----------------------------------------------------------------------
    // COULOMB COUPLING EXISTS: compute LIMOST using asymptotic formula
    // -----------------------------------------------------------------------

    //
    // FIRST WE FIND THE MAXIMUM LI TO BE USED FOR THE COULOMB INTEGRALS
    // WE USE AN ASYMPTOTIC FORM (Alder et al., Rev.Mod.Phys. 28, 432, 1956,
    // formula II E.83).
    // FOR CC THIS IS DONE SEPARATELY FOR EACH CHANNEL AND THE RESULT
    // STORED IN THE CHANNEL DEFINITION ARRAY.
    //

L200:
    if (PINFSW) std::printf("-CHN LX  LO-LI  MAXIMUM L NEEDED\n");

    DLCRIT = (double)LCRIT;
    LX = LXMIN;

    for (NCHN = 2; NCHN <= NCHNDF; NCHN++) {
        LI = LMAX;

        if (!COUPSW) goto L230;

        {
            // I = LCHNDF + (NCHN-1)*MCHNDF
            int I_loc = LCHNDF + (NCHN - 1) * MCHNDF;
            // IF ( ILLOC(I+14) .EQ. NOTDEF ) GO TO 249
            if (ILLOC(I_loc + 14) == static_cast<int>(INTRNL.NOTDEF)) goto L249;
            // ETAS(2) = ALLOC4(I+5)
            ETAS[2] = (double)ALLOC4(I_loc + 5);
            LX = ILLOC(I_loc + 14);
        }

L230:
        MLX = -LX;
        for (int LDEL = MLX; LDEL <= LX; LDEL += 2) {
            // ATERM = DABS( (ETAS(2)-ETAS(1))/ETAS(1) )
            ATERM = DABS((ETAS[2] - ETAS[1]) / ETAS[1]);
            // BTERM = (LX+LDEL+1.D0)/2
            BTERM = (LX + LDEL + 1.0) / 2.0;
            // CTERM = DLOG(DWCUT) - ATERM*(DLCRIT+.5) - BTERM*DLOG(DLCRIT+.5)
            CTERM = DLOG(DWCUT) - ATERM * (DLCRIT + 0.5) - BTERM * DLOG(DLCRIT + 0.5);
            // L = RTXLNX( ATERM, BTERM, CTERM, .1D0 ) + .5
            L = (int)(RTXLNX(ATERM, BTERM, CTERM, 0.1) + 0.5);
            LI = MAX0(LI, L);

            if (PINFSW) std::printf("%4d%3d%6d%13d\n", NCHN, LX, LDEL, L);
            if (PBUGSW) std::printf("+%49s%15.5G%15.5G%15.5G\n", "", ATERM, BTERM, CTERM);
        }

        //
        // LIMIT BY MAXLEX
        //
        LI = MIN0(LI, LMAX + MAXLEX);

        if (COUPSW) {
            int I_loc = LCHNDF + (NCHN - 1) * MCHNDF;
            ILLOC(I_loc + 14) = LI;
        }
        LIMOST = MAX0(LIMOST, LI);

    L249:
        ;
    }  // end NCHN loop

    // -----------------------------------------------------------------------

L300:
    LOMOST = LIMOST + LXMAX;

    //
    // NOW GENERATE POINTERS TO THE B AND FF2 ARRAYS
    //
    if (COUPSW) {
        //
        // CALL GENBNX ( 2, NBINDX, MBINDX, ILLOC(LBINDX), NBASCP,
        //   MBASCP, ILLOC(LBASCP), MBASDF, ILLOC(LBASDF), MCHNDF,
        //   ALLOC4(LCHNDF), NMBFAC, NMFFAC, MXLXGS, MXLXBF,
        //   COULSW, PBUGSW )
        //
        GENBNX_full(2, NBINDX, MBINDX, &ILLOC(LBINDX),
                    NBASCP, MBASCP, &ILLOC(LBASCP),
                    MBASDF, &ILLOC(LBASDF),
                    MCHNDF, &ALLOC4(LCHNDF),
                    NMBFAC, NMFFAC, MXLXGS, MXLXBF,
                    COULSW, PBUGSW);
    }

    //
    // EXPAND THE COULOMB PHASE SHIFT ARRAYS
    // WE ARE SLOPPY AND FIND SOME MORE THAN ARE REALLY NECESSARY
    //
    {
        int LMXMX_loc = LOMOST + 2 * LXMAX;
        IREDEF(LMXMX_loc + 1, ISIGS[1]);

        if (COUPSW) goto L370;

        //
        // NON-CC: expand both sigma arrays and fill in new entries
        //
        IREDEF(LMXMX_loc + 1, ISIGS[2]);

        LSIG1 = LOCPTRS.Z[ISIGS[1]];
        LSIG2 = LOCPTRS.Z[ISIGS[2]];

        CI = ALLOC(LSIG1 + LMAX);
        CO = ALLOC(LSIG2 + LMAX);

        LMP1 = LMAX + 1;
        for (L = LMP1; L <= LMXMX_loc; L++) {
            CI = CI + DATAN(ETAS[1] / (double)L);
            CO = CO + DATAN(ETAS[2] / (double)L);
            ALLOC(LSIG1 + L) = CI;
            ALLOC(LSIG2 + L) = CO;
        }

        //
        // WE SETUP INTEGRAL(R TO INF) FF, FG, GF AND GG FOR LMIN TO LMAX,
        // AND JUST FF FOR (0 TO INF) AND LMIN TO LIMOST.
        // THIS IS FOR NON-CC CALCULATIONS — FIRST THE (R, INF) STUFF
        //

        // Allocate the four integral arrays
        ICL1FF = NALLOC("FF1INTS ", NSMAT);
        ICL1FG = NALLOC("FG1INTS ", NSMAT);
        ICL1GF = NALLOC("GF1INTS ", NSMAT);
        ICL1GG = NALLOC("GG1INTS ", NSMAT);

        for (I = 1; I <= NSMAT; I++) {
            ALLOC(LOCPTRS.Z[ICL1FF] - 1 + I) = 0.0;
            ALLOC(LOCPTRS.Z[ICL1FG] - 1 + I) = 0.0;
            ALLOC(LOCPTRS.Z[ICL1GF] - 1 + I) = 0.0;
            ALLOC(LOCPTRS.Z[ICL1GG] - 1 + I) = 0.0;
        }
        goto L400;
    }

    // -----------------------------------------------------------------------
    // FOR COUPLED CHANNELS: store everything in one big array.
    // -----------------------------------------------------------------------

L370:
    {
        int LMXMX_loc = LOMOST + 2 * LXMAX;
        // ISIGS(2) = NALLOC( LMXMX+1, 'COULSIGB' )
        ISIGS[2] = NALLOC("COULSIGB", LMXMX_loc + 1);
        // ISMATS(2) = NALLOC( 2*LMXMX+2, 'SOPTICAL' )
        KANDM.ISMATS[2] = NALLOC("SOPTICAL", 2 * LMXMX_loc + 2);
        NUMJS[2] = NUMJS[1];

        TT = (double)0.0;   // TT = second() — timing stub
        ICL1FF = 1;
        ICL1FG = 1;

        if (NOBFAC == 1) goto L800;

        // ICL1FF = NALLOC( (2*NMBFAC*NUMLIS+1)/FACFR4, 'BFACTORS' )
        ICL1FF = NALLOC("BFACTORS", (2 * NMBFAC * NUMLIS + 1) / FACFR4);
        INIT4(ICL1FF, 0.0f);
        // ICL1FG = NALLOC( (4*NMBFAC+1)/FACFR4, 'BFACFORJ' )
        ICL1FG = NALLOC("BFACFORJ", (4 * NMBFAC + 1) / FACFR4);
    }

    // -----------------------------------------------------------------------
    // COMMON PATH FOR BOTH CC AND NON-CC: set up work array dimensions
    // -----------------------------------------------------------------------

L400:
    {
        // Compute work-array dimensions (valid for both CC and non-CC paths)
        MAXDL1  = MAX0(2, MXLXBF);
        LDLDIM  = MAXDL1 + 1;
        MAXDL1  = MAX0(MAXDL1, LXMAX);
        LMNMN1  = MAX0(LMIN - MAXDL1, 0);
        IWRKST  = NALLOC("COULSTRT", 16 * LDLDIM);

        TT = 0.0;  // TT = second()

        if (NOBFAC == 1) goto L500;

        // Non-CC branch: work arrays for FF/FG/GF/GG integrals
        LMXMX   = LMAX + LXMAX;
        LMXMX1  = LMAX + MAXDL1;
        IDIM1   = LDLDIM * (2 * (LMXMX1 - LMNMN1) + 1);
        IWRKFF  = NALLOC("FFWORK  ", 4 * IDIM1);
        IDIM2   = MAX0(LMXMX1 + MAXDL1, 4 * MXCOUL) + 1;
        IWRKFI  = NALLOC("FIWORK  ", 4 * IDIM2);
        IWRKWK  = NALLOC("COULWORK",
                          MAX0(2 * NPCOUL + 4 * MXCOUL,
                               2 * (LMXMX1 + MAXDL1) + 1,
                               2 * NPCOUL + LMXMX1 + MAXDL1) + 1);

        LWRKFF  = LOCPTRS.Z[IWRKFF];
        LWRKFI  = LOCPTRS.Z[IWRKFI];
        LWRKWK  = LOCPTRS.Z[IWRKWK];
        LWRKST  = LOCPTRS.Z[IWRKST];

        LLIS = FACFR4 * (LOCPTRS.Z[ILIS] - 1);

        LSIG1 = LOCPTRS.Z[ISIGS[1]];
        LSIG2 = LOCPTRS.Z[ISIGS[2]];

        if (COUPSW) goto L700;

        // NON-CC: set up FG, GF, GG, and FO/GI/GO work arrays
        // Declare all derived locals before any goto target that follows
        {
        LINDXS  = LOCPTRS.Z[IINDXS] * FACFR4 - FACFR4;
        LWRKFG  = LWRKFF + IDIM1;
        // LWRKFG = LWRKFF+IDIM1; LWRKGF = LWRKFG+IDIM1; LWRKGG = LWRKGF+IDIM1
        int LWRKGF_loc = LWRKFG + IDIM1;
        int LWRKGG_loc = LWRKGF_loc + IDIM1;
        int LWRKFO_loc = LWRKFI + IDIM2;
        int LWRKGI_loc = LWRKFO_loc + IDIM2;
        int LWRKGO_loc = LWRKGI_loc + IDIM2;

        //
        // DO 439 LX = LXMIN, LXMAX, 2
        //   LMNMN = MAX0( LMIN-(LX+1)/2, LX/2 )
        //   CALL COULIN ( LX+1, LX, LMNMN, LMXMX, ETAS(2), AKO,
        //     ALLOC(LSIG2),  ETAS(1), AKI, ALLOC(LSIG1), SUMMAX, .TRUE.,
        //     ALLOC(LWRKFF), ALLOC(LWRKFG), ALLOC(LWRKGF), ALLOC(LWRKGG),
        //     LDLDIM, ACCINE, COULML, MXCOUL, NPCOUL,
        //     ALLOC(LWRKWK), ALLOC(LWRKFI), ALLOC(LWRKFO),
        //     ALLOC(LWRKGI), ALLOC(LWRKGO), ALLOC(LWRKST),
        //     MOD(IPRINT/100, 10), IRET, CLTIME )
        //   IF ( IRET .NE. 0 )  GO TO 9960
        //
        for (LX = LXMIN; LX <= LXMAX; LX += 2) {
            LMNMN = MAX0(LMIN - (LX + 1) / 2, LX / 2);
            COULIN(
                LX + 1, LX, LMNMN, LMXMX,
                ETAS[2], KANDM.AKO, &ALLOC(LSIG2),
                ETAS[1], KANDM.AKI, &ALLOC(LSIG1),
                SUMMAX, 1 /*TRUE*/,
                &ALLOC(LWRKFF), &ALLOC(LWRKFG), &ALLOC(LWRKGF_loc), &ALLOC(LWRKGG_loc),
                LDLDIM,
                ACCINE, COULML, MXCOUL, NPCOUL,
                &ALLOC(LWRKWK), &ALLOC(LWRKFI), &ALLOC(LWRKFO_loc),
                &ALLOC(LWRKGI_loc), &ALLOC(LWRKGO_loc), &ALLOC(LOCPTRS.Z[IWRKST]),
                IPRINT / 100 % 10, IRET, CLTIME);

            if (IRET != 0) goto L9960;

            MAXDEL = MAX0(2, LX);

            //
            // STORE THE INTEGRALS IN THE SAME ORDER AS THE I(LX,LI,LO) ARRAYS.
            // MULTIPLY IN THE CHARGE AND RADIUS FACTORS (BUT NOT THE BETA'S)
            //
            // K = LX+1 + NLX*(LX-INDXJS(3,IEXCIT)/2)
            // INDXJS(3,IEXCIT) = (&INELCM.JPMIN)[3*IEXCIT + 2]   (0-based stride 4, col 3)
            // Fortran: INDXJS is equivalenced to JPMIN with dimension (4,2)
            //   INDXJS(3,IEXCIT) means row 3, col IEXCIT → offset (IEXCIT-1)*4 + (3-1) = 4*(IEXCIT-1)+2
            {
                int INDXJS_3_IEXCIT = (&INELCM.JPMIN)[4 * (IEXCIT - 1) + 2];
                K = LX + 1 + NLX * (LX - INDXJS_3_IEXCIT / 2);
                KBASE = ILLOC(LINDXS + 3 * K - 2);
                LDELMN = ILLOC(LINDXS + 3 * K - 1);
                NDEL  = ILLOC(LINDXS + 3 * K);

                for (ILI = 1; ILI <= NUMLIS; ILI++) {
                    LI = ILLOC(LLIS + ILI);
                    IDELMN = 1;
                    if (LI < LX) IDELMN = LX + 1 - LI;

                    for (IDEL = IDELMN; IDEL <= NDEL; IDEL++) {
                        LDEL = LDELMN - 2 + 2 * IDEL;
                        LO   = LI + LDEL;
                        // ID = iSHFT( MAXDEL-LDEL, -1 )  — arithmetic right shift by 1
                        ID   = ISHFT(MAXDEL - LDEL, -1);
                        IL   = LO + LI - 2 * LMNMN;
                        I1   = ID + IL * LDLDIM;
                        I    = NSPLI * (ILI - 1) + KBASE - 2 + IDEL;

                        // R2S(4) = RBNDS(4) = FORMF.RBNDS[4]
                        double R2S4 = FORMF.RBNDS[4];
                        ALLOC(LOCPTRS.Z[ICL1FF] + I) = -R2S4 * ALLOC(LWRKFF + I1);
                        ALLOC(LOCPTRS.Z[ICL1FG] + I) = -R2S4 * ALLOC(LWRKFG + I1);
                        ALLOC(LOCPTRS.Z[ICL1GF] + I) = -R2S4 * ALLOC(LWRKGF_loc + I1);
                        ALLOC(LOCPTRS.Z[ICL1GG] + I) = -R2S4 * ALLOC(LWRKGG_loc + I1);
                    }
                }
            }

        }  // end DO 439 LX loop

        }  // end non-CC inner block (LWRKGF_loc etc.)

        //
        //  490 TIMES(5) = CLTIME
        //      TIMES(7) = second() - TT
        //
L490:
        TIMES[5] = (float)CLTIME;
        TIMES[7] = 0.0f;  // second() - TT stub

        // -----------------------------------------------------------------------
        // NOW THE PURE COULOMB (FF) INTEGRALS OUT TO LOMOST
        // -----------------------------------------------------------------------

L500:
        IDIM1   = NMLOLX * (LIMOST - LMIN + 1);
        ICL2FF  = NALLOC("FF2INTS ", IDIM1);

        LMXMX   = LOMOST;
        LMXMX1  = LMXMX + MAXDL1 - LXMAX;

        for (I = 1; I <= IDIM1; I++) {
            ALLOC(LOCPTRS.Z[ICL2FF] - 1 + I) = 0.0;
        }

        if (!COULSW) goto L990;

        IDIM1   = LDLDIM * (2 * (LMXMX1 - LMNMN1) + 1);
        IWRKFF  = NALLOC("FFWORK  ", IDIM1);
        IDIM2   = MAX0(LMXMX1 + MAXDL1, 4 * MXCOUL) + 1;
        IWRKFI  = NALLOC("FIWORK  ", IDIM2);
        IWRKFO  = NALLOC("FOWORK  ", IDIM2);
        IWRKGI  = NALLOC("GIWORK  ", IDIM2);
        IWRKGO  = NALLOC("GOWORK  ", IDIM2);
        IWRKWK  = NALLOC("COULWORK",
                          MAX0(2 * NPCOUL + 4 * MXCOUL,
                               2 * (LMXMX1 + MAXDL1) + 1,
                               2 * NPCOUL + LMXMX1 + MAXDL1) + 1);

        LINDXS  = LOCPTRS.Z[IINDXS] * FACFR4 - FACFR4;

        TT = 0.0;  // TT = second()

        //
        // DO 599 LX = LXMIN, LXMAX, 2
        //
        for (LX = LXMIN; LX <= LXMAX; LX += 2) {
            LMNMN = MAX0(LMIN - (LX + 1) / 2, LX / 2);

            // FF-only integrals from 0 to INF (SUMMAX=0.D0, ALLSW=.FALSE.)
            // DUMMY8 placeholders for FG, GF, GG (not used)
            // Fortran passes DUMMY8 (a real*8 array of size 1) for unused args
            static double DUMMY8[2] = {0.0, 0.0};

            COULIN(
                LX + 1, LX, LMNMN, LMXMX,
                ETAS[2], KANDM.AKO, &ALLOC(LOCPTRS.Z[ISIGS[2]]),
                ETAS[1], KANDM.AKI, &ALLOC(LOCPTRS.Z[ISIGS[1]]),
                0.0 /*SUMMAX=0*/, 0 /*FALSE*/,
                &ALLOC(LOCPTRS.Z[IWRKFF]), DUMMY8, DUMMY8, DUMMY8,
                LDLDIM,
                ACCINE, COULML, MXCOUL, NPCOUL,
                &ALLOC(LOCPTRS.Z[IWRKWK]),
                &ALLOC(LOCPTRS.Z[IWRKFI]),
                &ALLOC(LOCPTRS.Z[IWRKFO]),
                &ALLOC(LOCPTRS.Z[IWRKGI]),
                &ALLOC(LOCPTRS.Z[IWRKGO]),
                &ALLOC(LOCPTRS.Z[IWRKST]),
                IPRINT / 100 % 10, IRET, CLTIME);

            if (IRET != 0) goto L9960;

            MAXDEL = MAX0(2, LX);

            {
                int INDXJS_3_IEXCIT = (&INELCM.JPMIN)[4 * (IEXCIT - 1) + 2];
                K = LX + 1 + NLX * (LX - INDXJS_3_IEXCIT / 2);
                KBASE  = ILLOC(LINDXS + 3 * K - 2);
                LDELMN = ILLOC(LINDXS + 3 * K - 1);
                NDEL   = ILLOC(LINDXS + 3 * K);

                for (LI = LMIN; LI <= LIMOST; LI++) {
                    IDELMN = 1;
                    if (LI < LX) IDELMN = LX + 1 - LI;

                    for (IDEL = IDELMN; IDEL <= NDEL; IDEL++) {
                        LDEL = LDELMN - 2 + 2 * IDEL;
                        LO   = LI + LDEL;
                        ID   = ISHFT(MAXDEL - LDEL, -1);
                        IL   = LO + LI - 2 * LMNMN;
                        I1   = ID + IL * LDLDIM;
                        I    = NSPLI * (LI - LMIN) + KBASE - 2 + IDEL;

                        double R2S4 = FORMF.RBNDS[4];
                        ALLOC(LOCPTRS.Z[ICL2FF] + I) = -R2S4 * ALLOC(LOCPTRS.Z[IWRKFF] + I1);
                    }
                }
            }

        }  // end DO 599 LX loop

        TIMES[6] = (float)CLTIME;
        TIMES[8] = 0.0f;  // second() - TT stub

        //
        // Free the work arrays by negating their LOCptrs entries
        // (Fortran idiom: Z(I) = -Z(I) releases the allocation)
        //
        LOCPTRS.Z[IWRKFO] = -LOCPTRS.Z[IWRKFO];
        LOCPTRS.Z[IWRKGI] = -LOCPTRS.Z[IWRKGI];
        LOCPTRS.Z[IWRKGO] = -LOCPTRS.Z[IWRKGO];
        LOCPTRS.Z[IWRKFF] = -LOCPTRS.Z[IWRKFF];
        LOCPTRS.Z[IWRKFI] = -LOCPTRS.Z[IWRKFI];
        LOCPTRS.Z[IWRKWK] = -LOCPTRS.Z[IWRKWK];
        LOCPTRS.Z[IWRKST] = -LOCPTRS.Z[IWRKST];

        goto L990;

        // -----------------------------------------------------------------------
        // FOR COUPLED CHANNELS: LOOP THROUGH BINDEX
        // -----------------------------------------------------------------------

L700:
        //
        //  700 LBINDX = FACFR4*(Z(IBINDX)-1) + 1
        //      LCL1FF = FACFR4*(Z(ICL1FF)-1) + 1
        //      LCHNDF = FACFR4*(Z(ICHNDF)-1) + 1
        //      LCHNDF = LCHNDF + ILLOC(LCHNDF+2)
        //
        LBINDX  = FACFR4 * (LOCPTRS.Z[IBINDX] - 1) + 1;
        LCL1FF  = FACFR4 * (LOCPTRS.Z[ICL1FF] - 1) + 1;
        LCHNDF2 = FACFR4 * (LOCPTRS.Z[ICHNDF_loc] - 1) + 1;
        LCHNDF2 = LCHNDF2 + ILLOC(LCHNDF2 + 2);

        //
        // CALL SETBFC
        //
        SETBFC_full(
            NBINDX, MBINDX, &ILLOC(LBINDX),
            MCHNDF, &ILLOC(LCHNDF2), &ALLOC4(LCHNDF2),
            NUMLIS, &ILLOC(LLIS + 1), &FORMF.RBNDS[1], LMXMX,
            &ALLOC(LSIG1), &ALLOC(LSIG2),
            SUMMAX, 1 /*TRUE*/, 2,
            NMBFAC, &ALLOC4(LCL1FF),
            &ALLOC(LWRKFF), IDIM1, LDLDIM,
            &ALLOC(LWRKWK), &ALLOC(LWRKFI), IDIM2,
            &ALLOC(LOCPTRS.Z[IWRKST]), IRET, CLTIME);

        if (IRET != 0) goto L9960;

        //
        // Free the work arrays
        //
        LOCPTRS.Z[IWRKFF] = -LOCPTRS.Z[IWRKFF];
        LOCPTRS.Z[IWRKFI] = -LOCPTRS.Z[IWRKFI];
        LOCPTRS.Z[IWRKWK] = -LOCPTRS.Z[IWRKWK];
        LOCPTRS.Z[IWRKST] = -LOCPTRS.Z[IWRKST];

        // -----------------------------------------------------------------------
        // REESTABLISH ELASTIC PHASE SHIFTS
        // -----------------------------------------------------------------------

L800:
        // CALL DSGMAL ( ETAS(1), LOMOST+LXMAX, ALLOC(Z(ISIGS(1))) )
        DSGMAL(ETAS[1], LOMOST + LXMAX, &ALLOC(LOCPTRS.Z[ISIGS[1]]));

        T_TT4 = 0.0;  // T4 = second()
        T2 = 0.0;     // T2 = T4 - TT

        //
        // NOW FIND ALL THE REQUIRED F'S AND G'S
        //
        // I = 4*NBASDF*NUMLIS/FACFR4
        {
            I = 4 * NBASDF * NUMLIS / FACFR4;
            IFG = NALLOC("FANDG    ", I);
            NFS[1] = IFG;
            INIT4(IFG, 0.0f);

            L2P1 = LMAX + LXMAX + 1;
            IWRKFG2 = NALLOC("F&GWORK ", 6 * L2P1);

            // Re-read LCHNDF and LBASDF pointers
            LCHNDF2 = FACFR4 * (LOCPTRS.Z[ICHNDF_loc] - 1) + 1;
            LCHNDF2 = LCHNDF2 + ILLOC(LCHNDF2 + 2);

            LBASDF2 = FACFR4 * (LOCPTRS.Z[IBASDF_loc] - 1) + 1;
            LBASDF2 = LBASDF2 + ILLOC(LBASDF2 + 2);

            //
            // CALL SETFG ( NCHNDF, MCHNDF, ILLOC(LCHNDF),
            //   ALLOC4(LCHNDF), NBASDF, MBASDF, ILLOC(LBASDF),
            //   MCHNVL, ALLOC(Z(ICHNVL)),
            //   NUMLIS, ALLOC(Z(ILIS)), ALLOC(Z(IFG)),
            //   MAX0(0, LMIN-LXMAX), L2P1, ALLOC(Z(IWRKFG)),
            //   SUMMAX-NBACK*RSTEPS(1), SUMMAX, IRTN )
            //
            int IRTN_loc = 0;
            SETFG_full(
                NCHNDF, MCHNDF, &ILLOC(LCHNDF2), &ALLOC4(LCHNDF2),
                NBASDF, MBASDF, &ILLOC(LBASDF2),
                MCHNVL_loc, &ALLOC(LOCPTRS.Z[ICHNVL_loc]),
                NUMLIS, &ALLOC(LOCPTRS.Z[ILIS]),
                &ALLOC(LOCPTRS.Z[IFG]),
                MAX0(0, LMIN - LXMAX), L2P1,
                &ALLOC(LOCPTRS.Z[IWRKFG2]),
                SUMMAX - NBACK * RSTEPS[1], SUMMAX, IRTN_loc);

            IRTN = IRTN_loc;

            LOCPTRS.Z[IWRKFG2] = -LOCPTRS.Z[IWRKFG2];

            // Timing output
            double T4_val = 0.0;   // second() - T4 stub
            double T3_val = T2 - CLTIME;
            double TALL_val = T4_val + T2;
            std::printf("-MATCHING CONDITION TIMES (SECONDS):\n");
            std::printf(" BELLINGS:                    %8.2f\n", CLTIME);
            std::printf(" RECURSIONS AND COUL. PHASES: %8.2f\n", T3_val);
            std::printf(" TOTAL FOR B FACTOR:          %8.2f\n", T2);  // FORMAT T40
            std::printf(" F AND G FUNCTIONS:           %8.2f\n", T4_val);
            std::printf(" TOTAL:                       %8.2f\n", TALL_val);

            //
            // ALLOCATIONS FOR COUPLN
            //

            // HOMOGENEOUS SOLUTION STORAGE
            NWRDS = (2 * (NUMPT - NPTMIN) * NBASDF + FACFR4 - 1) / FACFR4;
            NWRDS = FACFR4 * NWRDS;
            IHOMO = NALLOC("CCHOMO  ", NWRDS);

            // PARTICULAR INTEGRAL STORAGE
            NWRDS = ((NUMPT - NPTMIN) * NBASDF + FACFR4 - 1) / FACFR4;
            NWRDS = FACFR4 * NWRDS;
            IINHR = NALLOC("CCINHR  ", NWRDS);
            IINHI = NALLOC("CCINHI  ", NWRDS);

            // R.H.S. (ONLY 1 FOR SEQUENTIAL ITERATION)
            IRHSR = NALLOC("CCRHSR  ", NUMPT);
            IRHSI = NALLOC("CCRHSI  ", NUMPT);
            for (I = 1; I <= NUMPT; I++) {
                ALLOC(LOCPTRS.Z[IRHSR] + I - 1) = 0.0;
                ALLOC(LOCPTRS.Z[IRHSI] + I - 1) = 0.0;
            }

            // SEQUENCE OF STORED S-MATRIX ELEMENTS
            NWRDS = 4 * NBASDF;
            ISMSQ = NALLOC("CCSMSQ  ", NWRDS);

            // ELEMENTS OF THE PADE TABLE
            NWRDS = 6 * ((MAXITR + 1) / 2) * NBASDF;
            IPADE = NALLOC("PADE    ", NWRDS);

            //
            // FULL-PRECISION HOMO/INHR8/INHI8 STORAGE
            //
            IHOMOA_loc = NALLOC("CCHOMO8A", 4 * NBASDF);
            IHOMOB_loc = NALLOC("CCHOMO8B", 4 * NBASDF);
            IINHR8_loc = NALLOC("CCINHR8 ", 2 * NBASDF);
            IINHI8_loc = NALLOC("CCINHI8 ", 2 * NBASDF);

            if (ISAVSM == 0) goto L990;

            {
                int I_loc2 = 1;
                if (ISAVSM == 3) I_loc2 = 2;
                ISMTSV = NALLOC("SMATITER", 2 * MAXITR * I_loc2 * NBASDF);
            }
        }
    }  // end L400 block

L990:
    IRTN = 1;
    return;

    // -----------------------------------------------------------------------
    // ERROR HANDLER
    // -----------------------------------------------------------------------

L9960:
    std::printf("0**** ERROR RETURN FROM COULIN:%12d\n", IRET);
    return;
}

// ---------------------------------------------------------------------------
// Weak implementations for the extern routines that COULST calls but that
// have not yet been translated.  These forward to the stubs already in the
// code base where the name matches; for genuinely new signatures we provide
// minimal stubs here so that this translation unit compiles and links cleanly.
// The actual subroutines are expected to be translated in later phases.
// ---------------------------------------------------------------------------

// Full 4-argument RTXLNX — implementation (translated from source.f L30926-30958)
double RTXLNX(double A, double B, double C, double ACC) {
    if (A == 0.0) {
        if (B == 0.0) return 0.0;
        return std::exp(-C / B);
    }
    double BP = B / A;
    double CP = C / A;
    if (B == 0.0) return -CP;
    double X = -CP;
    if (X < 0.1) X = 0.1;
    for (int Ii = 1; Ii <= 100; Ii++) {
        double XOLD  = X;
        double DELTA = -(X + BP * std::log(X) + CP) / (1.0 + BP / X);
        X = X + DELTA;
        if (X <= 0.0) X = 0.1 * XOLD;
        if (std::fabs(DELTA) <= ACC) goto L300_rtx;
    }
    std::printf("0**** COULD NOT CONVERGE IN RTXLNX%20.10G%20.10G%20.10G%20.10G\n",
                A, B, C, X);
L300_rtx:
    return X;
}

// COULIN is now fully implemented in coulin_translated.cpp.

// GENBNX_full — stub
void GENBNX_full(
    int IPASS, int& NBINDX, int& MBINDX, int* BINDEX,
    int NBASCP, int MBASCP, int* BASCUP,
    int MBASDF, int* BASDEF,
    int MCHNDF, float* CHNDEF,
    int& NMBFAC, int& NMFFAC, int& MXLXGS, int& MXLXBF,
    int COULSW, int PBUGSW)
{
    (void)IPASS; (void)BINDEX;
    (void)NBASCP; (void)MBASCP; (void)BASCUP;
    (void)MBASDF; (void)BASDEF;
    (void)MCHNDF; (void)CHNDEF;
    (void)COULSW; (void)PBUGSW;
    NBINDX = 0; MBINDX = 9;
    NMBFAC = 0; NMFFAC = 0; MXLXGS = 0;
    // MXLXBF is an in/out — leave as is
    std::printf(" GENBNX_full: stub — translate GENBNX from source.f L17285\n");
}

// SETBFC_full — stub
void SETBFC_full(
    int NBINDX, int MBINDX, int* BINDEX,
    int MCHNDF, int* CHNDEF_i, float* RCHNDF,
    int NUMLIS, int* LIS, double* R2S_arg, int LMXMX,
    double* SIG1, double* SIG2,
    double RLOWER, int ALLSW, int IRORC,
    int NMBFAC, float* BFAC,
    double* FFWORK, int IDIM1, int LDLDIM,
    double* WRKWK, double* FIWORK, int IDIM2,
    double* WRKST, int& IRET, double& CLTIME)
{
    (void)NBINDX; (void)MBINDX; (void)BINDEX;
    (void)MCHNDF; (void)CHNDEF_i; (void)RCHNDF;
    (void)NUMLIS; (void)LIS; (void)R2S_arg; (void)LMXMX;
    (void)SIG1; (void)SIG2;
    (void)RLOWER; (void)ALLSW; (void)IRORC;
    (void)NMBFAC; (void)BFAC;
    (void)FFWORK; (void)IDIM1; (void)LDLDIM;
    (void)WRKWK; (void)FIWORK; (void)IDIM2;
    (void)WRKST;
    IRET   = 0;
    CLTIME = 0.0;
    std::printf(" SETBFC_full: stub — translate SETBFC from source.f L31512\n");
}

// SETFG_full — stub
void SETFG_full(
    int NCHNDF, int MCHNDF, int* CHNDEF, float* RCHNDF,
    int NBASDF, int MBASDF, int* BASDEF,
    int MCHNVL, double* CHNVAL,
    int NUMLIS, double* LIS,
    double* FG,
    int L1, int L2P1,
    double* FGWORK, double R1, double R2, int& IRET)
{
    (void)NCHNDF; (void)MCHNDF; (void)CHNDEF; (void)RCHNDF;
    (void)NBASDF; (void)MBASDF; (void)BASDEF;
    (void)MCHNVL; (void)CHNVAL;
    (void)NUMLIS; (void)LIS;
    (void)FG;
    (void)L1; (void)L2P1;
    (void)FGWORK; (void)R1; (void)R2;
    IRET = 1;
    std::printf(" SETFG_full: stub — translate SETFG from source.f L32481\n");
}
