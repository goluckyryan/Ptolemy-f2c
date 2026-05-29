// lintrp_translated.cpp — LINTRP
// Translated from source.f lines 23511-24796
// Interpolation and extrapolation in L of reaction S-matrices

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include "ptolemy_alloc.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_forward.h"
#include "ptolemy_io.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <complex>
#include <algorithm>

extern double second();

void LINTRP()
{
    // Local references to COMMON blocks
    auto* Z      = LOCPTRS.Z;
    auto& FACFR4 = ALLOCS.FACFR4;
    auto& PI     = CNSTNT.PI;
    auto& BIGNUM = CNSTNT.BIGNUM;
    auto& SMLNUM = CNSTNT.SMLNUM;
    auto& BIGLOG = CNSTNT.BIGLOG;
    auto& LMIN   = INTGER.LMIN;
    auto& LMAX   = INTGER.LMAX;
    auto& LSTEP  = INTGER.LSTEP;
    auto& MAXLEX = INTGER.MAXLEX;
    auto& IPRINT = INTGER.IPRINT;
    auto& LEBACK = INTGER.LEBACK;
    auto& PROBLM = SWITCH.PROBLM;
    auto& ISAVMX = SWITCH.ISAVMX;
    auto& ISAVB  = SWITCH.ISAVB;
    auto& IEXTYP = SWITCH.IEXTYP;
    auto& IBRNSB = SWITCH.IBRNSB;
    auto& ISTRIP = INTRNL.ISTRIP;
    auto& IHSAVE = INTRNL.IHSAVE;
    auto& LIMOST = INTRNL.LIMOST;
    auto& UNDEF  = INTRNL.UNDEF;
    auto& AKI    = KANDM.AKI;
    auto& AKO    = KANDM.AKO;
    auto& LOMOST = KANDM.LOMOST;
    auto& LXMIN  = INELCM.LXMIN;
    auto& LXMAX  = INELCM.LXMAX;
    auto& LXSTEP = INELCM.LXSTEP;
    auto& NSPLI  = INELCM.NSPLI;
    auto& ITOCS  = INELCM.ITOCS;
    auto& IRDINT = INELCM.IRDINT;
    auto& NUMLIS = INELCM.NUMLIS;
    auto& ILIS   = INELCM.ILIS;
    auto& NASPLI = INELCM.NASPLI;
    auto& IUNITR = INELCM.IUNITR;
    auto& LIFIT  = INELCM.LIFIT;
    auto& ILIFIT = INELCM.ILIFIT;
    auto& NOLFIT = INELCM.NOLFIT;
    auto& ISMATR = INELCM.ISMATR;
    auto& ISMATI = INELCM.ISMATI;
    auto& IDELSR = INELCM.IDELSR;
    auto& IDELSI = INELCM.IDELSI;
    auto& IDENSW = INELCM.IDENSW;

    // GRIDCM fields used by LINTRP (different overlay than GRDSET)
    // Access raw GRIDCM members via the existing struct
    // IBINDX maps to GRIDCM field at position 4 (after JACOB, STPINV, COSSTP)
    // For LINTRP's GRIDCM overlay, use reinterpret_cast
    int* GRIDCM_int = reinterpret_cast<int*>(&GRIDCM.IMSVAL);
    // LINTRP GRIDCM overlay: IBINDX=pos0, MBINDX=pos1, NBINDX=pos2, NMFFAC=pos3
    // These map to IMSVAL, IINTS, INTOFF, IA12VL in GRDSET's view
    int& IBINDX_l = GRIDCM_int[0];  // = GRIDCM.IMSVAL
    int& MBINDX_l = GRIDCM_int[1];  // = GRIDCM.IINTS
    int& NBINDX_l = GRIDCM_int[2];  // = GRIDCM.INTOFF
    int& NMFFAC_l = GRIDCM_int[3];  // = GRIDCM.IA12VL

    int ICL2FF_l = INELCM.ICL2FF;

    // Local variables
    int MWORK = 12;
    double TSTART, AP, BP, CP;
    bool TRANSW, XTRPSW, CCSW, PBUGSW;
    int IPRNT, LIPREV;
    int ICL, ICINTS, IWORK, ILVALS, IXINTS;
    int NCHNDF, KSMAT, KSMAT2, JMOST, MCHNDF, LMAXL;
    int MCHN, LNSMAT, LILOSZ;
    int LTOCS, LLIS, LWORK, LLVALS, LLVAL2, LXINTS, LXINT2;
    int LCL, LCINTS, LELIN, LELOUT, LUNITR_l, LRDINT_l;
    int LSMATR, LSMATI, LDELSR_l, LDELSI_l, LISTRT;
    int LX, LDEL, JP, JT, KOFFS, LPEAK;
    double BIGBOY, DLMAX, FACTOR, ELFIN;
    int LBRAT;
    bool BIGRSW;

    // Dummy arrays for INTRCF
    double DUMMY3[1];
    float DUMMY4[1];

    // Extrapolation type names
    const char* XWORDS[6][2] = {
        {"WOODS-SA", "XON"},
        {"POWER-LA", "W"},
        {"POWER-LA", "W 2"},
        {"WKB POWE", "R-LAW"},
        {"LEXTRAP4", " "},
        {"LEXTRAP5", " "}
    };

    TSTART = second();

    // Define things that might never get defined
    AP = 0; BP = 0; CP = 0;

    TRANSW = (ISTRIP != 0);
    IPRNT = IPRINT % 10;
    PBUGSW = ((IPRINT / 100) % 10) >= 4;
    CCSW = (PROBLM == 24);

    XTRPSW = (MAXLEX > 0);

    // Print header
    std::printf("1%58sP T O L E M Y\n"
                "%47sINTERPOLATION AND EXTRAPOLATION IN L\n"
                "0%.45sELAB =%7.2f MEV     %.65s\n\n",
                "", "",
                &HEDCOM.REACT[1], FLOAT_common.ELAB, &HEDCOM.HEADER[1]);

    int iex = IEXTYP - 1;
    if (iex < 0 || iex >= 6) iex = 0;
    std::printf("0USING %-8s%-8sEXTRAPOLATION FUNCTIONS.\n",
                XWORDS[iex][0], XWORDS[iex][1]);

    DLMAX = LMAX;
    LIPREV = LIMOST;
    if (TRANSW) LIMOST = LMAX;

    // These 2 work areas are used for continued-fraction interpolation.
    ICL = NALLOC("CMPLXLS ", 2 * NUMLIS);
    ICINTS = NALLOC("CMPLXIS ", 2 * NUMLIS);

    if (!XTRPSW) goto L90;

    // Extrapolation work arrays
    IWORK = NALLOC("LWORK   ", MWORK * NASPLI);
    ILVALS = NALLOC("FITLS   ", 2 * NOLFIT);
    IXINTS = NALLOC("FITIS   ", 2 * NOLFIT);
    INIT8(IWORK, 0.0);

    // Allocation is finished for now
    LTOCS = Z[ITOCS] * FACFR4 - FACFR4;
    LLIS = FACFR4 * Z[ILIS] - FACFR4;
    LWORK = Z[IWORK] - MWORK;
    LLVALS = Z[ILVALS] - 1;
    LLVAL2 = LLVALS + NOLFIT;
    LXINTS = Z[IXINTS] - 1;
    LXINT2 = LXINTS + NOLFIT;

    std::printf("\nEXTRAPOLATION AND CONTINUED FRACTION OVERLAP"
                " REGION:%4d < LI <%4d\n\n", LIFIT, LMAX);

    // ========================================================================
    // Possibly setup for multiple channels
    // ========================================================================
L90:
    NCHNDF = 2;
    KSMAT = 0;
    KSMAT2 = 0;
    JMOST = 0;
    int LCHNDF_cc = 0, LTOCSV = 0, LBINDX_cc = 0, LCHNSV = 0;
    if (CCSW) {
        LCHNDF_cc = FACFR4 * (Z[CCBLK.ICHNDF] - 1) + 1;
        NCHNDF = ILLOC(LCHNDF_cc);
        MCHNDF = ILLOC(LCHNDF_cc + 1);
        LCHNDF_cc = LCHNDF_cc + ILLOC(LCHNDF_cc + 2) + MCHNDF;
        LTOCSV = LTOCS;
        LBINDX_cc = FACFR4 * (Z[IBINDX_l] - 1);
    }

    // ========================================================================
    // PASS 1: Determine LIMOST and size requirements
    // ========================================================================
    for (int NCHN = 2; NCHN <= NCHNDF; NCHN++) {
        MCHN = NCHN;
        LMAXL = LMAX;
        if (CCSW) {
            MCHN = ILLOC(LCHNDF_cc);
            KANDM.ETAS[2] = (double)ALLOC4(LCHNDF_cc + 5);
            KSMAT = ILLOC(LCHNDF_cc + 9);
            NSPLI = ILLOC(LCHNDF_cc + 10);
            LTOCS = LTOCSV + ILLOC(LCHNDF_cc + 12);
            LIMOST = ILLOC(LCHNDF_cc + 14);
            if (LIMOST == NOTDEF_INT) LIMOST = LMAX;
        }

        // L110: if (!XTRPSW) goto L290
        if (!XTRPSW) goto L290;

        // Pass 1 extrapolation fitting loop (XTRPSW=true path)
        if (!CCSW) LTOCS = Z[ITOCS] * FACFR4 - FACFR4;
        LLIS = FACFR4 * Z[ILIS] - FACFR4;
        for (KOFFS = 1; KOFFS <= NSPLI; KOFFS++) {
            JT = ILLOC(LTOCS + 4 * KOFFS);
            if (JT < 0) continue;
            JP = ILLOC(LTOCS + 4 * KOFFS - 1);
            LX = ILLOC(LTOCS + 4 * KOFFS - 2);
            LDEL = ILLOC(LTOCS + 4 * KOFFS - 3);
            LDELSR_l = Z[IDELSR] + (KOFFS + KSMAT - 1);
            LDELSI_l = Z[IDELSI] + (KOFFS + KSMAT - 1);

            // Locate L-peak and check for decaying form
            int LPEAK_l = LMIN;
            BIGBOY = 0;
            bool BADSW = false;
            double SIZLAS = 0, SIGNRL = 0, SIGNIL = 0;

            for (int ILI = 1; ILI <= NUMLIS; ILI++) {
                int LI = ILLOC(LLIS + ILI);
                int LO = LI + LDEL;
                if (LO + LI < LX) continue;
                int IADD = (ILI - 1) * NSPLI;
                double RE = ALLOC(LDELSR_l + IADD);
                double HIM = ALLOC(LDELSI_l + IADD);
                double SIZE = RE * RE + HIM * HIM;
                if (LI >= LIFIT) {
                    double SIGNR = copysign(1.0, RE);
                    double SIGNI = copysign(1.0, HIM);
                    if (LI > LIFIT) {
                        if (SIZLAS > SIZE) {
                            // OK - decaying
                        } else {
                            BADSW = true;
                            double B_t = sqrt(SIZE);
                            double C_t = sqrt(SIZLAS);
                            std::printf("\n**** ERROR:  CANNOT EXTRAPOLATE FOR "
                                        "CHANNEL%3d  "
                                        "(JP,JT,LX,LO-LI) = (%3d/2%3d/2%2d%3d) ****\n"
                                        " ****%9sSINCE FOR LI =%4d"
                                        " |S| EXCEEDS PREVIOUS:%11.4G%11.4G ****\n",
                                        MCHN, JP, JT, LX, LDEL, "", LI, B_t, C_t);
                        }
                        if (SIGNR != SIGNRL || SIGNI != SIGNIL) {
                            std::printf("0**** WARNING:  PHASE OF S(L) FOR "
                                        "CHANNEL%3d  "
                                        "(JP,JT,LX,LO-LI,LI) = (%3d/2%3d/2%3d %3d %3d"
                                        ") HAS FLUCTUATING SIGN\n",
                                        MCHN, JP, JT, LX, LDEL, LI);
                        }
                    }
                    SIZLAS = SIZE;
                    SIGNRL = SIGNR;
                    SIGNIL = SIGNI;
                }
                SIZE = (2 * LI + 1) * SIZE;
                if (SIZE >= BIGBOY) {
                    LPEAK_l = LI;
                    BIGBOY = SIZE;
                }
            }

            if (BADSW) continue;

            // Find WEEBOY
            double WEEBOY = 1.0e+06;
            for (int ILI = 1; ILI <= NUMLIS; ILI++) {
                int LI = ILLOC(LLIS + ILI);
                int LO = LI + LDEL;
                if (LO + LI < LX) continue;
                if (LI > LPEAK_l) continue;
                int IADD = (ILI - 1) * NSPLI;
                double SIZE = ALLOC(LDELSR_l + IADD) * ALLOC(LDELSR_l + IADD)
                            + ALLOC(LDELSI_l + IADD) * ALLOC(LDELSI_l + IADD);
                SIZE = (2 * LI + 1) * SIZE;
                if (SIZE < WEEBOY) WEEBOY = SIZE;
            }

            WEEBOY = 0.2 * sqrt(WEEBOY);
            BIGBOY = FLOAT_common.DWCUT * sqrt(BIGBOY);
            WEEBOY = std::min(WEEBOY, BIGBOY);

            // Fit extrapolation function
            for (int ILI = ILIFIT; ILI <= NUMLIS; ILI++) {
                int I = ILI - ILIFIT + 1;
                ALLOC(LLVALS + I) = ILLOC(LLIS + ILI);
                ALLOC(LLVAL2 + I) = ILLOC(LLIS + ILI);
                int IADD = (ILI - 1) * NSPLI;
                double F = sqrt(ALLOC(LDELSR_l + IADD) * ALLOC(LDELSR_l + IADD)
                              + ALLOC(LDELSI_l + IADD) * ALLOC(LDELSI_l + IADD));
                ALLOC(LXINTS + I) = F;
                F = atan2(ALLOC(LDELSI_l + IADD), ALLOC(LDELSR_l + IADD));
                ALLOC(LXINT2 + I) = F;
            }

            double AA, WID, FLCRIT, BARL_l, BARA_l, B_l, BARC_l, CHISQ;
            int ISIG;
            LXTRP1(IEXTYP, NOLFIT, ISIG, IPRNT, ALLOC_base(LLVALS + 1),
                   ALLOC_base(LXINTS + 1), FLCRIT, AA, WID, 0.30, 1.0e-5,
                   BARL_l, BARA_l, B_l, BARC_l, DLMAX, CHISQ, LX, LDEL,
                   &KANDM.ETAS[0], MCHN, JP, JT);

            if (ISIG < 0) {
                std::printf("\n  **** ERROR RETURN FROM L-EXTRAP SUBROUTINE\n"
                            "   CHANNEL%3d  "
                            "       EXTRAPOLATION SUPPRESSED FOR (JP,JT,LX,LO-LI)=("
                            "%3d/2%3d/2%4d%4d) \n", MCHN, JP, JT, LX, LDEL);
                continue;
            }
            if (ISIG != 0) {
                std::printf("\n  SIMPLE EXP OR POWER  USED FOR |S(L)| "
                            "   CHANNEL%3d  "
                            "(JP,JT,LX,LO-LI)=(%3d/2%3d/2%4d%4d)\n",
                            MCHN, JP, JT, LX, LDEL);
            }

            // Phase extrapolation for power-law
            if (IEXTYP >= 2) {
                LXTRP1(3, NOLFIT, ISIG, IPRNT, ALLOC_base(LLVAL2 + 1),
                       ALLOC_base(LXINT2 + 1), FLCRIT, AA, WID, 0.30, 1.0e-5,
                       BARL_l, CP, BP, AP, DLMAX, CHISQ, LX, LDEL,
                       &KANDM.ETAS[0], MCHN, JP, JT);
                if (ISIG < 0) {
                    std::printf("\n**** FOLLOWING REFERS TO PHASE EXTRAP.\n");
                    continue;
                }
            }

            // Find LIMOST for this (LX, LDEL)
            int LI = (int)LXTRPM(IEXTYP, BARA_l, B_l, BARC_l, BARL_l, DLMAX,
                                  LX, LDEL, &KANDM.ETAS[0], WEEBOY);
            if (PBUGSW) {
                std::printf(" FOR CHANNEL%3d    FOR 2*JP, 2*JT, LX, LO-LI =%3d%3d%3d%3d   "
                            "ESTIMATED LARGEST SIGNIFICANT LI IS%4d\n",
                            MCHN, JP, JT, LX, LDEL, LI);
            }
            LMAXL = std::max(LMAXL, LI);

            // Save extrapolation parameters
            int LAREA = LWORK + MWORK * KOFFS;
            ALLOC(LAREA)     = AA;
            ALLOC(LAREA + 1) = B_l;
            ALLOC(LAREA + 2) = BARL_l;
            ALLOC(LAREA + 3) = FLCRIT;
            ALLOC(LAREA + 4) = AP;
            ALLOC(LAREA + 5) = CP;
            ALLOC(LAREA + 6) = BARA_l;
            ALLOC(LAREA + 7) = BARC_l;
            ALLOC(LAREA + 8) = LI;
            ALLOC(LAREA + 11) = BP;
        } // end KOFFS loop

L290:
        // Compute size requirements
        if (LIMOST <= LMAX) LIMOST = std::min(LMAXL, LMAX + MAXLEX);
        JMOST = std::max(JMOST, LIMOST);
        LNSMAT = NSPLI * (LIMOST - LMIN + 1);

        if (CCSW) {
            ILLOC(LCHNDF_cc + 13) = KSMAT2;
            ILLOC(LCHNDF_cc + 14) = LIMOST;
            LCHNDF_cc = LCHNDF_cc + MCHNDF;
            LWORK = Z[IWORK] - MWORK + MWORK * NSPLI;
            // Note: LWORK advancement handled below by not resetting it
        }

        KSMAT2 = KSMAT2 + LNSMAT;
    } // end NCHN loop (pass 1)

    LILOSZ = KSMAT2;
    INELCM.LILOSZ = LILOSZ;

    // Allocate space for reaction S-matrices to XSECTN
    int I1 = 2 * LILOSZ;
    IRDINT = NALLOC("S       ", I1 / FACFR4);
    IUNITR = NALLOC("UNITARIT", (JMOST - LMIN + 2) / FACFR4);
    INIT4(IUNITR, 0.0f);
    INIT4(IRDINT, 0.0f);

    BIGRSW = TRANSW || (JMOST > LIPREV);

    // Extend elastic S-matrix arrays if needed
    LOMOST = JMOST + LXMAX;
    if (JMOST > LMAX) {
        // Elastic S-matrix extension (Born approximation)
        if (CCSW) WAVCOM.NUMJS[2] = WAVCOM.NUMJS[1];

        for (int I = 1; I <= 2; I++) {
            if (BIGRSW) IREDEF(LOMOST + 1, KANDM.ISIGS[I]);
            int N = 2 * WAVCOM.NUMJS[I] * (LMAX + 1) + 1;
            int NN = 2 * WAVCOM.NUMJS[I] * (LOMOST + 1);
            IREDEF(NN, KANDM.ISMATS[I]);
            int LSMAT = Z[KANDM.ISMATS[I]];
            int LSIG = Z[KANDM.ISIGS[I]];
            for (int II = N; II <= NN; II++) {
                ALLOC(LSMAT + II - 1) = 0;
            }

            if (I == 2 && CCSW) continue;
            double CPH = ALLOC(LSIG + LMAX);
            double E2 = sqrt(KANDM.ETAS[I] * KANDM.ETAS[I] + DLMAX * (DLMAX + 1));
            double EXRFAC = 1.0e10;
            if (KANDM.ASCTS[I] != UNDEF) EXRFAC = 1.0 / (KANDM.AKI * KANDM.ASCTS[I]);
            // Use AKS(I) — but AKI=AKS(1), AKO=AKS(2)
            double AKS_I = (I == 1) ? AKI : AKO;
            if (KANDM.ASCTS[I] != UNDEF) EXRFAC = 1.0 / (AKS_I * KANDM.ASCTS[I]);
            double EXIFAC = 1.0e10;
            if (WAVCOM.AIS[I] != UNDEF) EXIFAC = 1.0 / (AKS_I * WAVCOM.AIS[I]);
            if (WAVCOM.ASIS[I] != UNDEF) EXIFAC = 1.0 / (AKS_I * WAVCOM.ASIS[I]);
            bool BADSW = (EXRFAC == 1.0e10 || EXIFAC == 1.0e10);
            if (BADSW) {
                std::printf("0**** WARNING:  ELASTIC S MATRIX FOR CHANNEL"
                            "%2d NOT EXTRAPOLATED BECAUSE A OR AI IS"
                            " NOT DEFINED.\n", I);
            }
            double DTRE = 0, DTIM = 0, SIZE_v;
            if (!BADSW) {
                int LOC = LSMAT + 2 * LMAX * WAVCOM.NUMJS[I];
                SIZE_v = sqrt(ALLOC(LOC) * ALLOC(LOC) + ALLOC(LOC + 1) * ALLOC(LOC + 1));
                DTIM = -0.5 * log(SIZE_v);
                DTRE = 0.5 * atan2(ALLOC(LOC + 1), ALLOC(LOC));
            }

            int LO_start = LIFIT - LXMAX;
            double RE_v = 1, HIM_v = 0;
            for (int L = LO_start; L <= LOMOST; L++) {
                int LOC = LSMAT + 2 * WAVCOM.NUMJS[I] * L;
                if (BADSW) goto L380;
                {
                    double E1 = sqrt(KANDM.ETAS[I] * KANDM.ETAS[I] + (double)L * (L + 1));
                    double DRO = E1 - E2;
                    double EXR = 0;
                    double X = DRO * EXRFAC;
                    if (X < BIGLOG) EXR = exp(-X);
                    double EXI = 0;
                    X = DRO * EXIFAC;
                    if (X < BIGLOG) EXI = exp(-X);
                    double DRE = EXR * DTRE;
                    double DIM = EXI * DTIM;
                    SIZE_v = exp(-2.0 * DIM);
                    RE_v = SIZE_v * cos(2.0 * DRE);
                    HIM_v = SIZE_v * sin(2.0 * DRE);
                    if (L <= LMAX) {
                        if (fabs(ALLOC(LOC)) + fabs(ALLOC(LOC + 1)) != 0) {
                            // Check extrapolation quality (diagnostic only)
                        }
                        goto L389;
                    }
                }
L380:           ALLOC(LOC) = RE_v;
                ALLOC(LOC + 1) = HIM_v;
                if (L > LMAX && BIGRSW) {
                    CPH = CPH + atan(KANDM.ETAS[I] / (double)L);
                    ALLOC(LSIG + L) = CPH;
                }
L389:           ;
            }
        } // end I=1,2 loop
    } // end if JMOST > LMAX

    // ========================================================================
    // PASS 2: Interpolation and extrapolation
    // ========================================================================
L500:
    SETLOG(2 * (LOMOST + LXMAX));
    LCL = Z[ICL];
    LCINTS = Z[ICINTS];
    LLIS = FACFR4 * Z[ILIS] - FACFR4;
    LTOCS = Z[ITOCS] * FACFR4 - FACFR4;
    LELIN = Z[KANDM.ISMATS[1]];
    LELOUT = Z[KANDM.ISMATS[2]];
    LUNITR_l = FACFR4 * (Z[IUNITR] - 1) + 1 - LMIN;

    if (TRANSW || CCSW) goto L680;
    FACTOR = sqrt(AKI * AKO / (WAVCOM.ES[1] * WAVCOM.ES[2] * PI));
    LBRAT = Z[INELCM.IBETAS[2]] - LXMIN / 2;

L680:
    if (!XTRPSW) goto L700;
    LLVALS = Z[ILVALS] - 1;
    LXINTS = Z[IXINTS] - 1;
    LWORK = Z[IWORK] - MWORK;

    // Phase extrapolation setup
    {
        int I = LELIN + 2 * WAVCOM.NUMJS[1] * LMAX;
        ELFIN = 0.5 * atan2(ALLOC(I + 1), ALLOC(I));
    }

L700:
    if (!TRANSW) DUMMY1();

    KSMAT2 = 0;
    if (CCSW) {
        LCHNDF_cc = FACFR4 * (Z[CCBLK.ICHNDF] - 1) + 1;
        LCHNDF_cc = LCHNDF_cc + ILLOC(LCHNDF_cc + 2) + MCHNDF;
        LCHNSV = LCHNDF_cc;
        LTOCSV = LTOCS;
        LBINDX_cc = FACFR4 * (Z[IBINDX_l] - 1);
        LELOUT = LELIN;
        goto L710;
    }

L710:
    for (int NCHN = 2; NCHN <= NCHNDF; NCHN++) {
        MCHN = NCHN;
        if (CCSW) {
            MCHN = ILLOC(LCHNDF_cc);
            KANDM.ETAS[2] = (double)ALLOC4(LCHNDF_cc + 5);
            KSMAT = ILLOC(LCHNDF_cc + 9);
            NSPLI = ILLOC(LCHNDF_cc + 10);
            LTOCS = LTOCSV + ILLOC(LCHNDF_cc + 12);
            KSMAT2 = ILLOC(LCHNDF_cc + 13);
            LIMOST = ILLOC(LCHNDF_cc + 14);
        }

        // L720: loop over (JT, JP, LX, LO-LI)
        for (KOFFS = 1; KOFFS <= NSPLI; KOFFS++) {
            JT = ILLOC(LTOCS + 4 * KOFFS);
            if (JT < 0) continue;
            JP = ILLOC(LTOCS + 4 * KOFFS - 1);
            LX = ILLOC(LTOCS + 4 * KOFFS - 2);
            LDEL = ILLOC(LTOCS + 4 * KOFFS - 3);
            LSMATR = Z[ISMATR] + (KOFFS + KSMAT - 1);
            LSMATI = Z[ISMATI] + (KOFFS + KSMAT - 1);
            LDELSR_l = Z[IDELSR] + (KOFFS + KSMAT - 1);
            LDELSI_l = Z[IDELSI] + (KOFFS + KSMAT - 1);
            LRDINT_l = FACFR4 * (Z[IRDINT] - 1) + 1 + 2 * (KOFFS - 1 + KSMAT2);
            LISTRT = std::max(LMIN, (LX - LDEL + 1) / 2);

            // Continued fraction interpolation
            INTRCF(MCHN, JP, JT, LX, LDEL,
                   NUMLIS, ILLOC_base(LLIS + 1), LSTEP, LXMAX,
                   &ALLOC(LSMATR), &ALLOC(LSMATI), 1, NSPLI, 0 /*false*/,
                   LISTRT, LMAX, DUMMY3, &ALLOC4(LRDINT_l),
                   1, NSPLI, LMIN - 1, 1 /*true*/,
                   BIGBOY, LPEAK, &ALLOC(LCL), &ALLOC(LCINTS));

            if (TRANSW && ((LPEAK == LMIN && LMIN > 0) || LPEAK >= LMAX)) {
                std::printf("0**** WARNING:  FOR (2*JP,2*JT,LX,LO-LI) = (%3d%3d%3d%3d"
                            ")    LPEAK = %3d IS AT END OF L-RANGE (FISHY) ***\n",
                            JP, JT, LX, LDEL, LPEAK);
            }

            if (!XTRPSW) goto L950;

            // Extrapolation (XTRPSW=true path)
            {
                int LAREA = LWORK + MWORK * KOFFS;
                double B_l = ALLOC(LAREA + 1);
                double BARL_l = ALLOC(LAREA + 2);
                double BARA_l = ALLOC(LAREA + 6);
                double BARC_l = ALLOC(LAREA + 7);
                ALLOC(LAREA + 9) = LPEAK;
                ALLOC(LAREA + 10) = sqrt(BIGBOY);

                if (BARA_l == 0) goto L950;

                AP = ALLOC(LAREA + 4);
                CP = ALLOC(LAREA + 5);
                BP = ALLOC(LAREA + 11);

                if (IEXTYP >= 2) goto L890;

                // Phase extrapolation using elastic phases
                {
                    int LO = LMAX + LDEL;
                    int IADD = (NUMLIS - 1) * NSPLI;
                    double PHASE = atan2(ALLOC(LDELSI_l + IADD), ALLOC(LDELSR_l + IADD));
                    CP = 0;
                    if (fabs(PHASE) > 0.25 * PI) CP = copysign(0.5 * PI, PHASE);
                    if (fabs(PHASE) > 0.75 * PI) CP = copysign(PI, PHASE);
                    PHASE = PHASE - CP;
                    int I = LELOUT + 2 * WAVCOM.NUMJS[2] * LO;
                    double ELFOUT = 0.5 * atan2(ALLOC(I + 1), ALLOC(I));
                    double ELPH = ELFOUT + ELFIN;
                    AP = PHASE / ELPH;
                    ALLOC(LAREA + 4) = AP;
                    ALLOC(LAREA + 5) = CP;
                }

L890:
                if (PBUGSW) {
                    std::printf("\n  * EXTRAPOLATED S-MATRICES WITH  CHANNEL =%3d   "
                                "(JP,JT,LX,LO-LI) = (%3d/2%3d/2%3d %3d) *\n",
                                MCHN, JP, JT, LX, LDEL);
                }

                for (int LI = LIFIT; LI <= LIMOST; LI++) {
                    int LO = LI + LDEL;
                    double SIZE_l;
                    LXTRP2(IEXTYP, BARA_l, B_l, BARC_l, BARL_l, DLMAX, LI,
                           SIZE_l, LX, LDEL, &KANDM.ETAS[0]);
                    if (SIZE_l < SMLNUM) goto L950;

                    double PHASE, RE_l, HIM_l;
                    int I = LELOUT + 2 * WAVCOM.NUMJS[2] * LO;
                    int II = LELIN + 2 * WAVCOM.NUMJS[1] * LI;
                    if (IEXTYP >= 2) {
                        double PHASE_tmp;
                        LXTRP2(3, AP, BP, CP, DLMAX, DLMAX, LI,
                               PHASE_tmp, LX, LDEL, &KANDM.ETAS[0]);
                        PHASE = PHASE_tmp;
                    } else {
                        double ELPH = 0.5 * (atan2(ALLOC(I + 1), ALLOC(I))
                                           + atan2(ALLOC(II + 1), ALLOC(II)));
                        PHASE = CP + AP * ELPH;
                    }
                    RE_l = SIZE_l * cos(PHASE);
                    HIM_l = SIZE_l * sin(PHASE);
                    if (fabs(RE_l) < 1.0e-30) RE_l = 0;
                    if (fabs(HIM_l) < 1.0e-30) HIM_l = 0;

                    // In overlap region, don't overwrite
                    if (LI <= LMAX) continue;

                    int IADD = 2 * (LI - LMIN) * NSPLI;
                    ALLOC4(LRDINT_l + IADD) = (float)RE_l;
                    ALLOC4(LRDINT_l + IADD + 1) = (float)HIM_l;
                }
            }

            // Add Coulomb excitation for inelastic scattering
L950:
            if (TRANSW) goto L980;
            if (CCSW) goto L960;

            {
                double FACTR = FACTOR;
                if (IDENSW) FACTR = sqrt(2.0) * FACTR;
                if ((LX % 4) < 2) FACTR = -FACTR;
                bool IMAGSW = (LX % 2) == 0;
                int LCL2FF_v = Z[ICL2FF_l] + KOFFS - 1;

                for (int LI = LISTRT; LI <= LIMOST; LI++) {
                    int LO = LI + LDEL;
                    int IADD = (LI - LMIN) * NSPLI;
                    double C = FACTR * ALLOC(LBRAT + LX / 2)
                             * fabs(CLEBSH(2 * LI, 2 * LX, 0, 0, 2 * LO, 0));
                    double FFI = C * ALLOC(LCL2FF_v + IADD);
                    if (!IMAGSW)
                        ALLOC4(LRDINT_l + 2 * IADD) = ALLOC4(LRDINT_l + 2 * IADD) + (float)FFI;
                    if (IMAGSW)
                        ALLOC4(LRDINT_l + 2 * IADD + 1) = ALLOC4(LRDINT_l + 2 * IADD + 1) + (float)FFI;
                }
            }
            goto L980;

L960:
            if (NMFFAC_l == 0 || !XTRPSW) goto L980;
            {
                int LMAXP1 = LMAX + 1;
                int LCL2FF = FACFR4 * (Z[ICL2FF_l] - 1);
                for (int IB = 1; IB <= NBINDX_l; IB++) {
                    int LBIN = LBINDX_cc + (IB - 1) * MBINDX_l;
                    int LB = ILLOC(LBIN + 9);
                    if (LB == 0) continue;
                    if (ILLOC(LBIN + 1) != LX || ILLOC(LBIN + 3) != NCHN) continue;
                    int NMBFRJ = (LX + 1) * ((ILLOC(LBIN + 4) - ILLOC(LBIN + 6)) / 4 + 1);
                    int INDX1 = LB + ((LDEL + LX) >> 1);
                    bool oddLX = (LX & 1) != 0;
                    for (int LI = LMAXP1; LI <= LIMOST; LI++) {
                        int I = LRDINT_l + 2 * NSPLI * (LI - LMIN);
                        float SBORN = ALLOC4(LCL2FF + INDX1 + NMBFRJ * (LI - LMIN));
                        if (oddLX)
                            ALLOC4(I) = ALLOC4(I) + SBORN;
                        else
                            ALLOC4(I + 1) = ALLOC4(I + 1) + SBORN;
                    }
                }
            }

            // Unitarity contributions
L980:
            for (int LI = LISTRT; LI <= LIMOST; LI++) {
                int I = LRDINT_l + 2 * NSPLI * (LI - LMIN);
                ALLOC4(LUNITR_l + LI) = ALLOC4(LUNITR_l + LI)
                    + ALLOC4(I) * ALLOC4(I) + ALLOC4(I + 1) * ALLOC4(I + 1);
            }
        } // end KOFFS loop (989)

        if (CCSW) {
            LCHNDF_cc = LCHNDF_cc + MCHNDF;
            LWORK = LWORK + MWORK * NSPLI;
        }
    } // end NCHN loop (999)

    // Skip pass 3 (CC only)
    if (!CCSW) goto L2800;
    // Pass 3 would go here for coupled channels
    // ...

L2800:
    if (!XTRPSW) goto L3000;

    // Print extrapolation summary
    if (IBRNSB == 0 && !TRANSW) {
        std::printf("\n\n%32sFOLLOWING REFER TO THE \"NUCLEAR\" PART OF THE INELASTIC "
                    "AMPLITUDES.\n", "");
    }
    std::printf("\n\n%48sSUMMARY OF EXTRAPOLATION PARAMETERS\n", "");
    std::printf("\n%19s|S| = A/( C + EXP( (LI-LCRIT)/LWIDTH ) )%23s"
                "PHASE = CP + AP*(SUM OF ELASTIC PHASE-SHIFTS)\n", "", "");

    LTOCS = Z[ITOCS] * FACFR4 - FACFR4;
    LWORK = Z[IWORK] - MWORK;
    if (CCSW) {
        LCHNDF_cc = LCHNSV;
        LTOCSV = LTOCS;
    }

    for (int NCHN = 2; NCHN <= NCHNDF; NCHN++) {
        MCHN = NCHN;
        if (CCSW) {
            MCHN = ILLOC(LCHNDF_cc);
            NSPLI = ILLOC(LCHNDF_cc + 10);
            LTOCS = LTOCSV + ILLOC(LCHNDF_cc + 12);
            LIMOST = ILLOC(LCHNDF_cc + 14);
        }

        for (KOFFS = 1; KOFFS <= NSPLI; KOFFS++) {
            JT = ILLOC(LTOCS + 4 * KOFFS);
            if (JT < 0) continue;
            JP = ILLOC(LTOCS + 4 * KOFFS - 1);
            LX = ILLOC(LTOCS + 4 * KOFFS - 2);
            LDEL = ILLOC(LTOCS + 4 * KOFFS - 3);
            int LAREA = LWORK + MWORK * KOFFS;
            double WID = ALLOC(LAREA + 1);
            if (IEXTYP == 1 && WID != 0) WID = 1.0 / WID;
            int I_v = (int)ALLOC(LAREA + 7);
            int LMAXL_v = (int)ALLOC(LAREA + 8);
            int LPEAK_v = (int)ALLOC(LAREA + 9);
            std::printf("%5d%6d/2%3d/2%3d%5d%6d%13.3G%6d"
                        "%17.5G%4d%10.3f%12.3G%6d%13.3f%10.3f%10.3f\n",
                        MCHN, JP, JT, LX, LDEL, LPEAK_v,
                        ALLOC(LAREA + 10), LIMOST,
                        ALLOC(LAREA), I_v, WID, ALLOC(LAREA + 3),
                        LMAXL_v, ALLOC(LAREA + 5), ALLOC(LAREA + 4),
                        ALLOC(LAREA + 11));
        }

        if (CCSW) {
            std::printf(" \n");
            LWORK = LWORK + NSPLI * MWORK;
            LCHNDF_cc = LCHNDF_cc + MCHNDF;
        }
    }

    Z[IWORK] = -Z[IWORK];
    Z[ILVALS] = -Z[ILVALS];
    Z[IXINTS] = -Z[IXINTS];

L3000:
    Z[ICL] = -Z[ICL];
    Z[ICINTS] = -Z[ICINTS];
    if (TRANSW) goto L3300;
    if (CCSW && NMFFAC_l == 0) goto L3300;
    Z[ICL2FF_l] = -Z[ICL2FF_l];

L3300:
    {
        double TIME = second() - TSTART;
        std::printf("0 TIME FOR L-INTERPOLATIONS : %7.3f SECONDS\n\n", TIME);
    }
    return;
}
