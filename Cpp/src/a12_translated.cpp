// a12_translated.cpp — A12
// Translated from source.f lines 1628-2162
// Computes the angular momentum transformation function A12

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include "ptolemy_alloc.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_forward.h"
#include "ptolemy_io.h"
#include <cstdio>
#include <cmath>
#include <algorithm>

void A12(int LI, int LXMIN, int LXMAX, int LMIN_a, int LMAX_a, double* XLAM_arg,
         int* NLAM_arg, double* A12VL_arg, double* DMSVAL_arg,
         int* JA12S_arg, int& JA12M, int& JA12N, int& JA12AN, int* IINDEX_flat,
         int& IHMAX, int& LOMNMN, int& LOMXMX, double* DINTS_arg, double* OUTTMP_arg,
         double* XLOTMP_arg, int* LXTMP_arg, int LHSM1, int LA12VL, int IPRINT, int INIT_arg,
         int* INDXDW_flat, int& NDW, int& NI, int* IDWFI_flat, int* IDWFO_flat)
{
    // All arrays are passed with Fortran 1-based convention.
    // Adjust base pointers so that arr[N] = Fortran arr(N).
    double* XLAM   = XLAM_arg - 1;
    int*    NLAM   = NLAM_arg - 1;
    double* A12VL  = A12VL_arg - 1;
    double* DMSVAL = DMSVAL_arg - 1;
    int*    JA12S  = JA12S_arg - 1;
    double* DINTS  = DINTS_arg - 1;
    double* OUTTMP = OUTTMP_arg - 1;
    double* XLOTMP = XLOTMP_arg - 1;
    int*    LXTMP  = LXTMP_arg - 1;

    // Access IINDEX, INDXDW, IDWFI, IDWFO as 2D arrays (column-major, 1-based)
    #define IINDEX(i,j) IINDEX_flat[((j)-1)*4 + ((i)-1)]
    #define INDXDW(i,j) INDXDW_flat[((j)-1)*4 + ((i)-1)]
    #define IDWFI(i,j)  IDWFI_flat[((j)-1)*3 + ((i)-1)]
    #define IDWFO(i,j)  IDWFO_flat[((j)-1)*4 + ((i)-1)]

    // COMMON block references
    auto& IHSAVE = INTRNL.IHSAVE;
    auto& LBP    = FORMF.LBP;
    auto& LBT    = FORMF.LBT;
    auto& JBT    = FORMF.JBT;
    auto& MSTOP  = INELCM.MSTOP;
    auto& IA12M  = INELCM.IA12M;
    auto& NMLOLX = INELCM.NMLOLX;
    auto& NPSUMI = GRIDCM.NPSUMI;
    auto& NWFI   = GRIDCM.NWFI;
    auto& NWFO   = GRIDCM.NWFO;
    auto* Z      = LOCPTRS.Z;

    // SAVE variables (static)
    static bool HALFSW;
    static int MTMIN, MPMIN, ICO, LOMOST_a;
    static int LLBT, LLBP;
    static double XN;
    static bool LBUGSW, MBUGSW, PMCHSW;

    // Local
    int LO, LX, MU, MX, MUPOS, LOSTRT, LXONE;
    int IH, IX, LL, MM, NUM, LLEVEN, MLOW;
    int IP1, IP2, IP3, IP4, ITES;
    int JA12ON, MULAST, MUONE;
    int IZERCT;
    double OUTTER, TEMP, TTT;
    double TEMPS, TEMPL, TEMPC;
    bool FRSTSW;

    // ========================================================================
    // INIT = 5: First call — compute lambda tables
    // ========================================================================
    if (INIT_arg != 5) goto L500;

    {
        int I = (IPRINT / 10000) % 10;
        LBUGSW = (I >= 4) || ((IPRINT % 10) >= 4);
        MBUGSW = (I >= 5);
        PMCHSW = (I >= 6);
    }

    // If LBT or LBP is odd, halve that interval
    HALFSW = false;
    MTMIN = -LBT;
    MPMIN = -LBP;
    if (LBT % 2 != 0) {
        // Halve MT
        HALFSW = true;
        MTMIN = 1;
        goto L50;
    }
    if (LBP % 2 != 0) {
        HALFSW = true;
        MPMIN = 1;
    }

L50:
    ICO = 1;
    LOMOST_a = LMAX_a + std::max({LXMAX, LBP, LBT});

    if (IHSAVE == 2) goto L500;

    NLAM[1] = 1;  // 1-based
    XLAM[1] = 1.0;
    {
        int M = 0;
        OUTTER = 1.0;
        if (LOMOST_a == 0) goto L200;

        if (MBUGSW) {
            std::printf("\nSTART%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d%8d%5s\n",
                        LBT, LBP, LO, LXMIN, LXMAX, LMIN_a, LMAX_a,
                        LOMOST_a, MTMIN, MPMIN, MSTOP, IA12M, NMLOLX,
                        HALFSW ? " T" : " F");
        }

        if (HALFSW) goto L160;

        // Generate LAMBDA(L, M)
        for (LL = 1; LL <= LOMOST_a; LL++) {
            M = 1 - M;
            OUTTER = OUTTER * sqrt(DINTS[LL + M - 1] / DINTS[LL + M]);
            if (M == 1) OUTTER = -OUTTER;
            ICO++;
            NLAM[LL + 1] = ICO;
            XLAM[ICO] = OUTTER;
            MU = std::min(MSTOP, LL);
            int MMTMIN = M + 2;
            if (MMTMIN > MU) continue;
            for (MM = MMTMIN; MM <= MU; MM += 2) {
                ICO++;
                XLAM[ICO] = -XLAM[ICO - 1] *
                    sqrt(DINTS[LL - MM + 2] * DINTS[LL + MM - 1] /
                        (DINTS[LL + MM] * DINTS[LL - MM + 1]));
            }
        }
        goto L200;

        // Generate LAMBDA(L, L-M) for halved interval
L160:
        for (LL = 1; LL <= LOMOST_a; LL++) {
            OUTTER = -sqrt(DINTS[2 * LL - 1] / DINTS[2 * LL]) * OUTTER;
            LLEVEN = LL - (LL % 2);
            MLOW = std::min(MSTOP, LLEVEN);
            NUM = MLOW / 2 + 1;
            ICO = ICO + NUM;
            XLAM[ICO] = OUTTER;
            NLAM[LL + 1] = ICO - LL / 2;
            if (NUM == 1) continue;
            MM = LL;
            int ICOO = ICO;
            for (int I = 2; I <= NUM; I++) {
                MM = MM - 2;
                ICOO--;
                XLAM[ICOO] = -XLAM[ICOO + 1] * sqrt(DINTS[LL - MM - 1]
                    * DINTS[LL + MM + 2] / (DINTS[LL - MM] * DINTS[LL + MM + 1]));
            }
        }
    }

L200:
    // INIT is modified in Fortran via SAVE — here we use static tracking
    // INIT = 0 means initialization done
    LLBT = 2 * LBT;
    LLBP = 2 * LBP;

    // ========================================================================
    // Main computation (INIT != 5 or after init)
    // ========================================================================
L500:
    {
        int LLI = 2 * LI;
        if (IHSAVE == 2) goto L510;

        XN = DINTS[LLI + 1] * DINTS[LLBT + 1] * DINTS[LLBP + 1];
        XN = 0.5 * sqrt(XN);

        // Set up LX table: for each LX store (IH_start, LOMIN, LOMAX)
L510:
        LOMNMN = LI + 500;
        IH = 0;
        IX = 0;
        for (LX = LXMIN; LX <= LXMAX; LX++) {
            int LOMIN = abs(LI - LX);
            int LOMAX = LI + LX;
            LOMIN = LOMIN + (LBT + LBP + LI + LOMIN) % 2;
            LOMAX = LOMAX - (LBT + LBP + LI + LOMAX) % 2;
            LXTMP[IX + 1] = IH + 1;   // 1-based
            LXTMP[IX + 2] = LOMIN;
            LXTMP[IX + 3] = LOMAX;
            IX += 3;
            LOMNMN = std::min(LOMIN, LOMNMN);
            if (LOMIN > LOMAX) continue;
            IH = IH + (LOMAX - LOMIN) / 2 + 1;
        }

        IHMAX = IH;
        LOMXMX = LX > LXMIN ? LXTMP[IX] : 0;  // LOMAX of last LX
        // Actually in Fortran, LOMXMX = LOMAX from last iteration
        // Let me use the last computed LOMAX
        {
            // Re-read from LXTMP for last LX
            int last_ix = 3 * (LXMAX - LXMIN);
            if (LXMAX >= LXMIN) LOMXMX = LXTMP[last_ix + 3];
        }

        // Build INDXDW array: pairs of (incident, outgoing) waves
        NDW = 0;
        int JPOMIN = std::max(0, 2 * LI - WAVCOM.JSPS[1] - JBT);  // JBT
        for (int KWI = 1; KWI <= NWFI; KWI++) {
            int LASI = IDWFI(1, KWI) + LI;
            if (LASI < LMIN_a || LASI > LMAX_a) continue;
            int JPI = IDWFI(2, KWI) + 2 * LI;
            if (JPI < abs(2 * LI - WAVCOM.JSPS[1])) continue;
            int JPOMAX = JPI + JBT;  // JBT
            if (WAVCOM.SOSWS[1]) JPOMIN = abs(JPI - JBT);
            for (int KWO = 1; KWO <= NWFO; KWO++) {
                LO = IDWFO(1, KWO) + LI;
                if (LO < LOMNMN) continue;
                if (!WAVCOM.SOSWS[2]) {
                    // No spin-orbit on channel 2
                    if (2 * LO + WAVCOM.JSPS[2] < JPOMIN
                        || 2 * LO - WAVCOM.JSPS[2] > JPOMAX) continue;
                } else {
                    int LASO = IDWFO(2, KWO) + LI;
                    if (LASO < 0) continue;
                    int JPO = IDWFO(3, KWO) + 2 * LI;
                    if (JPO < abs(2 * LO - WAVCOM.JSPS[2])) continue;
                    if (JPO < JPOMIN || JPO > JPOMAX) continue;
                }
                // Store pointers
                NDW++;
                INDXDW(1, NDW) = KWI;
                INDXDW(2, NDW) = KWO;
                INDXDW(3, NDW) = IDWFI(3, KWI);
                INDXDW(4, NDW) = IDWFO(4, KWO);
            }
        }

        // Build IINDEX array
        NI = 0;
        IX = 0;
        for (LX = LXMIN; LX <= LXMAX; LX++) {
            IH = LXTMP[IX + 1] - 1;
            int LOMIN = LXTMP[IX + 2];
            int LOMAX = LXTMP[IX + 3];
            IX += 3;
            for (int KDW = 1; KDW <= NDW; KDW++) {
                int KWO = INDXDW(2, KDW);
                LO = IDWFO(1, KWO) + LI;
                if (LO < LOMIN || LO > LOMAX) continue;
                NI++;
                IINDEX(1, NI) = Z[GRIDCM.ISMIVL] - 1 + (IH + (LO - LOMIN) / 2) * NPSUMI;
                IINDEX(2, NI) = Z[GRIDCM.IDW] - 1 + 2 * KDW;
                IINDEX(3, NI) = KDW;
                IINDEX(4, NI) = LX;
            }
        }

        if (IHSAVE == 2) return;

        IZERCT = 0;

        // Prepare for the triple sum on MP, MT, MU
        int MUSTRT = LI % 2;
        if (HALFSW) MUSTRT = -LI;

        JA12M = -4;
        JA12N = 0;
        JA12AN = 0;
        TEMPS = 10.0;
        TEMPL = 0;
        TEMPC = 0;

        if (LOMNMN > LOMXMX) return;

        // Main triple loop: MP, MT, MU
        for (int MP = MPMIN; MP <= LBP; MP += 2) {
            for (int MT = MTMIN; MT <= LBT; MT += 2) {
                MX = MT + MP;
                if (abs(MX) > LXMAX) continue;
                int ILAMLT = NLAM[LBT + 1] + abs(MT) / 2;
                int ILAMLP = NLAM[LBP + 1] + abs(MP) / 2;
                IP3 = 0;
                IP4 = 0;
                if (MT < 0) IP3 = LBT;
                if (MP < 0) IP4 = LBP;

                // Setup A12 factors for each LX
                FRSTSW = false;
                for (LX = LXMIN; LX <= LXMAX; LX++) {
                    OUTTER = XLAM[ILAMLT] * XLAM[ILAMLP] * XN *
                        THREEJ(LLBT, LLBP, 2 * LX, 2 * MT, 2 * MP, -2 * MX);
                    OUTTMP[LX - LXMIN + 1] = OUTTER;
                    if (OUTTER != 0) FRSTSW = true;
                }
                if (!FRSTSW) continue;  // skip this (MT, MP) pair

                FRSTSW = true;
                for (MU = MUSTRT; MU <= LI; MU += 2) {
                    MUPOS = abs(MX - MU);
                    if (MUPOS > LOMXMX) continue;

                    if (FRSTSW) {
                        // First valid MU: initialize XLOTMP
                        for (LO = LOMNMN; LO <= LOMXMX; LO++) {
                            XLOTMP[LO - LOMNMN + 1] = 0;
                        }
                    }

                    LOSTRT = std::max(LOMNMN, MUPOS);

                    // Build LO-dependent factors
                    for (LO = LOMNMN; LO <= LOMXMX; LO += 2) {
                        OUTTER = 0;
                        if (LO < LOSTRT) goto L355;
                        OUTTER = XLOTMP[LO - LOMNMN + 1];
                        if (OUTTER != 0) goto L350;

                        // First MU valid for this LO, LI
                        IP1 = 0;
                        IP2 = 0;
                        if (MU < 0) IP1 = LI;
                        if ((MX - MU) < 0) IP2 = LOMNMN;
                        ITES = IP4 + IP3 + IP2 + IP1;
                        OUTTER = XLAM[NLAM[LI + 1] + abs(MU) / 2] *
                            XLAM[NLAM[LO + 1] + MUPOS / 2] *
                            sqrt(DINTS[2 * LO + 1]);
                        if (ITES % 2 != 0) OUTTER = -OUTTER;
                        goto L355;

L350:                   // Iterate to new lambda's
                        OUTTER = OUTTER * sqrt(DINTS[LI - MU + 2]
                            * DINTS[LI + MU - 1] * DINTS[LO - MX + MU - 1]
                            * DINTS[LO + MX - MU + 2] / (DINTS[LI + MU]
                            * DINTS[LI - MU + 1] * DINTS[LO - MX + MU]
                            * DINTS[LO + MX - MU + 1]));
L355:                   XLOTMP[LO - LOMNMN + 1] = OUTTER;
                    }

                    // Compute A12 terms for all valid (LX, LO)
                    LXONE = std::max(abs(MX), LXMIN);
                    JA12ON = JA12N + 1;
                    for (LX = LXONE; LX <= LXMAX; LX++) {
                        int LOMIN = LXTMP[3 * (LX - LXMIN) + 2];
                        int LOMAX = LXTMP[3 * (LX - LXMIN) + 3];
                        if (LOMIN > LOMAX) continue;
                        for (LO = LOMIN; LO <= LOMAX; LO += 2) {
                            TTT = THREEJ(2 * LI, 2 * LO, 2 * LX,
                                         2 * MU, 2 * (MX - MU), -2 * MX);
                            TEMP = XLOTMP[LO - LOMNMN + 1] *
                                   OUTTMP[LX - LXMIN + 1] * TTT;

                            if (TEMP == 0) IZERCT++;
                            if (fabs(TEMP) > TEMPL) TEMPL = fabs(TEMP);
                            if (fabs(TEMP) < TEMPS && fabs(TEMP) != 0)
                                TEMPS = fabs(TEMP);
                            if (HALFSW || MU != 0) TEMP = TEMP + TEMP;
                            TEMPC += fabs(TEMP);
                            JA12N++;
                            A12VL[JA12N] = TEMP;
                        }
                    }
                    MULAST = MU;

                    if (!FRSTSW) continue;
                    FRSTSW = false;

                    // First MU of this (MP, MT) pair — store info
                    MUONE = MU;
                    JA12M += 5;
                    DMSVAL[JA12M]     = DINTS[MT];
                    DMSVAL[JA12M + 1] = DINTS[MP];
                    DMSVAL[JA12M + 2] = DINTS[MU] - 2;
                    JA12S[JA12M + 2]  = LXTMP[3 * (LXONE - LXMIN) + 1] + LHSM1;
                    JA12S[JA12M]      = LA12VL + JA12ON - JA12S[JA12M + 2];
                    JA12S[JA12M + 3]  = IHMAX - LXTMP[3 * (LXONE - LXMIN) + 1] + 1;
                } // end MU loop

                if (!FRSTSW) {
                    JA12S[JA12M + 1] = (MULAST - MUONE) / 2 + 1;
                    JA12AN += JA12S[JA12M + 1];
                }
            } // end MT loop
        } // end MP loop

        // Statistics
        if (JA12N > 0) TEMPC = TEMPC / JA12N;
        if (LBUGSW) {
            std::printf("\nLI, LOMNMN, LOMXMX, IHMAX, JA12M, JA12AN, JA12N:%4d%4d%4d%4d%4d%8d%8d\n"
                        " NUM OF ZEROS; SMALLEST, LARGEST AND AVERAGE A12VAL:%6d%15.5G%15.5G%15.5G\n",
                        LI, LOMNMN, LOMXMX, IHMAX, JA12M, JA12AN, JA12N,
                        IZERCT, TEMPS, TEMPL, TEMPC);
        }
    }

    #undef IINDEX
    #undef INDXDW
    #undef IDWFI
    #undef IDWFO

    return;
}
