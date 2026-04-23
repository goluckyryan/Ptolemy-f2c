// phsprt_translated.cpp — PHSPRT subroutine
// Translated from source.f lines 28242-28807
// Computes and prints elastic/reaction S-matrix magnitudes, phases, and unitarity.

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include "ptolemy_alloc.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_forward.h"
#include <cstdio>
#include <cmath>
#include <cstring>

// DMODPI: reduce X toward 0 (Fortran DINT = truncate toward zero, not floor)
static double DMODPI_f(double X) {
    const double SUBHAF = 0.4999999990;
    const double TWOPI  = 6.2831853071795864769;
    return X - TWOPI * std::trunc(X / TWOPI + (X >= 0 ? SUBHAF : -SUBHAF));
}
static double UNMOD_f(double PHASP, double PHASE) {
    return PHASP + DMODPI_f(PHASE - PHASP);
}

// Fortran G15.5: Gw.d uses F(w-4).k+4blanks when 0.1<=|val|<10^d, else Ew.d
static void print_G15_5(double val) {
    const double LIM = 1.0e5;
    double absv = std::fabs(val);
    if (absv == 0.0 || (absv >= 0.1 && absv < LIM)) {
        int e   = (absv >= 1.0) ? (int)std::floor(std::log10(absv)) + 1 : 0;
        int dec = 5 - e; if (dec < 0) dec = 0;
        std::printf("%11.*f    ", dec, val);
    } else {
        int exp_f = (int)std::floor(std::log10(absv)) + 1;
        double mant = val / std::pow(10.0, (double)exp_f);
        long long m = (long long)(std::fabs(mant) * 1e5 + 0.5);
        if (m >= 100000LL) { m /= 10; exp_f++; }
        std::printf("   %c0.%05lldE%+03d", (val < 0) ? '-' : ' ', m, exp_f);
    }
}

// Fortran F10.5 with overflow (Fortran prints "**********" when value won't fit)
static void print_F10_5(double val) {
    if (val >= 10000.0 || val < -1000.0)
        std::printf("**********");
    else
        std::printf("%10.5f", val);
}

// ============================================================================
void PHSPRT(bool PRNTSW)
{
    auto& FACFR4 = ALLOCS.FACFR4;
    auto* Z      = LOCPTRS.Z;

    int LMIN   = INTGER.LMIN;
    int LMAX   = INTGER.LMAX;
    int LSKIP  = INELCM.LSKIP;
    int NSPLI  = INELCM.NSPLI;
    int LIMOST = INTRNL.LIMOST;
    double ELAB = FLOAT_common.ELAB;

    bool CCSW = (SWITCH.PROBLM == 24);

    int IPROB = (SWITCH.PROBLM >= 23) ? 1 : 2;

    // Allocate elastic unitarity array (REAL*8, LMAX+1 elements, 0-indexed by LI)
    int IELUNI = NALLOC("ELUNITAR", LMAX + 1);
    INIT8(IELUNI, 0.0);
    int LELUNI = Z[IELUNI];  // 1-based REAL*8 index; ALLOC(LELUNI+LI) for LI=0..LMAX

    // Magnitude/phase pointer arrays (1-based: ICH=1..3)
    int ISMG[4] = {}, ISPH[4] = {};

    static const char* MAGNAM[4] = {"", "SMAG    ", "SINMAG  ", "SOUTMAG "};
    static const char* PHSNAM[4] = {"", "SPHASE  ", "SINPHASE", "SOUTPHAS"};

    // -----------------------------------------------------------------------
    // Phase 1: compute magnitudes and phases for ICH=1 (reaction), 2, 3 (elastic)
    // -----------------------------------------------------------------------
    for (int ICH = 1; ICH <= 3; ICH++) {
        bool ELSW = (ICH > 1);
        int  iprob_l = (SWITCH.PROBLM >= 23) ? 1 : (ICH - 1);

        int N = ELSW
              ? (WAVCOM.NUMJS[iprob_l] * (LMAX + 1) + FACFR4 - 1) / FACFR4
              : (INELCM.LILOSZ         + FACFR4 - 1) / FACFR4;

        ISMG[ICH] = NALLOC(MAGNAM[ICH], N);
        ISPH[ICH] = NALLOC(PHSNAM[ICH], N);
        INIT8(ISMG[ICH], 0.0);
        INIT8(ISPH[ICH], 0.0);

        int LMAG   = FACFR4 * (Z[ISMG[ICH]] - 1);
        int LPHASE = FACFR4 * (Z[ISPH[ICH]] - 1);

        // Channel parameters
        int    NPLI, LLAST, LBASE;
        double PHSMUL;
        int    LS;     // S-matrix data base (ALLOC4 for reaction, ALLOC for elastic)
        int    LTOCS_l;

        if (!ELSW) {
            NPLI   = NSPLI;   LLAST = LIMOST;  LBASE = LMIN;  PHSMUL = 1.0;
            LS     = FACFR4 * (Z[INELCM.IRDINT] - 1);
            LTOCS_l = Z[INELCM.ITOCS] * FACFR4 - FACFR4;
        } else {
            NPLI   = WAVCOM.NUMJS[iprob_l];
            LLAST  = LMAX;    LBASE = 0;       PHSMUL = 0.5;
            LS     = Z[KANDM.ISMATS[iprob_l]] - 1;
            LTOCS_l = Z[WAVCOM.ITOCE[iprob_l]] * FACFR4 - FACFR4;
        }

        for (int KOFFS = 1; KOFFS <= NPLI; KOFFS++) {
            if (ILLOC(LTOCS_l + 4*KOFFS) < 0) continue;  // skip invalid entries

            double PHASP = 0.0;
            // Backward sweep: LLI=LBASE..LLAST → LI=LLAST..LBASE
            for (int LLI = LBASE; LLI <= LLAST; LLI++) {
                int LI = LLAST + LBASE - LLI;
                int I  = (LI - LBASE) * NPLI + KOFFS;

                double SR, SI;
                if (!ELSW) {
                    SR = (double)ALLOC4(LS + 2*I - 1);
                    SI = (double)ALLOC4(LS + 2*I);
                } else {
                    SR = ALLOC(LS + 2*I - 1);
                    SI = ALLOC(LS + 2*I);
                }

                float SMAG_v = (float)std::sqrt(SR*SR + SI*SI);
                ALLOC4(LMAG + I) = SMAG_v;
                if (SMAG_v == 0.0f) continue;

                double PHASE = std::atan2(SI, SR);
                PHASE = UNMOD_f(PHASP, PHASE);
                PHASP = PHASE;
                ALLOC4(LPHASE + I) = (float)(PHSMUL * PHASE);

                if (ICH == 2)  // incoming elastic: accumulate unitarity
                    ALLOC(LELUNI + LI) += SR*SR + SI*SI;
            }
        }
    }

    if (!PRNTSW) return;

    // -----------------------------------------------------------------------
    // Phase 2: Print results
    // -----------------------------------------------------------------------
    {
        int LTOCI  = Z[WAVCOM.ITOCE[1]]     * FACFR4 - FACFR4;
        int LTOCO  = Z[WAVCOM.ITOCE[IPROB]] * FACFR4 - FACFR4;
        int LTOCS  = Z[INELCM.ITOCS]        * FACFR4 - FACFR4;

        int LMAGI  = FACFR4 * (Z[ISMG[2]] - 1);
        int LMAGO  = FACFR4 * (Z[ISMG[3]] - 1);
        int LPHASI = FACFR4 * (Z[ISPH[2]] - 1);
        int LPHASO = FACFR4 * (Z[ISPH[3]] - 1);

        int LSIGIN  = Z[KANDM.ISIGS[1]];
        int LSIGOT  = Z[KANDM.ISIGS[IPROB]];
        int LUNITR  = FACFR4 * (Z[INELCM.IUNITR] - 1) + 1 - LMIN;

        int LINES  = std::max(WAVCOM.NUMJS[1], WAVCOM.NUMJS[2]);
        int LOADD  = (!CCSW && (FORMF.LBP + FORMF.LBT) % 2 != 0) ? 1 : 0;

        // -----------------------------------------------------------
        // Elastic partial wave S-matrix table
        // -----------------------------------------------------------
        int  LINECT = 100;
        bool TOPSW  = false;

        for (int LI = LMIN; LI <= LMAX; LI += LSKIP) {

            if (LINECT + LINES > 58) {
                // FORMAT 533: elastic page header (4 output lines)
                std::printf("1%56sP T O L E M Y\n", "");
                std::printf("%38sELASTIC PARTIAL WAVE S-MATRIX ELEMENTS AND UNITARITY\n", "");
                std::printf("0%.45s%-6s%7.2f MEV     %.65s\n",
                            &HEDCOM.REACT[1], "ELAB =", ELAB, &HEDCOM.HEADER[1]);
                std::printf("\n");

                // FORMAT 543: column sub-headers (non-CC, IPROB=2)
                //   WORDS(1,2)='ELASTIC ', WORDS(2,2)='OUTGOING', WORDS(3,2)=' ELASTIC'
                std::printf("%15sINCOMING ELASTIC %25sUNITARITY%26sOUTGOING ELASTIC\n",
                            "", "", "");
                std::printf("    L   L' LX  MAGNITUDE   PHASE  COULOMB  "
                            "     ELASTIC  REACTION  RESIDUAL"
                            "        L   L' LX  MAGNITUDE   PHASE  COULOMB \n");
                std::printf("\n");

                LINECT = 9;
                TOPSW  = true;
            }

            // Unitarity values for this LI
            double UNITAR = (double)ALLOC4(LUNITR + LI);
            double ELUNIT = ALLOC(LELUNI + LI);
            double RESID  = 1.0 - ELUNIT - UNITAR;
            double CPHASI = ALLOC(LSIGIN + LI);
            int    LASO   = LI + LOADD;
            double CPHASO = ALLOC(LSIGOT + LASO);

            // Blank line before L block (except right after page header)
            if (!TOPSW) {
                if (LINES > 1 || LI % 10 == 0) {
                    std::printf(" \n");
                    LINECT++;
                }
            }

            // Inner loop: KOFFS=1..LINES, iterating incoming/outgoing together
            int KI = 0, KO = 0;
            for (int KOFFS = 1; KOFFS <= LINES; KOFFS++) {
                int JUMP = 0;
                int LXI = 0, LOI = 0;
                double DMAGI = 0.0, PHASI = 0.0;
                int LXO = 0, LOO = 0;
                double DMAGO = 0.0, PHASO = 0.0;

                // Incoming: advance KI, skip zeros
                {
                    int KI0 = KI;
                    for (;;) {
                        KI++;
                        if (KI > WAVCOM.NUMJS[1]) { KI = KI0; break; }
                        double mag = (double)ALLOC4(LMAGI + WAVCOM.NUMJS[1]*LI + KI);
                        if (mag != 0.0) {
                            DMAGI = mag;
                            PHASI = (double)ALLOC4(LPHASI + WAVCOM.NUMJS[1]*LI + KI);
                            LXI   = ILLOC(LTOCI + 4*KI - 2);
                            LOI   = LI + ILLOC(LTOCI + 4*KI - 3);
                            JUMP += 2;
                            break;
                        }
                    }
                }

                // Outgoing: advance KO, skip zeros
                {
                    int KO0 = KO;
                    for (;;) {
                        KO++;
                        if (KO > WAVCOM.NUMJS[IPROB]) { KO = KO0; break; }
                        double mag = (double)ALLOC4(LMAGO + WAVCOM.NUMJS[IPROB]*LASO + KO);
                        if (mag != 0.0) {
                            DMAGO = mag;
                            PHASO = (double)ALLOC4(LPHASO + WAVCOM.NUMJS[IPROB]*LASO + KO);
                            LXO   = ILLOC(LTOCO + 4*KO - 2);
                            LOO   = LASO + ILLOC(LTOCO + 4*KO - 3);
                            JUMP += 1;
                            break;
                        }
                    }
                }

                if (KOFFS == 1) {
                    // FORMAT 597: first row — includes unitarity columns
                    std::printf(" %4d%5d%3d%10.6f%9.3f%9.3f    ",
                                LI, LOI, LXI, DMAGI, PHASI, CPHASI);
                    print_F10_5(ELUNIT);
                    print_F10_5(UNITAR);
                    print_F10_5(RESID);
                    std::printf("     %4d%5d%3d%10.6f%9.3f%9.3f\n",
                                LASO, LOO, LXO, DMAGO, PHASO, CPHASO);
                } else {
                    if (JUMP == 0) break;  // no more data for this LI
                    if (JUMP == 1) {
                        // FORMAT 613: outgoing only (T81 = 80 spaces from col 1)
                        std::printf("%80s%4d%5d%3d%10.6f%9.3f\n",
                                    "", LASO, LOO, LXO, DMAGO, PHASO);
                    } else if (JUMP == 2) {
                        // FORMAT 623: incoming only (no outgoing data)
                        // 1X + I4+I5+I3+F10.6+F9.3, then record ends (no trailing 48X)
                        std::printf(" %4d%5d%3d%10.6f%9.3f\n",
                                    LI, LOI, LXI, DMAGI, PHASI);
                    } else {
                        // JUMP==3: both incoming and outgoing (FORMAT 623 full form)
                        // 1X + I4+I5+I3+F10.6+F9.3 + 48X + I4+I5+I3+F10.6+F9.3
                        std::printf(" %4d%5d%3d%10.6f%9.3f%48s%4d%5d%3d%10.6f%9.3f\n",
                                    LI, LOI, LXI, DMAGI, PHASI,
                                    "", LASO, LOO, LXO, DMAGO, PHASO);
                    }
                }
                LINECT++;
                TOPSW = false;
            }
        }

        // -----------------------------------------------------------
        // Reaction S-matrix table (non-CC: NCHNDF=2, NCHN=2 only)
        // -----------------------------------------------------------
        if (CCSW) return;

        {
            const int MCHN  = 2;  // channel number printed in header

            // Count valid KOFFS
            int LINES_r = 0;
            for (int K = 1; K <= NSPLI; K++)
                if (ILLOC(LTOCS + 4*K) >= 0) LINES_r++;
            int NSETS = (LINES_r + 4) / 5;

            int LMAG_r   = FACFR4 * (Z[ISMG[1]] - 1);
            int LPHASE_r = FACFR4 * (Z[ISPH[1]] - 1);

            int KF1 = 1;
            for (int NSET = 1; NSET <= NSETS; NSET++) {
                // Collect up to 5 valid partial waves
                int JTTS[6]={}, JTPS[6]={}, LXS_a[6]={}, LDELS[6]={}, KOFFSS[6]={};
                int NDO = 0;
                for (int K = KF1; K <= NSPLI && NDO < 5; K++) {
                    int JTT = ILLOC(LTOCS + 4*K);
                    if (JTT < 0) continue;
                    NDO++;
                    JTTS[NDO]   = JTT;
                    LDELS[NDO]  = ILLOC(LTOCS + 4*K - 3);
                    LXS_a[NDO]  = ILLOC(LTOCS + 4*K - 2);
                    JTPS[NDO]   = ILLOC(LTOCS + 4*K - 1);
                    KOFFSS[NDO] = K;
                    if (NDO == 5) { KF1 = K + 1; break; }
                    if (K == NSPLI) KF1 = K + 1;
                }

                // Force page break if this set won't fit on current page
                {
                    double est = 1.1 * (double)(LMAX - LMIN + LSKIP) / LSKIP;
                    if (LINECT + 8 + est > 58) LINECT = 100;
                }

                bool TOPSW_r = true;

                for (int LI = LMIN; LI <= LMAX; LI += LSKIP) {

                    bool need_hdr = (LINECT > 58) || TOPSW_r;
                    if (need_hdr) {
                        if (LINECT > 58) {
                            // FORMAT 803: new page
                            std::printf("1%47sP T O L E M Y\n", "");
                            std::printf(" %.45s%-6s%7.2f MEV     %.65s\n",
                                        &HEDCOM.REACT[1], "ELAB =", ELAB, &HEDCOM.HEADER[1]);
                        } else {
                            // FORMAT 813: blank separator between sets on same page
                            std::printf(" \n");
                        }
                        // FORMAT 804
                        std::printf("0 L IN%13sS-MATRICES FOR CHANNEL%3d"
                                    "  AND PARTIAL WAVES LABELED BY ( JP, JT, LX, L(OUT)-L(IN) )\n",
                                    "", MCHN);
                        // FORMAT 807: '0', 4X, NDO*( 3X, A2, '(', I2, '/2', I3, '/2', I3, I4, ' )' )
                        std::printf("0    ");
                        for (int I = 1; I <= NDO; I++)
                            std::printf("     (%2d/2%3d/2%3d%4d )",
                                        JTPS[I], JTTS[I], LXS_a[I], LDELS[I]);
                        std::printf("\n");
                        // FORMAT 808: 5X, NDO*( 7X, A1, '|S|', 8X, 'PHASE' )
                        std::printf("     ");
                        for (int I = 1; I <= NDO; I++)
                            std::printf("        |S|        PHASE");
                        std::printf("\n");
                        // FORMAT 813 (blank after header)
                        std::printf(" \n");
                        TOPSW_r = true;
                        LINECT += 8;
                        if (LINECT > 58) LINECT = 9;
                    }

                    // Blank line at multiples of 10 (not right after header)
                    if (LI % 10 == 0 && !TOPSW_r) {
                        std::printf(" \n");
                        LINECT++;
                    }

                    // FORMAT 863: I5, NDO*( G15.5, F9.3 )
                    std::printf("%5d", LI);
                    for (int I = 1; I <= NDO; I++) {
                        int IADD = NSPLI * (LI - LMIN) + KOFFSS[I];
                        print_G15_5((double)ALLOC4(LMAG_r   + IADD));
                        std::printf("%9.3f", (double)ALLOC4(LPHASE_r + IADD));
                    }
                    std::printf("\n");
                    TOPSW_r = false;
                    LINECT++;
                }
            }
        }
    }

}

