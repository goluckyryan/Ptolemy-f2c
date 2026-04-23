// ineld2_translated.cpp — Translated from Fortran source.f L21035-21285
// INELD2: finish INELDC — free work arrays and print timing statistics

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include "ptolemy_alloc.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_forward.h"
#include <cstdio>
#include <cmath>

extern double second();

void INELD2(float TSTART, int NUMHIN, int NUMTIM, int NUMTRM, int NUMANG,
            int NUM3, int NUM4, int NUM5, int NUM6)
{
    auto& IPRINT  = INTGER.IPRINT;
    auto& NPPHI   = INTGER.NPPHI;
    auto& NPDIF   = INTGER.NPDIF;
    auto& NPSUM   = INTGER.NPSUM;
    auto& IBSPAS  = SWITCH.IBSPAS;
    auto& IHSAVE  = INTRNL.IHSAVE;
    auto* Z       = LOCPTRS.Z;

    bool FREESW = ((IPRINT / 100) % 10) == 0;

    if (FREESW) {
        // Free areas not made for USEHS runs
        if (IHSAVE != 2) {
            Z[GRIDCM.IHS]    = -Z[GRIDCM.IHS];
            Z[GRIDCM.IHINT]  = -Z[GRIDCM.IHINT];
            Z[GRIDCM.IHABS]  = -Z[GRIDCM.IHABS];
            Z[GRIDCM.IPHI]   = -Z[GRIDCM.IPHI];
            Z[GRIDCM.IPHIP]  = -Z[GRIDCM.IPHIP];
            Z[GRIDCM.IPHIT]  = -Z[GRIDCM.IPHIT];
            Z[GRIDCM.LOCTRP] = -Z[GRIDCM.LOCTRP];
            Z[GRIDCM.IA12VL] = -Z[GRIDCM.IA12VL];
            Z[GRIDCM.IMSVAL] = -Z[GRIDCM.IMSVAL];
            Z[GRIDCM.IJA12S] = -Z[GRIDCM.IJA12S];
            Z[GRIDCM.IINTS]  = -Z[GRIDCM.IINTS];
            Z[GRIDCM.IA12TM] = -Z[GRIDCM.IA12TM];
            Z[GRIDCM.ICOSIN] = -Z[GRIDCM.ICOSIN];
            Z[GRIDCM.IXLAM]  = -Z[GRIDCM.IXLAM];
            Z[GRIDCM.INLAM]  = -Z[GRIDCM.INLAM];
            Z[GRIDCM.IRIOEX] = -Z[GRIDCM.IRIOEX];
        }
        // label 995: always free these areas
        Z[GRIDCM.IRI]    = -Z[GRIDCM.IRI];
        Z[GRIDCM.IRO]    = -Z[GRIDCM.IRO];
        Z[GRIDCM.ILOR]   = -Z[GRIDCM.ILOR];
        Z[GRIDCM.ILOI]   = -Z[GRIDCM.ILOI];
        Z[GRIDCM.ILIR]   = -Z[GRIDCM.ILIR];
        Z[GRIDCM.ILII]   = -Z[GRIDCM.ILII];
        Z[GRIDCM.IWIO]   = -Z[GRIDCM.IWIO];
        Z[GRIDCM.ISMHPT] = -Z[GRIDCM.ISMHPT];
        Z[GRIDCM.ISMIPT] = -Z[GRIDCM.ISMIPT];
        Z[GRIDCM.ISMHVL] = -Z[GRIDCM.ISMHVL];
        Z[GRIDCM.ISMIVL] = -Z[GRIDCM.ISMIVL];
        Z[GRIDCM.ISMHWK] = -Z[GRIDCM.ISMHWK];

        Z[GRIDCM.LOGIC]   = -Z[GRIDCM.LOGIC];
        Z[GRIDCM.IDW]     = -Z[GRIDCM.IDW];
        Z[GRIDCM.IABS1]   = -Z[GRIDCM.IABS1];
        Z[INELCM.ILILOR]  = -Z[INELCM.ILILOR];
        Z[INELCM.ILILOI]  = -Z[INELCM.ILILOI];
        Z[GRIDCM.IIINDX]  = -Z[GRIDCM.IIINDX];
        Z[GRIDCM.IDWI]    = -Z[GRIDCM.IDWI];
        Z[GRIDCM.IWFII]   = -Z[GRIDCM.IWFII];
        Z[GRIDCM.IWFIO]   = -Z[GRIDCM.IWFIO];

        // Keep wave-function areas for IBSPAS (multipass) calculations
        if (IBSPAS <= 0) {
            for (int I = 1; I <= 2; I++) {
                Z[WAVCOM.NFS[I]]  = -Z[WAVCOM.NFS[I]];
                Z[WAVCOM.NGS[I]]  = -Z[WAVCOM.NGS[I]];
                Z[WAVCOM.NF1S[I]] = -Z[WAVCOM.NF1S[I]];
                Z[WAVCOM.NG1S[I]] = -Z[WAVCOM.NG1S[I]];
                Z[WAVCOM.IRLVS[I]]  = -Z[WAVCOM.IRLVS[I]];
                Z[WAVCOM.IIMVS[I]]  = -Z[WAVCOM.IIMVS[I]];
                Z[WAVCOM.ICENTR[I]] = -Z[WAVCOM.ICENTR[I]];
                if (WAVCOM.V0SORS[I] != 0) Z[WAVCOM.ISORS[I]] = -Z[WAVCOM.ISORS[I]];
                if (WAVCOM.V0SOIS[I] != 0) Z[WAVCOM.ISOIS[I]] = -Z[WAVCOM.ISOIS[I]];
                for (int ITEN = 1; ITEN <= 6; ITEN++) {
                    if (WAVCOM.VTENS[I][ITEN] != 0)
                        Z[WAVCOM.ITENS[I][ITEN]] = -Z[WAVCOM.ITENS[I][ITEN]];
                }
                if (WAVCOM.TCSWS[I]) Z[WAVCOM.IWAVRS[I]] = -Z[WAVCOM.IWAVRS[I]];
                if (WAVCOM.TCSWS[I]) Z[WAVCOM.IWAVIS[I]] = -Z[WAVCOM.IWAVIS[I]];
            }
            Z[WAVCOM.IWAVR] = -Z[WAVCOM.IWAVR];
            Z[WAVCOM.IWAVI] = -Z[WAVCOM.IWAVI];
        }
    }

    // label 1500: print timing info
    // TIMES index mapping (from source.f comments):
    //   1 = WAVELJ initialization, 2 = WAVELJ loop, 3 = WAVELJ interp/norm
    //   4 = A12, 5 = phi/A12/first-U loops, 6 = H-interpolations
    //   7 = complete V-loop (5+6), 8 = complete LO-loop (1+2+3+4+7)
    float* T = FTIME.TIMES;  // 1-based: T[1..8]

    double TOTTIM = (double)second() - (double)TSTART;
    double OTHTIM = TOTTIM - (double)T[8];
    T[2] = T[2] + T[1];
    T[8] = T[8] - T[7] - T[4] - T[3] - T[2];
    T[7] = T[7] - T[6] - T[5];

    // FORMAT 2581: '0TIME IN A12...' T35 F7.3 ' SEC' (T35=col 35, so 34 chars label)
    std::printf("0%-33s%7.3f SEC\n"
                "%-34s%7.3f SEC\n"
                "%-34s%7.3f SEC\n"
                "%-34s%7.3f SEC\n"
                "%-34s%7.3f SEC\n"
                "%-34s%7.3f SEC\n"
                "%-34s%7.3f SEC\n"
                "%-34s%7.3f SEC\n"
                "%-34s%7.3f SEC\n"
                "%-34s%7d CALLS\n",
                "TIME IN A12",           (double)T[4],
                " IN PHI & A12 LOOP",   (double)T[5],
                " INTERPOLATING H'S",   (double)T[6],
                " COMPUTING SCATTERING WAVES",     (double)T[2],
                " INTERPOLATING SCATTERING WAVES", (double)T[3],
                " REST OF RI & RO LOOP",           (double)T[7],
                " REST OF LO LOOP",                (double)T[8],
                " ALL OTHER TIME",                 OTHTIM,
                " TOTAL TIME",                     TOTTIM,
                " NUMBER OF CALLS TO second",      NUMTIM);

    // Recompute TIMES as rates for FORMAT 2593/2594
    NUMTRM = NUMTRM * NPPHI;
    NUMANG = NUMANG * NPPHI;
    // T[7] = average H time (milliseconds per H computed)
    T[7] = (NUM6 > 0) ? 1.e3f * (T[4] + T[5] + T[7]) / NUM6 : 0.f;
    // T[4] = average term time (microseconds)
    T[4] = 1.e6f * T[5] / (float)(NUMTRM + 1);
    // T[5] = average cosine time (microseconds)
    T[5] = 1.e6f * T[5] / (float)(NUMANG + 1);
    // T[6] = average H-interpolation time (microseconds)
    T[6] = (NUMHIN > 0) ? 1.e6f * T[6] / (float)NUMHIN : 0.f;
    NUM3 = NPPHI * ((NUM3 + 4 * NUM4) / 5);
    NUM4 = NUM4 * NPPHI;
    double TERM_v = (double)INELCM.NUMLIS * (double)NPSUM * (double)NPDIF;
    TERM_v = (TERM_v > 0) ? 100.0 * (TERM_v - NUM5) / TERM_v : 0.0;

    // FORMAT 2587: blank line
    std::printf(" \n");

    // FORMAT 2593: T65, T63
    if (IHSAVE != 2) {
        // Col 65: I12(12) + label(45) = 57; T65 fills 7 spaces to col 65; ' AVERAGE' space at 65, A at 66
        std::printf("%12d PASSES WERE MADE THROUGH THE INNERMOST LOOP;       "
                    " AVERAGE TERM TIME =%9.3f MICROSECONDS\n",
                    NUMTRM, (double)T[4]);
        // Col 63: I12(12) + label(41) = 53; T63 fills 9 spaces to col 63; ' AVERAGE' space at 63, A at 64
        std::printf("%12d ITERATIONS OF COS(ARG-MU*PHI) WERE MADE;         "
                    " AVERAGE COSINE TIME =%9.3f MICROSECONDS\n",
                    NUMANG, (double)T[5]);
        std::printf("%12d PASSES WERE MADE THROUGH THE (M1, M2) LOOP\n", NUM3);
        std::printf("%12d PASSES WERE MADE THE PHI LOOP\n", NUM4);
    }

    // FORMAT 2594: T69, T55, then //
    // Col 69: I12(12) + label(26) = 38; T69 fills 30 spaces; "AVERAGE" (no lead space) at col 69
    std::printf("%12d H'S WERE COMPUTED OR READ                              "
                "AVERAGE H TIME =%9.3f MILLISECONDS\n",
                NUM6, (double)T[7]);
    // Col 55: I12(12) + label(32) = 44; T55 fills 10 spaces; "AVERAGE" (no lead space) at col 55
    std::printf("%12d H'S WERE FOUND BY INTERPOLATION          "
                "AVERAGE H-INTERPOLATION TIME =%9.3f MICROSECONDS\n",
                NUMHIN, (double)T[6]);
    std::printf("%12.2f%% OF THE (RI,RO) POINTS WERE SKIPPED DUE TO DWCUTOFF.\n\n\n",
                TERM_v);
}
