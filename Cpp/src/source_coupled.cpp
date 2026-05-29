// source_coupled.cpp — Coupled channels and inelastic
// Translated from source.f: A12, CCDUMP, CCMTCH, CCONV, COUPLN, DRIVE,
//   INELDC, INELD2, INGRST, INRDIN
// Phase 8 — stub bodies to be filled in from Fortran source

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include "ptolemy_alloc.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_forward.h"
#include "ptolemy_io.h"
#include <cstdio>
#include <cstring>
#include <cmath>

// A12 → a12_translated.cpp

// ============================================================================
// CCDUMP — Prints the various coupled channel arrays
// Translated from source.f lines 6424-6541
// ============================================================================
void CCDUMP()
{
    int LCHNDF, NCHNDF, MCHNDF, MCHN;
    int LBASDF, NBASDF, MBASDF, MBAS;
    int LBASCP, NBASCP, MBASCP;
    int L, NCHN, NBAS, LAST;

    if (CCBLK.ICHNDF == 0) goto L900;

    LCHNDF = ALLOCS.FACFR4 * (LOCPTRS.Z[CCBLK.ICHNDF] - 1) + 1;
    NCHNDF = ILLOC(LCHNDF);
    MCHNDF = ILLOC(LCHNDF + 1);
    LCHNDF = LCHNDF + ILLOC(LCHNDF + 2);

    std::printf("\f");
    for (int i = 1; i <= 65; i++) std::printf("%c", HEDCOM.HEADER[i]);
    std::printf("%20s", "");
    for (int i = 1; i <= 45; i++) std::printf("%c", HEDCOM.REACT[i]);
    std::printf("\n");
    std::printf("\n CHANNEL DEFINITION ARRAY:%5d CHANNELS,%5d ELEMENTS PER CHANNEL.\n",
        NCHNDF, MCHNDF);
    std::printf("\n  N  NAME   JT    KT    ET     JP    KP    EP"
        "%6sETA%8sK%7sBASIS"
        "%6sSMAT  NS%11sPOINTERS\n",
        "", "", "", "", "");
    std::printf("%68sPNT NUM%6sPNT PLI%6sINDX TOCS SMAT LIMOST\n\n",
        "", "", "");

    MCHN = MCHNDF - 1;
    for (NCHN = 1; NCHN <= NCHNDF; NCHN++) {
        L = LCHNDF + (NCHN - 1) * MCHNDF;
        std::printf("%5d%5d%4d/2%4d/2%7.3f%4d/2%4d/2%7.3f%9.2f%8.2f%6d%4d%6d%4d%5d%5d %6d%5d%5d%5d\n",
            NCHN, ILLOC(L), ILLOC(L + 1), ILLOC(L + 15),
            (double)ALLOC4(L + 2), ILLOC(L + 3), ILLOC(L + 16),
            (double)ALLOC4(L + 4),
            (double)ALLOC4(L + 5), (double)ALLOC4(L + 6),
            ILLOC(L + 7), ILLOC(L + 8),
            ILLOC(L + 9), ILLOC(L + 10),
            ILLOC(L + 11), ILLOC(L + 12),
            ILLOC(L + 13), ILLOC(L + 14),
            (MCHN >= 15 ? ILLOC(L + 15) : 0),
            (MCHN >= 16 ? ILLOC(L + 16) : 0));
    }

    if (CCBLK.IBASDF == 0) return;
    LBASDF = ALLOCS.FACFR4 * (LOCPTRS.Z[CCBLK.IBASDF] - 1) + 1;
    NBASDF = ILLOC(LBASDF);
    MBASDF = ILLOC(LBASDF + 1);
    LBASDF = LBASDF + ILLOC(LBASDF + 2);

    std::printf("\n-\n-BASIS DEFINITION ARRAY:%5d BASIS STATES,%5d ELEMENTS PER STATE.\n",
        NBASDF, MBASDF);
    std::printf("\n  N  CHAN   S    L-J%5sBASIS COUP%6sWAVE\n", "", "");
    std::printf("%24sPTR  NUM%6sSTRT\n", "", "");

    LAST = -987654;
    MBAS = MBASDF - 1;
    for (NBAS = 1; NBAS <= NBASDF; NBAS++) {
        L = LBASDF + (NBAS - 1) * MBASDF;
        if (ILLOC(L) != LAST) {
            LAST = ILLOC(L);
            std::printf(" \n");
        }
        std::printf("%5d%5d%4d/2%4d/2", NBAS, ILLOC(L), ILLOC(L + 1), ILLOC(L + 2));
        for (int i = 3; i <= MBAS; i++) std::printf("%5d", ILLOC(L + i));
        std::printf("\n");
    }

    LBASCP = ALLOCS.FACFR4 * (LOCPTRS.Z[CCBLK.IBASCP] - 1) + 1;
    NBASCP = ILLOC(LBASCP);
    MBASCP = ILLOC(LBASCP + 1);
    LBASCP = LBASCP + ILLOC(LBASCP + 2);

    std::printf("\n-\n-BASIS COUPLING ARRAY:%6d COUPLINGS,%5d ELEMENTS PER COUPLINGS.\n",
        NBASCP, MBASCP);
    std::printf("\n  N1-N2  ORD MULT%8sCOUP. CONS.%8sIFAC%5sIPT%12sPOINTERS\n",
        "", "", "", "");
    std::printf("%22sNUCLEAR%5sCOULOMB%11sBFAC%5sNUC%4sCOUL\n",
        "", "", "", "", "");

    LAST = -987654;
    MBAS = MBASCP - 1;
    for (NBAS = 1; NBAS <= NBASCP; NBAS++) {
        L = LBASCP + (NBAS - 1) * MBASCP;
        if (ILLOC(L) != LAST) {
            LAST = ILLOC(L);
            std::printf(" \n");
        }
        std::printf("%5d-%3d%4d%3d/2%13.4g%12.4g%4d%8d",
            ILLOC(L), ILLOC(L + 1), ILLOC(L + 2), ILLOC(L + 3),
            (double)ALLOC4(L + 4), (double)ALLOC4(L + 5),
            ILLOC(L + 6), ILLOC(L + 7));
        for (int i = 8; i <= MBAS; i++) std::printf("%5d", ILLOC(L + i));
        std::printf("\n");
    }

    return;

L900:
    std::printf("\n*** NO ARRAYS TO DUMP\n");
    return;
}

// ============================================================================
// CCMTCH — Extracts S-matrix element and homo C-factor from CC solution
// Translated from source.f lines 6761-7014
// ============================================================================
void CCMTCH(int ILI, int J, int IB, double* SMATS, double* CFAC_out,
            int NBASDF_arg, double* WAVR, double* WAVI, int NPTSAV,
            double* HOMO, int MBASDF_arg, int* BASDEF,
            int MBASCP_arg, int* BASCUP, float* RBASCP,
            int MCHNVL, double* CHNVAL,
            float* BFACS, int MBNDX, int* BINDEX,
            float* FG, double VSCALE)
{
    // GRIDCM CC overlay
    int* GI = reinterpret_cast<int*>(&GRIDCM.IMSVAL);
    int& NBACK_g  = GRIDCM.NOFLO;   // GI[17] = NOFLO, but CC overlay has NBACK...
    // Actually GRIDCM CC overlay has:
    //   NPTMIN=GI[16], NOFLO=GI[17], NUMPT=GI[18], NCRIT=GI[19]
    // NBACK is from INTGER common
    int& NUMPT_g  = GI[18];

    // Access BASCUP/RBASCP as 2D arrays: BASCUP(MBASCP, *), RBASCP(MBASCP, *)
    // BASCUP[col-1][row-1] = BASCUP((row), (col)) in Fortran
    // Using macros for clarity:
    #define BASCUP_F(row, col) BASCUP[(col-1)*MBASCP_arg + (row-1)]
    #define RBASCP_F(row, col) RBASCP[(col-1)*MBASCP_arg + (row-1)]
    #define BASDEF_F(row, col) BASDEF[(col-1)*MBASDF_arg + (row-1)]
    // FG is REAL*4 FG(2,2,NBASDF,*) — c4f code path
    // FG(i,j,k,l) = FG[ (l-1)*NBASDF_arg*4 + (k-1)*4 + (j-1)*2 + (i-1) ]
    #define FG_F(i, j, k, l) FG[ ((l)-1)*NBASDF_arg*4 + ((k)-1)*4 + ((j)-1)*2 + ((i)-1) ]
    // HOMO is COMPLEX*16 HOMO(2, NBASDF) — double array: HOMO(2*NBASDF*2)
    // HOMO(I, IB) = { real at [2*((IB-1)*2+(I-1))], imag at [2*((IB-1)*2+(I-1))+1] }
    #define HOMO_R(i, ib) HOMO[2*((ib-1)*2+(i-1))]
    #define HOMO_I(i, ib) HOMO[2*((ib-1)*2+(i-1))+1]
    // SMATS is COMPLEX*16 SMATS(NBASDF) — double pairs
    #define SMATS_R(ib) SMATS[2*(ib-1)]
    #define SMATS_I(ib) SMATS[2*(ib-1)+1]
    // CHNVAL is REAL*8 CHNVAL(MCHNVL, *)
    #define CHNVAL_F(row, col) CHNVAL[((col)-1)*MCHNVL + ((row)-1)]
    // BFACS is COMPLEX*8 BFACS(2,*) — float pairs: BFACS(2*2*N)
    // BFACS(II, I) = { real at [2*((I-1)*2+(II-1))], imag at [2*((I-1)*2+(II-1))+1] }
    #define BFACS_R(ii, i) ((double)BFACS[2*((i-1)*2+(ii-1))])
    #define BFACS_I(ii, i) ((double)BFACS[2*((i-1)*2+(ii-1))+1])
    // BINDEX is INTEGER BINDEX(MBNDX, *)
    #define BINDEX_F(row, col) BINDEX[((col)-1)*MBNDX + ((row)-1)]

    // ITOL array: stores -2, -1, 1, 2 to mean -I, -1, +1, +I
    // ITOL(-2..+2) with 0 -> (0,0)
    // Index: ITOL(K+3) where K is BASCUP(7,ICUP)
    double ITOL_R[5] = { 0.0, -1.0, 0.0, 1.0,  0.0};
    double ITOL_I[5] = {-1.0,  0.0, 0.0, 0.0,  1.0};

    double ROOT2 = 1.4142135623731e0;

    bool PBUGSW = (INTGER.IPRINT % 10) >= 8;

    // Local complex variables as real/imag pairs
    double S1_R[3], S1_I[3]; // 1-based
    double S2_R[3], S2_I[3];
    double S3_R[3], S3_I[3];
    double A_R[3][3], A_I[3][3]; // A(2,2) 1-based
    double D_R[3], D_I[3];
    double INCOU_R[3], INCOU_I[3];
    double OUTCOU_R[3], OUTCOU_I[3];

    for (int I = 1; I <= 2; I++) {
        S1_R[I] = 0; S1_I[I] = 0;
        S2_R[I] = 0; S2_I[I] = 0;
        S3_R[I] = 1; S3_I[I] = 0;

        // 6I*INCOMING and 6I*OUTGOING waves
        double fg1 = (double)FG_F(1, I, IB, ILI);
        double fg2 = (double)FG_F(2, I, IB, ILI);
        INCOU_R[I]  =  6.0 * fg1;
        INCOU_I[I]  =  6.0 * fg2;
        OUTCOU_R[I] = -6.0 * fg1;
        OUTCOU_I[I] =  6.0 * fg2;
    }

    double AK0 = CHNVAL_F(1, 1);
    int IC1 = BASDEF_F(1, IB);
    double AK1 = CHNVAL_F(1, IC1);
    int L1MJ = BASDEF_F(3, IB);

    if (SWITCH.NOBFAC == 1) goto L300;

    // Compute BFACTOR contributions
    {
        int ICPST  = BASDEF_F(4, IB);
        int ICPEND = ICPST + BASDEF_F(5, IB) - 1;

        for (int ICUP = ICPST; ICUP <= ICPEND; ICUP++) {
            int IBNDX = BASCUP_F(9, ICUP);
            if (IBNDX == 0) continue;
            int IBP = BASCUP_F(2, ICUP);
            int L2MJ = BASDEF_F(3, IBP);
            if (J < -L2MJ) continue;

            // Angular-momentum coupling terms and I**(L2-L1)
            int itol_idx = BASCUP_F(7, ICUP) + 3; // +3 to make 0-based from (-2..+2)
            double FACR = (double)RBASCP_F(6, ICUP) * (double)RBASCP_F(8, ICUP);
            double FACI = 0;
            // Multiply by ITOL
            double tr = FACR * ITOL_R[itol_idx-1] - FACI * ITOL_I[itol_idx-1];
            double ti = FACR * ITOL_I[itol_idx-1] + FACI * ITOL_R[itol_idx-1];
            FACR = tr; FACI = ti;
            // Multiply by VSCALE
            FACR *= VSCALE; FACI *= VSCALE;

            int origIBNDX = IBNDX;
            int LL1MJ = L1MJ;
            int LL2MJ = L2MJ;
            bool reversed = false;
            if (IBNDX < 0) {
                LL1MJ = L2MJ;
                LL2MJ = L1MJ;
                IBNDX = -IBNDX;
                reversed = true;
            }

            int LX = BINDEX_F(1, IBNDX);
            int Iii1 = ((LX+1)*(LL1MJ - BINDEX_F(6, IBNDX)) + LL2MJ - LL1MJ) / 2 + LX;
            int I_idx = ISHFT(Iii1, -1) + BINDEX_F(8, IBNDX);

            // Pick up 2 BFAC values
            double BFAC_R[3], BFAC_I[3]; // 1-based
            for (int II = 1; II <= 2; II++) {
                BFAC_R[II] = BFACS_R(II, I_idx);
                BFAC_I[II] = BFACS_I(II, I_idx);
                if (reversed) {
                    // Negate real part for reversed channels
                    BFAC_R[II] = -BFAC_R[II];
                }
                // Multiply by FAC (complex)
                double br = BFAC_R[II]*FACR - BFAC_I[II]*FACI;
                double bi = BFAC_R[II]*FACI + BFAC_I[II]*FACR;
                BFAC_R[II] = br;
                BFAC_I[II] = bi;
            }

            if (IBP == IB) goto L250;

            // Non-diagonal coupling
            {
                int ICP = BASDEF_F(1, IBP);
                double ratio = std::sqrt(AK0 / CHNVAL_F(1, ICP));
                // FAC2 = ratio * SMATS(IBP)
                double FAC2_R = ratio * SMATS_R(IBP);
                double FAC2_I = ratio * SMATS_I(IBP);
                if (BASDEF_F(7, IBP) == 2) {
                    FAC2_R /= ROOT2;
                    FAC2_I /= ROOT2;
                }
                for (int II = 1; II <= 2; II++) {
                    S2_R[II] += FAC2_R * BFAC_R[II] - FAC2_I * BFAC_I[II];
                    S2_I[II] += FAC2_R * BFAC_I[II] + FAC2_I * BFAC_R[II];
                }

                // If GS coupled to IB, need incoming term too
                if (ICP != 1) continue;
                if ((LX & 1) != 0) {
                    // Odd LX
                    S1_R[1] -= BFAC_R[1]; S1_I[1] -= BFAC_I[1];
                    S1_R[2] -= BFAC_R[2]; S1_I[2] -= BFAC_I[2];
                } else {
                    S1_R[1] += BFAC_R[1]; S1_I[1] += BFAC_I[1];
                    S1_R[2] += BFAC_R[2]; S1_I[2] += BFAC_I[2];
                }
            }
            continue;

        L250:
            // Diagonal coupling summed in S3
            S3_R[1] += BFAC_R[1]; S3_I[1] += BFAC_I[1];
            S3_R[2] += BFAC_R[2]; S3_I[2] += BFAC_I[2];

        } // end ICUP loop

        // Overall factors for B-factor stuff
        for (int I = 1; I <= 2; I++) {
            // S1(I) = -CONJG(S1(I)) * INCOU(I)
            double cr = -S1_R[I], ci = S1_I[I]; // -CONJG
            double rr = cr*INCOU_R[I] - ci*INCOU_I[I];
            double ri = cr*INCOU_I[I] + ci*INCOU_R[I];
            S1_R[I] = rr; S1_I[I] = ri;
            // S2(I) = S2(I) * OUTCOU(I)
            rr = S2_R[I]*OUTCOU_R[I] - S2_I[I]*OUTCOU_I[I];
            ri = S2_R[I]*OUTCOU_I[I] + S2_I[I]*OUTCOU_R[I];
            S2_R[I] = rr; S2_I[I] = ri;
        }

        if (PBUGSW) {
            std::printf(" B-CONTRIBS:\n");
            std::printf("%14s%15.5G%13.5G%15.5G%13.5G%15.5G%13.5G%15.5G%13.5G\n", "",
                S1_R[1], S1_I[1], S1_R[2], S1_I[2],
                S2_R[1], S2_I[1], S2_R[2], S2_I[2]);
            std::printf("%14s%15.5G%13.5G%15.5G%13.5G\n", "",
                S3_R[1], S3_I[1], S3_R[2], S3_I[2]);
        }
    }

L300:
    {
        int N = NUMPT_g - INTGER.NBACK;
        double sqrtk = std::sqrt(AK0 / AK1);

        for (int I = 1; I <= 2; I++) {
            // A(I,1) = -SQRT(AK0/AK1) * S3(I)*OUTCOU(I)
            double s3out_r = S3_R[I]*OUTCOU_R[I] - S3_I[I]*OUTCOU_I[I];
            double s3out_i = S3_R[I]*OUTCOU_I[I] + S3_I[I]*OUTCOU_R[I];
            A_R[I][1] = -sqrtk * s3out_r;
            A_I[I][1] = -sqrtk * s3out_i;
            // A(I,2) = -HOMO(I,IB)
            A_R[I][2] = -HOMO_R(I, IB);
            A_I[I][2] = -HOMO_I(I, IB);
            // D(I) = CMPLX(WAVR(N), WAVI(N)) + S1(I) + S2(I)
            D_R[I] = WAVR[N] + S1_R[I] + S2_R[I];
            D_I[I] = WAVI[N] + S1_I[I] + S2_I[I];
            if (IC1 == 1) {
                D_R[I] -= INCOU_R[I];
                D_I[I] -= INCOU_I[I];
            }
            N = NUMPT_g;
        }

        // DET = 1/(A(1,1)*A(2,2) - A(1,2)*A(2,1))
        double det_r = (A_R[1][1]*A_R[2][2] - A_I[1][1]*A_I[2][2])
                      -(A_R[1][2]*A_R[2][1] - A_I[1][2]*A_I[2][1]);
        double det_i = (A_R[1][1]*A_I[2][2] + A_I[1][1]*A_R[2][2])
                      -(A_R[1][2]*A_I[2][1] + A_I[1][2]*A_R[2][1]);
        // Invert det
        double denom = det_r*det_r + det_i*det_i;
        double DET_R = det_r / denom;
        double DET_I = -det_i / denom;

        if (PBUGSW) {
            std::printf(" MATCH:\n");
            std::printf("%14s%15.5G%13.5G%15.5G%13.5G%15.5G%13.5G%15.5G%13.5G\n", "",
                A_R[1][1], A_I[1][1], A_R[1][2], A_I[1][2],
                A_R[2][1], A_I[2][1], A_R[2][2], A_I[2][2]);
            std::printf("%14s%15.5G%13.5G%15.5G%13.5G%15.5G%13.5G\n", "",
                D_R[1], D_I[1], D_R[2], D_I[2], DET_R, DET_I);
        }

        // SMATS(IB) = (D(1)*A(2,2) - D(2)*A(1,2)) * DET
        double nr = (D_R[1]*A_R[2][2] - D_I[1]*A_I[2][2])
                   -(D_R[2]*A_R[1][2] - D_I[2]*A_I[1][2]);
        double ni = (D_R[1]*A_I[2][2] + D_I[1]*A_R[2][2])
                   -(D_R[2]*A_I[1][2] + D_I[2]*A_R[1][2]);
        SMATS_R(IB) = nr*DET_R - ni*DET_I;
        SMATS_I(IB) = nr*DET_I + ni*DET_R;

        // CFAC = (D(2)*A(1,1) - D(1)*A(2,1)) * DET
        nr = (D_R[2]*A_R[1][1] - D_I[2]*A_I[1][1])
            -(D_R[1]*A_R[2][1] - D_I[1]*A_I[2][1]);
        ni = (D_R[2]*A_I[1][1] + D_I[2]*A_R[1][1])
            -(D_R[1]*A_I[2][1] + D_I[1]*A_R[2][1]);
        CFAC_out[0] = nr*DET_R - ni*DET_I;  // real part
        CFAC_out[1] = nr*DET_I + ni*DET_R;  // imag part

        if (BASDEF_F(7, IB) == 2) {
            SMATS_R(IB) *= ROOT2;
            SMATS_I(IB) *= ROOT2;
        }
    }

    #undef BASCUP_F
    #undef RBASCP_F
    #undef BASDEF_F
    #undef FG_F
    #undef HOMO_R
    #undef HOMO_I
    #undef SMATS_R
    #undef SMATS_I
    #undef CHNVAL_F
    #undef BFACS_R
    #undef BFACS_I
    #undef BINDEX_F
    return;
}

// ============================================================================
// CCONV — Check convergence of S-matrix elements, Wynn/Pade acceleration
// Translated from source.f lines 7016-7298
// ============================================================================
void CCONV(double* CSMATS, int ITR, int& ITEST, int NBASDF_arg,
           double* PADE, int MAXIT1, double* ELAST, double* EXR,
           double* EXI, int* INDXS, int* CHNDEF, int MCHNDF_arg,
           int* BASDEF, int MBASDF_arg, int IJ, int NCHNDF_arg,
           int JT, int& NLINE, int& IHEAD, int& NFLOP)
{
    // CSMATS is COMPLEX*16 CSMATS(NBASDF, 2) — double pairs, 2 columns
    // CSMATS(I,col) real = CSMATS[2*((col-1)*NBASDF_arg + (I-1))]
    #define CSMR(i, col) CSMATS[2*((col-1)*NBASDF_arg + (i-1))]
    #define CSMI(i, col) CSMATS[2*((col-1)*NBASDF_arg + (i-1))+1]
    // PADE is COMPLEX*16 PADE(3, MAXIT1, NBASDF) — double pairs
    // PADE(k, j, i) real = PADE[2*((i-1)*MAXIT1*3 + (j-1)*3 + (k-1))]
    #define PADE_R(k, j, i) PADE[2*((i-1)*MAXIT1*3 + (j-1)*3 + (k-1))]
    #define PADE_I(k, j, i) PADE[2*((i-1)*MAXIT1*3 + (j-1)*3 + (k-1))+1]
    // CHNDEF is INTEGER CHNDEF(MCHNDF, *)
    #define CHNDEF_F(r, c) CHNDEF[((c)-1)*MCHNDF_arg + ((r)-1)]
    // BASDEF is INTEGER BASDEF(MBASDF, *)
    #define BASDEF_F2(r, c) BASDEF[((c)-1)*MBASDF_arg + ((r)-1)]

    auto& BIGEST = CNSTNT.BIGEST;
    auto& ACCURA = FLOAT_common.ACCURA;
    auto& ELAB   = FLOAT_common.ELAB;
    auto& CONVRG = SWITCH.CONVRG;
    auto& MAXITR = INTGER.MAXITR;
    auto& IPRINT = INTGER.IPRINT;

    // QWIKAB(CA) = ABS(REAL(CA)) + ABS(IMAG(CA))
    auto QWIKAB = [](double r, double i) -> double { return std::fabs(r) + std::fabs(i); };

    NFLOP = 0;
    bool PRSW = (IPRINT % 10) > 0;
    bool ITRSW = (IPRINT % 10) >= 3;
    int NADD = 1;
    if (PRSW) NADD = NBASDF_arg;

    ITEST = 0;
    int II2 = 0;
    int ITR1 = (ITR - 1) / 2;
    int ITR1P1 = ITR1 + 1;
    double DD1 = 0;
    if (ITR < 3) goto L700;

    // Use Pade table (Wynn algorithm) to accelerate convergence
    {
        for (int I = 1; I <= NBASDF_arg; I++) {
            PADE_R(3, 1, I) = CSMR(I, 1);
            PADE_I(3, 1, I) = CSMI(I, 1);
            if (QWIKAB(CSMR(I,1), CSMI(I,1)) > 1.0e-10 * BIGEST) ITEST = 2;

            for (int JJ = 1; JJ <= ITR1; JJ++) {
                double CW1_R = PADE_R(2, JJ, I), CW1_I = PADE_I(2, JJ, I);
                double CW2_R = PADE_R(1, JJ, I) - CW1_R, CW2_I = PADE_I(1, JJ, I) - CW1_I;
                double CW3_R = PADE_R(3, JJ, I), CW3_I = PADE_I(3, JJ, I);

                if (QWIKAB(CW3_R, CW3_I) == 0) {
                    PADE_R(3, JJ+1, I) = CW3_R;
                    PADE_I(3, JJ+1, I) = CW3_I;
                    continue;
                }

                CW3_R -= CW1_R; CW3_I -= CW1_I;
                double CW4_R = 0, CW4_I = 0;
                if (JJ > 1) {
                    double t_r = PADE_R(1, JJ-1, I) - CW1_R;
                    double t_i = PADE_I(1, JJ-1, I) - CW1_I;
                    double d = t_r*t_r + t_i*t_i;
                    CW4_R = t_r / d; CW4_I = -t_i / d;
                }
                // CW5 = 1/CW2 + 1/CW3 - CW4
                double d2 = CW2_R*CW2_R + CW2_I*CW2_I;
                double d3 = CW3_R*CW3_R + CW3_I*CW3_I;
                double CW5_R = CW2_R/d2 + CW3_R/d3 - CW4_R;
                double CW5_I = -CW2_I/d2 + (-CW3_I/d3) - CW4_I;
                // CW5 = 1/CW5
                double d5 = CW5_R*CW5_R + CW5_I*CW5_I;
                CW5_R = CW5_R/d5; CW5_I = -CW5_I/d5;
                // CW3 = CW5 + CW1
                PADE_R(3, JJ+1, I) = CW5_R + CW1_R;
                PADE_I(3, JJ+1, I) = CW5_I + CW1_I;
            }
        }

        NFLOP = (12 * ITR1 + 1) * NBASDF_arg;
        DD1 = 0;
        if (ITR <= 3) goto L700;

        // Check convergence of Pade S-matrix
        for (int ICHAN = 1; ICHAN <= NCHNDF_arg; ICHAN++) {
            double D1 = 0, D2 = 0;
            int I1 = CHNDEF_F(8, ICHAN);
            II2 = CHNDEF_F(9, ICHAN) + I1 - 1;

            for (int I3 = I1; I3 <= II2; I3++) {
                double CS1_R = PADE_R(3, ITR1P1, I3);
                double CS1_I = PADE_I(3, ITR1P1, I3);
                double CS2_R, CS2_I;
                if ((ITR & 1) != 0) {
                    CS2_R = PADE_R(1, ITR1, I3);
                    CS2_I = PADE_I(1, ITR1, I3);
                } else {
                    CS2_R = PADE_R(2, ITR1P1, I3);
                    CS2_I = PADE_I(2, ITR1P1, I3);
                }
                double SMOD = QWIKAB(CS1_R, CS1_I);
                if (SMOD == 0) continue;
                D2 = std::max(D2, SMOD);
                double DIFMOD = QWIKAB(CS1_R - CS2_R, CS1_I - CS2_I);
                if (CONVRG == 2) DIFMOD = DIFMOD / SMOD;
                D1 = std::max(D1, DIFMOD);
            }

            if (CONVRG == 0) D1 = D1 / D2;
            DD1 = std::max(DD1, D1);
            if (D1 <= ACCURA) continue;
            if (ITR < MAXITR && !ITRSW) goto L700;
        }
        if (DD1 < ACCURA) ITEST = 1;

    L700:
        if (ITEST == 0 && !ITRSW && ITR < MAXITR) goto L950;
        if (IHEAD == 1) goto L720;
        if (NLINE + NADD + 3 <= 59) goto L750;

        std::printf("1%59sP T O L E M Y\n"
                    " COMPUTATION OF COUPLED INELASTIC CHANNELS S-MATRICES"
                    "%41sWIN A FEW, LOOSE A FEW\n"
                    "0%.45sELAB =%7.2f MEV     %.65s\n\n",
                    "", "", &HEDCOM.REACT[1], ELAB, &HEDCOM.HEADER[1]);
        NLINE = 5;

    L720:
        std::printf("0CHANNEL  SOUT  LOUT%8sS-MATRIX(ITER)%8sS-MATRIX(PADE)"
                    "%7s|S(PADE)|%13sRELATIVE DIFFERENCES\n"
                    "%87sPADE-ITER%11sITER CONVRG.%9sPADE CONVRG.\n",
                    "", "", "", "", "", "", "");
        NLINE += 3;
        IHEAD = 2;

    L750:
        {
            const char* stop_str = (ITEST != 0) ? " STOPPING FOR OVERFLOW" : "";
            std::printf("0J =%5d,%5d ITERATIONS     CONVERGENCE =%9.2G%s\n",
                        JT, ITR, DD1, stop_str);
        }

        // Store converged S-matrix, print convergence
        for (int ICHAN = 1; ICHAN <= NCHNDF_arg; ICHAN++) {
            int NAME = CHNDEF_F(1, ICHAN);
            int I1 = CHNDEF_F(8, ICHAN);
            int I2 = CHNDEF_F(9, ICHAN) + I1 - 1;
            int IPOINT = CHNDEF_F(10, ICHAN);
            int NSPLI_v = CHNDEF_F(11, ICHAN);
            int IPINDX = CHNDEF_F(12, ICHAN);

            for (int I3 = I1; I3 <= I2; I3++) {
                double CS2_R = CSMR(I3, 1), CS2_I = CSMI(I3, 1);
                double CS1_R = CS2_R, CS1_I = CS2_I;
                if (ITR >= 3) {
                    CS1_R = PADE_R(3, ITR1P1, I3);
                    CS1_I = PADE_I(3, ITR1P1, I3);
                }
                double CS3_R = 0, CS3_I = 0;
                if (ITR >= 2 || I3 == 1) {
                    CS3_R = CSMR(I3, 2); CS3_I = CSMI(I3, 2);
                }

                int JOUT = BASDEF_F2(2, I3) / 2;
                int I4 = JOUT + 1;
                int LOUT = BASDEF_F2(3, I3) / 2;
                double SMOD = QWIKAB(CS1_R, CS1_I) + 1.0e-30;
                double D2_v = 0;
                double SITMOD = QWIKAB(CS2_R, CS2_I) + 1.0e-30;
                double D1_v = QWIKAB(CS3_R - CS2_R, CS3_I - CS2_I) / SITMOD;

                if (ITR > 3) {
                    double CS4_R, CS4_I;
                    if ((ITR & 1) != 0) {
                        CS4_R = PADE_R(1, ITR1, I3);
                        CS4_I = PADE_I(1, ITR1, I3);
                    } else {
                        CS4_R = PADE_R(2, ITR1P1, I3);
                        CS4_I = PADE_I(2, ITR1P1, I3);
                    }
                    D2_v = QWIKAB(CS4_R - CS1_R, CS4_I - CS1_I) / SMOD;
                }
                double D3_v = QWIKAB(CS1_R - CS2_R, CS1_I - CS2_I) / SMOD;
                SMOD = std::sqrt(CS1_R*CS1_R + CS1_I*CS1_I);

                if (I3 == 1) goto L850;
                {
                    // Fortran: INDXS(3*I4-1+IPINDX) is 1-based
                    int KBASE = INDXS[3*I4 - 2 + IPINDX];
                    int LDELMN = INDXS[3*I4 - 1 + IPINDX];
                    int KOFFS_v = KBASE + (LOUT - LDELMN) / 2;
                    int LABEL = (IJ - 1) * NSPLI_v + IPOINT + KOFFS_v + 1;
                    EXR[LABEL - 1] = CS1_R;
                    EXI[LABEL - 1] = CS1_I;
                }
                goto L860;
            L850:
                ELAST[0] = CS1_R;
                ELAST[1] = CS1_I;
            L860:
                LOUT = LOUT + JT;
                if (LOUT < std::abs(JT - JOUT)) continue;
                if (ICHAN != 1 && !PRSW) continue;
                NLINE++;

                if (SITMOD <= 10) {
                    std::printf(" %5d%6d%6d  %12.6f+%9.6fI%12.6f+%9.6fI%16.5G%14.3G%14.3G%14.3G\n",
                        NAME, JOUT, LOUT, CS2_R, CS2_I, CS1_R, CS1_I, SMOD, D3_v, D1_v, D2_v);
                } else {
                    std::printf(" %5d%6d%6d %12.4G%11.4GI%12.6f+%9.6fI%16.5G%14.3G%14.3G%14.3G\n",
                        NAME, JOUT, LOUT, CS2_R, CS2_I, CS1_R, CS1_I, SMOD, D3_v, D1_v, D2_v);
                }
            }
        }
        NLINE += 3;

        if (ITR >= MAXITR && ITEST != 1) {
            std::printf("0*** ERROR - FOR J =%5d, S-MATRICES NOT CONVERGED"
                        " TO SPECIFIED ACCURACY OF%12.3G\n", JT, ACCURA);
            NLINE += 2;
        }

    L950:
        // Store previous iteration
        for (int I = 1; I <= NBASDF_arg; I++) {
            CSMR(I, 2) = CSMR(I, 1);
            CSMI(I, 2) = CSMI(I, 1);
        }
        NFLOP += 7 * II2;

        if (ITR >= 3) {
            // Store Pade approximants
            for (int I = 1; I <= NBASDF_arg; I++) {
                for (int JJ = 1; JJ <= ITR1P1; JJ++) {
                    PADE_R(1, JJ, I) = PADE_R(2, JJ, I);
                    PADE_I(1, JJ, I) = PADE_I(2, JJ, I);
                    PADE_R(2, JJ, I) = PADE_R(3, JJ, I);
                    PADE_I(2, JJ, I) = PADE_I(3, JJ, I);
                }
            }
        } else {
            for (int I = 1; I <= NBASDF_arg; I++) {
                PADE_R(ITR, 1, I) = CSMR(I, 1);
                PADE_I(ITR, 1, I) = CSMI(I, 1);
            }
        }
    }

    #undef CSMR
    #undef CSMI
    #undef PADE_R
    #undef PADE_I
    #undef CHNDEF_F
    #undef BASDEF_F2
    return;
}

// ============================================================================
// COUPLN — Iterative solution of multi-channel coupled equations
// Translated from source.f lines 11875-12803
// ============================================================================
void COUPLN(int& IRET)
{
    extern double second();

    auto* Z      = LOCPTRS.Z;
    auto& FACFR4 = ALLOCS.FACFR4;

    // GRIDCM CC overlay
    int* GI = reinterpret_cast<int*>(&GRIDCM.IMSVAL);
    double& RCUEFF_g = GRIDCM.JACOB;
    int& IBINDX_g = GI[0];
    int& MBINDX_g = GI[1];
    int& NBINDX_g = GI[2];
    int& NMFFAC_g = GI[3];
    int& NMBFAC_g = GI[5];
    int& IINHR_g  = GI[6];
    int& IINHI_g  = GI[7];
    int& IRHSR_g  = GI[8];
    int& IRHSI_g  = GI[9];
    int& INUCH_g  = GI[10];
    int& ICOULH_g = GI[11];
    int& NPTMIN_g = GI[16];
    int& NOFLO_g  = GI[17];
    int& NUMPT_g  = GI[18];
    int& NCRIT_g  = GI[19];
    int& ISMX_g   = GI[20];
    int& ISMSQ_g  = GI[21];
    int& IPADE_g  = GI[22];
    int& IHOMO_g  = GI[23];
    int& ISMTSV_g = GI[24];
    int& ICL1FF_g = GI[45];
    int& ICL1FG_g = GI[46];

    auto& LXMAX   = INELCM.LXMAX;
    auto& NUMLIS  = INELCM.NUMLIS;
    auto& ILIS    = INELCM.ILIS;
    auto& NSPLI   = INELCM.NSPLI;
    auto& ISMATR  = INELCM.ISMATR;
    auto& ISMATI  = INELCM.ISMATI;
    auto& IINDXS  = INELCM.IINDXS;
    auto& IELCUP  = INELCM.IELCUP;
    auto& NFIROF  = INTGER.NFIROF;
    auto& MAXITR  = INTGER.MAXITR;
    auto& IPRINT  = INTGER.IPRINT;
    auto& NBACK   = INTGER.NBACK;
    auto& LMIN    = INTGER.LMIN;
    auto& ELAB    = FLOAT_common.ELAB;
    auto& ACCURA  = FLOAT_common.ACCURA;
    auto& STEPSZ  = FLOAT_common.STEPSZ;
    auto& UNDEF   = INTRNL.UNDEF;
    auto& AKI     = KANDM.AKI;
    auto& LOMOST  = KANDM.LOMOST;
    auto& NOBFAC  = SWITCH.NOBFAC;
    auto& CONVRG  = SWITCH.CONVRG;
    auto& IWRTWV  = SWITCH.IWRTWV;
    auto& ISAVSM  = SWITCH.ISAVSM;
    auto& ISW3    = SWITCH.ISW3;
    auto& JP      = JBLOCK.JP;
    auto& ICHNDF  = CCBLK.ICHNDF;
    auto& IBASDF  = CCBLK.IBASDF;
    auto& IBASCP  = CCBLK.IBASCP;
    auto& MCHNVL  = CCBLK.MCHNVL;
    auto& ICHNVL  = CCBLK.ICHNVL;
    auto& IHOMOA  = CCBLK.IHOMOA;
    auto& IHOMOB  = CCBLK.IHOMOB;
    auto& IINHR8  = CCBLK.IINHR8;
    auto& IINHI8  = CCBLK.IINHI8;

    // COMPLEX*16 CFAC
    double CFAC[2]; // {real, imag}
    double CFACR, CFACI;
    double CFACR8, CFACI8;
    bool PBUGSW, CONVSW, BLOKSW;
    float DUMMY4[1] = {0};
    char ACCWRD[3][9] = {"BLOCKREL", "ABSOLUTE", "RELATIVE"};

    double TSTART = second();
    double TRIGHT = 0, TUOPT = 0, TINHO = 0, TWAVE = 0;
    double TNORM = 0, TITER = 0, TPADE = 0;
    int NMINHO = 0, NMINPT = 0;
    double FLPINH = 0, FLPRHS = 0, FLPUOP = 0, FLPADE = 0, FLPNRM = 0;
    int IPRNT = IPRINT % 10;

    std::printf("1%59sP T O L E M Y\n"
                " COMPUTATION OF COUPLED INELASTIC CHANNELS S-MATRICES"
                "%41sWIN A FEW, LOOSE A FEW\n"
                "0%.45sELAB =%7.2f MEV     %.65s\n\n",
                "", "", &HEDCOM.REACT[1], ELAB, &HEDCOM.HEADER[1]);

    // Temporary coupling parameters
    double BLEND = FLOAT_arr(131);
    if (BLEND == UNDEF) BLEND = 1;
    double VSCALE = FLOAT_arr(130);
    if (VSCALE == UNDEF) VSCALE = 0;
    BLOKSW = (ISW3 == 1);

    std::printf(" %.8s ACCURACY REQUIREMENT IS%12.3G\n"
                " BLENDING FACTOR IS%15.5G\n"
                " COUPLING STRENGTH SCALED BY 1 - (%14.5G)**ITERATION.\n",
                ACCWRD[CONVRG], ACCURA, BLEND, VSCALE);
    if (!BLOKSW) std::printf(" USING SEQUENTIAL ITERATION\n\n");
    if (BLOKSW)  std::printf(" USING BLOCK (BORN) ITERATION\n\n");

    int NLINE = 10;
    int IHEAD = 1;
    PBUGSW = (IPRNT >= 4);

    int ISPUN = 1; // no spin-orbit

    // Setup variables
    int NSTEP = NUMPT_g - 1;
    int NPTSAV = NUMPT_g - NPTMIN_g;
    int NPSTRT = NPTMIN_g + 1;

    // Get LBASDF to read NBASDF
    int LBASDF_t = FACFR4 * Z[IBASDF] - FACFR4 + 1;
    int NBASDF_v = ILLOC(LBASDF_t);
    int NUMINH = NPTSAV * NBASDF_v;

    int IBLOCK = 0, ISBLOK = 0;
    if (BLOKSW) IBLOCK = NALLOC("BLOCKSAV", (2*NUMINH+1)/FACFR4);
    if (BLOKSW) ISBLOK = NALLOC("SBLOCKSA", 2*NBASDF_v);

    SETLOG(2 * (LOMOST + LXMAX));

    // Get addresses (no more ALLOC calls after this)
    // c8h code path: LHOMO = Z(IHOMO) - 2*NPTMIN
    int LHOMO = Z[IHOMO_g] - 2*NPTMIN_g;
    // c8h: LINHR = Z(IINHR) - NPTMIN
    int LINHR = Z[IINHR_g] - NPTMIN_g;
    int LINHI = Z[IINHI_g] - NPTMIN_g;
    int LHOMOA = Z[IHOMOA];
    int LHOMOB = Z[IHOMOB];
    int LINHR8 = Z[IINHR8] - 1;
    int LINHI8 = Z[IINHI8] - 1;
    int LBLOCK = 0, LSBLOK = 0;
    if (BLOKSW) LBLOCK = FACFR4 * Z[IBLOCK] - FACFR4;
    if (BLOKSW) LSBLOK = Z[ISBLOK];
    int LRHSR = Z[IRHSR_g];
    int LRHSI = Z[IRHSI_g];
    int LSMSQ = Z[ISMSQ_g];
    int LPADE = Z[IPADE_g];
    int LLIS = FACFR4 * Z[ILIS] - FACFR4;
    int LFG = FACFR4 * Z[WAVCOM.NFS[1]] - FACFR4 + 1;
    int LWAVR = Z[WAVCOM.IWAVR];
    int LWAVI = Z[WAVCOM.IWAVI];
    int LWRK1 = LWAVR - 1;
    int LWRK2 = LWAVI - 1;
    int LVRL1 = Z[WAVCOM.IRLVS[1]];
    int LVRL2 = Z[WAVCOM.IRLVS[2]];
    int LVIM1 = Z[WAVCOM.IIMVS[1]];
    int LVIM2 = Z[WAVCOM.IIMVS[2]];
    int LNUCH = FACFR4 * Z[INUCH_g] - FACFR4 + 1;
    int LCOULH = FACFR4 * Z[ICOULH_g] - FACFR4 + 1;
    int LCENT1 = Z[WAVCOM.ICENTR[1]];
    int LCENT2 = Z[WAVCOM.ICENTR[2]];
    int LCL1FF = FACFR4 * Z[ICL1FF_g] - FACFR4 + 1;
    int LCL1FG = FACFR4 * Z[ICL1FG_g] - FACFR4 + 1;
    int LBINDX = FACFR4 * Z[IBINDX_g] - FACFR4 + 1;
    int LSMATS = Z[KANDM.ISMATS[1]];
    int LCHNDF = FACFR4 * Z[ICHNDF] - FACFR4 + 1;
    int LBASDF = FACFR4 * Z[IBASDF] - FACFR4 + 1;
    int LBASCP = FACFR4 * Z[IBASCP] - FACFR4 + 1;
    int NCHNDF_v = ILLOC(LCHNDF);
    NBASDF_v = ILLOC(LBASDF);
    int NBASCP_v = ILLOC(LBASCP);
    int MCHNDF_v = ILLOC(LCHNDF + 1);
    int MBASDF_v = ILLOC(LBASDF + 1);
    int MBASCP_v = ILLOC(LBASCP + 1);
    LCHNDF += ILLOC(LCHNDF + 2);
    LBASDF += ILLOC(LBASDF + 2);
    LBASCP += ILLOC(LBASCP + 2);
    int LCHNVL = Z[ICHNVL];
    int LSMATR = Z[ISMATR] - 1;
    int LSMATI = Z[ISMATI] - 1;
    int LINDXS = FACFR4 * (Z[IINDXS] - 1);
    int LELCUP = Z[IELCUP] - 2;

    if (IWRTWV != 0) {
        char thedat[10] = {};
        get_date(thedat);
        // Write to unit 1 (file output — not implemented in DWBA stub)
    }

    // ========================================================================
    // MAJOR LOOP OVER TOTAL ANGULAR MOMENTUM J
    // ========================================================================
    for (int IJ = 1; IJ <= NUMLIS; IJ++) {
        CONVSW = (IWRTWV % 10) == 2;

        // Initialize arrays
        for (int I = 1; I <= NUMINH; I++) {
            ALLOC(LINHR + NPTMIN_g + I - 1) = 0.0;
            ALLOC(LINHI + NPTMIN_g + I - 1) = 0.0;
        }
        int II = 2 * NBASDF_v;
        for (int I = 1; I <= II; I++) {
            ALLOC(LINHR8 + I) = 0;
            ALLOC(LINHI8 + I) = 0;
        }
        int NWRDS = 4 * NBASDF_v;
        for (int I = 1; I <= NWRDS; I++) {
            ALLOC(LSMSQ + I - 1) = 0.0;
        }

        int JT = ILLOC(LLIS + IJ);
        int J = 2 * JT;
        int NELCUP = LELCUP + 2 * IJ;
        int NADD = NBASDF_v;
        if (IPRNT == 0) NADD = 1;
        if (NLINE + NADD + 3 > 59) {
            std::printf("1%59sP T O L E M Y\n"
                        " COMPUTATION OF COUPLED INELASTIC CHANNELS S-MATRICES"
                        "%41sWIN A FEW, LOOSE A FEW\n"
                        "0%.45sELAB =%7.2f MEV     %.65s\n\n",
                        "", "", &HEDCOM.REACT[1], ELAB, &HEDCOM.HEADER[1]);
            NLINE = 5;
        }

        // Calculate B factors
        if (NOBFAC == 0) {
            GETBFC(IJ, J, NBACK * STEPSZ, FLOAT_common.SUMMAX,
                   NMBFAC_g, &ALLOC4(LCL1FF), &ALLOC4(LCL1FG),
                   MBINDX_g, &ILLOC(LBINDX),
                   MCHNDF_v, &ILLOC(LCHNDF), MCHNVL, &ALLOC(LCHNVL),
                   MBASDF_v, &ILLOC(LBASDF), NBASDF_v, &ALLOC4(LFG));
        }

        // Loop over channel array to setup
        for (int I1 = 1; I1 <= NCHNDF_v; I1++) {
            int LABEL = (I1 - 1) * MCHNDF_v + LCHNDF;
            double CKC = ALLOC(LCHNVL + MCHNVL * (I1 - 1));

            int I2 = ILLOC(LABEL + 7);
            int I3 = ILLOC(LABEL + 8);
            I3 = I3 + I2 - 1;

            // Loop over basis states
            for (int I4 = I2; I4 <= I3; I4++) {
                int LABEL1 = (I4 - 1) * MBASDF_v + LBASDF;
                int J1 = ILLOC(LABEL1 + 1);
                int L1X2 = ILLOC(LABEL1 + 2);
                L1X2 = L1X2 + J;
                if (L1X2 < 0) continue;
                int L1 = L1X2 / 2;

                int I5 = ILLOC(LABEL1 + 3);
                int I6 = ILLOC(LABEL1 + 4);
                I6 = I6 + I5 - 1;

                // Calculate coupling coefficients
                for (int I7 = I5; I7 <= I6; I7++) {
                    int LABEL2 = (I7 - 1) * MBASCP_v + LBASCP;
                    int IMULT = ILLOC(LABEL2 + 3);
                    int IPT = ILLOC(LABEL2 + 1);
                    int J2 = ILLOC((IPT - 1) * MBASDF_v + LBASDF + 1);
                    int L2X2 = ILLOC((IPT - 1) * MBASDF_v + LBASDF + 2);
                    L2X2 = L2X2 + J;
                    {
                        double akrest_tmp;
                        AKFCR(L1X2, J1, L2X2, J2, IMULT, J, akrest_tmp);
                        ALLOC4(LABEL2 + 7) = (float)akrest_tmp;
                    }
                }

                // Generate homogeneous wavefunction
                int NWP = 2;
                if (I2 == 1) NWP = 1;

                AKS()[NWP] = CKC;
                WAVCOM.HS[NWP] = CKC * WAVCOM.RSTEPS[NWP];

                // Label for Coulomb functions
                int LABEL2 = 4 * (I4 - 1) + 4 * NBASDF_v * (IJ - 1);
                LABEL2 = LFG + LABEL2;
                double TT = second();

                if (NWP == 1)
                    WAVELJ(L1, ISPUN, 1, 0, DUMMY4, ALLOC4_base(LABEL2), DUMMY4,
                           ALLOC_base(LWAVR), ALLOC_base(LWAVI),
                           ALLOC_base(LVRL1), ALLOC_base(LVIM1), ALLOC_base(LCENT1));
                if (NWP == 2)
                    WAVELJ(L1, ISPUN, 2, 0, DUMMY4, ALLOC4_base(LABEL2), DUMMY4,
                           ALLOC_base(LWAVR), ALLOC_base(LWAVI),
                           ALLOC_base(LVRL2), ALLOC_base(LVIM2), ALLOC_base(LCENT2));

                double TT1 = second();
                TWAVE += TT1 - TT;

                // Starting value for this basis state
                int NFIRST = MAX0(WAVCOM.LASTNF[NWP] + NFIROF, NPTMIN_g + 2);
                ILLOC(LABEL1 + 5) = NFIRST;
                int ISTRT = NFIRST - 2;

                // Store homogeneous wavefunction (c8h path)
                II = 2 * (I4 - 1) * NPTSAV;
                if (PBUGSW) {
                    std::printf(" J, BASIS, NFIRST =%5d/2%9d%9d%9d%9d%9d%9d%9d%9d%9d\n",
                                J, I4, NFIRST, ISTRT, NSTEP, II,
                                NPTSAV, NPTMIN_g, LHOMO, LWAVR);
                }
                for (int I7 = ISTRT; I7 <= NSTEP; I7++) {
                    ALLOC(LHOMO + II + 2*I7)     = ALLOC(LWAVR + I7);
                    ALLOC(LHOMO + II + 2*I7 + 1) = ALLOC(LWAVI + I7);
                }

                // Full precision for first and last values
                int ITEMP = NSTEP - NBACK;
                for (int I = 1; I <= 2; I++) {
                    ALLOC(LHOMOA + 4*(I4-1) + 2*I - 2) = ALLOC(LWAVR + ISTRT + I - 1);
                    ALLOC(LHOMOA + 4*(I4-1) + 2*I - 1) = ALLOC(LWAVI + ISTRT + I - 1);
                    ALLOC(LHOMOB + 4*(I4-1) + 2*I - 2) = ALLOC(LWAVR + ITEMP);
                    ALLOC(LHOMOB + 4*(I4-1) + 2*I - 1) = ALLOC(LWAVI + ITEMP);
                    ITEMP = NSTEP;
                }
                if (I2 != 1) goto L180;

                // Store elastic wavefunction for 0th guess
                ALLOC(LSMSQ)     = ALLOC(LSMATS + 2*L1);
                ALLOC(LSMSQ + 2*NBASDF_v) = ALLOC(LSMATS + 2*L1);
                ALLOC(LSMSQ + 1) = ALLOC(LSMATS + 2*L1 + 1);
                ALLOC(LSMSQ + 1 + 2*NBASDF_v) = ALLOC(LSMATS + 2*L1 + 1);

                for (int I = 1; I <= 2; I++) {
                    ALLOC(LINHR8 + I) = ALLOC(LWAVR + ISTRT + I - 1);
                    ALLOC(LINHI8 + I) = ALLOC(LWAVI + ISTRT + I - 1);
                }
                II = II / 2;
                for (int I7 = ISTRT; I7 <= NSTEP; I7++) {
                    ALLOC(LINHR + II + I7) = ALLOC(LWAVR + I7);
                    ALLOC(LINHI + II + I7) = ALLOC(LWAVI + I7);
                }
            L180:
            } // end I4 loop (194)
        } // end I1 loop (199)

        // Block iteration: save initial solution
        if (BLOKSW) {
            for (int I = 1; I <= NUMINH; I++) {
                ALLOC4(LBLOCK + I) = (float)ALLOC(LINHR + NPTMIN_g + I - 1);
                ALLOC4(LBLOCK + NUMINH + I) = (float)ALLOC(LINHI + NPTMIN_g + I - 1);
            }
            for (int I = 1; I <= NBASDF_v; I++) {
                ALLOC(LSBLOK + 2*I - 2) = ALLOC(LSMSQ + 2*I - 2);
                ALLOC(LSBLOK + 2*I - 1) = ALLOC(LSMSQ + 2*I - 1);
            }
        }

        // ====================================================================
        // START ITERATIVE SOLUTION OF INHOMOGENEOUS EQUATIONS
        // ====================================================================
        double TT = second();
        double VFACTR = 1;
        int ITR;
        int ITEST = 0;

        for (ITR = 1; ITR <= MAXITR; ITR++) {
            VFACTR *= VSCALE;
            double VMULT = 1 - VFACTR;

            // Loop over channels (elastic last)
            for (int I1 = 1; I1 <= NCHNDF_v; I1++) {
                int ICHAN = I1 + 1;
                if (I1 == NCHNDF_v) ICHAN = 1;
                int LABEL = (ICHAN - 1) * MCHNDF_v + LCHNDF;

                int NWP = 2;
                if (ICHAN == 1) NWP = 1;

                WAVCOM.HS[NWP] = ALLOC(LCHNVL + MCHNVL * (ICHAN - 1)) * WAVCOM.RSTEPS[NWP];

                int I2 = ILLOC(LABEL + 7);
                int I3 = ILLOC(LABEL + 8);
                I3 = I3 + I2 - 1;

                // Loop over basis states
                for (int I4 = I2; I4 <= I3; I4++) {
                    int LABEL1 = (I4 - 1) * MBASDF_v + LBASDF;
                    int L1X2 = ILLOC(LABEL1 + 2);
                    L1X2 = L1X2 + J;
                    if (L1X2 < 0) continue;
                    int L1 = L1X2 / 2;
                    int NFIRST = ILLOC(LABEL1 + 5);

                    int I5 = ILLOC(LABEL1 + 3);
                    int I6 = ILLOC(LABEL1 + 4);
                    I6 = I6 + I5 - 1;

                    // Construct driving term (c8h path)
                    double TT2 = second();
                    double FLP;
                    DRIVE(ALLOC_base(LRHSR), ALLOC_base(LRHSI),
                          I5, I6, NFIRST, ITR, MBASCP_v,
                          &ILLOC(LBASCP), &ALLOC4(LBASCP), &ALLOC4(LNUCH),
                          &ALLOC4(LCOULH), NUMPT_g, NPTSAV,
                          &ALLOC(LINHR), &ALLOC(LINHI),
                          FLP, VMULT);

                    if (FLP == 0) continue;

                    double TT3 = second();
                    TRIGHT += TT3 - TT2;
                    FLPRHS += FLP;

                    // Setup pointers for this basis state
                    int ISTRT_v = (I4 - 1) * NPTSAV;
                    int NHOMO = LHOMO + 2 * ISTRT_v;
                    int NINHR = LINHR + ISTRT_v;
                    int NINHI = LINHI + ISTRT_v;
                    int NINHR8 = LINHR8 + 2*(I4 - 1);
                    int NINHI8 = LINHI8 + 2*(I4 - 1);
                    int NHOMOA = LHOMOA + 4*(I4 - 1) - 2;
                    if (BLOKSW) NINHR = LBLOCK - NPTMIN_g + ISTRT_v + 1;
                    if (BLOKSW) NINHI = NINHR + NUMINH;
                    ISTRT_v = NFIRST - 2;

                    // Get first two 12*U's (c8h path)
                    II = ISTRT_v + 1;
                    for (int I = 1; I <= 2; I++) {
                        ALLOC(LWRK1 + ISTRT_v + I) = ALLOC(NINHR8 + I);
                        ALLOC(LWRK2 + ISTRT_v + I) = ALLOC(NINHI8 + I);
                    }

                    // Construct diagonal potential and integrate
                    WAVELJ(L1, (int)JP, NWP, -NFIRST, DUMMY4, DUMMY4, DUMMY4,
                           ALLOC_base(LWRK1 + 1), ALLOC_base(LWRK2 + 1),
                           ALLOC_base(Z[WAVCOM.IRLVS[NWP]]),
                           ALLOC_base(Z[WAVCOM.IIMVS[NWP]]),
                           ALLOC_base(Z[WAVCOM.ICENTR[NWP]]));

                    if (IPRNT >= 5) {
                        std::printf(" ++%3d%5d", I4, ISTRT_v);
                        std::printf("%15.6G%15.6G%15.6G%15.6G%15.6G%15.6G%15.6G%15.6G\n",
                            ALLOC(LWRK1+II+1), ALLOC(LWRK2+II+1),
                            ALLOC(LWRK1+II+4), ALLOC(LWRK2+II+4),
                            ALLOC(LRHSR+NFIRST-1), ALLOC(LRHSI+NFIRST-1),
                            0.0, 0.0);
                    }

                    // Integrate inhomogeneous equation
                    TT2 = second();
                    TUOPT += TT2 - TT3;
                    FLPUOP += 3.0 * (NSTEP - NFIRST);

                    WAVINH(NFIRST, NSTEP, ALLOC_base(LWRK1 + 1), ALLOC_base(LWRK2 + 1),
                           ALLOC_base(LRHSR), ALLOC_base(LRHSI), FLP, NSTEP + 5);

                    TINHO += second() - TT2;
                    FLPINH += FLP;
                    NMINHO++;
                    NMINPT += NSTEP - NFIRST + 1;

                    // Calculate matching condition
                    CCMTCH(IJ, J, I4, &ALLOC(LSMSQ), CFAC, NBASDF_v,
                           ALLOC_base(LWRK1 + 1), ALLOC_base(LWRK2 + 1), NPTSAV,
                           &ALLOC(LHOMOB), MBASDF_v, &ILLOC(LBASDF),
                           MBASCP_v, &ILLOC(LBASCP), &ALLOC4(LBASCP),
                           MCHNVL, &ALLOC(LCHNVL), &ALLOC4(LCL1FG),
                           MBINDX_g, &ILLOC(LBINDX), &ALLOC4(LFG), VMULT);

                    if (PBUGSW) {
                        std::printf(" CCMTCH,ITER,J,I4=%5d%5d/2%5d   S=%15.5G%15.5G   C=%15.5G%15.5G\n",
                                    ITR, J, I4, ALLOC(LSMSQ + 2*I4 - 2),
                                    ALLOC(LSMSQ + 2*I4 - 1), CFAC[0], CFAC[1]);
                    }

                    // Block iteration S-matrix patch
                    if (BLOKSW) {
                        for (int I = 1; I <= 2; I++) {
                            double TEMP = ALLOC(LSMSQ + 2*I4 + I - 3);
                            ALLOC(LSMSQ + 2*I4 + I - 3) = ALLOC(LSBLOK + 2*I4 + I - 3);
                            ALLOC(LSBLOK + 2*I4 + I - 3) = TEMP;
                        }
                    }

                    // Calculate normalised wavefunction and blend
                    CFACR8 = CFAC[0];
                    CFACR = CFACR8;
                    CFACI8 = CFAC[1];
                    CFACI = CFACI8;
                    TT2 = second();
                    double TEMP;

                    if (BLEND != 1) goto L640;
                    // Direct replacement (BLEND==1)
                    for (int I = 1; I <= 2; I++) {
                        ALLOC(NINHR8 + I) = ALLOC(LWRK1 + ISTRT_v + I)
                            + (CFACR8 * ALLOC(NHOMOA + 2*I) - CFACI8 * ALLOC(NHOMOA + 2*I + 1));
                        ALLOC(NINHI8 + I) = ALLOC(LWRK2 + ISTRT_v + I)
                            + (CFACR8 * ALLOC(NHOMOA + 2*I + 1) + CFACI8 * ALLOC(NHOMOA + 2*I));
                    }
                    for (int I = ISTRT_v; I <= NSTEP; I++) {
                        ALLOC(NINHR + I) = ALLOC(LWRK1 + I + 1)
                            + (CFACR * ALLOC(NHOMO + 2*I) - CFACI * ALLOC(NHOMO + 2*I + 1));
                        ALLOC(NINHI + I) = ALLOC(LWRK2 + I + 1)
                            + (CFACR * ALLOC(NHOMO + 2*I + 1) + CFACI * ALLOC(NHOMO + 2*I));
                    }
                    TEMP = 8;
                    goto L650;

                L640:
                    // Blending
                    TEMP = 1 - BLEND;
                    for (int I = 1; I <= 2; I++) {
                        ALLOC(NINHR8 + I) = TEMP * ALLOC(NINHR8 + I)
                            + BLEND * (ALLOC(LWRK1 + ISTRT_v + I)
                            + (CFACR8 * ALLOC(NHOMOA + 2*I) - CFACI8 * ALLOC(NHOMOA + 2*I + 1)));
                        ALLOC(NINHI8 + I) = TEMP * ALLOC(NINHI8 + I)
                            + BLEND * (ALLOC(LWRK2 + ISTRT_v + I)
                            + (CFACR8 * ALLOC(NHOMOA + 2*I + 1) + CFACI8 * ALLOC(NHOMOA + 2*I)));
                    }
                    for (int I = ISTRT_v; I <= NSTEP; I++) {
                        ALLOC(NINHR + I) = TEMP * ALLOC(NINHR + I)
                            + BLEND * (ALLOC(LWRK1 + I + 1)
                            + (CFACR * ALLOC(NHOMO + 2*I) - CFACI * ALLOC(NHOMO + 2*I + 1)));
                        ALLOC(NINHI + I) = TEMP * ALLOC(NINHI + I)
                            + BLEND * (ALLOC(LWRK2 + I + 1)
                            + (CFACR * ALLOC(NHOMO + 2*I + 1) + CFACI * ALLOC(NHOMO + 2*I)));
                    }
                    TEMP = 14;
                L650:
                    {
                        double TT4 = second();
                        TNORM += TT4 - TT2;
                        FLPNRM += TEMP * (NSTEP - ISTRT_v + 1);
                    }

                    if (!CONVSW) continue;
                    // Write wavefunctions (stub)
                    if (IWRTWV >= 10 && ICHAN != 1) continue;
                } // end I4 loop (689)
            } // end I1 loop (699)

            // Move S-matrices for block iteration
            NWRDS = 2 * NBASDF_v;
            if (BLOKSW) {
                for (int I = 1; I <= NWRDS; I++) {
                    ALLOC(LSMSQ + I - 1) = ALLOC(LSBLOK + I - 1);
                }
            }

            // Test convergence
            {
                double TT4 = second();
                int NFLOP_v;
                CCONV(&ALLOC(LSMSQ), ITR, ITEST, NBASDF_v, &ALLOC(LPADE),
                      (MAXITR + 1) / 2,
                      &ALLOC(NELCUP), &ALLOC(LSMATR), &ALLOC(LSMATI), &ILLOC(LINDXS),
                      &ILLOC(LCHNDF), MCHNDF_v, &ILLOC(LBASDF), MBASDF_v,
                      IJ, NCHNDF_v, JT, NLINE, IHEAD, NFLOP_v);
                FLPADE += NFLOP_v;
                TPADE += second() - TT4;
            }

            // Save S-matrices if requested
            if (ISAVSM != 0) {
                int IP = ISAVSM / 2;
                int IS = ISAVSM - 2 * IP;
                int LSMTSV = Z[ISMTSV_g] + 2 * (IS + IP) * NBASDF_v * (ITR - 1) - 1;
                if (IS != 0) {
                    for (int I = 1; I <= NWRDS; I++) {
                        ALLOC(LSMTSV + I) = ALLOC(LSMSQ + I - 1);
                    }
                    LSMTSV += 2 * NBASDF_v;
                }
                if (IP != 0) {
                    IS = 2 * 3 * ((MAXITR + 1) / 2);
                    int IP2 = LPADE + 2 * (3 * ((ITR - 1) / 2) + 3 - 1);
                    for (int I = 1; I <= NBASDF_v; I++) {
                        ALLOC(LSMTSV + 2*I - 1) = ALLOC(IP2 + (I - 1) * IS);
                        ALLOC(LSMTSV + 2*I) = ALLOC(IP2 + (I - 1) * IS + 1);
                    }
                }
            }

            if (ITEST != 0) {
                // Converged — must we do one more?
                if (IWRTWV == 0 || CONVSW) goto L950;
                CONVSW = true;
            }

            // Block: transfer new solutions
            if (BLOKSW) {
                for (int I = 1; I <= NUMINH; I++) {
                    ALLOC(LINHR + NPTMIN_g + I - 1) = (double)ALLOC4(LBLOCK + I);
                    ALLOC(LINHI + NPTMIN_g + I - 1) = (double)ALLOC4(LBLOCK + NUMINH + I);
                }
            }
        } // end ITR loop (899)
        ITR = MAXITR;

    L950:
        {
            double TT1_v = second();
            TITER += TT1_v - TT;
        }
    } // end IJ loop (999)

    // Timing information
    {
        double TFIN = second();
        double TTOT = TFIN - TSTART;
        double TREST = TTOT - TITER - TWAVE;
        double TITER2 = TITER - TRIGHT - TINHO - TNORM - TUOPT - TPADE;
        if (NLINE > 38) {
            std::printf("1%59sP T O L E M Y\n"
                        " COMPUTATION OF COUPLED INELASTIC CHANNELS S-MATRICES"
                        "%41sWIN A FEW, LOOSE A FEW\n"
                        "0%.45sELAB =%7.2f MEV     %.65s\n\n",
                        "", "", &HEDCOM.REACT[1], ELAB, &HEDCOM.HEADER[1]);
        }
        std::printf("0COUPLED EQUATIONS COMPUTATION TIME (SECONDS):\n"
                    "0HOMOGENEOUS SOLUTIONS:%12s%9.3f\n"
                    " DRIVING TERM:%21s%9.3f\n"
                    " CONSTRUCTING OPT. POT.%11s%9.3f\n"
                    " INHOMOGENEOUS SOLUTIONS:%9s%9.3f\n"
                    " NORMALIZATION:%19s%9.3f\n"
                    " PADE AND CONVERG. CHECK:%9s%9.3f\n"
                    " OTHER ITERATION TIME:%12s%9.3f\n"
                    " REST OF COUPLN TIME:%13s%9.3f\n"
                    " TOTAL COUPLN TIME:%23s%9.3f\n\n",
                    "", TWAVE, "", TRIGHT, "", TUOPT, "", TINHO,
                    "", TNORM, "", TPADE, "", TITER2, "", TREST,
                    "", TTOT);

        int NUMITS = NMINHO / NBASDF_v;
        double TPRJ = TTOT / NUMLIS;
        double TPRIT = TITER / (NUMITS > 0 ? NUMITS : 1);
        double TPRRHS = 1000.0 * TRIGHT / (NMINHO > 0 ? NMINHO : 1);
        FLPRHS = 1.0e-6 * FLPRHS;
        double TFLRHS = FLPRHS / (TRIGHT + 1.0e-20);
        double TPRINH = 1000.0 * TINHO / (NMINHO > 0 ? NMINHO : 1);
        FLPINH = 1.0e-6 * FLPINH;
        double TFLINH = FLPINH / (TINHO + 1.0e-20);
        double TPRSTP = 1.0e+6 * TINHO / (NMINPT > 0 ? NMINPT : 1);
        double TPRNRM = 1000.0 * TNORM / (NMINHO > 0 ? NMINHO : 1);
        FLPNRM = 1.0e-6 * FLPNRM;
        double TFLNRM = FLPNRM / (TNORM + 1.0e-20);
        double TPRPAD = 1000.0 * TPADE / (NMINHO > 0 ? NMINHO : 1);
        FLPADE = 1.0e-6 * FLPADE;
        double TFLPAD = FLPADE / (TPADE + 1.0e-20);

        std::printf("0TOTAL TIME PER J-VALUE:%30s%13.3f SECONDS\n"
                    " TIME PER ITERATION:%33s%13.3f SECONDS\n"
                    "%8d RIGHT-HAND SIDES; TIME PER RHS:%9s%13.3f MILLISECONDS"
                    "%18s(%12.3f MFLOP,%10.3f MFLOP/SEC.)\n"
                    "%9sTIME PER INHOMO. SOLUTION:%16s%13.3f MILLISECONDS"
                    "%18s(%12.3f MFLOP,%10.3f MFLOP/SEC.)\n"
                    "%9sTIME PER NORMALIZATION:%19s%13.3f MILLISECONDS"
                    "%18s(%12.3f MFLOP,%10.3f MFLOP/SEC.)\n"
                    "%9sTIME PER PADE & CONVRG:%19s%13.3f MILLISECONDS"
                    "%18s(%12.3f MFLOP,%10.3f MFLOP/SEC.)\n"
                    "%10d NUMEROV STEPS, TIME PER STEP:%11s%13.3f MICROSECONDS\n",
                    "", TPRJ, "", TPRIT,
                    NMINHO, "", TPRRHS, "", FLPRHS, TFLRHS,
                    "", "", TPRINH, "", FLPINH, TFLINH,
                    "", "", TPRNRM, "", FLPNRM, TFLNRM,
                    "", "", TPRPAD, "", FLPADE, TFLPAD,
                    NMINPT, "", TPRSTP);
    }

    // Free arrays
    if ((IPRINT / 100) % 10 >= 1) goto L2000;

    if (BLOKSW) Z[IBLOCK] = -Z[IBLOCK];
    Z[IHOMO_g] = -Z[IHOMO_g];
    Z[IINHR_g] = -Z[IINHR_g];
    Z[IINHI_g] = -Z[IINHI_g];
    Z[IRHSR_g] = -Z[IRHSR_g];
    Z[IRHSI_g] = -Z[IRHSI_g];
    Z[ISMSQ_g] = -Z[ISMSQ_g];
    Z[WAVCOM.IWAVR] = -Z[WAVCOM.IWAVR];
    Z[WAVCOM.IWAVI] = -Z[WAVCOM.IWAVI];
    Z[IPADE_g] = -Z[IPADE_g];
    Z[INUCH_g] = -Z[INUCH_g];
    Z[ICOULH_g] = -Z[ICOULH_g];
    Z[ICL1FF_g] = -Z[ICL1FF_g];
    Z[ICL1FG_g] = -Z[ICL1FG_g];
    Z[WAVCOM.NFS[1]] = -Z[WAVCOM.NFS[1]];
    for (int I = 1; I <= 2; I++) {
        Z[WAVCOM.IRLVS[I]] = -Z[WAVCOM.IRLVS[I]];
        Z[WAVCOM.IIMVS[I]] = -Z[WAVCOM.IIMVS[I]];
        Z[WAVCOM.ICENTR[I]] = -Z[WAVCOM.ICENTR[I]];
    }

L2000:
    {
        int ITR_v = MAXITR; // Last used ITR
        int I = 1;
        if (ISAVSM == 3) I = 2;
        if (ISAVSM != 0) IREDEF(2 * I * NBASDF_v * ITR_v, ISMTSV_g);
    }
    return;
}

// ============================================================================
// DRIVE — Construct driving term of inhomogeneous equation
// Translated from source.f lines 14635-14830 (c8h code path)
// ============================================================================
void DRIVE(double* RHSR, double* RHSI, int I5, int I6, int NFIRST,
           int ITR, int MBASCP_arg, int* BASCUP, float* RBASCP,
           float* NUCH, float* COULH, int NPTIN, int NPTSAV,
           double* INHR, double* INHI, double& FLP, double VMULT)
{
    // GRIDCM CC overlay
    int* GI = reinterpret_cast<int*>(&GRIDCM.IMSVAL);
    int& NPTMIN_g = GI[16];
    int& NUMPT_g  = GI[18];

    // Array access macros (all 1-based, Fortran-style)
    // BASCUP(MBASCP, *) — integer
    #define BASCUP_D(r, c) BASCUP[((c)-1)*MBASCP_arg + ((r)-1)]
    // RBASCP(MBASCP, *) — float
    #define RBASCP_D(r, c) RBASCP[((c)-1)*MBASCP_arg + ((r)-1)]
    // NUCH(2, NPTIN, *) — float, NUCH(i, j, k)
    #define NUCH_D(i, j, k) NUCH[((k)-1)*NPTIN*2 + ((j)-1)*2 + ((i)-1)]
    // COULH(NPTIN, *) — float
    #define COULH_D(i, k) COULH[((k)-1)*NPTIN + ((i)-1)]
    // INHR(NPTSAV, *) and INHI(NPTSAV, *) — c8h: REAL*8
    #define INHR_D(i, k) INHR[((k)-1)*NPTSAV + ((i)-1)]
    #define INHI_D(i, k) INHI[((k)-1)*NPTSAV + ((i)-1)]

    int ISTRT = NFIRST - 1;
    int IFLOP = 0;

    // Zero out RHS
    int NPSTRT = NPTMIN_g + 1;
    int NDO2 = NUMPT_g - (NUMPT_g - NPTMIN_g) % 4;
    for (int I = NPSTRT; I <= NDO2; I += 4) {
        RHSR[I] = 0; RHSI[I] = 0;
        RHSR[I+1] = 0; RHSI[I+1] = 0;
        RHSR[I+2] = 0; RHSI[I+2] = 0;
        RHSR[I+3] = 0; RHSI[I+3] = 0;
    }
    if (NDO2 != NUMPT_g) {
        NDO2++;
        for (int I = NDO2; I <= NUMPT_g; I++) {
            RHSR[I] = 0; RHSI[I] = 0;
        }
    }

    // Loop over all RHS components
    for (int I7 = I5; I7 <= I6; I7++) {
        int IPTWF = BASCUP_D(2, I7);
        double XN = (double)RBASCP_D(5, I7);
        double XC = (double)RBASCP_D(6, I7);
        int ITRI = BASCUP_D(7, I7);
        int ITRAB = std::abs(ITRI);
        double XKF = (double)RBASCP_D(8, I7);

        // IMPORTANT: use fabs check for denormalized floats
        if (std::fabs(XKF) < 1.0e-30) continue;

        int IPTN = BASCUP_D(10, I7);
        if (IPTN == NOTDEF_INT || std::abs(IPTN) > 100000) continue;
        IPTN = (IPTN + 1) / 2;
        int IPTC = BASCUP_D(11, I7);
        if (IPTC == NOTDEF_INT || std::abs(IPTC) > 100000) continue;

        // Factor of 12 corrects simplified Numerov algorithm
        double FAC = XKF * VMULT;
        if (ITRI < 0) FAC = -FAC;
        if (ITR != 1 || IPTWF != 1) FAC = FAC / 12.0;
        XN = FAC * XN;
        XC = FAC * XC;

        if (ITRAB > 1) goto L500;

        // REAL TIME REVERSAL INVARIANCE
        if (XN == 0) goto L200;
        if (XC != 0) goto L150;

        // Nuclear potential only
        IFLOP += 10;
        for (int I = ISTRT; I <= NUMPT_g; I++) {
            RHSR[I] += ((double)NUCH_D(1,I,IPTN) * INHR_D(I,IPTWF)
                      - (double)NUCH_D(2,I,IPTN) * INHI_D(I,IPTWF)) * XN;
            RHSI[I] += ((double)NUCH_D(1,I,IPTN) * INHI_D(I,IPTWF)
                      + (double)NUCH_D(2,I,IPTN) * INHR_D(I,IPTWF)) * XN;
        }
        continue;

    L150:
        // Both Coulomb and nuclear
        IFLOP += 12;
        {
            double XCBYXN = XC / XN;
            for (int I = ISTRT; I <= NUMPT_g; I++) {
                double TERM = (double)NUCH_D(1,I,IPTN) + XCBYXN * (double)COULH_D(I,IPTC);
                RHSR[I] += XN * (TERM * INHR_D(I,IPTWF) - (double)NUCH_D(2,I,IPTN) * INHI_D(I,IPTWF));
                RHSI[I] += XN * (TERM * INHI_D(I,IPTWF) + (double)NUCH_D(2,I,IPTN) * INHR_D(I,IPTWF));
            }
        }
        continue;

    L200:
        // Coulomb only
        if (XC == 0) continue;
        IFLOP += 5;
        for (int I = ISTRT; I <= NUMPT_g; I++) {
            RHSR[I] += (XC * (double)COULH_D(I,IPTC)) * INHR_D(I,IPTWF);
            RHSI[I] += (XC * (double)COULH_D(I,IPTC)) * INHI_D(I,IPTWF);
        }
        continue;

    L500:
        // IMAGINARY TIME REVERSAL INVARIANCE
        if (XN == 0) goto L600;
        if (XC != 0) goto L550;

        // Nuclear only (imaginary TRI)
        IFLOP += 10;
        for (int I = ISTRT; I <= NUMPT_g; I++) {
            RHSR[I] += -((double)NUCH_D(1,I,IPTN) * INHI_D(I,IPTWF)
                       + (double)NUCH_D(2,I,IPTN) * INHR_D(I,IPTWF)) * XN;
            RHSI[I] += ((double)NUCH_D(1,I,IPTN) * INHR_D(I,IPTWF)
                      - (double)NUCH_D(2,I,IPTN) * INHI_D(I,IPTWF)) * XN;
        }
        continue;

    L550:
        // Both Coulomb and nuclear (imaginary TRI)
        IFLOP += 12;
        {
            double XCBYXN = XC / XN;
            for (int I = ISTRT; I <= NUMPT_g; I++) {
                double TERM = (double)NUCH_D(1,I,IPTN) + XCBYXN * (double)COULH_D(I,IPTC);
                RHSR[I] += -XN * (TERM * INHI_D(I,IPTWF) + (double)NUCH_D(2,I,IPTN) * INHR_D(I,IPTWF));
                RHSI[I] += XN * (TERM * INHR_D(I,IPTWF) - (double)NUCH_D(2,I,IPTN) * INHI_D(I,IPTWF));
            }
        }
        continue;

    L600:
        // Coulomb only (imaginary TRI)
        if (XC == 0) continue;
        IFLOP += 5;
        for (int I = ISTRT; I <= NUMPT_g; I++) {
            RHSR[I] += -(XC * (double)COULH_D(I,IPTC)) * INHI_D(I,IPTWF);
            RHSI[I] += (XC * (double)COULH_D(I,IPTC)) * INHR_D(I,IPTWF);
        }
    } // end I7 loop (999)

    FLP = (double)IFLOP * (NUMPT_g - NFIRST + 1);
    #undef BASCUP_D
    #undef RBASCP_D
    #undef NUCH_D
    #undef COULH_D
    #undef INHR_D
    #undef INHI_D
    return;
}

// INELDC → ineldc_translated.cpp

// INELD2 → ineld2_translated.cpp

// INGRST → ingrst_translated.cpp

// INRDIN → inrdin_translated.cpp
