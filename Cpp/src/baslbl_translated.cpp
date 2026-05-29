// BASLBL — translated from source.f lines 3216-3514
// Constructs basis states definitions for coupled channels

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include "ptolemy_alloc.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_forward.h"
#include "ptolemy_io.h"
#include <cstdio>
#include <cstring>

void BASLBL()
{
    // Local variables
    extern int GIVEAL(int IWRDS);
    extern int CHOPIT(int IWRDS);

    char8 NAME9("BASISSTA");
    int LENDAT = 7;
    int LPL = 3;
    int JLTP = 10000;
    char SIGNS[3] = {' ', '+', '-'}; // 1-based

    int LCHNDF = ALLOCS.FACFR4 * (LOCPTRS.Z[CCBLK.ICHNDF] - 1) + 1;
    int NCHNDF = ILLOC(LCHNDF);
    int MCHNDF = ILLOC(LCHNDF + 1);
    int MCHNVL = 2;
    CCBLK.MCHNVL = MCHNVL;
    CCBLK.ICHNVL = NALLOC("CHANVALS", MCHNVL * NCHNDF);
    int I = NAMLOC(NAME9.data);
    if (I != 0) LOCPTRS.Z[I] = -LOCPTRS.Z[I];
    int IWRDS;
    CCBLK.IBASDF = GIVEAL(IWRDS);
    NAMCOM.NAMES[CCBLK.IBASDF] = NAME9;
    int LBASDF = ALLOCS.FACFR4 * (LOCPTRS.Z[CCBLK.IBASDF] - 1) + 1 + LPL;
    LCHNDF = ALLOCS.FACFR4 * (LOCPTRS.Z[CCBLK.ICHNDF] - 1) + 1;
    double AWN = 2.0 * KANDM.REDMI / (CNSTNT.HBARC * CNSTNT.HBARC);
    LCHNDF = ILLOC(LCHNDF + 2) + LCHNDF;
    int LCHNVL = LOCPTRS.Z[CCBLK.ICHNVL];

    // Accumulate max|LIN-LOUT|
    INELCM.LXMAX = 0;
    INELCM.LXMIN = 0;
    KANDM.QMIN = 0;

    for (int II = 1; II <= NCHNDF; II++) {
        // Labels for ETA and K
        int I1 = LCHNDF + 5 + (II - 1) * MCHNDF;
        int I2 = I1 + 1;

        // Calculate ETA and K in each channel
        double AKC = AWN * (ALLOC4(I1 - 1) + ALLOC4(I1 - 3));
        AKC = DSQRT(KANDM.AKI * KANDM.AKI - AKC);
        double ETAC = KANDM.ETAS[1] * KANDM.AKI / AKC;
        ALLOC(LCHNVL + MCHNVL * (II - 1)) = AKC;
        ALLOC(LCHNVL + MCHNVL * (II - 1) + 1) = ETAC;
        ALLOC4(I1) = ETAC;
        ALLOC4(I2) = AKC;
    }

    // In this version we assume G.S. spins=0+
    // and hence elastic channel parities = (-1)**L(EL)
    int I1 = 0;
    int I3 = 0;
    ILLOC(LCHNDF + 7) = 1;

    int LABEL2 = 0;

    // Loop over channels
    for (int II = 1; II <= NCHNDF; II++) {
        // Don't forget NXHINT gives 2*JP(JT)
        int JP = IABS(ILLOC(LCHNDF + (II - 1) * MCHNDF + 3));
        int JT = IABS(ILLOC(LCHNDF + (II - 1) * MCHNDF + 1));
        double EP = ALLOC4(LCHNDF + (II - 1) * MCHNDF + 4);
        double ET = ALLOC4(LCHNDF + (II - 1) * MCHNDF + 2);
        I1 = I1 + 1;
        int I2 = 0;
        int J1 = IABS(JP - JT);
        int J2 = JP + JT;
        int I5 = 0;
        int J4 = J2 / 2;
        if (BTEST(J4, 0)) I5 = 1;
        int I4 = LCHNDF + (II - 1) * MCHNDF + 8;
        INELCM.LXMAX = MAX0(INELCM.LXMAX, J2 / 2);

        // Loop over channel spins
        for (int J = J1; J <= J2; J += 2) {
            int J3 = J / 2;

            // For identical spin-0 bosons, omit odd channel spins
            if (WAVCOM.LSKIPS[2] == 2) {
                if (JP == JT && EP == ET) {
                    if (BTEST(J3, 0)) continue;
                }
            }

            // Set up symmetry indicator
            int ISYM = 1;
            if (WAVCOM.ISTATS[2] != 3) {
                if (JP != JT || EP != ET) ISYM = 2;
            }

            int L1 = IABS(JLTP - J);
            int L2 = JLTP + J;
            if (I5 == 1) goto L27;
            if (BTEST(J3, 0)) L1 = L1 + 2;
            if (BTEST(J3, 0)) L2 = L2 - 2;
            goto L28;
        L27:
            J3 = J3 + 1;
            if (BTEST(J3, 0)) L1 = L1 + 2;
            if (BTEST(J3, 0)) L2 = L2 - 2;

            // Loop over L-J
        L28:
            for (int L = L1; L <= L2; L += 4) {
                I2 = I2 + 1;
                I3 = I3 + 1;
                LABEL2 = (I3 - 1) * LENDAT + LBASDF;
                ILLOC(LABEL2) = I1;
                ILLOC(LABEL2 + 1) = J;
                ILLOC(LABEL2 + 2) = L - JLTP;
                ILLOC(LABEL2 + 6) = ISYM;
            }
            ILLOC(I4) = I2;
            if (II == NCHNDF) continue;
            ILLOC(I4 + MCHNDF - 1) = I2 + ILLOC(I4 - 1);
        }
    }

    KANDM.LOMOST = INTGER.LMAX + INELCM.LXMAX;

    // Print channel information
    std::printf("\n-  CHANNEL%*sPROJECTILE%*sTARGET"
                "%*sQ%*sETA%*sK%*sNUMBER OF\n",
                10, "", 21, "", 13, "", 2, "", 6, "", 2, "");
    std::printf("%8s%6sE(MEV)%4sJ%5sK %6sE(MEV)%4sJ%5sK "
                "%*sMEV%*s1/FM%*sEQUATIONS\n\n",
                "", "", "", "", "", "", "",
                5, "", 14, "", 3, "");

    int L = 0;
    for (int II = 1; II <= NCHNDF; II++) {
        int J = LCHNDF + (II - 1) * MCHNDF;
        int J1 = ILLOC(J);
        int JP = ILLOC(J + 3);
        int IP = 1;
        if (JP < 0) IP = 2;
        JP = IABS(JP);
        int JT = ILLOC(J + 1);
        int IT = 1;
        if (JT < 0) IT = 2;
        JT = IABS(JT);
        double EP = ALLOC4(J + 4);
        double ET = ALLOC4(J + 2);
        double Q = -(EP + ET);
        KANDM.QMIN = DMIN1(KANDM.QMIN, Q);
        int N = ILLOC(J + 8);
        L = L + N;
        std::printf(" %7d%12.4f%4d/2%c%3d/2%12.4f%4d/2%c%3d/2%12.3f%10.3f%10.3f%10d\n",
                    J1, EP, JP, SIGNS[IP], ILLOC(J + 16),
                    ET, JT, SIGNS[IT], ILLOC(J + 15),
                    Q, (double)ALLOC4(J + 5), (double)ALLOC4(J + 6), N);
    }
    std::printf("\nWE MUST SOLVE%5d COUPLED EQUATIONS\n\n", L);
    ILLOC(LBASDF - LPL) = L;
    ILLOC(LBASDF - LPL + 1) = LENDAT;
    ILLOC(LBASDF - LPL + 2) = LPL;
    LABEL2 = LABEL2 + LENDAT;
    IWRDS = (LABEL2 + LPL - LBASDF) / ALLOCS.FACFR4 + 1;
    int IL = CHOPIT(IWRDS);
}
