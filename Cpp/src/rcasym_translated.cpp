// rcasym_translated.cpp — RCASYM subroutine
// Translated from source.f lines 30438-30608
// Asymptotic series for Coulomb functions

#include <cstdio>
#include <cmath>

void RCASYM(int L, double ETA, double RHO, int IPRN, double SIGL,
            double* ZETA_p, double* PHI_p, double* DZETA_p,
            double* F_p, double* FP_p, double* G_p, double* GP_p,
            double* Z, double* DZSQ, double* S, double* ZINV,
            double EPS, int& NMAX, int& NTZ, int& ISIG)
{
    static const double PI = 3.14159265358979300;

    double ZETA, DZETA, PHI, FLL, ETASQ, RHOSQ, BF1, BF2, RHOT;
    double CF, BF, TF, D2F, FACTOR, RTZ;
    int I, J, JTOP, JHI, JJTOP, ILOC, KSER;

    ISIG = 0;
    FLL = (double)L * (L + 1);
    ETASQ = ETA * ETA;
    RHOSQ = RHO * RHO;
    BF1 = 2.0 * (ETA / RHO);
    BF2 = FLL / RHOSQ;
    if (ETA > 0) RHOT = ETA * (1.0 + std::sqrt(1.0 + FLL / ETASQ));
    if (ETA == 0) RHOT = std::sqrt(FLL);
    if (ETA == 0 && L == 0) RHOT = 0.5;
    if (ETA < 0) RHOT = std::fabs(ETA) * (std::sqrt(1.0 + FLL / ETASQ) - 1.0);
    if (RHO <= RHOT) {
        if (IPRN > -4) std::printf("\n -- RCASYM. ASYMPTOTIC EXPANSION INSIDE TP --\n");
        ISIG = -5;
        return;
    }

    Z[1] = 1.0;
    S[1] = 1.0;
    DZSQ[1] = 0.0;
    ZINV[1] = 1.0;
    ZETA = Z[1];
    DZETA = 0.0;
    PHI = RHO - ETA * std::log(2.0 * RHO) + SIGL - 0.5 * L * PI;
    KSER = 2;
    NTZ = 1;
    I = 2;

L500:
    Z[I] = 0.0;
    JTOP = (I + 1) / 2;
    S[I] = 0.0;
    DZSQ[I] = 0.0;
    CF = 0.0;
    if (I == 2) goto L515;
    for (J = 2; J <= JTOP; J++) {
        FACTOR = 2.0;
        if (J == (I - J + 1)) FACTOR = 1.0;
        S[I] = S[I] + FACTOR * Z[J] * Z[I - J + 1];
    }
L515:
    JHI = I - 1;
    for (J = 1; J <= JHI; J++)
        CF = CF + Z[J] * S[I - J + 1];
    BF = -BF1 * Z[I - 1];
    if (I > 2) BF = BF - BF2 * Z[I - 2];
    if (I <= 5) goto L545;
    ILOC = I - 4;
    JJTOP = ILOC - 1;
    ZINV[ILOC] = 0.0;
    for (J = 1; J <= JJTOP; J++)
        ZINV[ILOC] = ZINV[ILOC] - ZINV[J] * Z[ILOC - J + 1];
L545:
    TF = 0.0;
    if (I <= 4) goto L560;
    for (J = 3; J <= JTOP; J++) {
        FACTOR = 2.0 / RHOSQ;
        if (J == (I - J + 1)) FACTOR = 1.0 / RHOSQ;
        DZSQ[I] = DZSQ[I] + FACTOR * (J - 2) * (I - J - 1) * Z[J - 1] * Z[I - J];
    }
    JJTOP = I - 4;
    for (J = 1; J <= JJTOP; J++)
        TF = TF + ZINV[J] * DZSQ[I - J + 1];
L560:
    D2F = 0.0;
    if (I > 3) D2F = (double)(I - 2) * (I - 3) * (Z[I - 2] / RHOSQ);
    Z[I] = (-CF + BF + 0.75 * TF - 0.5 * D2F) * 0.5;
    if (IPRN >= 5) std::printf("%5d Z = %13.6G%13.6G%13.6G%13.6G%13.6G\n",
        I, Z[I], CF, BF, TF, D2F);
    NTZ = NTZ + 1;
    S[I] = S[I] + 2.0 * Z[I];
    ZETA = ZETA + Z[I];
    if (std::fabs(Z[I]) >= EPS) goto L570;
    KSER = KSER - 1;
L570:
    DZETA = DZETA - (double)(I - 2) * Z[I - 1] / RHO;
    if (I == 2 || std::fabs((double)(I - 2) * Z[I - 1] / RHO) >= EPS) goto L580;
    KSER = KSER - 1;
L580:
    if (I == 2) goto L590;
    PHI = PHI - (Z[I] * RHO) / (double)(I - 2);
L590:
    if (KSER <= 0) goto L600;
    I = I + 1;
    if (I > NMAX) goto L1000;
    if (I < 8) goto L500;
    if (std::fabs(Z[I - 1]) > std::fabs(Z[I - 7])) goto L1020;
    goto L500;

L600:
    if (IPRN > 0) std::printf("\n WITH %2d TERMS,  ZETA = %22.14G\n DZETA = %22.14G PHI = %22.14G\n",
        NTZ, ZETA, DZETA, PHI);

    RTZ = std::sqrt(ZETA);
    *F_p = std::sin(PHI) / RTZ;
    *G_p = std::cos(PHI) / RTZ;
    FACTOR = (0.5 * DZETA) / (ZETA * ZETA);
    *FP_p = ZETA * (*G_p) - FACTOR * (*F_p);
    *GP_p = -ZETA * (*F_p) - FACTOR * (*G_p);

    *ZETA_p = ZETA;
    *PHI_p = PHI;
    *DZETA_p = DZETA;
    return;

L1000:
    if (IPRN > -4) std::printf(" **** MORE THAN %4d TERMS NEEDED IN SERIES FOR ZETA ****\n", NMAX);
    goto L1050;

L1020:
    if (IPRN > -4) std::printf(" **** COULOMB SERIES TURNS AT%6d TERMS ****\n", I);

L1050:
    if (IPRN > -4) std::printf(" **** L, ETA, RHO, ZETA, LAST = %4d%15.5G%15.5G%15.5G%15.5G\n",
        L, ETA, RHO, ZETA, Z[I - 1]);
    ISIG = -5;
    return;
}
