// rcwfn.cpp -- translated from rcwfn.f
// Computes regular and irregular Coulomb wavefunctions.

#include "ptolemy_types.h"
#include "ptolemy_intrinsics.h"
#include <cstdio>

void RCWFN(double RHO, double ETA, int MINL, int MAXL,
           double* FC_arg, double* FCP_arg, double* GC_arg, double* GCP_arg,
           double ACCUR, int& IRET)
{
    // Fortran arrays are 1-based: FC(1)..FC(LMAX+1). The caller passes
    // a 0-based pointer where [0]=L_min position. Adjust to 1-based
    // so that FC[L+1] in C++ = Fortran FC(L+1).
    double* FC  = FC_arg  - 1;
    double* FCP = FCP_arg - 1;
    double* GC  = GC_arg  - 1;
    double* GCP = GCP_arg - 1;
    // implicit real*8 (a-h, o-z)
    double ACC, AR, AI, B1, B2, B3, BR, BI, C, D, DEL, DELINV;
    double DFACTR, DI, DP, DQ, DR, ETA2, F, FMAG, G, GP, H;
    double P, PL, PLSAVE, Q, R, RHOUSE, S, STOP, SUM, T, TURN;
    double W, WI, X, XLL1, ACCR;
    int I, IGOTO, L, LMAX, LMIN, LMIN1, LL, LP;
    bool FRSTSW;

    // Machine-dependent constants (IEEE / RS6000)
    double VRYBIG = 1.79769e+308;
    double BIG    = 1.0e+300;
    double SMALL  = 1.0e-300;
    double SMALLN = -690.775527898214;
    double PRERT3 = 2.7e-5;
    double PRECIS = 2.0e-14;
    double PRECLN = 32.0e0;

    double PI = 3.14159265358979323800e+00;

    // Here we limit accuracy to reasonable values for the machine

    ACC = ACCUR;
    ACC = DMAX1(ACC, PRECIS);
    ACC = DMIN1(ACC, PRERT3);

    LMAX = MAXL;
    LMIN = MINL;
    LMIN1 = LMIN + 1;
    XLL1 = (double)LMIN * (double)LMIN1;
    ETA2 = ETA * ETA;

    // Determine which region we are in
    // For RHO < .45, Q of P+IQ is poorly determined so we don't use it.
    // Except that for large negative ETA the Maclauren series also has problems

    if (RHO > 0.45e0) goto label_20;
    if (ETA >= 0) goto label_10;
    if (-ETA * RHO > 7) goto label_20;

    // For RHO < .005, we only return F and F' since the P+IQ recursion
    // is very slowly convergent

label_10:
    if (RHO > 0.005e0) goto label_60;
    IGOTO = 5;
    goto label_70;

label_20:
    TURN = ETA + DSQRT(ETA2 + XLL1);
    IGOTO = 1;
    if (RHO >= TURN - 1.0e-4) goto label_100;

    // We are inside the turning point for MINL, can we get outside
    // of it by reducing MINL. (This is always possible for ETA < 0).

    if (RHO < ETA + DABS(ETA)) goto label_60;

    // Yes, do so

    LMIN = (int)(0.5 * (DSQRT(1 + 4 * ((RHO - ETA) * (RHO - ETA) - ETA2)) - 1));
    LMIN1 = LMIN + 1;
    goto label_80;

    // Must use a different method to supplement the bad I Q value.
    // Always start with LMIN = 0 for simplicity.
    // Note only ETA > 0 gets to here (except when RHO < .45)

label_60:
    IGOTO = 2;
label_70:
    LMIN = 0;
    LMIN1 = 1;
    if (ETA < 10 || RHO <= ETA) goto label_80;
    IGOTO = 3;

label_80:
    XLL1 = (double)LMIN * (double)LMIN1;

    // Here we compute F'/F for L = MAXL
    // We then recurse down to LMIN to generate the unnormalized F's
    // This section is used for all RHO.

label_100:
    PL = LMAX + 1;
    RHOUSE = RHO;
label_105:
    PLSAVE = PL;
label_110:
    FRSTSW = true;
    // Continued fraction for R = FP(MAXL)/F(MAXL)
    R  = ETA / PL + PL / RHOUSE;
    DQ = (ETA * RHOUSE) * 2.0 + 6 * PL * PL;
    DR = 12 * PL + 6;
    DEL = 0.0;
    D   = 0.0;
    F   = 1.0;
    X   = (PL * PL - PL + (ETA * RHOUSE)) * (2.0 * PL - 1.0);
    AI  = RHOUSE * PL * PL;
    DI  = (2 * PL + 1) * RHOUSE;

    // Loop and converge on R

    for (I = 1; I <= 100000; I++) {
        H = (AI + RHOUSE * ETA2) * (RHOUSE - AI);
        X = X + DQ;
        D = D * H + X;

        // If we pass near a zero of the divisor, start over at larger LMAX

        if (DABS(D) > PRERT3 * DABS(DR)) goto label_130;
        PL = PL + 1;
        if (PL < PLSAVE + 10) goto label_110;
        IRET = 5;
        goto label_990;

    label_130:
        D = 1 / D;
        DQ = DQ + DR;
        DR = DR + 12;
        AI = AI + DI;
        DI = DI + 2 * RHOUSE;
        DEL = DEL * (D * X - 1.0);
        if (FRSTSW) DEL = -RHOUSE * (PL * PL + ETA2) * (PL + 1.0) * D / PL;
        FRSTSW = false;
        R = R + DEL;
        if (D < 0.0) F = -F;
        if (DABS(DEL) < DABS(R * ACC)) goto label_140;
    } // end for I (label 139)
    IRET = 6;
    goto label_990;

    // R has converged; did we increase LMAX

label_140:
    if (PL == PLSAVE) goto label_160;

    // Recurse down on R to LMAX
    // Here the only part of F that is of interest is the sign

    PL = PL - 1;
label_150:
    D = ETA / PL + PL / RHOUSE;
    F = (R + D) * F;
    R = D - (1 + ETA2 / (PL * PL)) / (R + D);
    PL = PL - 1;
    if (PL > PLSAVE) goto label_150;

    // Now have R(LMAX, RHO) or if IGOTO=4, R(LMIN, 2*ETA)

label_160:
    if (IGOTO == 4) goto label_210;
    FC [LMAX + 1] = F;
    FCP[LMAX + 1] = F * R;
    if (LMAX == LMIN) goto label_200;
    // Downward recursion to LMIN for F and FP, arrays GC,GCP are storage
    L  = LMAX;
    PL = LMAX;
    AR = 1 / RHO;
    for (LP = LMIN1; LP <= LMAX; LP++) {
        GC [L + 1] = ETA / PL + PL * AR;
        GCP[L + 1] = DSQRT((ETA / PL) * (ETA / PL) + 1);
        FC [L]     = (GC[L + 1] * FC[L + 1] + FCP[L + 1]) / GCP[L + 1];
        FCP[L]     =  GC[L + 1] * FC[L]     - GCP[L + 1] * FC[L + 1];
        PL = PL - 1;
        L  = L - 1;

        // If we are getting near an overflow, renormalize everything down

        if (DABS(FC[L + 1]) < BIG) goto label_189;
        for (LL = L; LL <= LMAX; LL++) {
            FC [LL + 1] = SMALL * FC [LL + 1];
            FCP[LL + 1] = SMALL * FCP[LL + 1];
        } // end for LL (label 179)
    label_189:;
    } // end for LP (label 189)
    F = FC [LMIN1];
    R = FCP[LMIN1] / F;

    // Here we find  P + IQ = (G'+IF')/(G+IF)
    // This section is used in all cases except when
    //    15 < ETA < RHO < 2*ETA

label_200:
    if (IGOTO == 3) goto label_500;
    if (IGOTO == 5) goto label_400;

    // Now obtain P + I.Q for LMIN from continued fraction (32)
    // Real arithmetic to facilitate conversion to IBM using REAL*8
label_210:
    P  = 0.0;
    Q  = RHOUSE - ETA;
    PL = 0.0;
    AR = -(ETA2 + XLL1);
    AI = ETA;
    BR = Q + Q;
    BI = 2.0;
    WI = ETA + ETA;
    DR =  BR / (BR * BR + BI * BI);
    DI = -BI / (BR * BR + BI * BI);
    DP = -(AR * DI + AI * DR);
    DQ =  (AR * DR - AI * DI);

    // Loop and converge on P + IQ

label_230:
    P  = P + DP;
    Q  = Q + DQ;
    PL = PL + 2.0;
    AR = AR + PL;
    AI = AI + WI;
    BI = BI + 2.0;
    D  = AR * DR - AI * DI + BR;
    DI = AI * DR + AR * DI + BI;
    T  = 1.0 / (D * D + DI * DI);
    DR =  T * D;
    DI = -T * DI;
    H  = BR * DR - BI * DI - 1.0;
    X  = BI * DR + BR * DI;
    T  = DP * H  - DQ * X;
    DQ = DP * X  + DQ * H;
    DP = T;
    if (PL > 46000.0) goto label_920;
    if (DABS(DP) + DABS(DQ) >= (DABS(P) + DABS(Q)) * ACC) goto label_230;
    P = P / RHOUSE;
    Q = Q / RHOUSE;

    // We now have R and P+IQ, is this enough

    if (IGOTO == 2) goto label_400;

    // Solve for FP,G,GP and normalise F at L=LMIN
    // Since this is for RHO > RHO(TURN), F and G are reasonable numbers

    X = (R - P) / Q;
    FMAG = DSQRT(1 / (Q * (1 + X * X)));
    W = FMAG / DABS(F);
    F = W * F;
    G = F * X;
    GP = R * G - 1 / F;
    if (IGOTO == 4) goto label_600;
    goto label_800;

    // Here RHO < ETA or RHO < 2*ETA < 20 or RHO < .45
    // We use the Maclauren series to get F(L=0, ETA, RHO)
    // First compute RHO*C(L=0, ETA)

label_400:
    C = 2 * PI * ETA;
    if (DABS(C) > 0.5) goto label_410;

    // Use Maclaurin expansion of X / (EXP(X)-1)

    X = 0;
    T = 1;
    AR = 1;
    BR = C;
    AI = 1;
    C = 1;
label_405:
    AI = AI + 1;
    AR = AR * BR / AI;
    C = C + AR;
    if (DABS(AR) >= ACC * C) goto label_405;
    C = 1 / C;
    goto label_430;

    // Here ETA is not tiny.

label_410:
    if (ETA > 0) goto label_420;
    C = -C;
    X = 0;
    T = 1;
    goto label_425;
label_420:
    X = -SMALLN - PI * ETA;
    T = SMALL;
label_425:
    if (C < PRECLN) C = C / (1 - DEXP(-C));
label_430:
    C = RHO * DSQRT(C);
    B1 = 1;
    B2 = ETA * RHO;
    SUM = B1 + B2;
    AI = 6;
    DI = 6;
    for (I = 1; I <= 10000; I++) {
        B3 = ((2 * ETA * RHO) * B2 - (RHO * RHO) * B1) / AI;
        AI = AI + DI;
        DI = DI + 2;
        SUM = SUM + B3;
        STOP = DABS(B1) + DABS(B2) + DABS(B3);
        B1 = B2;
        B2 = B3;
        if (DABS(SUM) < BIG) goto label_445;
        X = X - SMALLN;
        SUM = SUM * SMALL;
        B1 = B1 * SMALL;
        B2 = B2 * SMALL;
    label_445:
        if (STOP < ACC * DABS(SUM)) goto label_450;
    } // end for I (label 449)
    IRET = 8;
    goto label_990;

label_450:
    SUM = (C * DEXP(X) * SUM) * T;

    // Did it underflow

    if (SUM == 0) goto label_900;

    // We now have F (=SUM), R, and P (P only if RHO > .005)
    // Use the Wronskian as the 4th condition

    W = SUM / F;
    F = SUM;
    if (IGOTO == 5) goto label_850;
    X = (R - P) * F;
    if (DABS(X) > PRERT3) goto label_480;

    // Here F**3 and F**4 terms are less than machine precision

    G = 1 / X;
    GP = P * G;
    goto label_800;

    // Here we must include F**3, F**4; we must also worry about
    // which sign of the root to use.
    // The positive root applies for G > F; else the negative root
    // We use Q in determining which is correct

label_480:
    B1 = 0.5 / X;
    B2 = B1 * DSQRT(1 - 4 * (X * F) * (X * F));
    G = B1 + B2;
    // G > F in all of region 2 for ETA > 0
    if (ETA >= 0) goto label_490;
    SUM = 1 / Q - F * F;
    GP = B1 - B2;
    if (DABS(G * G - SUM) > DABS(GP * GP - SUM)) G = GP;
label_490:
    GP = P * G - X * F / G;
    goto label_800;

    // ETA > 15 and ETA < RHO < 2*ETA
    // We find G and G' for LMIN, RHO=2*ETA using the above method
    // consisting of R, P+IQ, and W.

label_500:
    RHOUSE = ETA + ETA;
    PL = LMIN + 1;
    IGOTO = 4;
    goto label_105;

    // Now we have G, G' at the turning point, go in using Taylor

label_600:
    DEL = RHOUSE - RHO;
    B1 = G;
    B2 = -DEL * GP;
    B3 = 0;
    G = B1 + B2;
    ACCR = ACC / 2;
    DELINV = -1 / DEL;
    DFACTR = 3 * DELINV;
    X = DEL / RHOUSE;
    AI = X + X;
    DI = AI + AI;
    AR = 6;
    DR = 6;
    for (I = 1; I <= 10000; I++) {
        S = (AI * B3 + (X * DEL * DEL) * B1) / AR;
        AR = AR + DR;
        DR = DR + 2;
        AI = AI + DI;
        DI = DI + 2 * X;
        G = G + S;
        GP = GP + DFACTR * S;
        if (G >= VRYBIG) goto label_900;
        DFACTR = DFACTR + DELINV;
        B1 = B2;
        B2 = B3;
        B3 = S;
        if (S < ACCR * G) goto label_650;
    } // end for I (label 639)
    IRET = 9;
    goto label_990;

    // Here we have R = F'/F, G, G'
    // Use Wronskian as the 4th condition

label_650:
    F = FC[LMIN1];
    R = FCP[LMIN1] / F;
    SUM = 1 / (R * G - GP);
    W = SUM / F;
    F = SUM;

    // We now have F, R = F'/F, G, G' at LMIN
    // Upward recursion from GC(LMIN) and GCP(LMIN), stored values are RHO
    // Renormalise FC,FCP for each L-value

label_800:
    GC [LMIN1] = G;
    GCP[LMIN1] = GP;
    FC [LMIN1] = F;
    FCP[LMIN1] = R * F;
    IRET = 0;
    if (LMAX == LMIN) return;
    for (L = LMIN1; L <= LMAX; L++) {
        T         = GC[L + 1];
        GC [L + 1] = (GC[L] * GC [L + 1] - GCP[L]) / GCP[L + 1];
        GCP[L + 1] =  GC[L] * GCP[L + 1] - GC[L + 1] * T;
        FC [L + 1] = W * FC [L + 1];
        FCP[L + 1] = W * FCP[L + 1];
    } // end for L (label 829)
label_840:
    if (DABS(FC[LMAX + 1]) + DABS(FCP[LMAX + 1]) == 0) IRET = IRET + 1;
    return;

    // RHO < .005; we cannot find P or Q and so return only F, F'.

label_850:
    FC [LMIN1] = F;
    FCP[LMIN1] = R * F;
    IRET = 2;
    if (LMAX == LMIN) return;
    for (L = LMIN1; L <= LMAX; L++) {
        FC [L + 1] = W * FC [L + 1];
        FCP[L + 1] = W * FCP[L + 1];
    } // end for L (label 859)
    goto label_840;

    // F and G are out of the machine exponent range for LMIN.
    // It will be even worse for L > LMIN so give up and return

label_900:
    IRET = 4;
    goto label_990;

    // P + IQ failed to converge
label_920:
    IRET = 7;

    // Error detected, print the input and some stuff

label_990:
    printf("0***RCWFIN IRET = %5d      INPUT = %20.10G%20.10G%10d%10d%15.5G\n",
           IRET, RHO, ETA, MINL, MAXL, ACCUR);
    printf(" ***%18.8G%18.8G%18.8G%18.8G%18.8G%18.8G\n",
           RHOUSE, P, Q, R, T, X);
    printf(" %7d%18.8G\n", I, SUM);
    return;
}
