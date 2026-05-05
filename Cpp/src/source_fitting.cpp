// source_fitting.cpp — Fitting, chi-square, L-extrapolation, derivatives
// Translated from source.f: CALFUN+CALGRA, FITEL, FITEL2, FITINP, LINTRP, LMCFUN,
//   LXTRPM, LXTRP1, LXTRP2, MAKDER, SDERIV, SECDER, DERCHK, GENBNX, GETBFC,
//   SETFIT, SETBFC, SETBRN
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

void FITEL(int& IRET)
{ IRET = 0; std::printf(" FITEL: stub — translate from source.f L15647-16315\n"); }

void FITEL2(int& IRET)
{ IRET = 0; std::printf(" FITEL2: stub — translate from source.f L16317-16915\n"); }

void FITINP(int& IRET, int NUMFLT, char8* FLTWRD, int NUMALI, char8* ALIAI)
{ IRET = 0; std::printf(" FITINP: stub — translate from source.f L16917-17283\n"); }

// LINTRP → lintrp_translated.cpp

void LMCFUN(int ITYPE, int MPARA, int MPTS, double* X, double* FS,
            double* FJACOB, int& IERR)
{ IERR = 0; std::printf(" LMCFUN: stub — translate from source.f L24798-24850\n"); }

// 4-argument RTXLNX: finds real solution of  A*X + B*LN(X) + C = 0
extern double RTXLNX(double A, double B, double C, double ACC);

double LXTRPM(int IEXTYP, double BARA, double B, double BARC, double BARL,
              double DLMAX, int LX, int LDEL, double* ETAS, double WEEBOY)
{
    // Returns L such that WEEBOY = SQRT(2L) * EXTRAP FUNC(L)
    // The extrapolation functions are defined in LXTRP1
    // Translated from source.f L25149-25192

    double DL = DLMAX;
    switch (IEXTYP) {
    case 1: goto L100;
    case 2: goto L200;
    case 3: goto L10;
    case 4: goto L300;
    default: goto L10;
    }

L10:
    std::printf("\n*** CANNOT FIND LIMOST FOR FUNCTION TYPE%3d\n", IEXTYP);
    return 2 * DLMAX;

    // Woods Saxon
    //    WEEBOY = SQRT(2L)*BARA / (1 + EXP(B*(L-FLCRIT)) )
L100:
    DL = RTXLNX(B, -0.50, DLOG(WEEBOY / BARA) - .5 * DLOG(2.0) - B * BARL, 0.010);
    goto L900;

    // Power law
    //    WEEBOY = SQRT(2L)*BARA * (DLMAX/L)**B
L200:
    if (B > .5) DL = DEXP((DLOG(1.41 * DABS(BARA) / WEEBOY)
                     + B * DLOG(DLMAX)) / (B - .5));
    goto L900;

    // WKB form
    // WEEBOY = SQRT(2L)*BARA*ETA**LX*DLMAX**B
    //        * EXP((ETA2-ETA1)/ETA * L) / L**(LX+B)
L300:
    {
        double ETA = .5 * (ETAS[1] + ETAS[2]);
        if (LX + B - .5 > 0) DL = RTXLNX((ETAS[2] - ETAS[1]) / ETA,
            LX + B - .5, DLOG(DABS(WEEBOY / (1.4142 * BARA * std::pow(ETA, LX))))
            - B * DLOG(DLMAX), 0.010);
    }

L900:
    return (int)(DL + .5);
}

void LXTRP1(int IEXTYP, int N, int& ISIG, int IPRINT, double* XVALS, double* SVALS,
            double& FLCRIT, double& AVAL, double& WID, double STEPIN, double DELIN,
            double& BARL, double& BARA, double& B, double& BARC, double DLMAX,
            double& CHISQ, int LX, int LDEL, double* ETAS, int MCHN, int JP, int JT)
{
    // Least-squares fit of L-extrapolating functions
    // Translated from source.f L25194-25554

    double AEXP, BEXP, CEXP, CHIEXP;
    double A, C, A1, C1, A2, C2;
    double F, F1, F2, F3, FPL, FMI, FM;
    double X, X1, X2, X3, XPL, XMI, XM, XNEW;
    double S, SM, DX, DF, DTES, TES, F2P, FP;
    double XBAR, ETA, DL, CVAL, POWER;
    int IT, IXTYP2;
    int PBUGSW, PLINSW;

    PBUGSW = IPRINT >= 4;
    PLINSW = IPRINT >= 6;
    double DEL = DELIN;
    double STEP = STEPIN;
    IXTYP2 = 2 * IEXTYP;
    if (IEXTYP == 4) IXTYP2 = 4;
    ISIG = 0;
    double STPMUL = 1;
    double BMIN = .001;
    if (IEXTYP == 4) BMIN = -LX;
    double BMAX;

    // Bias the L's to 0 to avoid overflow problems
    ETA = .5 * (ETAS[1] + ETAS[2]);
    XBAR = (XVALS[1] + XVALS[2]) / 2;
    if (IEXTYP >= 2) XBAR = DLMAX;

    for (int I = 1; I <= N; I++) {
        switch (IEXTYP) {
        case 1:
            XVALS[I] = XVALS[I] - XBAR;
            break;
        case 2:
        case 3:
            XVALS[I] = XBAR / XVALS[I];
            break;
        case 4:
            DL = XVALS[I] + .5 * LDEL;
            SVALS[I] = SVALS[I] / (std::pow(DSIN(DATAN(ETA / DL)), LX)
                * DEXP(-(ETAS[2] - ETAS[1]) * (DATAN(ETA / DL) + DL / ETA)));
            XVALS[I] = XBAR / DL;
            break;
        }
    }

    if (PBUGSW) std::printf("\nEXTRAPOLATION FOR MCHN =%5d%5d%5d%5d%5d   IEXTYP =%5d%5d\n",
                             MCHN, JP, JT, LX, LDEL, IEXTYP, N);
    if (PLINSW) {
        std::printf("\nX AND S:\n");
        for (int I = 1; I <= N; I++)
            std::printf("%4d%17.5G%17.5G\n", I, XVALS[I], SVALS[I]);
    }

    // Get initial guess of the non-linear parameter
    // For IEXTYP=1 use exp fit; for 2 use simple power law.
    // For IEXTYP=3 we just use B = 3
    CVAL = 1;
    if (IEXTYP == 3) CVAL = 3;
    LINLSQ(IXTYP2 - 1, N, XVALS, SVALS, BEXP, AEXP, CVAL, F, PLINSW);
    if (IEXTYP != 3) goto L105;
    CEXP = BEXP;
    IXTYP2 = 5;
    BEXP = 3;
    goto L110;

L105:
    if (BEXP > BMIN) goto L110;
    std::printf("\n***  INITIAL B = %12.6G (OUT OF BOUNDS)--\n", BEXP);
    goto L940;

    // Evaluate chisq and derivative at initial width
L110:
    if (PBUGSW) std::printf("\nB FOR EXP OR POWER FORM = %12.6G\n", BEXP);
    STEP = STEP * DABS(BEXP / 10);
    if (IEXTYP == 3) STEP = .1;
    BMAX = 10 * DABS(BEXP);
    if (IEXTYP == 1) BMAX = 1.2;

    X = BEXP;
    LINLSQ(IXTYP2, N, XVALS, SVALS, C1, A1, X, F, PLINSW);
    CHIEXP = F;

L130:
    XPL = X * (1 + DEL);
    XMI = X * (1 - DEL);
    LINLSQ(IXTYP2, N, XVALS, SVALS, CVAL, AVAL, XPL, FPL, PLINSW);
    LINLSQ(IXTYP2, N, XVALS, SVALS, CVAL, AVAL, XMI, FMI, PLINSW);
    DF = (FPL - FMI) / (2 * X * DEL);
    DTES = DABS(DF / F);
    if (DTES > 1.0e-10) goto L200;
    if (DEL >= .5) goto L700;
    DEL = 2 * DEL;
    goto L130;

L200:
    S = DABS(STEP);
    if (DF > 0) S = -S;

    // Step downhill to box in minimum
    IT = 1;
    F1 = F;
    X1 = X;

L220:
    X2 = X1 + S;

    // A Woods-Saxon with width less than 1 doesn't make sense in this context
    if (IEXTYP == 1 && X2 > 1.2) goto L950;
    if (IEXTYP >= 2 && X2 < BMIN) goto L960;
    if (IEXTYP == 1 && X2 < BMIN) goto L700;
    if (IEXTYP >= 2 && X2 > BMAX) goto L700;

    LINLSQ(IXTYP2, N, XVALS, SVALS, C, A, X2, F2, PLINSW);
    if (F2 <= F1) goto L250;

    // Passed a minimum; try to get three points with F2 the least
    S = .5 * S;
    X3 = X2;
    F3 = F2;
    X2 = X3 - S;
    LINLSQ(IXTYP2, N, XVALS, SVALS, C, A, X2, F, PLINSW);
    if (F < F1) goto L240;
    X3 = X2;
    F3 = F;
    X2 = X1;
    F2 = F1;
    A2 = A1;
    C2 = C1;
    X1 = X2 - S;
    LINLSQ(IXTYP2, N, XVALS, SVALS, C, A, X1, F1, PLINSW);
    goto L300;

L240:
    F2 = F;
    C2 = C;
    A2 = A;
    goto L300;

    // Still looking for the minimum
L250:
    X1 = X2;
    F1 = F2;
    A1 = A;
    C1 = C;
    S = S * STPMUL;
    IT = IT + 1;
    if (IT <= 1000) goto L220;
    std::printf("\n*** MORE THAN 1000 STEPS TO BOX MINIMUM IN B-- \n");
    goto L700;

    // Find minimum of fitted parabola
L300:
    F2P = F1 - 2 * F2 + F3;
    if (F2P == 0) goto L310;
    SM = (.5 * S * (F1 - F3)) / F2P;
    XM = X2 + SM;
    if ((XM <= BMIN) || (XM > BMAX)) goto L700;
    LINLSQ(IXTYP2, N, XVALS, SVALS, CVAL, AVAL, XM, F, PLINSW);
    if (F <= F2) goto L320;

L310:
    XM = X2;
    AVAL = A2;
    CVAL = C2;
    F = F2;

L320:
    B = XM;
    CHISQ = F;
    goto L400;

L350:
    if (PBUGSW) std::printf("\n--- BEST-FIT IS AN EXP OR SIMPLE POWER --\n");
    goto L700;

    // Use local Taylor expansion to refine estimate of minimum
L400:
    DX = DEL * XM;
    FM = F;
    XPL = XM + DX;
    LINLSQ(IXTYP2, N, XVALS, SVALS, C, A, XPL, FPL, PLINSW);
    XMI = XM - DX;
    LINLSQ(IXTYP2, N, XVALS, SVALS, C, A, XMI, FMI, PLINSW);
    FP = (FPL - FMI) / (2 * DX);
    F2P = (FPL - 2 * FM + FMI) / (DX * DX);
    if (F2P == 0) goto L430;
    XNEW = XM - FP / F2P;
    LINLSQ(IXTYP2, N, XVALS, SVALS, C, A, XNEW, FM, PLINSW);
    if (FM >= F) goto L430;
    AVAL = A;
    CVAL = C;
    B = XNEW;
    CHISQ = FM;

L430:
    TES = CVAL + AVAL;
    if (PBUGSW) std::printf("\n-*  (A,C) = (%12.6G %12.6G)\n"
                             "     (B,B2,CH)=(%12.6G %12.6G %12.6G)\n",
                             AVAL, CVAL, B, X2, CHISQ);
    if (TES <= 0 && IEXTYP != 3) goto L350;
    if (PBUGSW) std::printf("\n-COMPLETE IN %6d ITERATIONS\n", IT);

    switch (IEXTYP) {
    case 1: goto L540;
    case 2: goto L550;
    case 3: goto L570;
    case 4: goto L550;
    }

    // Return values for Woods-Saxon
L540:
    BARA = 1 / DABS(AVAL);
    WID = 1 / B;
    BARC = DSIGN(1.0, AVAL);
    FLCRIT = XBAR - WID * DLOG(BARA * CVAL);
    BARL = FLCRIT;
    AVAL = BARA;
    goto L800;

    // Return values for power fit
L550:
    POWER = 1;
    WID = B;
    if (AVAL < 0) goto L700;
    BARA = AVAL;
    BARC = 1;
    BARL = CVAL * XBAR / AVAL;
    FLCRIT = BARL;
    goto L800;

    // Values for phase power (IEXTYP=3)
L570:
    BARA = AVAL;
    BARC = CVAL;
    FLCRIT = CVAL;
    WID = B;
    goto L800;

    // WID-minimum not located within bounds
    // Revert to exponential form
L700:
    if (IEXTYP == 3) goto L730;
    BARA = DEXP(AEXP);
    B = BEXP;
    CHISQ = CHIEXP;
    BARC = 0;
    ISIG = 5;
    if (IEXTYP != 1) goto L720;
    BARL = XBAR;
    AVAL = BARA;
    WID = 1 / BEXP;
    FLCRIT = BARL;
    goto L800;

L720:
    WID = BEXP;
    FLCRIT = 0;
    BARL = 0;
    AVAL = BARA;
    goto L800;

L730:
    AVAL = AEXP;
    CVAL = CEXP;
    B = BEXP;
    goto L570;

    // Make a final test of the correctness of everything
L800:
    if (B > BMIN) goto L810;
    std::printf("\n*** ERROR IN EXTRAP: B OUT OF BOUNDS:%15.5G%15.5G%15.5G%15.5G\n",
                 B, BARA, BARC, BARL);
    goto L940;

L810:
    if (PBUGSW) std::printf("\nLXTRP1 END%3d%13.5G%13.5G%13.5G%13.5G%13.5G%13.5G%13.5G%13.5G%13.5G\n\n",
                             IEXTYP, BARA, B, BARC, BARL, AVAL, WID, CVAL, FLCRIT, CHISQ);

    // Check for slowly decaying extrapolation
    if (MCHN == 1) goto L900;
    switch (IEXTYP) {
    case 1: goto L820;
    case 2: goto L830;
    case 3: goto L900;
    case 4: goto L830;
    }

L820:
    if (B > .01) goto L900;
    goto L890;

L830:
    if (B > BMIN + .6) goto L900;
    goto L890;

L890:
    std::printf("\n*** WARNING:  EXTRAPOLATION DECAYS VERY SLOWLY;"
                "  A, B, C, LCRIT = %15.5G%15.5G%15.5G%15.5G\n",
                 AVAL, B, CVAL, FLCRIT);
    std::printf("      FOR CHANNEL%3d, JP, JT, LX, LO-LI =%3d/2%3d/2%4d%4d\n",
                 MCHN, JP, JT, LX, LDEL);

L900:
    return;

    // One-parameter minimization w.r.t. width has failed
    // Error return; suppress extrapolation
L940:
    std::printf("\n*** MINIMIZATION W.R.T. WIDTH HAS FAILED \n"
                "       EXTRAPOLATION MUST BE SUPPRESSED --\n");
    ISIG = -5;
    return;

L950:
    std::printf("\n**** ATTEMPTING TO FIT VERY SHARP WOODS-SAXON -"
                " ATTEMPT ABANDONED\n");
L955:
    std::printf("      FOR CHANNEL%3d, JP, JT, LX, LO-LI =%3d/2%3d/2%4d%4d\n",
                 MCHN, JP, JT, LX, LDEL);
    goto L700;

L960:
    std::printf("\n*** ATTEMPTING TO FIT WITH POWER OF 0 --"
                " DEFAULT STARTING VALUE USED.\n");
    goto L955;
}

void LXTRP2(int IEXTYP, double BARA, double B, double BARC, double BARL,
            double DLMAX, int LI, double& SIZE, int LX, int LDEL, double* ETAS)
{
    // Computes the L-extrapolation functions
    // See LXTRP1 for the parameters.
    // Except note that for IEXTYP=3, BARA and BARC have been reversed.
    // SIZE is set to the computed value at LI.
    // Translated from source.f L25556-25594

    double DLI = LI;

    switch (IEXTYP) {
    case 1: goto L110;
    case 2: goto L120;
    case 3: goto L130;
    case 4: goto L140;
    }

L110:
    {
        double X = B * (DLI - BARL);
        if (X < -46) X = -46;
        SIZE = BARA / (BARC + DEXP(X));
        return;
    }

L120:
    SIZE = BARA * std::pow(DLMAX / DLI, B) * (1 + BARL / DLI);
    return;

L130:
    SIZE = BARC + BARA * std::pow(DLMAX / DLI, B);
    return;

L140:
    {
        double ETA = .5 * (ETAS[1] + ETAS[2]);
        DLI = DLI + .5 * LDEL;
        SIZE = BARA * std::pow(DLMAX / DLI, B) * (1 + BARL / DLI);
        SIZE = SIZE * std::pow(DSIN(DATAN(ETA / DLI)), LX)
            * DEXP(-(ETAS[2] - ETAS[1]) * (DATAN(ETA / DLI) + DLI / ETA));
        return;
    }
}

void MAKDER(int NSTEPS, int NUMVS, int MXVSIN, int MPARA,
            double* DERIV, double* POTRL, double* POTIM, double STEPSZ)
{ std::printf(" MAKDER: stub — translate from source.f L25596-25748\n"); }

void SDERIV(int& IRTN)
{ IRTN = 0; std::printf(" SDERIV: stub — translate from source.f L31178-31291\n"); }

void SECDER(int& IRTN)
{ IRTN = 0; std::printf(" SECDER: stub — translate from source.f L31293-31510\n"); }

void DERCHK(int NPARA, int NPTS, double STEP, double* PARAM_arr, double* FS,
            double* FGRAD, double* TGRAD, double* FSAVE)
{ std::printf(" DERCHK: stub — translate from source.f L14564-14633\n"); }

void GENBNX(int IPASS, int* NBINDX, int* MBINDX, int* BINDEX, int NBASCP,
            int* BASCP, int MBASCP, int NBASDF, int* BASDF, int MBASDF)
{ std::printf(" GENBNX: stub — translate from source.f L17285-17556\n"); }

void GETBFC(int ILI, int J, double* DELTA, double RASYMP, int NMBFC,
            double* BFACIN, double* BFACOT, int NBASDF, int* BASDF,
            int MBASDF)
{ std::printf(" GETBFC: stub — translate from source.f L17558-17753\n"); }

void SETFIT(int& IRET)
{ IRET = 0; std::printf(" SETFIT: stub — translate from source.f L32576-32878\n"); }

void SETBFC(int NBINDX, int MBINDX, int* BINDEX, int MCHNDF, int* CHNDF_i,
            float* CHNDF_f, int* NUMJS_arr, int LMIN2, double* R2S, int LMXMX,
            double* SIG1, double* SIG2, double DUM, int FLGSW, int IPASS,
            int NMFFAC, float* CL2FF, double* WRKFF, int IDIM1, int LDLDIM,
            double* WRKWK, double* WRKFI, int IDIM2, double* WRKST,
            int& IRET, double& CLTIME)
{ IRET = 0; CLTIME = 0.0; std::printf(" SETBFC: stub — translate from source.f L31512-31798\n"); }

void SETBRN(int NBINDX, int MBINDX, int* BINDEX, int MCHNDF, int* CHNDF_i,
            int MBASDF, int* BASDF_i, int MBASCP, int* BASCP_i, float* BASCP_f,
            int NUMLIS, double* LIS, int LMIN2, float* CL2FF,
            double* INDXS_arr, double* DELSR, double* DELSI)
{ std::printf(" SETBRN: stub — translate from source.f L31800-32000\n"); }
