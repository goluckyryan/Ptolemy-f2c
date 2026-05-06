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

// FUNCTION LXTRPM — source.f L25149-25192
double LXTRPM(int IEXTYP, double BARA, double B, double BARC, double BARL,
              double DLMAX, int LX, int LDEL, double* ETAS, double WEEBOY)
{
    double DL = DLMAX;
    double ETA;
    switch (IEXTYP) {
    case 1: // Woods-Saxon
        DL = RTXLNX(B, -0.50, DLOG(WEEBOY/BARA) - 0.5*DLOG(2.0) - B*BARL, 0.010);
        goto L900;
    case 2: // Power law
        if (B > 0.5) DL = DEXP((DLOG(1.41*DABS(BARA)/WEEBOY) + B*DLOG(DLMAX)) / (B - 0.5));
        goto L900;
    case 3: // Error — cannot find LIMOST for phase type
        std::printf("\n*** CANNOT FIND LIMOST FOR FUNCTION TYPE%3d\n", IEXTYP);
        return 2*DLMAX;
    case 4: // WKB form — falls through from case 2 label to case 4 code
        ETA = 0.5*(ETAS[1]+ETAS[2]);
        if (LX+B-0.5 > 0) DL = RTXLNX((ETAS[2]-ETAS[1])/ETA, LX+B-0.5,
            DLOG(DABS(WEEBOY/(1.4142*BARA*std::pow(ETA,(double)LX)))) - B*DLOG(DLMAX), 0.010);
        goto L900;
    }
L900:
    return DL + 0.5;
}

// SUBROUTINE LXTRP1 — source.f L25194-25554
void LXTRP1(int IEXTYP, int N, int& ISIG, int IPRINT, double* XVALS, double* SVALS,
            double& FLCRIT, double& AVAL, double& WID, double STEPIN, double DELIN,
            double& BARL, double& BARA, double& B, double& BARC, double DLMAX,
            double& CHISQ, int LX, int LDEL, double* ETAS, int MCHN, int JP, int JT)
{
    int PBUGSW = IPRINT >= 4;
    int PLINSW = IPRINT >= 6;
    double DEL = DELIN;
    double STEP = STEPIN;
    int IXTYP2 = 2*IEXTYP;
    if (IEXTYP == 4) IXTYP2 = 4;
    ISIG = 0;
    double STPMUL = 1;
    double BMIN = 0.001;
    if (IEXTYP == 4) BMIN = -LX;

    double ETA = 0.5*(ETAS[1]+ETAS[2]);
    double XBAR = (XVALS[1]+XVALS[2]) / 2.0;
    if (IEXTYP >= 2) XBAR = DLMAX;

    for (int I = 1; I <= N; I++) {
        switch (IEXTYP) {
        case 1: XVALS[I] = XVALS[I] - XBAR; break;
        case 2: case 3: XVALS[I] = XBAR / XVALS[I]; break;
        case 4: {
            double DL = XVALS[I] + 0.5*LDEL;
            SVALS[I] = SVALS[I] / (std::pow(DSIN(DATAN(ETA/DL)), (double)LX)
                * DEXP(-(ETAS[2]-ETAS[1])*(DATAN(ETA/DL)+DL/ETA)));
            XVALS[I] = XBAR / DL;
            } break;
        }
    }

    if (PBUGSW) std::printf("\nEXTRAPOLATION FOR MCHN =%5d%5d%5d%5d%5d   IEXTYP =%5d%5d\n",
        MCHN, JP, JT, LX, LDEL, IEXTYP, N);
    if (PLINSW) { std::printf("\nX AND S:\n"); for(int I=1;I<=N;I++) std::printf("%4d%17.5G%17.5G\n",I,XVALS[I],SVALS[I]); }

    double CVAL = 1;
    if (IEXTYP == 3) CVAL = 3;
    double BEXP=0, AEXP=0, CEXP=0, F=0, CHIEXP=0;
    LINLSQ(IXTYP2-1, N, XVALS, SVALS, BEXP, AEXP, CVAL, F, PLINSW);
    if (IEXTYP == 3) { CEXP = BEXP; IXTYP2 = 5; BEXP = 3; goto L110; }
    if (BEXP <= BMIN) {
        std::printf("\n***  INITIAL B = %12.6G (OUT OF BOUNDS)--\n", BEXP);
        goto L940;
    }

L110:
    { // scope for variables initialized after goto target
    if (PBUGSW) std::printf("\nB FOR EXP OR POWER FORM = %12.6G\n", BEXP);
    STEP = STEP * DABS(BEXP/10.0);
    if (IEXTYP == 3) STEP = 0.1;
    double BMAX = 10*DABS(BEXP);
    if (IEXTYP == 1) BMAX = 1.2;

    double X = BEXP;
    double C1=0, A1=0;
    LINLSQ(IXTYP2, N, XVALS, SVALS, C1, A1, X, F, PLINSW);
    CHIEXP = F;

    { // Block for L130 loop
    double XPL, XMI, FPL, FMI, DF, DTES;
L130:
    XPL = X*(1+DEL);
    XMI = X*(1-DEL);
    LINLSQ(IXTYP2, N, XVALS, SVALS, CVAL, AVAL, XPL, FPL, PLINSW);
    LINLSQ(IXTYP2, N, XVALS, SVALS, CVAL, AVAL, XMI, FMI, PLINSW);
    DF = (FPL-FMI) / (2*X*DEL);
    DTES = DABS(DF/F);
    if (DTES > 1.0e-10) goto L200;
    if (DEL >= 0.5) goto L700;
    DEL = 2*DEL;
    goto L130;

L200:
    { double S_v = DABS(STEP);
    if (DF > 0) S_v = -S_v;

    int IT = 1;
    double F1 = F, X1 = X;
    double X2, F2, X3, F3, C, A, C2=0, A2=0;

L220:
    X2 = X1 + S_v;
    if (IEXTYP == 1 && X2 > 1.2) goto L950;
    if (IEXTYP >= 2 && X2 < BMIN) goto L960;
    if (IEXTYP == 1 && X2 < BMIN) goto L700;
    if (IEXTYP >= 2 && X2 > BMAX) goto L700;

    LINLSQ(IXTYP2, N, XVALS, SVALS, C, A, X2, F2, PLINSW);
    if (F2 <= F1) goto L250;

    // Passed minimum — bracket it
    S_v = 0.5*S_v;
    X3 = X2; F3 = F2;
    X2 = X3 - S_v;
    LINLSQ(IXTYP2, N, XVALS, SVALS, C, A, X2, F, PLINSW);
    if (F < F1) { F2=F; C2=C; A2=A; goto L300; }
    X3 = X2; F3 = F;
    X2 = X1; F2 = F1; A2 = A1; C2 = C1;
    X1 = X2 - S_v;
    LINLSQ(IXTYP2, N, XVALS, SVALS, C, A, X1, F1, PLINSW);
    goto L300;

L250:
    X1=X2; F1=F2; A1=A; C1=C;
    S_v = S_v*STPMUL;
    IT = IT+1;
    if (IT <= 1000) goto L220;
    std::printf("\n*** MORE THAN 1000 STEPS TO BOX MINIMUM IN B-- \n");
    goto L700;

L300:
    { double F2P_v = F1 - 2*F2 + F3;
    double SM, XM;
    if (F2P_v == 0) goto L310;
    SM = (0.5*S_v*(F1-F3)) / F2P_v;
    XM = X2 + SM;
    if (XM <= BMIN || XM > BMAX) goto L700;
    LINLSQ(IXTYP2, N, XVALS, SVALS, CVAL, AVAL, XM, F, PLINSW);
    if (F <= F2) goto L320;
L310:
    XM = X2; AVAL = A2; CVAL = C2; F = F2;
L320:
    B = XM; CHISQ = F;

    // L400: Taylor refinement
    { double DX = DEL*XM;
    double FM = F;
    XPL = XM + DX;
    LINLSQ(IXTYP2, N, XVALS, SVALS, C, A, XPL, FPL, PLINSW);
    XMI = XM - DX;
    LINLSQ(IXTYP2, N, XVALS, SVALS, C, A, XMI, FMI, PLINSW);
    double FP_v = (FPL-FMI) / (2*DX);
    F2P_v = (FPL - 2*FM + FMI) / (DX*DX);
    if (F2P_v != 0) {
        double XNEW = XM - FP_v/F2P_v;
        LINLSQ(IXTYP2, N, XVALS, SVALS, C, A, XNEW, FM, PLINSW);
        if (FM < F) { AVAL=A; CVAL=C; B=XNEW; CHISQ=FM; }
    }
    } // end Taylor block

    // L430
    { double TES = CVAL + AVAL;
    if (PBUGSW) std::printf("\n-*  (A,C) = (%12.6G %12.6G)\n     (B,B2,CH)=(%12.6G %12.6G %12.6G)\n",
        AVAL, CVAL, B, X2, CHISQ);
    if (TES <= 0 && IEXTYP != 3) {
        if (PBUGSW) std::printf("\n--- BEST-FIT IS AN EXP OR SIMPLE POWER --\n");
        goto L700;
    }
    if (PBUGSW) std::printf("\n-COMPLETE IN %6d ITERATIONS\n", IT);
    }

    // Return values per IEXTYP
    switch (IEXTYP) {
    case 1: // Woods-Saxon
        BARA = 1.0/DABS(AVAL); WID = 1.0/B;
        BARC = DSIGN(1.0, AVAL);
        FLCRIT = XBAR - WID*DLOG(BARA*CVAL);
        BARL = FLCRIT; AVAL = BARA;
        goto L800;
    case 2: case 4: // Power fit
        WID = B;
        if (AVAL < 0) goto L700;
        BARA = AVAL; BARC = 1;
        BARL = CVAL*XBAR/AVAL; FLCRIT = BARL;
        goto L800;
    case 3: // Phase power
        BARA = AVAL; BARC = CVAL; FLCRIT = CVAL; WID = B;
        goto L800;
    }
    } // end L300 block
    } // end L200 block
    } // end L130 block
    } // end L110 scope

L700:
    if (IEXTYP == 3) { AVAL=AEXP; CVAL=CEXP; B=BEXP; goto L570; }
    BARA = DEXP(AEXP); B = BEXP; CHISQ = CHIEXP; BARC = 0; ISIG = 5;
    if (IEXTYP != 1) { WID=BEXP; FLCRIT=0; BARL=0; AVAL=BARA; goto L800; }
    BARL = XBAR; AVAL = BARA; WID = 1.0/BEXP; FLCRIT = BARL;
    goto L800;

L570:
    BARA = AVAL; BARC = CVAL; FLCRIT = CVAL; WID = B;

L800:
    if (B <= BMIN) {
        std::printf("\n*** ERROR IN EXTRAP: B OUT OF BOUNDS:%15.5G%15.5G%15.5G%15.5G\n", B, BARA, BARC, BARL);
        goto L940;
    }
    if (PBUGSW) std::printf("\nLXTRP1 END%3d%13.5G%13.5G%13.5G%13.5G%13.5G%13.5G%13.5G%13.5G%13.5G\n\n",
        IEXTYP, BARA, B, BARC, BARL, AVAL, WID, CVAL, FLCRIT, CHISQ);
    if (MCHN == 1) return;
    switch (IEXTYP) {
    case 1: if (B > 0.01) return; break;
    case 2: case 4: if (B > BMIN+0.6) return; break;
    case 3: return;
    }
    std::printf("\n*** WARNING:  EXTRAPOLATION DECAYS VERY SLOWLY;  A, B, C, LCRIT = %15.5G%15.5G%15.5G%15.5G\n",
        AVAL, B, CVAL, FLCRIT);
    std::printf("      FOR CHANNEL%3d, JP, JT, LX, LO-LI =%3d/2%3d/2%4d%4d\n", MCHN, JP, JT, LX, LDEL);
    return;

L940:
    std::printf("\n*** MINIMIZATION W.R.T. WIDTH HAS FAILED \n       EXTRAPOLATION MUST BE SUPPRESSED --\n");
    ISIG = -5;
    return;

L950:
    std::printf("\n**** ATTEMPTING TO FIT VERY SHARP WOODS-SAXON - ATTEMPT ABANDONED\n");
    std::printf("      FOR CHANNEL%3d, JP, JT, LX, LO-LI =%3d/2%3d/2%4d%4d\n", MCHN, JP, JT, LX, LDEL);
    goto L700;

L960:
    std::printf("\n*** ATTEMPTING TO FIT WITH POWER OF 0 -- DEFAULT STARTING VALUE USED.\n");
    std::printf("      FOR CHANNEL%3d, JP, JT, LX, LO-LI =%3d/2%3d/2%4d%4d\n", MCHN, JP, JT, LX, LDEL);
    goto L700;
}

// SUBROUTINE LXTRP2 — source.f L25556-25594
void LXTRP2(int IEXTYP, double BARA, double B, double BARC, double BARL,
            double DLMAX, int LI, double& SIZE, int LX, int LDEL, double* ETAS)
{
    double DLI = LI;
    double ETA, X;
    switch (IEXTYP) {
    case 1:
        X = B*(DLI-BARL);
        if (X < -46) X = -46;
        SIZE = BARA / (BARC + DEXP(X));
        return;
    case 2:
        SIZE = BARA * std::pow(DLMAX/DLI, B) * (1 + BARL/DLI);
        return;
    case 3:
        SIZE = BARC + BARA * std::pow(DLMAX/DLI, B);
        return;
    case 4:
        ETA = 0.5*(ETAS[1]+ETAS[2]);
        DLI = DLI + 0.5*LDEL;
        SIZE = BARA * std::pow(DLMAX/DLI, B) * (1 + BARL/DLI);
        SIZE = SIZE * std::pow(DSIN(DATAN(ETA/DLI)), (double)LX)
             * DEXP(-(ETAS[2]-ETAS[1]) * (DATAN(ETA/DLI) + DLI/ETA));
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
