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

void GETBFC(int ILI, int J, double DELTA, double RASYMP, int NMBFC,
            float* BFACIN, float* BFACOT, int MBNDX, int* BINDEX,
            int MCHNDF, int* CHNDEF, int MCHNVL, double* CHNVAL,
            int MBASDF, int* BASDEF, int NBASDF, float* FG)
{
    // Computes the B-factors at RASYMP-DELTA from those at RASYMP
    // source.f L17558-17753
    // BFACIN is REAL*4(2, NMBFC, *) — 1-based on ILI dimension
    // BFACOT is REAL*4(2, 2, NMBFC) — 1-based
    // BINDEX is INTEGER(MBNDX, *) — 1-based
    // CHNDEF is INTEGER(MCHNDF, *) — 1-based
    // CHNVAL is REAL*8(MCHNVL, *) — 1-based
    // BASDEF is INTEGER(MBASDF, NBASDF) — 1-based
    // FG is REAL*4(4, NBASDF, *) — 1-based (c4f path)

    // Access macros for column-major Fortran arrays (1-based)
    #define BFACIN_F(ic, ib, ili)  BFACIN[((ili)-1)*NMBFC*2 + ((ib)-1)*2 + ((ic)-1)]
    #define BFACOT_F(ic, ir, ib)   BFACOT[((ib)-1)*2*2 + ((ir)-1)*2 + ((ic)-1)]
    #define BINDEX_F(i, j)         BINDEX[((j)-1)*MBNDX + ((i)-1)]
    #define CHNDEF_F(i, j)         CHNDEF[((j)-1)*MCHNDF + ((i)-1)]
    #define CHNVAL_F(i, j)         CHNVAL[((j)-1)*MCHNVL + ((i)-1)]
    #define BASDEF_F(i, j)         BASDEF[((j)-1)*MBASDF + ((i)-1)]
    #define FG_F(i, ib, ili)       FG[((ili)-1)*NBASDF*4 + ((ib)-1)*4 + ((i)-1)]

    int DBUGSW = (INTGER.IPRINT % 10) >= 5;

    int& NBINDX = GRIDCM.INTOFF;  // GRIDCM overlay for CC

    double FAC = 3.0 * INTGER.IZS[1] * INTGER.IZS[3] * CNSTNT.HBARC * DELTA
               / (4.0 * CNSTNT.AFINE) * 2.0 * FLOAT_common.AM / (CNSTNT.HBARC * CNSTNT.HBARC);

    if (DBUGSW)
        std::printf(" BFACTORS:%5d%5d/2%16.6G%16.6G%16.6G\n", ILI, J, DELTA, RASYMP, FAC);

    // Loop over all sets of BFACTORS
    for (int IBSET = 1; IBSET <= NBINDX; IBSET++) {
        int LX     = BINDEX_F(1, IBSET);
        int IC1    = BINDEX_F(2, IBSET);
        int IC2    = BINDEX_F(3, IBSET);
        int MXLMJ1 = BINDEX_F(4, IBSET);
        int MNLMJ1 = BINDEX_F(6, IBSET);
        int LB     = BINDEX_F(8, IBSET);

        int MINL1 = (MNLMJ1 + J) / 2;

        double AK1 = CHNVAL_F(1, IC1);
        int IB1ST   = CHNDEF_F(8, IC1);
        int IB1END  = IB1ST + CHNDEF_F(9, IC1) - 1;
        double AK2 = CHNVAL_F(1, IC2);
        int IB2ST   = CHNDEF_F(8, IC2);
        int IB2END  = IB2ST + CHNDEF_F(9, IC2) - 1;

        if (DBUGSW)
            std::printf(" CHANNELS%4d%4d     LB =%6d\n", IC1, IC2, LB);

        double FACTOR = FAC / std::sqrt(AK1 * AK2);

        // Loop through both sets of basis states
        for (int IB1 = IB1ST; IB1 <= IB1END; IB1++) {
            int L1 = J + BASDEF_F(3, IB1);
            if (L1 < 0) continue;
            L1 = L1 >> 1;  // ishft(L1, -1)

            for (int IB2 = IB2ST; IB2 <= IB2END; IB2++) {
                int L2 = J + BASDEF_F(3, IB2);
                if (L2 < 0) continue;
                L2 = L2 >> 1;
                if (std::abs(L1 - L2) > LX) continue;
                int I = (LX + 1) * (L1 - MINL1) + LX - L1 + L2;
                I = (I >> 1) + LB;

                // Use one trapezoid to integrate from RASYMP to RASYMP-DELTA
                BFACOT_F(1, 1, I) = BFACIN_F(1, I, ILI) + (float)(FACTOR * (
                    (FG_F(1,IB1,ILI) * FG_F(2,IB2,ILI) -
                     FG_F(2,IB1,ILI) * FG_F(1,IB2,ILI)) / std::pow(RASYMP - DELTA, LX + 1)
                  + (FG_F(3,IB1,ILI) * FG_F(4,IB2,ILI) -
                     FG_F(4,IB1,ILI) * FG_F(3,IB2,ILI)) / std::pow(RASYMP, LX + 1)));

                BFACOT_F(1, 2, I) = BFACIN_F(1, I, ILI);

                BFACOT_F(2, 1, I) = BFACIN_F(2, I, ILI) + (float)(FACTOR * (
                    (FG_F(1,IB2,ILI) * FG_F(1,IB1,ILI) +
                     FG_F(2,IB2,ILI) * FG_F(2,IB1,ILI)) / std::pow(RASYMP - DELTA, LX + 1)
                  + (FG_F(3,IB2,ILI) * FG_F(3,IB1,ILI) +
                     FG_F(4,IB2,ILI) * FG_F(4,IB1,ILI)) / std::pow(RASYMP, LX + 1)));

                BFACOT_F(2, 2, I) = BFACIN_F(2, I, ILI);

                if (DBUGSW)
                    std::printf("%4d%4d%4d%20.10G%18.10G%20.10G%18.10G\n",
                        L1, L2, I,
                        (double)BFACOT_F(1,1,I), (double)BFACOT_F(2,1,I),
                        (double)BFACOT_F(1,2,I), (double)BFACOT_F(2,2,I));
            }
        }
    }

    #undef BFACIN_F
    #undef BFACOT_F
    #undef BINDEX_F
    #undef CHNDEF_F
    #undef CHNVAL_F
    #undef BASDEF_F
    #undef FG_F
}

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
            int NUMLIS, int* LIS, int JMIN, float* FF2INT,
            int* INDXS, double* SMATR, double* SMATI)
{
    #define BINDEX_B(i, j) BINDEX[((j)-1)*MBINDX + ((i)-1)]
    #define CHNDEF_B(i, j) CHNDF_i[((j)-1)*MCHNDF + ((i)-1)]
    #define BASDEF_B(i, j) BASDF_i[((j)-1)*MBASDF + ((i)-1)]
    #define BASCUP_B(i, j) BASCP_i[((j)-1)*MBASCP + ((i)-1)]
    #define RBASCP_B(i, j) BASCP_f[((j)-1)*MBASCP + ((i)-1)]

    const double ROOT2 = 1.4142135623731e0;

    bool PRNTSW = (INTGER.IPRINT % 10) >= 2;
    int L1MJ = 0;

    for (int IBF = 1; IBF <= NBINDX; IBF++) {
        int LB = BINDEX_B(9, IBF);
        if (LB == 0) continue;
        int LX = BINDEX_B(1, IBF);
        int IC1 = BINDEX_B(2, IBF);
        int IC2 = BINDEX_B(3, IBF);
        int MXLMJ1 = BINDEX_B(4, IBF);
        int MNLMJ1 = BINDEX_B(6, IBF);

        int NM1FRJ = (MXLMJ1 - MNLMJ1) / 4 + 1;
        int NMBFRJ = (LX + 1) * NM1FRJ;

        int MCHN = CHNDEF_B(1, IC2);
        int IB2ST = CHNDEF_B(8, IC2);
        int IB2ND = IB2ST + CHNDEF_B(9, IC2) - 1;
        int KSMAT = CHNDEF_B(10, IC2);
        int NSPLI = CHNDEF_B(11, IC2);
        int KINDXS = CHNDEF_B(12, IC2) + 3 + 3 * LX;
        int KBASE = KSMAT + (int)INDXS[KINDXS - 3];
        int LDELMN = (int)INDXS[KINDXS - 2];
        int LIMOST = CHNDEF_B(15, IC2);
        int JMOST = 2 * LIMOST;
        int NUMJ = (JMOST - JMIN) / 2 + 1;

        bool BADSW = true;

        for (int IB2 = IB2ST; IB2 <= IB2ND; IB2++) {
            int L2MJ = BASDEF_B(3, IB2);
            int ICUPST = BASDEF_B(4, IB2);
            int ICUPND = ICUPST + BASDEF_B(5, IB2) - 1;
            int LDEL = (L2MJ - L1MJ) / 2;
            int KOFFS = KBASE + ((LDEL - LDELMN) >> 1);

            double FAC1 = 0;
            for (int ICUP = ICUPST; ICUP <= ICUPND; ICUP++) {
                if (abs(BASCUP_B(9, ICUP)) == IBF)
                    FAC1 = FAC1 + (double)RBASCP_B(6, ICUP);
            }
            if (FAC1 == 0) continue;
            BADSW = false;

            FAC1 = FAC1 / CNSTNT.RT4PI;
            if (BASDEF_B(7, IB2) == 2) FAC1 = ROOT2 * FAC1;
            int LXP1 = LX + 1;
            if (!(LXP1 & 2)) FAC1 = -FAC1;

            int INDX1 = LB + ((LDEL + LX) >> 1);

            for (int NJ = 1; NJ <= NUMJ; NJ++) {
                int J = 2 * NJ + (JMIN - 2);
                int L1X2 = J + L1MJ;
                int L2X2 = J + L2MJ;
                if (L2X2 < 2 * LX - L1X2) continue;
                int I = INDX1 + NMBFRJ * (NJ - 1);
                FF2INT[I - 1] = (float)(FAC1 * (double)FF2INT[I - 1]
                    * fabs(CLEBSH(L1X2, 2 * LX, 0, 0, L2X2, 0)));
            }

            bool IMAGSW = !(LX & 1);

            if (PRNTSW)
                std::printf("\n COMPARASON OF C.C. AND COULOMB BORN"
                    " S MATRICES FOR CHANNEL%3d, BASIS STATE%4d"
                    "     LX =%3d     LOUT-LIN =%3d\n"
                    "\n    J   L IN  L OUT%13sS(C.C.)%17sS(BORN)"
                    "%13sDIFFERENCE%12s|DIFFERENCE|%9s%%\n",
                    MCHN, IB2, LX, LDEL, "", "", "", "", "");

            for (int ILI = 1; ILI <= NUMLIS; ILI++) {
                int J = 2 * LIS[ILI - 1];
                int L1X2 = J + L1MJ;
                int L2X2 = J + L2MJ;
                if (L2X2 < 2 * LX - L1X2) continue;
                int L1 = L1X2 >> 1;
                int L2 = L2X2 >> 1;
                int I1 = (ILI - 1) * NSPLI + KOFFS;
                int I2 = INDX1 + NMBFRJ * ((J - JMIN) >> 1);
                double SR = SMATR[I1 - 1];
                double SI = SMATI[I1 - 1];
                double SBORN = (double)FF2INT[I2 - 1];
                if (!IMAGSW)
                    SMATR[I1 - 1] = SR - SBORN;
                else
                    SMATI[I1 - 1] = SI - SBORN;

                if (PRNTSW) {
                    double SMAG = sqrt(SMATR[I1 - 1] * SMATR[I1 - 1] + SMATI[I1 - 1] * SMATI[I1 - 1]);
                    double PERCNT = 100.0 * SMAG / (fabs(SBORN) + 1.0e-20);
                    std::printf("%5d/2%5d%6d%15.5G%13.5G%16.5G%16.5G%13.5G%14.4G%13.2f\n",
                        J, L1, L2, SR, SI, SBORN, SMATR[I1 - 1], SMATI[I1 - 1], SMAG, PERCNT);
                }
            }
        }
        if (BADSW)
            std::printf("\n *** INTERNAL ERROR IN SETBRN:%8d%8d%8d%8d%8d%8d%8d\n",
                IBF, LB, IC1, IC2, LX, IB2ST, IB2ND);
    }

    #undef BINDEX_B
    #undef CHNDEF_B
    #undef BASDEF_B
    #undef BASCUP_B
    #undef RBASCP_B
}
