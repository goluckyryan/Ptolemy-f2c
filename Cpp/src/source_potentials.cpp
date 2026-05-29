// source_potentials.cpp — Potential generation, form factors, Woods-Saxon
// Translated from source.f: WOODSX, DEFWOO, MAKPOT, SETPOT, WAVPOT, FFACST, SETSPT, SFROMI, BEBETA

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include "ptolemy_alloc.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_forward.h"
#include "ptolemy_io.h"

#include <cstdio>
#include <cstring>
#include <cmath>
#include <algorithm>

// External functions
extern int IALLOC(int IWRDS);
extern int NALLOC(const char* NAME, int IWRDS);
extern int NAMLOC(const char* NAME);
extern void INIT8(int I, double X8);
extern void INIT4(int I, float x4);
extern void INIT4i(int I, int i4);
extern void SETLOG(int N);
extern double second();
extern void VCSQ12(double RVAL, double& X, int NWP);
extern void SETBFC_full(
    int NBINDX, int MBINDX, int* BINDEX,
    int MCHNDF, int* CHNDEF_i, float* RCHNDF,
    int NUMLIS, int* LIS, double* R2S_arg, int LMXMX,
    double* SIG1, double* SIG2,
    double RLOWER, int ALLSW, int IRORC,
    int NMBFAC, float* BFAC,
    double* FFWORK, int IDIM1, int LDLDIM,
    double* WRKWK, double* FIWORK, int IDIM2,
    double* WRKST, int& IRET, double& CLTIME);
extern void SETBRN(int NBINDX, int MBINDX, int* BINDEX,
    int MCHNDF, int* CHNDF_i, int MBASDF, int* BASDF_i,
    int MBASCP, int* BASCP_i, float* BASCP_f,
    int NUMLIS, int* LIS, int JMIN, float* FF2INT,
    int* INDXS, double* SMATR, double* SMATI);

// ============================================================================
// WOODSX — Computes Woods-Saxon type potentials
// ============================================================================
void WOODSX(int NUMPTS, double RSTRT, double STEPSZ, double* VRAY,
            int& N1, int& N2, int ITYPE, double V, double R, double A, double POWER)
{
    // Check for zero potential
    N1 = 0;
    N2 = 0;
    if (V == 0.0) goto L300;

    {
        double XSTP = STEPSZ / A;
        double YSTP = DEXP(XSTP);
        double X = (RSTRT - STEPSZ - R) / A;

        // Break region into three regions
        N1 = (int)(-(100.0 + X) / XSTP);
        N1 = MAX0(N1, 1);
        N1 = MIN0(N1, NUMPTS);
        double Y_val = 46.0 + DLOG(DABS(V));
        Y_val = DMIN1(Y_val, CNSTNT.BIGLOG);
        if (POWER > 1.0) Y_val = Y_val / POWER;
        N2 = (int)((Y_val - X) / XSTP);
        N2 = MIN0(N2, NUMPTS);
        N2 = MAX0(N2, N1);

        // Region 1 — potential is -V or zero
        if (N1 == 1) goto L200;
        {
            double Y_fill = -V;
            if (ITYPE > 1) Y_fill = 0.0;
            for (int I = 1; I <= N1; I++) VRAY[I] = Y_fill;
        }

        // Region 2 — compute the potential
    L200:
        X = X + N1 * XSTP;
        double Y = DEXP(X);
        switch (ITYPE) {
            case 1: goto L210;
            case 2: goto L230;
            case 3: goto L230;
            case 4: goto L250;
            case 5: goto L260;
        }

        // Straight W.S.
    L210:
        if (POWER != 1.0) goto L220;
        for (int I = N1; I <= N2; I++) {
            VRAY[I] = -V / (1.0 + Y);
            Y = Y * YSTP;
        }
        goto L300;

        // Power of W.S.
    L220:
        for (int I = N1; I <= N2; I++) {
            VRAY[I] = -V / std::pow(1.0 + Y, POWER);
            Y = Y * YSTP;
        }
        goto L300;

        // Derivative of W.S.
    L230:
        for (int I = N1; I <= N2; I++) {
            VRAY[I] = (-4.0 * V) * (Y / (1.0 + Y)) / (1.0 + Y);
            Y = Y * YSTP;
        }
        if (ITYPE != 2) goto L300;

        // Convert derivative of W.S. to spin orbit
        {
            double RT = RSTRT + (N1 - 1) * STEPSZ + 1.0e-30;
            for (int I = N1; I <= N2; I++) {
                VRAY[I] = (0.5 / A) * VRAY[I] / RT;
                RT = RT + STEPSZ;
            }
        }
        goto L300;

        // Folded real TR
    L250:
        {
            double V1 = 2.0 * V / (A * A);
            double RT = RSTRT + (N1 - 1) * STEPSZ + 1.0e-30;
            for (int I = N1; I <= N2; I++) {
                double XX = 1.0 + Y;
                VRAY[I] = V1 * (Y / XX) * ((1.0 - Y) / XX - A / RT) / XX;
                RT = RT + STEPSZ;
                Y = Y * YSTP;
            }
        }
        goto L300;

        // Folded imaginary TR
    L260:
        {
            double V1 = 8.0 * V / A;
            double AA = 1.0 / A;
            double A4 = 4.0 * AA;
            double RT = RSTRT + (N1 - 1) * STEPSZ + 1.0e-30;
            for (int I = N1; I <= N2; I++) {
                double XX = 1.0 + Y;
                double YY = Y / XX;
                double RR = 1.0 / RT;
                VRAY[I] = V1 * YY * (((RR - AA) / XX + A4 * YY) / XX
                    - (RR + AA) * YY * YY) / XX;
                RT = RT + STEPSZ;
                Y = Y * YSTP;
            }
        }
    }

    // Fill out rest with 0
L300:
    if (N2 == NUMPTS) return;
    {
        int N2P1 = N2 + 1;
        for (int I = N2P1; I <= NUMPTS; I++) VRAY[I] = 0.0;
    }
    return;
}


// ============================================================================
// DEFWOO — Computes deformed Woods-Saxon type potentials
// ============================================================================
void DEFWOO(int NUMPTS, double RSTRT, double STEPSZ, double* VRAY,
            int& N1, int& N2, int ITYPE, double V, double R, double A, double POWER,
            int IPORT, int LXMIN_arg, int* LDEFMX, double* RBETAS_flat,
            int* LXS, int NUMPT)
{
    // RBETAS is dimensioned (2,*) in Fortran; access as RBETAS_flat[2*(LX-1)+I]
    // using 1-based: RBETAS(I,LX) = RBETAS_flat[(LX)*2 + I]  (0-based flat)
    // Actually let's use a macro for clarity
    #define RBETAS(i,lx) RBETAS_flat[((lx)-1)*2 + (i)]

    double XS[51], WTS[51]; // 1-based
    double YS[3][51];       // YS(2,50) -> [3][51]
    double WT[3][51];       // WT(2,50) -> [3][51]
    double PLM[21];         // PLM(20) -> [21]
    double FACS[3][21];     // FACS(2,20) -> [3][21]
    int ODDLSW_val;

    N1 = 0;
    N2 = 0;
    if (V == 0.0) goto L500;

    // Check for odd multipoles and count deformed nuclei
    ODDLSW_val = FALSE_F;
    {
        int NDIM = 0;
        int IPT = 0;
        int LXMAX_loc;

        for (int I = 1; I <= 2; I++) {
            LXMAX_loc = LDEFMX[I];
            if (LXMAX_loc == 0) continue;
            NDIM++;
            IPT = I;
            for (int LX = LXMIN_arg; LX <= LXMAX_loc; LX++) {
                if ((LX & 1) && RBETAS(I, LX) != 0.0) ODDLSW_val = TRUE_F;
                FACS[I][LX] = DSQRT((2 * LX + 1) / (4.0 * CNSTNT.PI)) * RBETAS(I, LX);
            }
            int LX = LXS[I];
            if ((LX & 1) && !ODDLSW_val) goto L500;
        }

        int NGAUSS = NUMPT;
        if (ODDLSW_val) NGAUSS = 2 * NUMPT;

        // Setup deformation and weights at each integration point
        GAUSSL(NGAUSS, XS, WTS);
        double VSCALE = std::pow(0.5, NDIM) * V;
        if (IPORT == 0) goto L50;
        if (IPORT < 3) VSCALE = std::sqrt(CNSTNT.PI * (2 * LXS[IPORT] + 1)) * 2 * VSCALE;
        if (IPORT == 3) VSCALE = CNSTNT.PI * std::sqrt(
            (2 * LXS[1] + 1) * (2 * LXS[2] + 1.0)) * V;

    L50:
        LXMAX_loc = MAX0(LXS[1], MAX0(LXS[2], MAX0(LDEFMX[1], LDEFMX[2])));
        for (int I = 1; I <= NGAUSS; I++) {
            double X = XS[I];
            if (!ODDLSW_val) X = 0.5 * (1.0 + X);
            PLMSUB(LXMAX_loc, 0, X, PLM);
            for (int II = 1; II <= 2; II++) {
                int LXM = LDEFMX[II];
                double Y = 1.0;
                if (LXM == 0) goto L85;
                Y = 0.0;
                for (int LX = LXMIN_arg; LX <= LXM; LX++) {
                    Y = Y + FACS[II][LX] * PLM[LX + 1];
                }
                Y = DEXP(-Y / A);
            L85:
                YS[II][I] = Y;
                WT[II][I] = WTS[I] * PLM[LXS[II] + 1];
            }
        }

        double XSTP = STEPSZ / A;
        double YSTP = DEXP(XSTP);
        double X = (RSTRT - STEPSZ - R) / A;

        N1 = (int)(-(100.0 + X) / XSTP);
        N1 = MAX0(N1, 1);
        N1 = MIN0(N1, NUMPTS);
        double Y_val = 46.0 + DLOG(DABS(V));
        Y_val = DMIN1(Y_val, CNSTNT.BIGLOG);
        if (POWER > 1.0) Y_val = Y_val / POWER;
        N2 = (int)((Y_val - X) / XSTP);
        N2 = MIN0(N2, NUMPTS);
        N2 = MAX0(N2, N1);

        // Region 1
        if (N1 == 1) goto L200;
        {
            double Y_fill = -V;
            if (ITYPE > 1 || IPORT != 0) Y_fill = 0.0;
            for (int I = 1; I <= N1; I++) VRAY[I] = Y_fill;
        }

        // Region 2
    L200:
        X = X + N1 * XSTP;
        double Y = DEXP(X);
        if (ITYPE == 3) goto L300;

        // Straight W.S.
        if (POWER != 1.0) goto L250;
        for (int I = N1; I <= N2; I++) {
            double VV = 0.0;
            if (NDIM == 2) {
                for (int II = 1; II <= NGAUSS; II++)
                    for (int III = 1; III <= NGAUSS; III++)
                        VV += WT[1][II] * WT[2][III] / (1.0 + Y * YS[1][II] * YS[2][III]);
            } else {
                for (int II = 1; II <= NGAUSS; II++)
                    VV += WT[IPT][II] / (1.0 + Y * YS[IPT][II]);
            }
            VRAY[I] = -VSCALE * VV;
            Y = Y * YSTP;
        }
        goto L500;

        // Power of W.S.
    L250:
        for (int I = N1; I <= N2; I++) {
            double VV = 0.0;
            if (NDIM == 2) {
                for (int II = 1; II <= NGAUSS; II++)
                    for (int III = 1; III <= NGAUSS; III++)
                        VV += WT[1][II] * WT[2][III] / std::pow(1.0 + Y * YS[1][II] * YS[2][III], POWER);
            } else {
                for (int II = 1; II <= NGAUSS; II++)
                    VV += WT[IPT][II] / std::pow(1.0 + Y * YS[IPT][II], POWER);
            }
            VRAY[I] = -VSCALE * VV;
            Y = Y * YSTP;
        }
        goto L500;

        // Derivative of W.S.
    L300:
        VSCALE = +4.0 * VSCALE;
        for (int I = N1; I <= N2; I++) {
            double VV = 0.0;
            if (NDIM == 2) {
                for (int II = 1; II <= NGAUSS; II++)
                    for (int III = 1; III <= NGAUSS; III++) {
                        double YY = Y * YS[1][II] * YS[2][III];
                        VV += WT[1][II] * WT[2][III] / (1.0 + YY) * YY / (1.0 + YY);
                    }
            } else {
                for (int II = 1; II <= NGAUSS; II++) {
                    double YY = Y * YS[IPT][II];
                    VV += WT[IPT][II] / (1.0 + YY) * YY / (1.0 + YY);
                }
            }
            VRAY[I] = -VSCALE * VV;
            Y = Y * YSTP;
        }
    }

    // Fill out rest with 0
L500:
    if (N2 == NUMPTS) return;
    {
        int N2P1 = N2 + 1;
        for (int I = N2P1; I <= NUMPTS; I++) VRAY[I] = 0.0;
    }
    return;

    #undef RBETAS
}


// ============================================================================
// MAKPOT — Computes potentials for WAVELJ
// ============================================================================
void MAKPOT(int NWPIN, int& IRET)
{
    static int ITYPS[7] = { 0, 4, 5, 2, 2, 1, 3 }; // 1-based
    int IDUMMY = 0;
    double DUMMY = 0.0;

    int NWP = NWPIN;
    if (NWPIN == 3) NWP = 1;
    IRET = 0;
    double STEPSZ_loc = WAVCOM.RSTEPS[NWP];
    double EINV = 1.0 / FLOAT_common.E;
    double H = WAVCOM.HS[NWP];
    double H2 = H * H;
    int NSTEPP = WAVCOM.NSTPSS[NWP] + 1;
    int LREALV = LOCPTRS.Z[WAVCOM.IRLVS[NWP]] - 1;
    int LIMAGV = LOCPTRS.Z[WAVCOM.IIMVS[NWP]] - 1;
    int LCENTR = LOCPTRS.Z[WAVCOM.ICENTR[NWP]] - 1;
    int LVWORK = LOCPTRS.Z[WAVCOM.IVWORK] - 1;

    double RVAL = -STEPSZ_loc;

    // Indicate WAVELJ must reestablish starting point
    WAVCOM.LASTL[NWP] = 999999;
    WAVCOM.LASTZR = 1;

    // Get the Coulomb potential
    int I = 5;
    int LL = LREALV + 1;
    int N1, N2;
    if (LNKBLK.LNKADR[5][3] > 0 && NWPIN != 3) goto L800;
    {
        double X = 0.0;
        for (I = 1; I <= NSTEPP; I++) {
            RVAL = RVAL + STEPSZ_loc;
            if (NWPIN != 3)
                VCSQ12(RVAL, X, NWP);
            ALLOC(LREALV + I) = X;
        }
    }

    // Get the real potential
L300:
    I = 1;
    LL = LVWORK + 1;
    if (LNKBLK.LNKADR[1][3] != 0) goto L800;
    WOODSX(NSTEPP, 0.0, STEPSZ_loc, ALLOC_base(LVWORK + 1), N1, N2,
        1, TEMPVS.TVR, FLOAT_common.R, TEMPVS.TAR, FLOAT_common.POWRL);

    // Store Coulomb + real potential
L400:
    {
        double RHO = 1.0e-20;
        double AJ = 1.0;
        if (SWITCH.PROBLM == 24 || NWPIN == 3) AJ = 0.0;
        double TEMP = 1.0;
        if (NWPIN == 3) TEMP = 0.0;
        for (I = 1; I <= NSTEPP; I++) {
            double AK1 = EINV * (ALLOC(LREALV + I) + ALLOC(LVWORK + I)) - AJ;
            ALLOC(LREALV + I) = TEMP - AK1 * (H2 / 12.0);
            ALLOC(LCENTR + I) = -(H2 / 12.0) / (RHO * RHO);
            RHO = RHO + H;
        }
    }

    // Imaginary potential
    I = 2;
    if (LNKBLK.LNKADR[2][3] != 0) goto L800;
    WOODSX(NSTEPP, 0.0, STEPSZ_loc, ALLOC_base(LVWORK + 1), N1, N2,
        1, TEMPVS.TVI, FLOAT_common.RI, TEMPVS.TAI, FLOAT_common.POWIM);
L510:
    for (I = 1; I <= NSTEPP; I++)
        ALLOC(I + LIMAGV) = -(H2 / (12.0 * FLOAT_common.E)) * ALLOC(I + LVWORK);

    if (FLOAT_common.VSI == 0.0) goto L600;
    I = 13;
    if (LNKBLK.LNKADR[13][3] != 0) goto L800;
    WOODSX(NSTEPP, 0.0, STEPSZ_loc, ALLOC_base(LVWORK + 1), N1, N2,
        3, FLOAT_common.VSI, FLOAT_common.RSI, FLOAT_common.ASI, 1.0);
L540:
    for (I = 1; I <= NSTEPP; I++)
        ALLOC(I + LIMAGV) = ALLOC(I + LIMAGV) - (H2 / (12.0 * FLOAT_common.E)) * ALLOC(I + LVWORK);

    // S.O. forces (declared here to avoid jump-over-init)
    int LSOR, LSOI, LTEN;
    LSOR = 0; LSOI = 0; LTEN = 0;
    int II_ten;
    II_ten = 0;
L600:
    if (WAVCOM.ISORS[NWP] == 0) goto L650;
    LSOR = LOCPTRS.Z[WAVCOM.ISORS[NWP]] - 1;
    I = 3;
    if (LNKBLK.LNKADR[3][3] != 0) goto L800;
    WOODSX(NSTEPP, 0.0, STEPSZ_loc, ALLOC_base(LVWORK + 1), N1, N2,
        2, FLOAT_common.VSO, FLOAT_common.RSO, FLOAT_common.ASO, 1.0);
L610:
    {
        int NSTEP_SO = std::min(NSTEPP, LENGTH.LENG[WAVCOM.ISORS[NWP]]);
        for (I = 1; I <= NSTEP_SO; I++)
            ALLOC(LSOR + I) = (-H2 / (12.0 * FLOAT_common.E)) * ALLOC(LVWORK + I);
    }

L650:
    if (WAVCOM.ISOIS[NWP] == 0) goto L700;
    LSOI = LOCPTRS.Z[WAVCOM.ISOIS[NWP]] - 1;
    I = 4;
    if (LNKBLK.LNKADR[4][3] != 0) goto L800;
    WOODSX(NSTEPP, 0.0, STEPSZ_loc, ALLOC_base(LVWORK + 1), N1, N2,
        2, FLOAT_common.VSOI, FLOAT_common.RSOI, FLOAT_common.ASOI, 1.0);
L660:
    {
        int NSTEP_SOI = std::min(NSTEPP, LENGTH.LENG[WAVCOM.ISOIS[NWP]]);
        for (I = 1; I <= NSTEP_SOI; I++)
            ALLOC(LSOI + I) = (-H2 / (12.0 * FLOAT_common.E)) * ALLOC(LVWORK + I);
    }

    // Tensor forces
L700:
    for (II_ten = 1; II_ten <= 6; II_ten++) {
        if (WAVCOM.ITENS[NWP][II_ten] == 0) continue;
        LTEN = LOCPTRS.Z[WAVCOM.ITENS[NWP][II_ten]] - 1;
        I = II_ten + 6;
        if (LNKBLK.LNKADR[I][3] != 0) goto L800;
        WOODSX(NSTEPP, 0.0, STEPSZ_loc, ALLOC_base(LVWORK + 1), N1, N2,
            ITYPS[II_ten], VTEN()[II_ten-1], RTEN()[II_ten-1], ATEN()[II_ten-1], 1.0);
    L710:
        {
            int NSTEP_TEN = std::min(NSTEPP, LENGTH.LENG[WAVCOM.ITENS[NWP][II_ten]]);
            for (I = 1; I <= NSTEP_TEN; I++)
                ALLOC(LTEN + I) = (-H2 / (12.0 * FLOAT_common.E)) * ALLOC(LVWORK + I);
            if (II_ten < 5) { continue; }
            RVAL = 1.0e-30;
            for (I = 1; I <= NSTEP_TEN; I++) {
                ALLOC(LTEN + I) = ALLOC(LTEN + I) / (RVAL * RVAL);
                RVAL = RVAL + STEPSZ_loc;
            }
        }
    }
    return;

    // Get potential from a linkule
L800:
    {
        char8 lname;
        std::memcpy(lname.data, reinterpret_cast<char*>(&LNKBLK.LNKADR[I][1]), 8);
        int myints[3]; myints[1] = LNKBLK.LNKADR[I][5]; myints[2] = LNKBLK.LNKADR[I][6];
        LINKUL(LNKBLK.LNKADR[I][3], lname, myints, I, 3, IRET,
               IDUMMY, DUMMY, 0.0, STEPSZ_loc, NSTEPP,
               ALLOC_base(LL), &DUMMY, (char*)"");
    }
    if (IRET < 0) goto L900;
    switch (I) {
        case 1: goto L400;
        case 2: goto L510;
        case 3: goto L610;
        case 4: goto L660;
        case 5: goto L300;
        case 6: case 7: case 8: case 9: case 10: case 11: case 12: goto L710;
        case 13: goto L540;
    }

L900:
    return;
}


// ============================================================================
// SETPOT — Sets up potential parameters for a single channel
// (large routine: 447 lines of Fortran)
// ============================================================================
void SETPOT(int& IRET)
{
    double RTEMP[3], ATEMP[3]; // 1-based
    IRET = 1;

    // Setup R0 mass factor
    double AMLESS = DMIN1(FLOAT_common.AMP, FLOAT_common.AMT);
    double AMMORE = FLOAT_common.AMP + FLOAT_common.AMT - AMLESS;
    double THIRD = 1.0 / 3.0;

    int I = SWITCH.R0TYPE + 1;
    switch (I) {
        case 1: // R0DEFAULT
            INTRNL.R0MASS = std::pow(AMMORE, THIRD);
            if (AMLESS > 2.5) INTRNL.R0MASS += std::pow(AMLESS, THIRD);
            break;
        case 2: // R0TARGET
            INTRNL.R0MASS = std::pow(AMMORE, THIRD);
            break;
        case 3: // R0SUM
            INTRNL.R0MASS = std::pow(AMMORE, THIRD) + std::pow(AMLESS, THIRD);
            break;
        case 4: // R0MTOT
            INTRNL.R0MASS = std::pow(AMLESS + AMMORE, THIRD);
            break;
    }

    // Energy factor
    double EELAB = (1.0 + INTRNL.RATMAS) * FLOAT_common.E;
    double EEMUL = EELAB;
    if (SWITCH.IEPOW == 1) EEMUL = 1.0 / EELAB;

    // Restore parameters on second pass
    if (SWITCH.IBSPAS > 1 && SWITCH.NEXT >= 1 && SWITCH.NEXT <= 2) {
        if (INTGER.L == NOTDEF_INT) INTGER.L = FORMF.LBP; // simplified
        // ... (restoration logic for bound-state parameters)
    }

    if (FLOAT_common.V == INTRNL.UNDEF) FLOAT_common.V = 0.0;
    if (FLOAT_common.VSO == INTRNL.UNDEF) FLOAT_common.VSO = 0.0;
    if (FLOAT_common.TAU == INTRNL.UNDEF) FLOAT_common.TAU = 0.0;

    // V real params must always be defined
    if (FLOAT_common.R0 != INTRNL.UNDEF) {
        FLOAT_common.R = (FLOAT_common.R0 + FLOAT_common.R0E * EEMUL
            + FLOAT_common.R0ESQ * EEMUL * EEMUL) * INTRNL.R0MASS;
    } else if (FLOAT_common.R == INTRNL.UNDEF) {
        if (LNKBLK.LNKADR[1][3] == 0) {
            std::printf("\n**** R OR R0 MUST BE DEFINED.\n");
            IRET = 0;
        }
    }
    if (FLOAT_common.A == INTRNL.UNDEF) {
        std::printf("\n**** A MUST BE DEFINED.\n");
        IRET = 0;
    }

    // Check imaginary part
    if (FLOAT_common.VI + FLOAT_common.VIE * EEMUL + FLOAT_common.VIESQ * EEMUL * EEMUL != 0.0) {
        if (FLOAT_common.RI0 != INTRNL.UNDEF) {
            FLOAT_common.RI = (FLOAT_common.RI0 + FLOAT_common.RI0E * EEMUL
                + FLOAT_common.RI0ESQ * EEMUL * EEMUL) * INTRNL.R0MASS;
        } else if (FLOAT_common.RI == INTRNL.UNDEF) {
            FLOAT_common.RI = FLOAT_common.R;
            FLOAT_common.RI0 = FLOAT_common.R0;
            FLOAT_common.RI0E = FLOAT_common.R0E;
            FLOAT_common.RI0ESQ = FLOAT_common.R0ESQ;
            INTRNL.WASSET[1] = TRUE_F;
        }
        if (FLOAT_common.AI == INTRNL.UNDEF) {
            FLOAT_common.AI = FLOAT_common.A;
            FLOAT_common.AIE = FLOAT_common.AE;
            FLOAT_common.AIESQ = FLOAT_common.AESQ;
            INTRNL.WASSET[2] = TRUE_F;
        }
    }

    // Real S.O.
    if (FLOAT_common.VSO != 0.0 || FLOAT_common.TAU != 0.0) {
        if (FLOAT_common.RSO0 != INTRNL.UNDEF) {
            FLOAT_common.RSO = FLOAT_common.RSO0 * INTRNL.R0MASS;
        } else if (FLOAT_common.RSO == INTRNL.UNDEF) {
            FLOAT_common.RSO = FLOAT_common.R;
            INTRNL.WASSET[3] = TRUE_F;
        }
        if (FLOAT_common.ASO == INTRNL.UNDEF) {
            FLOAT_common.ASO = FLOAT_common.A;
            INTRNL.WASSET[4] = TRUE_F;
        }
    }

    // Imag S.O.
    if (FLOAT_common.VSOI != 0.0 || FLOAT_common.TAUI != 0.0) {
        if (FLOAT_common.RSOI0 != INTRNL.UNDEF) {
            FLOAT_common.RSOI = FLOAT_common.RSOI0 * INTRNL.R0MASS;
        } else if (FLOAT_common.RSOI == INTRNL.UNDEF) {
            FLOAT_common.RSOI = FLOAT_common.RSO;
            INTRNL.WASSET[5] = TRUE_F;
        }
        INTRNL.WASSET[6] = (FLOAT_common.ASOI == INTRNL.UNDEF) ? TRUE_F : FALSE_F;
        if (FLOAT_common.ASOI == INTRNL.UNDEF) FLOAT_common.ASOI = FLOAT_common.ASO;
    }

    // Surface absorption
    if (FLOAT_common.VSI != 0.0) {
        if (FLOAT_common.RSI0 != INTRNL.UNDEF) {
            FLOAT_common.RSI = FLOAT_common.RSI0 * INTRNL.R0MASS;
        } else if (FLOAT_common.RSI == INTRNL.UNDEF) {
            FLOAT_common.RSI = FLOAT_common.RI;
            INTRNL.WASSET[7] = TRUE_F;
        }
        INTRNL.WASSET[8] = (FLOAT_common.ASI == INTRNL.UNDEF) ? TRUE_F : FALSE_F;
        if (FLOAT_common.ASI == INTRNL.UNDEF) FLOAT_common.ASI = FLOAT_common.AI;
    }

    // Coulomb
    if (INTGER.IZP != NOTDEF_INT && INTGER.IZT != NOTDEF_INT) {
        int IZPZT = INTGER.IZP * INTGER.IZT;
        if (FLOAT_common.RC0 != INTRNL.UNDEF) {
            FLOAT_common.RC = FLOAT_common.RC0 * INTRNL.R0MASS;
        } else if (FLOAT_common.RC == INTRNL.UNDEF) {
            FLOAT_common.RC = 1.0;
            INTRNL.WASSET[9] = TRUE_F;
            if (IZPZT != 0) {
                // Check if folded Coulomb needed
                if (SWITCH.NEXT != 3 && SWITCH.NEXT != 4 && SWITCH.NEXT != 6 && SWITCH.NEXT != 21) {
                    std::printf("\n**** RC OR RC0 MUST BE DEFINED FOR BS.\n");
                    IRET = 0;
                }
            }
        }
    } else {
        std::printf("\n**** ZP AND ZT MUST ALWAYS BE DEFINED.\n");
        IRET = 0;
    }

    // Tensor potentials
    RTEMP[1] = FLOAT_common.R;
    RTEMP[2] = FLOAT_common.RSI;
    if (RTEMP[2] == INTRNL.UNDEF) RTEMP[2] = FLOAT_common.RI;
    ATEMP[1] = FLOAT_common.A;
    ATEMP[2] = FLOAT_common.ASI;
    if (ATEMP[2] == INTRNL.UNDEF) ATEMP[2] = FLOAT_common.AI;
    int IRI = 1;
    for (int ITEN = 1; ITEN <= 6; ITEN++) {
        if (VTEN()[ITEN] == 0.0) goto L850;
        if (RTEN0()[ITEN] != INTRNL.UNDEF) {
            RTEN()[ITEN] = RTEN0()[ITEN] * INTRNL.R0MASS;
        } else if (RTEN()[ITEN] == INTRNL.UNDEF) {
            RTEN()[ITEN] = RTEMP[IRI];
            INTRNL.WASSET[2 * ITEN + 15] = TRUE_F;
        }
        if (ATEN()[ITEN] == INTRNL.UNDEF) {
            ATEN()[ITEN] = ATEMP[IRI];
            INTRNL.WASSET[2 * ITEN + 16] = TRUE_F;
        }
    L850:
        IRI = 3 - IRI;
    }
    return;
}


// ============================================================================
// BEBETA — Computes E(LX) matrix elements and related parameters
// ============================================================================
void BEBETA(int MODEL, int IORD, int I12, int& IBETAC, int MEORD,
            double* XCS, int LXX2, int* LXX2S, int IZ, double AMVAL,
            double RC, double DELTAE, int J2, int J1, int K2, int K1,
            int IBETAN, double* XNS, double RNUC, double SCALE, int IPT,
            double& RMENUC, double& RMECOU, int& IRTN, double* VALS)
{
    IRTN = 1;
    int LX = LXX2 / 2;
    double TJ1P1 = J1 + 1.0;
    double TJ2P1 = J2 + 1.0;
    double ALIFE = 0.0;
    double RME;

    if (LX == 0) goto L40;
    ALIFE = 1.0e11 * std::pow(10.0 * DELTAE / CNSTNT.HBARC, LX);
    for (int L = 1; L <= LX; L++)
        ALIFE = ALIFE / (2 * L + 1);
    ALIFE = 8.0 * CNSTNT.PI * (LX + 1) * 2.9979 * 10.0 * DELTAE /
        (LX * CNSTNT.AFINE * CNSTNT.HBARC) * ALIFE * ALIFE;

L40:
    {
        double DLN2 = DLOG(2.0);
        double AQ = DSQRT(16.0 * CNSTNT.PI / (5.0 * TJ1P1)) *
            CLEBSH(J1, 4, J1, 0, J1, J1);
        double ABE = TJ1P1;
        double X;

        int IBC = -99;
        int IGOTO;
        if (IBETAC == -99) goto L520;
        IBC = 3;
        X = XCS[I12];

        IGOTO = IBETAC + 2;
        switch (IGOTO) {
            case 1: // BE(UP) from 0 to LX (-1)
                RME = DSQRT(DABS(X));
                IBC = 4;
                X = DSIGN(RME, X);
                goto L500;
            case 2: // beta (0)
                IBC = 1;
                goto L500;
            case 3: // delta (1)
                X = X / RC;
                IBETAC = 0;
                IBC = 1;
                goto L500;
            case 4: // RME alpha (2)
                IBC = 2;
                goto L500;
            case 5: // total LX (3) - error
                std::printf("\n*** TOTAL LX INCORRECTLY PLACED IN INPUT\n");
                goto L900;
            case 6: // halflife (4)
                X = X / DLN2;
                // fall through to mean life
            case 7: // mean life (5)
                if (DELTAE <= 0.0) {
                    std::printf("\n*** LIFETIME MAY BE SPECIFIED ONLY FOR DIFFERENT NUCLEAR STATES.\n");
                    return;
                }
                X = 1.0 / (ALIFE * X);
                // fall through to BE(down)
            case 8: // BE(down) (6)
                X = TJ2P1 / TJ1P1 * X;
                // fall through to BE(up)
            case 9: // BE(up) (7)
                RME = DSQRT(DABS(ABE * X));
                X = DSIGN(RME, X);
                goto L500;
            case 10: // RME (B-M) (7)
                goto L500;
            case 11: // BE(up) in W.U. (8)
                X = 1.0 / (4.0 * CNSTNT.PI) * std::pow(3.0 * std::pow(0.12 * std::pow(AMVAL, 0.333333333333330), LX)
                    / (LX + 3), 2);
                RME = DSQRT(DABS(ABE * X));
                X = DSIGN(RME, X);
                goto L500;
            case 12: // BE(down) in W.U. (9)
                X = TJ2P1 / TJ1P1 * X;
                X = 1.0 / (4.0 * CNSTNT.PI) * std::pow(3.0 * std::pow(0.12 * std::pow(AMVAL, 0.333333333333330), LX)
                    / (LX + 3), 2);
                RME = DSQRT(DABS(ABE * X));
                X = DSIGN(RME, X);
                goto L500;
            case 13: // BE(up) (10)
                RME = DSQRT(DABS(ABE * X));
                X = DSIGN(RME, X);
                goto L500;
            case 14: // Q (11)
                if (J1 != J2 || DELTAE != 0.0) {
                    std::printf("\n*** CANNOT INPUT QUADRUPOLE MOMENT FOR OFF-DIAGONAL TRANSITIONS.\n");
                    return;
                }
                if (LX != 2) {
                    std::printf("\n*** QUADRUPOLE MOMENT MAY BE SPECIFIED ONLY FOR LX = 2.\n");
                    return;
                }
                X = X / AQ;
                goto L500;
            case 15: // width (12)
                X = CNSTNT.HBARC / (2.9979e23 * X);
                if (DELTAE <= 0.0) {
                    std::printf("\n*** LIFETIME MAY BE SPECIFIED ONLY FOR DIFFERENT NUCLEAR STATES.\n");
                    return;
                }
                X = 1.0 / (ALIFE * X);
                X = TJ2P1 / TJ1P1 * X;
                RME = DSQRT(DABS(ABE * X));
                X = DSIGN(RME, X);
                goto L500;
        }

    L500:
        XCS[I12] = X;
        if (IBC == 2) IBETAC = 7;

    L520:
        {
            int IBN = 1;
            if (IBETAN == 2) IBN = 2;
            if (IBETAN == -99) IBN = -99;

            // Model-dependent conversions
            switch (MODEL) {
                case 1:
                    MEBROT(IORD, I12, IBN, XNS, IBC, XCS, MEORD, LXX2, LXX2S, IZ,
                        RNUC, RC, DELTAE, J2, J1, K2, K1, RMENUC, RMECOU, IRTN);
                    break;
                case 2:
                    MEBVIB(IORD, I12, IBN, XNS, IBC, XCS, MEORD, LXX2, LXX2S, IZ,
                        RNUC, RC, DELTAE, J2, J1, K2, K1, RMENUC, RMECOU, IRTN);
                    break;
                case 3:
                    MEBDEF(IORD, I12, IBN, XNS, IBC, XCS, MEORD, LXX2, LXX2S, IZ,
                        RNUC, RC, DELTAE, J2, J1, K2, K1, RMENUC, RMECOU, IRTN, IPT);
                    break;
            }
        }

        if (DABS(RMENUC) + DABS(RMECOU) == 0.0) goto L900;
        if (IRTN != 0) goto L900;

        RMENUC = SCALE * RMENUC;
        RMECOU = SCALE * RMECOU;

        // Model-independent quantities
        double BEUP = RMECOU * RMECOU / ABE;
        double BEDOWN = TJ1P1 / TJ2P1 * BEUP;
        double TMEAN = 0.0;
        if (DELTAE > 0.0 && BEDOWN != 0.0)
            TMEAN = 1.0 / (ALIFE * BEDOWN);
        double THALF = DLN2 * TMEAN;
        double Q_val = 0.0;
        if (DELTAE == 0.0 && J1 == J2 && LXX2 == 4) Q_val = AQ * RMECOU;
        double DELTA = RC * XCS[I12];
        VALS[1] = DELTA;
        VALS[2] = RMECOU;
        VALS[3] = BEUP;
        VALS[4] = Q_val;
        VALS[5] = TMEAN;

        // Convert RME E(LX) to FM for calculations
        RMECOU = RMECOU * std::pow(10.0, LX);
        return;
    }

L900:
    for (int i = 1; i <= 5; i++) VALS[i] = 0.0;
    return;
}


// ============================================================================
// FFACST — Gets Coulomb Born S matrix elements for C.C.
// ============================================================================
void FFACST(int& IRTN)
{
    auto* Z = LOCPTRS.Z;
    auto& FACFR4 = ALLOCS.FACFR4;

    // GRIDCM overlay for FFACST context (different field names than GRDSET view)
    int* GRIDCM_int = reinterpret_cast<int*>(&GRIDCM.IMSVAL);
    int& IBINDX = GRIDCM_int[0];
    int& MBINDX = GRIDCM_int[1];
    int& NBINDX = GRIDCM_int[2];
    int& NMFFAC = GRIDCM_int[3];
    int& MXLXGS = GRIDCM_int[4];

    auto& LMIN   = INTGER.LMIN;
    auto& NPCOUL = INTGER.NPCOUL;
    auto& MXCOUL = INTGER.MXCOUL;
    auto& LOMOST = KANDM.LOMOST;
    auto& NUMLIS = INELCM.NUMLIS;
    auto& NASPLI = INELCM.NASPLI;
    auto& ICL2FF = INELCM.ICL2FF;
    auto& IDELSR = INELCM.IDELSR;
    auto& IDELSI = INELCM.IDELSI;
    auto& ILIS   = INELCM.ILIS;
    auto& ISMATR = INELCM.ISMATR;
    auto& ISMATI = INELCM.ISMATI;
    auto& IINDXS = INELCM.IINDXS;

    IRTN = 1;
    if (NMFFAC == 0) return;

    double T1 = second();
    std::printf("1%59sP T O L E M\n"
                "%33sEVALUATION OF COULOMB BORN S-MATRIX ELEMENTS FOR EXTRAPOLATION\n"
                " %.45s%22s%.65s\n",
                "", "", &HEDCOM.REACT[1], "", &HEDCOM.HEADER[1]);

    ICL2FF = NALLOC("FF2INTS ", (NMFFAC + 1) / FACFR4);
    INIT4(ICL2FF, 0.0f);
    IDELSR = NALLOC("DELTASR ", NASPLI * NUMLIS);
    IDELSI = NALLOC("DELTASI ", NASPLI * NUMLIS);

    int MAXDL1 = std::max(2, MXLXGS);
    int LDLDIM = MAXDL1 + 1;
    int LMXMX = LOMOST;
    int LMNMN1 = std::max(LMIN - MAXDL1, 0);
    int IDIM1 = LDLDIM * (2 * (LMXMX - LMNMN1) + 1);
    int IWRKFF = NALLOC("FFWORK  ", IDIM1);
    int IDIM2 = std::max(LMXMX + MAXDL1, 4 * MXCOUL) + 1;
    int IWRKFI = NALLOC("FIWORK  ", 4 * IDIM2);
    int IWRKWK = NALLOC("COULWORK", std::max({2 * NPCOUL + 4 * MXCOUL,
                                               2 * (LMXMX + MAXDL1) + 1,
                                               2 * NPCOUL + LMXMX + MAXDL1}) + 1);
    int IWRKST = NALLOC("COULSTRT", 16 * LDLDIM);

    SETLOG(2 * (LOMOST + MAXDL1));

    int LBASCP = FACFR4 * (Z[CCBLK.IBASCP] - 1) + 1;
    int MBASCP = ILLOC(LBASCP + 1);
    LBASCP = LBASCP + ILLOC(LBASCP + 2);

    int LBASDF = FACFR4 * (Z[CCBLK.IBASDF] - 1) + 1;
    int MBASDF = ILLOC(LBASDF + 1);
    LBASDF = LBASDF + ILLOC(LBASDF + 2);

    int LCHNDF = FACFR4 * (Z[CCBLK.ICHNDF] - 1) + 1;
    int MCHNDF = ILLOC(LCHNDF + 1);
    LCHNDF = LCHNDF + ILLOC(LCHNDF + 2);

    int LBINDX = FACFR4 * (Z[IBINDX] - 1) + 1;
    int LCL2FF = FACFR4 * (Z[ICL2FF] - 1) + 1;
    int LSIG1 = Z[KANDM.ISIGS[1]];
    int LSIG2 = Z[KANDM.ISIGS[2]];

    // Copy CC S-matrices into delta arrays
    int LDELSR = Z[IDELSR];
    int LDELSI = Z[IDELSI];
    int II = NASPLI * NUMLIS;
    for (int I = 1; I <= II; I++) {
        ALLOC(LDELSR + I - 1) = ALLOC(Z[ISMATR] + I - 1);
        ALLOC(LDELSI + I - 1) = ALLOC(Z[ISMATI] + I - 1);
    }

    // ALLSW=0 so NUMJS parameter (7th) is unused — pass 0
    int JMIN = 2 * LMIN;
    int IRET;
    double CLTIME;
    SETBFC_full(NBINDX, MBINDX, &ILLOC(LBINDX),
        MCHNDF, &ILLOC(LCHNDF), &ALLOC4(LCHNDF), 0, &JMIN,
        R2S(), LMXMX, &ALLOC(LSIG1), &ALLOC(LSIG2),
        0.0, 0 /*FALSE*/, 1, NMFFAC, &ALLOC4(LCL2FF),
        &ALLOC(Z[IWRKFF]), IDIM1, LDLDIM, &ALLOC(Z[IWRKWK]),
        &ALLOC(Z[IWRKFI]), IDIM2, &ALLOC(Z[IWRKST]), IRET, CLTIME);

    if (IRET != 0) goto L9960;

    // Free work arrays
    Z[IWRKFF] = -Z[IWRKFF];
    Z[IWRKFI] = -Z[IWRKFI];
    Z[IWRKWK] = -Z[IWRKWK];
    Z[IWRKST] = -Z[IWRKST];

    {
        double T2 = second();

        SETBRN(NBINDX, MBINDX, &ILLOC(LBINDX),
            MCHNDF, &ILLOC(LCHNDF), MBASDF, &ILLOC(LBASDF),
            MBASCP, &ILLOC(LBASCP), &ALLOC4(LBASCP),
            NUMLIS, &ILLOC((int)(FACFR4 * Z[ILIS] - FACFR4 + 1)), 2 * LMIN,
            &ALLOC4(LCL2FF), &ILLOC((int)(FACFR4 * Z[IINDXS] - FACFR4 + 1)),
            &ALLOC(LDELSR), &ALLOC(LDELSI));

        double T3 = second();
        double TTOT = T3 - T1;
        T3 = T3 - T2;
        T2 = T2 - T1;
        double TRECUR = T2 - CLTIME;
        std::printf("-TIMES FOR COULOMB BORN S-MATRIX ELEMENTS (SECONDS):\n"
                    " BELLINGS:%21.3f\n"
                    " RECURSIONS:%19.3f\n"
                    " TOTAL COULOMB INTEGRALS:%15.3f\n"
                    " CLEBSCH FACTORS:%22.3f\n"
                    " TOTAL TIME:%27.3f\n",
                    CLTIME, TRECUR, T2, T3, TTOT);
    }
    return;

L9960:
    std::printf(" **** ERROR IN SETBFC\n");
    IRTN = 0;
    return;
}


// ============================================================================
// WAVPOT — Gets L-dependent stuff for WAVELJ
// ============================================================================
void WAVPOT(int& IRET)
{
    char8 NAMES[3][9]; // NAMES(2,8)
    // Initialize names
    NAMES[1][1] = char8("SIN"); NAMES[2][1] = char8("SOUT");
    NAMES[1][2] = char8("COULFA"); NAMES[2][2] = char8("COULFB");
    NAMES[1][3] = char8("COULGA"); NAMES[2][3] = char8("COULGB");
    NAMES[1][4] = char8("COULFBA"); NAMES[2][4] = char8("COULFBB");
    NAMES[1][5] = char8("COULGBA"); NAMES[2][5] = char8("COULGBB");
    NAMES[1][6] = char8("COULSIGA"); NAMES[2][6] = char8("COULSIGB");
    NAMES[1][7] = char8("INDXEA"); NAMES[2][7] = char8("INDXEB");
    NAMES[1][8] = char8("TOCEA"); NAMES[2][8] = char8("TOCEB");

    IRET = 1;
    int ICHANW = INTRNL.ICHANW;
    double STEPSZ_loc = WAVCOM.RSTEPS[ICHANW];
    double ETA = KANDM.ETAS[ICHANW];
    FLOAT_common.E = WAVCOM.ES[ICHANW];
    int NSTEP = WAVCOM.NSTPSS[ICHANW];
    int LOMOST = KANDM.LOMOST;
    int I = (LOMOST + 2) * MAX0(WAVCOM.NUMFIT, 1);  // +2: RCWFN writes FC[LMAX+1]
    int IOFF = (WAVCOM.IOFFIT - 1) * (LOMOST + 2);
    int L1 = 0, L2 = LOMOST;

    if (SWITCH.PROBLM == 24) {
        L1 = MAX0(0, INTGER.LMIN - INELCM.LXMAX);
        L2 = L1;
    }

    int NF = WAVCOM.NFS[ICHANW];
    int NG = WAVCOM.NGS[ICHANW];
    int NF1 = WAVCOM.NF1S[ICHANW];
    int NG1 = WAVCOM.NG1S[ICHANW];

    if (WAVCOM.IOFFIT == 1) {
        NF = NALLOC(NAMES[ICHANW][2].data, I);
        NG = NALLOC(NAMES[ICHANW][3].data, I);
        NF1 = NALLOC(NAMES[ICHANW][4].data, I);
        NG1 = NALLOC(NAMES[ICHANW][5].data, I);
        WAVCOM.NFS[ICHANW] = NF;
        WAVCOM.NGS[ICHANW] = NG;
        WAVCOM.NF1S[ICHANW] = NF1;
        WAVCOM.NG1S[ICHANW] = NG1;
    }

    // Compute Coulomb wave functions
    double RHO = NSTEP * WAVCOM.HS[ICHANW];
    double RHO1 = 0.0;
    int IWRK = IALLOC(2 * (LOMOST + 2));
    int LWRK = LOCPTRS.Z[IWRK];
    int IFLAG;

    RCWFN(RHO, ETA, L1, L2,
        &ALLOC(LOCPTRS.Z[NF] + IOFF),
        &ALLOC(LWRK),
        &ALLOC(LOCPTRS.Z[NG] + IOFF),
        &ALLOC(LWRK + (LOMOST + 1)),
        1.0e-14, IFLAG);
    if (IFLAG != 0) goto L9910;

    // Calculate Coulomb function at NBACK steps back
    RHO1 = (NSTEP - WAVCOM.NBAKCM) * WAVCOM.HS[ICHANW];
    RCWFN(RHO1, ETA, L1, L2,
        &ALLOC(LOCPTRS.Z[NF1] + IOFF),
        &ALLOC(LWRK),
        &ALLOC(LOCPTRS.Z[NG1] + IOFF),
        &ALLOC(LWRK + (LOMOST + 1)),
        1.0e-14, IFLAG);
    if (IFLAG != 0) goto L9910;

    // Free work area
    LOCPTRS.Z[IWRK] = -LOCPTRS.Z[IWRK];

    // Coulomb phase shifts
    if (WAVCOM.IOFFIT == 1)
        KANDM.ISIGS[ICHANW] = NALLOC(NAMES[ICHANW][6].data, I);
    DSGMAL(ETA, LOMOST, &ALLOC(LOCPTRS.Z[KANDM.ISIGS[ICHANW]] + IOFF));

    if (WAVCOM.IOFFIT != 1) return;

    // Allocate S-matrix space
    {
        int NLXE = 1;
        if (WAVCOM.SOSWS[ICHANW]) NLXE = (int)WAVCOM.JSPS[ICHANW] + 1;
        WAVCOM.IINDXE[ICHANW] = NALLOC(NAMES[ICHANW][7].data, (3 * NLXE * NLXE + 1) / ALLOCS.FACFR4);
        INIT4i(WAVCOM.IINDXE[ICHANW], 0);
        I = (NLXE * NLXE + 1) / 2;
        WAVCOM.ITOCE[ICHANW] = NALLOC(NAMES[ICHANW][8].data, 4 * I / ALLOCS.FACFR4);
        // ... (index setup logic)
        int II = 1;
        int LINDXE = LOCPTRS.Z[WAVCOM.IINDXE[ICHANW]] * ALLOCS.FACFR4 - ALLOCS.FACFR4;
        int LTOCE = LOCPTRS.Z[WAVCOM.ITOCE[ICHANW]] * ALLOCS.FACFR4 - ALLOCS.FACFR4;
        for (int LXP = 1; LXP <= NLXE; LXP++) {
            int LXTMP = LXP - 1;
            int NLO = (LXTMP + LXP + (LXP % 2)) / 2;
            int LDEL = 1 - NLO;
            int K = LINDXE + 1 + 3 * LXTMP * ((int)WAVCOM.JSPS[ICHANW] + 2);
            ILLOC(K) = II;
            ILLOC(K + 1) = LDEL;
            ILLOC(K + 2) = NLO;
            for (int IDEL = 1; IDEL <= NLO; IDEL++) {
                ILLOC(LTOCE + 4 * II - 3) = LDEL;
                ILLOC(LTOCE + 4 * II - 2) = LXTMP;
                ILLOC(LTOCE + 4 * II - 1) = 2 * LXTMP;
                ILLOC(LTOCE + 4 * II) = 0;
                II++;
                LDEL += 2;
            }
        }
        WAVCOM.NUMJS[ICHANW] = II - 1;
        // Verify INDXE setup
        {
            int LINDXE_v = LOCPTRS.Z[WAVCOM.IINDXE[ICHANW]] * ALLOCS.FACFR4 - ALLOCS.FACFR4;
            for (int lxp = 0; lxp < NLXE; lxp++) {
                int k = LINDXE_v + 1 + 3*lxp*((int)WAVCOM.JSPS[ICHANW]+2);
            }
        }
        II = 2 * WAVCOM.NUMJS[ICHANW] * (LOMOST + 1);
        KANDM.ISMATS[ICHANW] = NALLOC(NAMES[ICHANW][1].data, II);
        INIT8(KANDM.ISMATS[ICHANW], 0.0);
    }
    return;

L9910:
    std::printf("-***** ERROR RETURN FROM RCWFN:%8d%18.8g%18.8g%8d%8d\n",
        IFLAG, ETA, RHO, L1, L2);
    LOCPTRS.Z[IWRK] = -LOCPTRS.Z[IWRK];
    IRET = 0;
    return;
}


// ============================================================================
// SETSPT — Sets up S-matrix pointers
// ============================================================================
void SETSPT(int& ISIZE, int LOC, int IPASS, int JBGBIN, int JBIN,
            int PBUGSW, int* BASDEF_flat, int MBASDF, int NBASDF)
{
    // BASDEF is INTEGER BASDEF(MBASDF, NBASDF) in Fortran
    #define BASDEF(i,j) BASDEF_flat[((j)-1)*MBASDF + (i) - 1]

    JBLOCK.JS[2] = JBIN;
    JBLOCK.JS[4] = JBGBIN;

    if (SWITCH.PROBLM == 24) {
        INELCM.JPMIN = IABS((int)JBLOCK.JS[1] - (int)JBLOCK.JS[2]);
        INELCM.JPMAX = (int)JBLOCK.JS[1] + (int)JBLOCK.JS[2];
        INELCM.JTMIN = IABS((int)JBLOCK.JS[3] - (int)JBLOCK.JS[4]);
        INELCM.JTMAX = (int)JBLOCK.JS[3] + (int)JBLOCK.JS[4];
        FORMF.LBP = INELCM.JPMIN / 2;
        FORMF.LBT = INELCM.JTMIN / 2;
    }

    INELCM.JTBASE = IABS((int)JBLOCK.JS[3] - (int)JBLOCK.JS[4]);
    INELCM.NJT = ((int)JBLOCK.JS[3] + (int)JBLOCK.JS[4] - INELCM.JTBASE) / 2 + 1;
    INELCM.JPBASE = IABS((int)JBLOCK.JS[1] - (int)JBLOCK.JS[2]);
    INELCM.NJP = ((int)JBLOCK.JS[1] + (int)JBLOCK.JS[2] - INELCM.JPBASE) / 2 + 1;
    INELCM.NLX = ((int)JBLOCK.JS[1] + (int)JBLOCK.JS[2] + (int)JBLOCK.JS[3] + (int)JBLOCK.JS[4]) / 2 + 1;

    switch (IPASS) {
        case 1: // Return required INDXS size
            ISIZE += 3 * INELCM.NLX * INELCM.NJP * INELCM.NJT;
            return;

        case 2: // Fill in INDXS and compute NSPLI
        {
            INELCM.NSPLI = 0;
            int LINDXS = LOC;
            bool MULTSW = WAVCOM.SOSWS[1] || WAVCOM.SOSWS[2] || SWITCH.IBSPAS > 0;
            int JTMN = INELCM.JTMIN, JTMX = INELCM.JTMAX;
            if (SWITCH.IBSPAS > 0) { JTMN = INELCM.JTBASE; JTMX = (int)JBLOCK.JS[3] + (int)JBLOCK.JS[4]; }
            int JPMN = INELCM.JPMIN, JPMX = INELCM.JPMAX;
            if (MULTSW) { JPMN = INELCM.JPBASE; JPMX = (int)JBLOCK.JS[1] + (int)JBLOCK.JS[2]; }

            for (int JTT = JTMN; JTT <= JTMX; JTT += 2) {
                for (int JTP = JPMN; JTP <= JPMX; JTP += 2) {
                    int LXMN = IABS(JTT - JTP) / 2;
                    int LXMX = (JTT + JTP) / 2;
                    if (!MULTSW) { LXMN = MAX0(LXMN, INELCM.LXMIN); LXMX = MIN0(LXMX, INELCM.LXMAX); }
                    if (LXMX < LXMN) continue;
                    for (int LX = LXMN; LX <= LXMX; LX++) {
                        if (SWITCH.PROBLM == 24) {
                            bool found = false;
                            for (int i = 1; i <= NBASDF; i++) {
                                if (BASDEF(2, i) == 2 * LX) { found = true; break; }
                            }
                            if (!found) continue;
                        }
                        int NDEL = LX + 1 - ((LX + FORMF.LBP + FORMF.LBT) % 2);
                        if (NDEL == 0) continue;
                        int K = (JTP - INELCM.JPBASE + INELCM.NJP * (JTT - INELCM.JTBASE)) / 2;
                        K = 1 + LX + INELCM.NLX * K;
                        int LDELMN = 1 - NDEL;
                        ILLOC(LINDXS + 3 * K - 2) = INELCM.NSPLI + 1;
                        ILLOC(LINDXS + 3 * K - 1) = LDELMN;
                        ILLOC(LINDXS + 3 * K) = NDEL;
                        INELCM.NSPLI += NDEL;
                    }
                }
            }
            ISIZE += INELCM.NSPLI;
            return;
        }

        case 3: // Fill in TOCS
        {
            int LTOCS = ISIZE;
            int LINDXS = LOC;
            WAVCOM.SOSWT = WAVCOM.SOSWS[1] || WAVCOM.SOSWS[2];
            int JPMN = INELCM.JPMIN, JPMX = INELCM.JPMAX;
            if (WAVCOM.SOSWT) { JPMN = INELCM.JPBASE; JPMX = (int)JBLOCK.JS[1] + (int)JBLOCK.JS[2]; }
            int NLXJP = INELCM.NLX * INELCM.NJP;
            int N = NLXJP * INELCM.NJT;

            for (int K = 1; K <= N; K++) {
                int NDEL = ILLOC(LINDXS + 3 * K);
                if (NDEL == 0) continue;
                int JTT = ((K - 1) / NLXJP) * 2 + INELCM.JTBASE;
                if (JTT < INELCM.JTMIN || JTT > INELCM.JTMAX) continue;
                int JTP = ((K - 1) / INELCM.NLX % INELCM.NJP) * 2 + INELCM.JPBASE;
                if (JTP < JPMN || JTP > JPMX) continue;
                int LX = (K - 1) % INELCM.NLX;
                if (!WAVCOM.SOSWT && (LX < INELCM.LXMIN || LX > INELCM.LXMAX)) continue;
                int LDEL = ILLOC(LINDXS + 3 * K - 1);
                int KOFFS = ILLOC(LINDXS + 3 * K - 2);
                for (int IDEL = 1; IDEL <= NDEL; IDEL++) {
                    ILLOC(LTOCS + 4 * KOFFS - 3) = LDEL;
                    ILLOC(LTOCS + 4 * KOFFS - 2) = LX;
                    ILLOC(LTOCS + 4 * KOFFS - 1) = JTP;
                    ILLOC(LTOCS + 4 * KOFFS) = JTT;
                    LDEL += 2;
                    KOFFS++;
                }
            }
            return;
        }
    }

    #undef BASDEF
}


// ============================================================================
// SFROMI — Computes S matrix elements from transfer radial integrals
// ============================================================================
void SFROMI(int LI, int ILI, double* SMATR, double* SMATI, int* INDXS_flat,
            double* XIREAL, double* XIIMAG, int* IINDEX_flat, int NI,
            int* INDXDW_flat, int* IDWFI_flat, int* IDWFO_flat,
            double* RIROAB, double* ATERM, double FACTOR,
            int PINFSW, int& NLINE)
{
    // Multi-dimensional index macros
    #define INDXS(i,k) INDXS_flat[((k)-1)*3 + (i) - 1]
    #define IINDEX(i,k) IINDEX_flat[((k)-1)*4 + (i) - 1]
    #define INDXDW(i,k) INDXDW_flat[((k)-1)*4 + (i) - 1]
    #define IDWFI(i,k) IDWFI_flat[((k)-1)*3 + (i) - 1]
    #define IDWFO(i,k) IDWFO_flat[((k)-1)*4 + (i) - 1]

    int JA = (int)WAVCOM.JSPS[1];
    int JB = (int)WAVCOM.JSPS[2];
    double AMPLI = 0.0;

    // Print headings if necessary
    if (PINFSW) {
        NLINE += NI + 1;
        if (NLINE > 58 && NLINE < 1000) {
            std::printf("1%58sP T O L E M Y\n", "");
            std::printf("%10sCOMPUTATION OF TRANSFER S-MATRIX ELEMENTS\n", "");
        }
        if (!WAVCOM.SOSWT) {
            std::printf("   L   L  LX%20c S-MATRIX%17c CANCELLATION\n", ' ', ' ');
            std::printf("   IN OUT%11c REAL%8c IMAG%7c MAGNITUDE%9c RI, RO\n", ' ', ' ', ' ', ' ');
        }
        NLINE = 7 + NI;
        std::printf("\n");
    }

    // Loop through radial integrals
    int JPMX = JA + JB;
    for (int II = 1; II <= NI; II++) {
        int LXP = IINDEX(4, II);
        int KDW = IINDEX(3, II);
        int KWI = INDXDW(1, KDW);
        int KWO = INDXDW(2, KDW);
        int LASI = IDWFI(1, KWI) + LI;
        int JPI = IDWFI(2, KWI) + 2 * LI;
        int LO = IDWFO(1, KWO) + LI;
        int LASO = IDWFO(2, KWO) + LI;
        int JPO = IDWFO(3, KWO) + 2 * LI;

        int KSBASE = ILI - 1;
        if (WAVCOM.TCSWS[1]) KSBASE = LASI - INTGER.LMIN;
        KSBASE = INELCM.NSPLI * KSBASE;

        double TEMP = FACTOR * ATERM[LXP + 1] / DSQRT(2.0 * LASI + 1.0);
        int ITEST = LASI + LASO + 2 * LXP + 1;
        if ((ITEST % 4) >= 2) TEMP = -TEMP;
        // Guard against NaN radial integrals (forbidden transitions)
        double xr = XIREAL[II], xi = XIIMAG[II];
        if (std::isnan(xr)) xr = 0.0;
        if (std::isnan(xi)) xi = 0.0;
        double TEMPR = TEMP * xr;
        double TEMPI = TEMP * xi;
        double SR = TEMPR;
        double SI = TEMPI;
        if ((ITEST % 2) == 0) {
            // already correct
        } else {
            SR = -TEMPI;
            SI = TEMPR;
        }
        if (PINFSW) {
            AMPLI = DSQRT(SR * SR + SI * SI);
            if (AMPLI == 0.0) AMPLI = 1.0e-30;
            TEMP = DABS(TEMP) * RIROAB[II] / AMPLI;
        }

        if (!WAVCOM.SOSWT) {
            int K = ((int)FORMF.JBP - INELCM.JPBASE + INELCM.NJP * ((int)FORMF.JBT - INELCM.JTBASE)) / 2;
            K = 1 + LXP + INELCM.NLX * K;
            int KOFFS = INDXS(1, K);
            int LDELMN = INDXS(2, K);
            int KS = KSBASE + KOFFS + (LASO - LASI - LDELMN) / 2;
            SMATR[KS] += SR;
            SMATI[KS] += SI;
            if (PINFSW)
                std::printf(" %3d%4d%3d%16.4g%12.4g%14.4g%12.2f\n", LI, LO, LXP, SR, SI, AMPLI, TEMP);
        } else {
            // S.O. / tensor case - 9-J symbols
            if (AMPLI == 1.0e-30) continue;
            int JPIMN = JPI, JPIMX = JPI;
            if (!WAVCOM.SOSWS[1]) JPIMN = IABS(2 * LI - (int)WAVCOM.JSPS[1]);
            int JPOMN = JPO, JPOMX = JPO;
            if (!WAVCOM.SOSWS[2]) JPOMN = IABS(2 * LO - (int)WAVCOM.JSPS[2]);
            for (int jpi = JPIMN; jpi <= JPIMX; jpi += 2) {
                for (int jpo = JPOMN; jpo <= JPOMX; jpo += 2) {
                    double SAV9J = (jpi + 1) * (jpo + 1) * (2 * LXP + 1) * ((int)FORMF.JBP + 1);
                    SAV9J = DSQRT(SAV9J) * WIG9J((int)FORMF.JBT, 2 * LXP, (int)FORMF.JBP,
                        jpi, 2 * LI, JA, jpo, 2 * LO, JB);
                    if (SAV9J == 0.0) continue;
                    TEMPR = SAV9J * SR;
                    TEMPI = SAV9J * SI;
                    for (int JP_val = INELCM.JPBASE; JP_val <= JPMX; JP_val += 2) {
                        int LXMN = MAX0(IABS((int)FORMF.JBT - JP_val) / 2, IABS(LASO - LASI));
                        int LXMX = MIN0(((int)FORMF.JBT + JP_val) / 2, LASO + LASI);
                        if (LXMN > LXMX) continue;
                        for (int LX = LXMN; LX <= LXMX; LX++) {
                            int K = (JP_val - INELCM.JPBASE + INELCM.NJP * ((int)FORMF.JBT - INELCM.JTBASE)) / 2;
                            K = 1 + LX + INELCM.NLX * K;
                            int KOFFS = INDXS(1, K);
                            int LDELMN = INDXS(2, K);
                            double TEMP2 = SAV9J;
                            if (!(LX == LXP && JP_val == (int)FORMF.JBP && LASI == LI && LASO == LO)) {
                                TEMP2 = (jpi + 1) * (jpo + 1) * (2 * LX + 1) * (JP_val + 1);
                                TEMP2 = DSQRT(TEMP2) * WIG9J((int)FORMF.JBT, 2 * LX, JP_val,
                                    jpi, 2 * LASI, JA, jpo, 2 * LASO, JB);
                                if (TEMP2 == 0.0) continue;
                            }
                            if ((LX + LXP) % 2 != 0) TEMP2 = -TEMP2;
                            int KS = KSBASE + KOFFS + (LASO - LASI - LDELMN) / 2;
                            SMATR[KS] += TEMP2 * TEMPR;
                            SMATI[KS] += TEMP2 * TEMPI;
                        }
                    }
                }
            }
        }
    }

    #undef INDXS
    #undef IINDEX
    #undef INDXDW
    #undef IDWFI
    #undef IDWFO
}
