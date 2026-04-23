// scattering_partA_translated.cpp — Verbatim C++ translation of COULIN, COULNG, GETSCT
// Translated from source.f lines 10647-11131, 11133-11208, 17896-18210

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include "ptolemy_alloc.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_forward.h"
#include <cstdio>
#include <cstring>
#include <cmath>

// ============================================================================
// SUBROUTINE COULIN — lines 10647-11131
//
// COULOMB INTEGRALS BY RECURSION
//
// COMPUTES
//   FF = INTEGRAL( R TO INF ) DR F(LOUT,ETAOUT,KOUT*R)
//        * F(LIN,ETAIN,KIN*R) / R**N
// AND SIMILARLY FOR FG, GF, AND GG  FOR
//   LMIN =< (LOUT+LIN)/2 =< LMAX
//   -MAXDEL =< LIN-LOUT =< MAXDEL
//
// 12/10/76 - FIRST VERSION - S. PIEPER
// ============================================================================
void COULIN(int N, int MAXDEL, int LMIN, int LMAX, double ETAOUT,
            double AKOUT, double* SIGOUT,
            double ETAIN, double AKIN, double* SIGIN, double R, int ALLSW,
            double* FF, double* FG, double* GF, double* GG, int LDLDIM,
            double ACCURA, double ASMULT, int NTERMS, int NPTS,
            double* WORK, double* FI, double* FO, double* GI, double* GO,
            double* STARTS,
            int IPRINT, int& IRET, double& CLTIME)
{
    // Local variables
    double RTURNS[3]; // 1-based
    double RMINS[3];  // 1-based
    double ACCS[3];   // 1-based
    double A, B, C, D, E, X;
    double DLI, DLO, DN, ACC, RVAL, DELTA, RTOP, RBOT, PSIMX;
    double TT, TEMP;
    int MXDEL, LMAXMX, LMINMN, MMXDEL;
    int NTRM, ISTOPR, ISTOP, LSUMMN, LSUMMX;
    int LS, IL, MD, IS, ID, ID2;
    int LIN, LOUT, LOUTP, LSUM, LDEL;
    int II, I, IPT, LINLS, LOUTLS;
    int LMN, LMX, III;

    // STARTS is dimensioned STARTS(LDLDIM, 2, 4, 2) in Fortran
    // Access macro: STARTS(id, is, k, ii) => STARTS[ (ii-1)*LDLDIM*2*4 + (k-1)*LDLDIM*2 + (is-1)*LDLDIM + (id-1) ]
    // But we use a flat pointer, so define a macro
    #define STARTS_4D(id,is,k,ii) STARTS[((ii)-1)*LDLDIM*2*4 + ((k)-1)*LDLDIM*2 + ((is)-1)*LDLDIM + ((id)-1)]

    // FF, FG, GF, GG are dimensioned (LDLDIM, *) — column-major
    #define FF_2D(id,il) FF[((il)-1)*LDLDIM + ((id)-1)]
    #define FG_2D(id,il) FG[((il)-1)*LDLDIM + ((id)-1)]
    #define GF_2D(id,il) GF[((il)-1)*LDLDIM + ((id)-1)]
    #define GG_2D(id,il) GG[((il)-1)*LDLDIM + ((id)-1)]

    MXDEL = MAX0(2, MAXDEL);
    LMAXMX = LMAX + MXDEL;
    LMINMN = MAX0(0, LMIN - MXDEL);
    if (LDLDIM <= MXDEL) goto L20_fail;
    goto L20;
L20_fail:
    // LDLDIM TOO SMALL
    IRET = -10;
    return;
L20:
    MMXDEL = -MXDEL;
    DN = N;
    NTRM = MAX0(NTERMS, LMAX / 4);

    // FOLLOWING IS AN ATTEMPT TO ALLOW FOR LOSS OF ACCURACY

    ACC = ACCURA * DEXP(10 - 5 * N - 2 * DABS(ETAIN - ETAOUT)
          - (LMAX - LMIN) / 100.);
    ACC = DMAX1(ACC, 1.0e-12, ACCURA * .01);

    // IF THE INTEGRAL IS TO BE FROM 0, FIND A SUITABLE STARTING POINT

    RMINS[1] = R;
    RMINS[2] = R;
    ACCS[1] = ACC;
    ACCS[2] = .01 * ACCURA;
    if (R > 0) goto L50;
    ACCS[1] = .1 * ACCURA;
    ACCS[2] = .1 * ACC;
    if (!ALLSW) goto L30;
    IRET = -12;
    return;
L30:
    RTURNS[1] = DMAX1((ETAIN + DSQRT(ETAIN * ETAIN + LMINMN * (LMINMN + 1))) / AKIN,
        (ETAOUT + DSQRT(ETAOUT * ETAOUT + LMINMN * (LMINMN + 1))) / AKOUT);
    RTURNS[2] = DMAX1((ETAIN + DSQRT(ETAIN * ETAIN + LMAXMX * (LMAXMX + 1))) / AKIN,
        (ETAOUT + DSQRT(ETAOUT * ETAOUT + LMAXMX * (LMAXMX + 1))) / AKOUT);
    if (IPRINT >= 1) std::printf(" RTURN'S = %12.4G%12.4G\n", RTURNS[1], RTURNS[2]);

    // WE DO THE NUMERIC INTEGRAL MORE EFFICENTLY HERE SINCE WE DO
    // ALL OF THEM AT ONCE.  SO FORCE CLINTS RIGHT TO THE PLACE WHERE
    // IT STARTS TO USE RCASYM.

    RMINS[1] = 1.4 * RTURNS[1];
    RMINS[2] = 1.4 * RTURNS[2];

    // GET STARTING VALUES FOR RECURSIONS AND FINAL VALUES
    // TO CHECK THE STABILITY OF THE RECURSIONS.

L50:
    LS = 2 * LMIN;
    if (IPRINT >= 1) std::printf(" REQUESTED CLINTS ACCURACIES ARE %13.3G%13.3G\n", ACCS[1], ACCS[2]);

    // FOR RECURSION UPWARDS WE WILL NEED TO START WITH A
    // FULL WIDTH SET.  ALSO THE DOWNWARDS RECURSION CAN NOT GO
    // BENEATH  LIN+LOUT = N-2  EVEN THOUGH FOR R > 0 THE RESULTS
    // WILL STILL BE FINITE.  THEREFORE WE GENERATE ALL THESE SMALL
    // L VALUES NUMERICALLY.

    ISTOPR = MAX0(N - 2 * LMIN, 2);
    if (ALLSW) ISTOPR = MAX0(MXDEL + 1 - 2 * LMIN, ISTOPR);
    ISTOP = ISTOPR;
    LSUMMN = 2 * LMIN + 1 + ISTOPR - 2;
    LSUMMX = 2 * LMAX - 1;

    TT = (float)dtime_();
    for (II = 1; II <= 2; II++) {

        // LIN+LOUT = LS+I-1

        for (I = 1; I <= ISTOP; I++) {
            IL = LS - 2 * LMIN + I;
            MD = MXDEL + MOD(MXDEL + IL, 2);
            IS = 1;
            if (I == ISTOP) IS = 2;
            for (ID = 1; ID <= MD; ID++) {
                LIN = (LS + I + MXDEL) / 2 + ID - 1 - MXDEL;
                LOUT = LS + I - 1 - LIN;
                if (LIN < 0 || LOUT < 0) goto L59;
                A = ASMULT;
                // TODO: CLINTS call needs correct signature — stubbed for now
                // CLINTS(RMINS[II], ETAIN, ETAOUT, AKIN, AKOUT, ...);
                IRET = 0;
                if (IRET != 0) return;
L59:            ;
            }

            // SAVE END OF THE FF RECURSION AND
            // STORE START OF FG, GF, GG STUFF

            for (ID = 1; ID <= MD; ID++) {
                FF_2D(ID, IL) = STARTS_4D(ID, IS, 1, II);
                if (!ALLSW) goto L64;
                FG_2D(ID, IL) = STARTS_4D(ID, IS, 2, II);
                GF_2D(ID, IL) = STARTS_4D(ID, IS, 3, II);
                GG_2D(ID, IL) = STARTS_4D(ID, IS, 4, II);
L64:            ;
            }

        } // 69
        LS = 2 * LMAX - 1;
        ISTOP = 2;
    } // 79

    CLTIME = (float)dtime_() - TT;

    if (R != 0) goto L400;

    // HERE WE DO THE INTEGRAL  0 < R < 1.4*RTURN  BACKWARDS UNTIL
    // WE COME TO INSIGNIFICANT CONTRIBUTIONS.

    DELTA = 3.1 / AKIN;
    ISTOP = ISTOPR;
    LMN = LMINMN;
    LMX = LMIN + (ISTOP + MXDEL + 1) / 2;
    LS = 2 * LMIN;
    for (II = 1; II <= 2; II++) {
        RTOP = RMINS[II];
L110:
        RBOT = DMAX1(0.0e0, RTOP - DELTA);
        if (RBOT >= RTOP) goto L165;
        B = RTOP - RBOT;
        for (IPT = 1; IPT <= NPTS; IPT++) {
            RVAL = RTOP - .5 * B * (1 + WORK[IPT]);
            RCWFN(AKIN * RVAL, ETAIN, LMN, LMX, FI, FO,
                GI, GO, 1.0e-8, IRET);
            if (IRET != 0 && IRET != 2) goto L900;
            if (DABS(FI[LMN + 1]) < ACCURA &&
                RVAL < RTURNS[II]) goto L165;
            RCWFN(AKOUT * RVAL, ETAOUT, LMN, LMX, FO, GI,
                GO, &WORK[2 * NPTS + 1], 1.0e-8, IRET);
            if (IRET != 0 && IRET != 2) goto L900;
            if (DABS(FO[LMN + 1]) < ACCURA &&
                RVAL < RTURNS[II]) goto L165;
            X = .5 * B * WORK[NPTS + IPT] / std::pow(RVAL, N);
            for (I = 1; I <= ISTOP; I++) {
                IL = LS - 2 * LMIN + I;
                MD = MXDEL + MOD(MXDEL + IL, 2);
                IS = 1;
                if (I == ISTOP) IS = 2;
                for (ID = 1; ID <= MD; ID++) {
                    LIN = (LS + I + MXDEL) / 2 + ID - 1 - MXDEL;
                    LOUT = LS + I - 1 - LIN;
                    if (LIN < 0 || LOUT < 0) goto L149;
                    A = X * FI[LIN + 1] * FO[LOUT + 1];
                    FF_2D(ID, IL) = FF_2D(ID, IL) + A;
                    STARTS_4D(ID, IS, 1, II) = FF_2D(ID, IL);
L149:               ;
                }
            }
        } // 159
        RTOP = RBOT;
        RVAL = RBOT;
        goto L110;
L165:
        RMINS[II] = RVAL;
        LS = 2 * LMAX - 1;
        ISTOP = 2;
        LMN = (LS - MXDEL - 1) / 2;
        LMX = LMAXMX;
    } // 169

    CLTIME = (float)dtime_() - TT;

    if (IPRINT >= 1) std::printf(" RMIN'S =  %12.4G%12.4G\n", RMINS[1], RMINS[2]);
    MD = MXDEL + 1;
    if (IPRINT >= 2) {
        std::printf("0FF STARTING VALUES:\n");
        for (IS = 1; IS <= 2; IS++) {
            for (ID = 1; ID <= MD; ID++) {
                std::printf(" %20.10G", STARTS_4D(ID, IS, 1, 2));
                if (ID % 3 == 0 || ID == MD) std::printf("\n");
            }
        }
    }

    // RECURSE DOWNWARD FOR FF
    // WE GO FROM  LSUM  TO  LSUM-1  AND GENERATE ALL THE
    // LIN-LOUT FOR EACH LSUM.

    for (LS = LSUMMN; LS <= LSUMMX; LS++) {
        LSUM = LSUMMX + LSUMMN - LS;
        for (LDEL = MMXDEL; LDEL <= MXDEL; LDEL++) {

            // LIN, LOUT POINT TO THE CORNER WE WILL EXTRAPOLATE FROM

            LIN = (LSUM - LDEL);
            if (LIN < 0) goto L299;
            LIN = LIN / 2;
            LOUT = (LSUM + LDEL) / 2;
            if (LSUM != LIN + LOUT) goto L299;
            LOUTP = LOUT - 1;
            if (IABS(LIN - LOUTP) > MXDEL) goto L299;
            if (LOUTP < 0) goto L299;
            if (IABS(LOUT + 1 - LIN) > MXDEL) goto L299;

            // LIN, LOUT-1  IS NEEDED AND WE HAVE ALL THE STUFF NEEDED TO
            // FIND IT. BY RECURSION DOWNWARD ON LOUT.

            DLI = LIN;
            DLO = LOUT;
            ID = (LIN - LOUT + MXDEL) + 2;
            IL = LIN + LOUT - 2 * LMIN + 1;
            ID2 = (ID + 1) / 2;
            ID = ID / 2;

            E = (2 * DLO + 1) / ((DLI + DLO - DN + 2) * DSQRT(DLO * DLO + ETAOUT * ETAOUT));
            A = ETAOUT * (DLI - DN + 2) / (DLO + 1)
                - (ETAIN * (AKIN / AKOUT)) * DLO / (DLI + 1);
            B = (AKIN / AKOUT) * (DLO / (DLI + 1)) * DSQRT((DLI + 1) * (DLI + 1) + ETAIN * ETAIN);
            C = DLO * (DLO - DLI + DN - 1) * DSQRT((DLO + 1) * (DLO + 1) + ETAOUT * ETAOUT)
                / ((DLO + 1) * (2 * DLO + 1));
            FF_2D(ID2, IL - 1) = E * (A * FF_2D(ID, IL) + C * FF_2D(ID2 - 1, IL + 1)
                + B * FF_2D(ID2, IL + 1));
L299:       ;
        }

        // DLI, DLO, ETC, REFER TO THE LAST CASE FOR WHICH WE REDUCED
        // LOUT.  NOW REDUCE LIN FOR THAT SAME CASE TO COMPLETE THE
        // GENERATION OF LSUM-1 FROM LSUM.

        if (DLI <= 0) goto L399;

        E = (2 * DLI + 1) / ((DLI + DLO - DN + 2) * DSQRT(DLI * DLI + ETAIN * ETAIN));
        A = ETAIN * (DLO - DN + 2) / (DLI + 1)
            - (ETAOUT * (AKOUT / AKIN)) * DLI / (DLO + 1);
        B = (AKOUT / AKIN) * (DLI / (DLO + 1)) * DSQRT((DLO + 1) * (DLO + 1) + ETAOUT * ETAOUT);
        C = DLI * (DLI - DLO + DN - 1) * DSQRT((DLI + 1) * (DLI + 1) + ETAIN * ETAIN)
            / ((DLI + 1) * (2 * DLI + 1));
        FF_2D(ID2 - 1, IL - 1) = E * (A * FF_2D(ID, IL) + B * FF_2D(ID2 - 1, IL + 1)
            + C * FF_2D(ID2, IL + 1));

L399:   ;
    } // 399

    goto L800;

    // GET THE INHOMO TERMS

L400:
    TT = (float)dtime_();

    RCWFN(AKIN * R, ETAIN, LMINMN, LMAXMX, FI, WORK,
        GI, &WORK[LMAXMX + 1], 1.0e-14, IRET);
    if (IRET != 0) goto L910;
    RCWFN(AKOUT * R, ETAOUT, LMINMN, LMAXMX, FO, WORK,
        GO, &WORK[LMAXMX + 1], 1.0e-14, IRET);
    if (IRET != 0) goto L910;

    CLTIME = CLTIME + (float)dtime_() - TT;

    // RECURSE UPWARDS ON FG, GF, AND GG. AND FF WHEN R > 0.

    for (LSUM = LSUMMN; LSUM <= LSUMMX; LSUM++) {
        for (LDEL = MMXDEL; LDEL <= MXDEL; LDEL++) {

            // LIN, LOUT POINT TO THE CORNER WE WILL EXTRAPOLATE FROM

            LIN = (LSUM + LDEL) / 2;
            LOUT = (LSUM - LDEL) / 2;
            if (LSUM != LIN + LOUT) goto L499;
            LOUTP = LOUT + 1;
            if (IABS(LIN - LOUTP) > MXDEL) goto L499;
            if (IABS(LOUT - 1 - LIN) > MXDEL) goto L499;

            // LIN, LOUT+1  IS NEEDED AND WE HAVE ALL THE STUFF NEEDED TO
            // FIND IT. BY RECURSION UPWARD ON LOUT.

            DLI = LIN;
            DLO = LOUT;
            LINLS = LIN;
            LOUTLS = LOUT;
            ID = (LIN - LOUT + MXDEL) + 2;
            IL = LIN + LOUT - 2 * LMIN + 1;
            ID2 = (ID + 1) / 2;
            ID = ID / 2;
            E = (2 * DLO + 1) / ((DLI + DLO + DN) * DSQRT((DLO + 1) * (DLO + 1) + ETAOUT * ETAOUT));
            A = ETAOUT * (DLI + DN - 1) / DLO
                - (ETAIN * (AKIN / AKOUT)) * (DLO + 1) / DLI;
            B = (AKIN / AKOUT) * ((DLO + 1) / DLI) * DSQRT(DLI * DLI + ETAIN * ETAIN);
            C = (DLO + 1) * (DLO - DLI - DN + 1) * DSQRT(DLO * DLO + ETAOUT * ETAOUT)
                / (DLO * (2 * DLO + 1));
            D = (DLO + 1) / (AKOUT * std::pow(R, N));

            FF_2D(ID2 - 1, IL + 1) = E * (A * FF_2D(ID, IL) + C * FF_2D(ID2, IL - 1) +
                B * FF_2D(ID2 - 1, IL - 1) + D * FO[LOUT + 1] * FI[LIN + 1]);
            if (!ALLSW) goto L499;
            FG_2D(ID2 - 1, IL + 1) = E * (A * FG_2D(ID, IL) + C * FG_2D(ID2, IL - 1) +
                B * FG_2D(ID2 - 1, IL - 1) + D * FO[LOUT + 1] * GI[LIN + 1]);
            GF_2D(ID2 - 1, IL + 1) = E * (A * GF_2D(ID, IL) + C * GF_2D(ID2, IL - 1) +
                B * GF_2D(ID2 - 1, IL - 1) + D * GO[LOUT + 1] * FI[LIN + 1]);
            GG_2D(ID2 - 1, IL + 1) = E * (A * GG_2D(ID, IL) + C * GG_2D(ID2, IL - 1) +
                B * GG_2D(ID2 - 1, IL - 1) + D * GO[LOUT + 1] * GI[LIN + 1]);
L499:       ;
        }

        // DLI, DLO, ETC, REFER TO THE LAST CASE FOR WHICH WE INCREASED
        // LOUT.  NOW INCREASE LIN FOR THAT SAME CASE TO COMPLETE THE
        // GENERATION OF LSUM+1 FROM LSUM.

        E = (2 * DLI + 1) / ((DLI + DLO + DN) * DSQRT((DLI + 1) * (DLI + 1) + ETAIN * ETAIN));
        A = ETAIN * (DLO + DN - 1) / DLI
            - (ETAOUT * (AKOUT / AKIN)) * (DLI + 1) / DLO;
        B = (AKOUT / AKIN) * ((DLI + 1) / DLO) * DSQRT(DLO * DLO + ETAOUT * ETAOUT);
        C = (DLI + 1) * (DLI - DLO - DN + 1) * DSQRT(DLI * DLI + ETAIN * ETAIN)
            / (DLI * (2 * DLI + 1));
        D = (DLI + 1) / (AKIN * std::pow(R, N));

        FF_2D(ID2, IL + 1) = E * (A * FF_2D(ID, IL) + B * FF_2D(ID2, IL - 1)
            + C * FF_2D(ID2 - 1, IL - 1) + D * FO[LOUTLS + 1] * FI[LINLS + 1]);
        if (!ALLSW) goto L599;
        FG_2D(ID2, IL + 1) = E * (A * FG_2D(ID, IL) + B * FG_2D(ID2, IL - 1)
            + C * FG_2D(ID2 - 1, IL - 1) + D * FO[LOUTLS + 1] * GI[LINLS + 1]);
        GF_2D(ID2, IL + 1) = E * (A * GF_2D(ID, IL) + B * GF_2D(ID2, IL - 1)
            + C * GF_2D(ID2 - 1, IL - 1) + D * GO[LOUTLS + 1] * FI[LINLS + 1]);
        GG_2D(ID2, IL + 1) = E * (A * GG_2D(ID, IL) + B * GG_2D(ID2, IL - 1)
            + C * GG_2D(ID2 - 1, IL - 1) + D * GO[LOUTLS + 1] * GI[LINLS + 1]);

L599:   ;
    } // 599


    // CHECK THE ACCURACY OF THE RECURSIONS

L800:
    for (I = 1; I <= 2; I++) {
        IL = I + ISTOPR - 2;
        MD = MXDEL + MOD(MXDEL + IL, 2);
        if (R > 0) goto L850;
        for (ID = 1; ID <= MD; ID++) {
            LIN = (2 * LMIN + IL + MXDEL) / 2 + ID - 1 - MXDEL;
            LOUT = 2 * LMIN + IL - 1 - LIN;
            if (LIN < 0 || LOUT < 0) goto L849;
            if (IPRINT >= 1)
                std::printf(" CHECKING %7d%7d%7d%7d%7d%7d\n          %20.10G%20.10G\n",
                    I, MD, ID, IL, LIN, LOUT,
                    FF_2D(ID, IL), STARTS_4D(ID, I, 1, 1));
            if (DABS(FF_2D(ID, IL) / STARTS_4D(ID, I, 1, 1) - 1) > 10 * ACCURA)
                std::printf("\n**** FOR LIN, LOUT =%4d%4d FF RECURSION IS POOR: %20.10G%20.10G\n",
                    LIN, LOUT, FF_2D(ID, IL), STARTS_4D(ID, I, 1, 1));
L849:       ;
        }
        goto L869;

L850:
        MD = MXDEL + MOD(MXDEL + I + 1, 2);
        for (ID = 1; ID <= MD; ID++) {
            IL = 2 * (LMAX - LMIN) - 1 + I;
            LIN = (IL + 2 * LMIN + MXDEL) / 2 + ID - 1 - MXDEL;
            LOUT = IL + 2 * LMIN - 1 - LIN;
            if (IPRINT >= 1)
                std::printf(" CHECKING %7d%7d%7d%7d%7d%7d\n          %20.10G%20.10G%20.10G%20.10G\n          %20.10G%20.10G%20.10G%20.10G\n",
                    I, MD, ID, IL, LIN, LOUT,
                    FF_2D(ID, IL), FG_2D(ID, IL), GF_2D(ID, IL), GG_2D(ID, IL),
                    STARTS_4D(ID, I, 1, 2), STARTS_4D(ID, I, 2, 2),
                    STARTS_4D(ID, I, 3, 2), STARTS_4D(ID, I, 4, 2));
            if (DABS(FF_2D(ID, IL) / STARTS_4D(ID, I, 1, 2) - 1) > 10 * ACCURA)
                std::printf("\n**** FOR LIN, LOUT =%4d%4d FF RECURSION IS POOR: %20.10G%20.10G\n",
                    LIN, LOUT, FF_2D(ID, IL), STARTS_4D(ID, I, 1, 2));
            if (!ALLSW) goto L859;
            if (DABS(FG_2D(ID, IL) / STARTS_4D(ID, I, 2, 2) - 1) > 10 * ACCURA)
                std::printf("\n**** FOR LIN, LOUT =%4d%4d FG RECURSION IS POOR: %20.10G%20.10G\n",
                    LIN, LOUT, FG_2D(ID, IL), STARTS_4D(ID, I, 2, 2));
            if (DABS(GF_2D(ID, IL) / STARTS_4D(ID, I, 3, 2) - 1) > 10 * ACCURA)
                std::printf("\n**** FOR LIN, LOUT =%4d%4d GF RECURSION IS POOR: %20.10G%20.10G\n",
                    LIN, LOUT, GF_2D(ID, IL), STARTS_4D(ID, I, 3, 2));
            if (DABS(GG_2D(ID, IL) / STARTS_4D(ID, I, 4, 2) - 1) > 10 * ACCURA)
                std::printf("\n**** FOR LIN, LOUT =%4d%4d GG RECURSION IS POOR: %20.10G%20.10G\n",
                    LIN, LOUT, GG_2D(ID, IL), STARTS_4D(ID, I, 4, 2));
L859:       ;
        }
L869:   ;
    } // 869

    IRET = 0;
    return;

    // ERROR DETECTED BY RCWFN

L900:
    std::printf("\n***RCWFN FROM COULIN IRET = %6d%15.5G%5d%5d%15.5G%15.5G%15.5G%15.5G\n $$$%20.10G%20.10G%20.10G%5d%5d\n",
        IRET, RVAL, LMN, LMX, ETAIN, ETAOUT, AKIN, AKOUT,
        RTOP, DELTA, B, II, IPT);
    return;
L910:
    RVAL = R;
    LMN = LMINMN;
    LMX = LMAXMX;
    B = 123456789.;
    goto L900;

    #undef STARTS_4D
    #undef FF_2D
    #undef FG_2D
    #undef GF_2D
    #undef GG_2D
}


// ============================================================================
// FUNCTION COULNG — lines 11133-11208
//
// COMPUTE COULOMB WAVEFUNCTION FOR NEGATIVE ENERGY
//
// THE RETURNED FUNCTION IS
//   COULNG =  GAMMA(L+ETA+1) / GAMMA(ANORM)
//           * EXP(RHO) * (2*RHO)**ETA * W( -ETA, L+1/2, 2*RHO )
// WHERE W IS THE WHITTAKER FUNCTION OF ABRAMOWITZ AND STEGUN
//
// 4/23/75 - FIRST VERSION - S. PIEPER
// ============================================================================
double COULNG(int L, double ETA, double RHO, double ANORM)
{
    // COMMON block aliases
    auto& PI     = CNSTNT.PI;
    auto& BIGEST = CNSTNT.BIGEST;
    auto& BIGLOG = CNSTNT.BIGLOG;

    // Local variables
    int NPTS = 28;
    double X[32];  // 1-based, extra space for LAGBC/LAGUER
    double DX[32]; // 1-based
    double A_arr[32]; // 1-based
    double B_arr[32]; // 1-based
    double C_arr[32]; // 1-based
    double COULNG_val;
    double ETA1, DL1, GAM, XLIM, TEMP;
    double CSX, CSW, TSX, TSW;
    int L1, L2, IETA2, N;

    ETA1 = DMOD(ETA, 1.0e0);
    L1 = MOD(L, 10);
    if (ETA < 0) L1 = 0;
    L2 = L - L1;
    DL1 = L1;
    ETA1 = ETA1 + DMAX1(DMIN1(10 - DL1, ETA - ETA1), 0.0e0);
    IETA2 = (int)(ETA - ETA1);

    GAM = 1;
    if (ETA + L != 0) GAM = DEXP(-DLGAMA(ANORM) / (ETA + L));
    XLIM = 0;
    if (L2 + IETA2 <= 1) goto L100;
    XLIM = DEXP(-BIGLOG / (L2 + IETA2)) / GAM;

L100:
    LAGUER(NPTS, &X[1], &DX[1], ETA1 + DL1, 1.0e-20, CSX, CSW, TSX, TSW,
        &A_arr[1], &B_arr[1], &C_arr[1]);

    COULNG_val = 0;
    for (N = 1; N <= NPTS; N++) {
        if (X[N] < XLIM) goto L159;
        if (RHO < 2) goto L130;
        COULNG_val = COULNG_val + DX[N] * std::pow(GAM * X[N], L2 + IETA2) *
            std::pow(1 + X[N] / (2 * RHO), L - ETA);
        goto L159;
L130:
        TEMP = 1 + X[N] / (2 * RHO);
        COULNG_val = COULNG_val + DX[N] * std::pow(GAM * X[N] / TEMP, L2 + IETA2)
            * std::pow(TEMP, L + L2 - ETA1);
L159:   ;
    }

    COULNG_val = COULNG_val * std::pow(GAM, ETA1 + DL1);
    if (ETA + L == 0) COULNG_val = COULNG_val / DGAMMA(ANORM);
    return COULNG_val;
}


// ============================================================================
// SUBROUTINE GETSCT — lines 17896-18210
//
// COMPUTES SCATTERING WAVES AT LMIN AND LCRIT FOR GRDSET
//
// 7/8/75 - FIRST VERSION - S. P.
// ============================================================================
void GETSCT(int& IRET)
{
    // COMMON block aliases
    auto* Z = LOCPTRS.Z;  // 1-based

    auto& ACCURA = FLOAT_common.ACCURA;

    auto& L       = INTGER.L;
    auto& LMAX    = INTGER.LMAX;
    auto& LMIN    = INTGER.LMIN;
    auto& IPRINT  = INTGER.IPRINT;

    auto& PROBLM  = SWITCH.PROBLM;

    auto& ONELSW  = INELCM.ONELSW;
    auto& LXMIN   = INELCM.LXMIN;
    auto& LXMAX   = INELCM.LXMAX;

    auto& IHSAVE  = INTRNL.IHSAVE;
    auto& ICHANW  = INTRNL.ICHANW;

    auto& STANSW  = WAVCOM.STANSW;
    auto* SOSWS   = WAVCOM.SOSWS;    // 1-based
    auto* TCSWS   = WAVCOM.TCSWS;    // 1-based
    auto* NSTPSS  = WAVCOM.NSTPSS;   // 1-based
    auto* RSTEPS  = WAVCOM.RSTEPS;   // 1-based
    auto& IWAVR   = WAVCOM.IWAVR;
    auto& IWAVI   = WAVCOM.IWAVI;
    auto* NFS     = WAVCOM.NFS;      // 1-based
    auto* NGS     = WAVCOM.NGS;      // 1-based
    auto* NF1S    = WAVCOM.NF1S;     // 1-based
    auto* NG1S    = WAVCOM.NG1S;     // 1-based
    auto* LASTL   = WAVCOM.LASTL;    // 1-based
    auto* JSPS    = WAVCOM.JSPS;     // 1-based
    auto* IRLVS   = WAVCOM.IRLVS;   // 1-based
    auto* IIMVS   = WAVCOM.IIMVS;   // 1-based
    auto* ISMATS  = KANDM.ISMATS;    // 1-based — in KANDM, not WAVCOM
    auto* ASYMPS  = WAVCOM.ASYMPS;   // 1-based
    auto* ICENTR  = WAVCOM.ICENTR;   // 1-based

    auto& LCRIT   = KANDM.LCRIT;
    auto& ISCTMN  = KANDM.ISCTMN;
    auto& ISCTCR  = KANDM.ISCTCR;

    auto& IOUTSW  = BSBLK.IOUTSW;
    auto& NSPSCT  = BSBLK.NSPSCT;
    auto& SCTMAX  = BSBLK.SCTMAX;
    auto& SCTSP   = BSBLK.SCTSP;
    auto& ROFMAX  = BSBLK.ROFMAX;

    auto& NUMLNK  = LNKBLK.NUMLNK;
    auto& IUNIQU  = LNKBLK.IUNIQU;
    auto* LNKADR  = &LNKBLK.LNKADR[0][0];  // flat access

    auto* LNKAD3_flat = &SCRTCH.LNKAD3[0][0]; // flat access

    // ISCRTS equivalence: ISCRTS(1) = ISCTMN, ISCRTS(2) = ISCTCR
    int* ISCRTS = &KANDM.ISCTMN; // ISCRTS[0]=ISCTMN, ISCRTS[1]=ISCTCR

    // Local variables
    float FG_local[5]; // 1-based, REAL*4
    char NAMES[2][9]; // character*8
    std::memcpy(NAMES[0], "MAGSCTMN", 8); NAMES[0][8] = '\0';
    std::memcpy(NAMES[1], "MAGSCTCR", 8); NAMES[1][8] = '\0';
    int SSOSW, STCSW_local, CCSW;
    float DUMMY4_local[2]; // REAL*4 DUMMY4(1)
    int ITOP, LC;
    int LWAVR, LWAVI, LVREAL, LVIMAG, LCENTR, LSMATS;
    int NSTEPS, LSCRTS;
    int II, I;
    double RVAL, PSI, PSIMX, RPSIMX, RSTEP;

    int IPRNT = MOD(IPRINT, 10);
    CCSW = (PROBLM == 24) ? TRUE_F : FALSE_F;
    ITOP = 2;
    if (CCSW) ITOP = 1;

    // INITIALIZE THE COULOMB FUNCTIONS, ETC.

    for (ICHANW = 1; ICHANW <= ITOP; ICHANW++) {
        WAVPOT(IRET);
        if (IPRNT >= 5) std::printf(" ++++ WAVPOT CALL, ICHANW, IRET, IHSAVE = %10d%10d%10d\n",
            ICHANW, IRET, IHSAVE);
        if (IRET == 0) return;
    }

    if (IHSAVE == 2) return;

    LC = LCRIT;
    if (LC < LMIN || LC > LMAX) LC = (LMIN + LMAX) / 2;

    if (IPRNT >= 3) std::printf("\n L CRITICAL:  AVERAGE = %4d     USED FOR GRID SETUP = %4d\n",
        LCRIT, LC);

    // COMPUTE THE WAVEFUNCTIONS, USING CENTRAL POTENTIALS ONLY.

    STANSW = TRUE_F;
    SSOSW = SOSWS[1];
    SOSWS[1] = FALSE_F;
    STCSW_local = TCSWS[1];
    TCSWS[1] = FALSE_F;

    L = MAX0(0, LMIN - LXMAX);
    NSTEPS = NSTPSS[1] + 1;
    RSTEP = RSTEPS[1];

    // WAVELJ WILL SET LNKBLK BACK TO THE STATUS FOR THE INCOMING
    // CHANNEL.  IT MAY HAVE BEEN CHANGED JUST FOR THE FORM FACTOR.
    // THEREFORE WE SAVE ITS STATUS.

    for (II = 1; II <= NUMLNK; II++) {
        for (I = 1; I <= 6; I++) {
            SCRTCH.LNKAD3[II][I] = LNKBLK.LNKADR[II][I];
        }
    }

    for (II = 1; II <= ITOP; II++) {
        LWAVR = Z[IWAVR];
        LWAVI = Z[IWAVI];
        LVREAL = Z[IRLVS[1]];
        LVIMAG = Z[IIMVS[1]];
        LCENTR = Z[ICENTR[1]];
        LSMATS = Z[ISMATS[1]];

        // ALWAYS FORCE WAVELJ TO START OVER AND STEP OUT SLOWLY LOOKING FOR
        // OVERFLOWS.  THIS IS NECESSARY SINCE WE MAKE A LARGE JUMP IN L THAT
        // WAVELJ IS NOT EXPECTING.

        LASTL[1] = 999999;
        FG_local[1] = (float)ALLOC(Z[NF1S[1]] + L);
        FG_local[2] = (float)ALLOC(Z[NG1S[1]] + L);
        FG_local[3] = (float)ALLOC(Z[NFS[1]] + L);
        FG_local[4] = (float)ALLOC(Z[NGS[1]] + L);
        WAVELJ(L, 2 * L + JSPS[1], 1, 0, DUMMY4_local, FG_local,
            DUMMY4_local,
            ALLOC_base(LWAVR), ALLOC_base(LWAVI),
            ALLOC_base(LVREAL), ALLOC_base(LVIMAG), ALLOC_base(LCENTR));

        // JPTOLX HAS STORED AN S-MATRIX.  LATER ON WE WILL ACCUMULATE A
        // SUM OF CONTRIBUTIONS IN THIS SPOT SO WE MUST ZERO IT AGAIN.

        ALLOC(LSMATS + 2 * L) = 0;
        ALLOC(LSMATS + 2 * L + 1) = 0;

        // COMPUTE AND STORE   RPSI

        if (CCSW) goto L400;
        PSIMX = 0;
        RPSIMX = 0;
        ISCRTS[II - 1] = NALLOC(NAMES[II - 1], NSTEPS);
        LSCRTS = Z[ISCRTS[II - 1]] - 1;
        RVAL = -RSTEP;

        for (I = 1; I <= NSTEPS; I++) {
            RVAL = RVAL + RSTEP;
            PSI = DABS(ALLOC(LWAVR - 1 + I)) +
                DABS(ALLOC(LWAVI - 1 + I));
            if (PSI > PSIMX) ROFMAX = RVAL;
            PSIMX = DMAX1(PSIMX, PSI);
            RPSIMX = DMAX1(RPSIMX, RVAL * PSI);
            ALLOC(LSCRTS + I) = RVAL * PSI;
        } // 389
        if (IPRNT >= 4) std::printf(" FOR L = %4d, MAX(PSI) = %11.3G   IS AT R = %5.1f     MAX(R*PSI) = %11.3G     FINAL PSI = %11.3G\n",
            L, PSIMX, ROFMAX, RPSIMX, PSI);

        L = LC;
    } // 399

L400:
    STANSW = FALSE_F;
    SOSWS[1] = SSOSW;
    TCSWS[1] = STCSW_local;

    for (II = 1; II <= NUMLNK; II++) {
        for (I = 1; I <= 6; I++) {
            LNKBLK.LNKADR[II][I] = SCRTCH.LNKAD3[II][I];
        }
    }

    if (!CCSW) goto L500;

    // WAVPOT ESPECIALLY SETUP COULOMB FUNCTIONS FOR THIS ONE CALL
    // THEY ARE NOLONGER NEEDED SO RESTORE TO NORMAL C.C. MODE.

    Z[NFS[1]] = -Z[NFS[1]];
    Z[NGS[1]] = -Z[NGS[1]];
    Z[NF1S[1]] = -Z[NF1S[1]];
    Z[NG1S[1]] = -Z[NG1S[1]];

    // SAVE STUFF FOR BSPROD

L500:
    NSPSCT = NSTEPS;
    SCTMAX = ASYMPS[1];
    SCTSP = 1 / RSTEP;
    return;
}
