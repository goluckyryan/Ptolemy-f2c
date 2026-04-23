// source_misc.cpp — Miscellaneous utilities, printing, setup, sort
// Translated from source.f: AITLAG, AKFCR, ALSORT, CUBMAP, CUPAB, CUPSPN,
//   DUMMY1/2, EPSLON, GRDSET, INIT8+entries, INTRCF, IXSORT, LINLSQ,
//   MEBDEF, MEBROT, MEBVIB, MUELCO, PARPRT, PHSPRT, PRBPRT, PRTCHI, PRTDER,
//   RCASYM, REQUIV, RTXLNX, SAVEHS+entries, SETINT, SETLOG, TCOEF,
//   TEMP1/2/3, VCSQ12+SETVSQ
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
#include <complex>
#include <vector>
#include <algorithm>

double AITLAG(double X, double STPINV, double* TABLE, int NTABLE, int NAIT)
{
    // Aitken-Lagrange interpolation for DWBA
    // TABLE(I+1) = F(I/STPINV), TABLE(1) = F(0)
    // NAIT = order of interpolation (NAIT+1 points used)

    static const double INVERS[17] = { 0,
        1.0, 0.5, 1.0/3.0, 0.25, 0.2,
        1.0/6.0, 1.0/7.0, 0.125, 1.0/9.0, 0.1,
        1.0/11.0, 1.0/12.0, 1.0/13.0, 1.0/14.0, 1.0/15.0, 0.0625
    };

    double FS[17], DELS[17];
    double F, DEL1, XBYH;
    int ITABLE, NAIT2, NSTART, NEND;

    if (X < 0) {
        std::printf("\n****** AITLAG - X < 0:%15.5G\n", X);
        return TABLE[1];
    }

    XBYH = X * STPINV;
    ITABLE = (int)XBYH;

    // Allow going one stepsize beyond table
    if (ITABLE >= NTABLE) {
        std::printf("\n***** AITLAG - X TOO LARGE:%15.5G%10d%10d\n", X, ITABLE, NTABLE);
        return TABLE[NTABLE];
    }

    // Get range of table to use
    NAIT2 = NAIT / 2;
    NSTART = ITABLE - NAIT2;
    if (NSTART < 0) NSTART = 0;
    NEND = NSTART + NAIT;
    if (NEND >= NTABLE) NEND = NTABLE - 1;
    NSTART = NEND - NAIT;

    // Get first value
    FS[1] = TABLE[NSTART + 1];
    DEL1 = XBYH - NSTART;
    DELS[1] = DEL1;

    // Aitken loop
    for (int I = 1; I <= NAIT; I++) {
        F = TABLE[NSTART + 1 + I];
        DEL1 = DEL1 - 1;
        for (int J = 1; J <= I; J++) {
            F = (F * DELS[J] - FS[J] * DEL1) * INVERS[I + 1 - J];
        }
        DELS[I + 1] = DEL1;
        FS[I + 1] = F;
    }

    return F;
}

void AKFCR(int L1, int J1, int L2, int J2, int LXX2, int J, double& AKREST)
{ AKREST = 0.0; /* translate from source.f L143-165 */ }

void ALSORT(char8* A, int* IX, int N)
{
    // Simple insertion sort on CHARACTER*8 names, permuting IX
    for (int I = 2; I <= N; I++) {
        char8 KEY = A[IX[I]];
        int JVAL = IX[I];
        int J = I - 1;
        while (J >= 1 && std::memcmp(A[IX[J]].data, KEY.data, 8) > 0) {
            IX[J + 1] = IX[J];
            J--;
        }
        IX[J + 1] = JVAL;
    }
}

void CUBMAP(int MAPTYP, double XLO, double XMID_in, double XHI, double GAMMA,
            double* ARGS, double* WTS, int NPTS)
{
    // Cubic-sinh mapping of [-1,1] into [XLO,XHI]
    // MAPTYP: 0=linear, 1=cubic-sinh, 2=rational-sinh, 3=linear-sinh
    // Arrays ARGS, WTS are 1-based

    double TAU, XLEN, XADD, A, B, C, D, TU, XI, SH, DENOM;
    double XMID = XMID_in;

    // Compute arcsinh(GAMMA)
    if (GAMMA > 1.0e-6)
        TAU = std::log(GAMMA + std::sqrt(GAMMA*GAMMA + 1.0));
    else
        TAU = GAMMA * (1.0 - GAMMA*GAMMA / 6.0);

    GAUSSL(NPTS, ARGS, WTS);
    XLEN = XHI - XLO;
    XADD = XLO + XHI;

    switch (MAPTYP) {
    case 0: // Linear map
        for (int I = 1; I <= NPTS; I++) {
            ARGS[I] = XLO + 0.5*XLEN * (ARGS[I] + 1.0);
            WTS[I] = 0.5*XLEN * WTS[I];
        }
        return;

    case 1: // Cubic-sinh mapping
        XMID = std::max(XMID, XLO + XLEN/7.0);
        XMID = std::min(XMID, 0.5*XADD);
        A = 0.5*XADD - XMID;
        B = 0.5*XLEN;
        C = 0.5*XADD;
        for (int I = 1; I <= NPTS; I++) {
            TU = TAU * ARGS[I];
            XI = std::sinh(TU) / GAMMA;
            ARGS[I] = A*(XI*XI - 1.0)*(XI + 1.0) + B*XI + C;
            WTS[I] = WTS[I] * (TAU/GAMMA) * std::cosh(TU)
                   * ((3.0*XI - 1.0)*(XI + 1.0)*A + B);
        }
        return;

    case 2: // Rational-sinh mapping
        A = -XMID * XLEN;
        B = XLEN;
        C = XMID*XADD - 2.0*XLO*XHI;
        D = XADD - 2.0*XMID;
        for (int I = 1; I <= NPTS; I++) {
            TU = TAU * ARGS[I];
            SH = std::sinh(TU);
            DENOM = B - (D/GAMMA)*SH;
            ARGS[I] = (-A + (C/GAMMA)*SH) / DENOM;
            WTS[I] = WTS[I] * (TAU/GAMMA) * std::cosh(TU) * (B*C - A*D) / (DENOM*DENOM);
        }
        return;

    case 3: // Linear-sinh mapping
        for (int I = 1; I <= NPTS; I++) {
            TU = TAU * ARGS[I];
            ARGS[I] = XLO + 0.5*XLEN * (std::sinh(TU)/GAMMA + 1.0);
            WTS[I] = 0.5*XLEN*TAU/GAMMA * std::cosh(TU) * WTS[I];
        }
        return;
    }
}

void CUPAB(int IORDER, int IPT, int* LXX2S, int LXX2,
           int JIN, int JOUT, int JBIG, int LX, double& COUP)
{ COUP = 0.0; /* translate from source.f L13083-13177 */ }

void CUPSPN(int& IRTN, int LXX2, int IPT, double* SCHN,
            int JIN, int JOUT, int JBIGIN, int JBIGOUT, double& COUP)
{ IRTN = 0; COUP = 0.0; /* translate from source.f L13179-13292 */ }

void DUMMY1() { }
void DUMMY2() { }

// EPSLON — Wynn epsilon algorithm for complex sequences.
// XIN(2,N): input partial sums as real*8 pairs (real, imag).
// FRET(2): output accelerated estimate (real, imag).
// Translated from source.f lines 15269-15464.
void EPSLON(double* XIN, int N, double* FRET, double& DEPS, double& DERR, int& IER)
{
    using C16 = std::complex<double>;
    const double& SMLNUM = CNSTNT.SMLNUM;
    const double BIG = 1.0e38;

    if (N > 2000) {
        std::fprintf(stderr, "EPSLON: N > 2000 NOT ALLOWED %d\n", N);
        std::exit(9876);
    }

    std::vector<C16> X(N + 1);
    for (int I = 1; I <= N; I++)
        X[I] = C16(XIN[2*(I-1)], XIN[2*(I-1)+1]);

    auto APXABS = [](C16 z) { return std::abs(std::real(z)) + std::abs(std::imag(z)); };
    auto NONZER = [&](C16 z) { return APXABS(z) != 0.0; };
    auto ZERO_F = [&](C16 z) { return APXABS(z) == 0.0; };

    IER = -1;
    if (N < 5) {
        if (N > 0) { FRET[0] = std::real(X[N]); FRET[1] = std::imag(X[N]); }
        return;
    }

    double ACC = std::max(1.0e-8, DEPS * DEPS);
    int n = N;
    C16 FIN, T, W1, W2, W3, W4, W5, W6, W7;
    int IMIN = 4, IAUS = 0;
    bool ISW1 = false, ISW2 = false;

L_OUTER:
    if (n <= 0) return;
    FIN = X[n];
    FRET[0] = std::real(FIN);
    FRET[1] = std::imag(FIN);
    if (n == 1) return;
    DERR = 0;
    for (int I = 1; I <= n; I++) {
        double d = APXABS(FIN - X[I]);
        if (d > DERR) DERR = d;
    }
    DERR = DERR / (APXABS(FIN) + SMLNUM);
    IER = 0;
    if (DERR <= ACC) return;
    ACC = DEPS;
    if (n < 6) { IER = 1; return; }

    // L2: initial epsilon values
    ISW1 = false; ISW2 = false;
    W1 = C16(BIG, 0);
    W7 = X[4] - X[3];
    if (NONZER(W7)) W1 = C16(1.0, 0) / W7;
    W5 = C16(BIG, 0);
    W7 = X[2] - X[1];
    if (NONZER(W7)) W5 = C16(1.0, 0) / W7;
    W4 = X[3] - X[2];
    if (NONZER(W4)) goto L_STEP9;
    W4 = C16(BIG, 0);
    T  = X[2];
    W2 = X[3];
    W3 = C16(BIG, 0); // L8
    goto L17;

L_STEP9:
    W4 = C16(1.0, 0) / W4;
    T  = C16(BIG, 0);
    W7 = W4 - W5;
    if (NONZER(W7)) T = X[2] + C16(1.0, 0) / W7;
    W2 = W1 - W4;
    if (NONZER(W2)) goto L15;
    W2 = C16(BIG, 0);
    ISW2 = (std::real(T) != BIG);
    W3 = W4;
    goto L17;

L15:
    W2 = X[3] + C16(1.0, 0) / W2;
    W7 = W2 - T;
    if (ZERO_F(W7)) { W3 = C16(BIG, 0); goto L17; }
    W3 = W4 + C16(1.0, 0) / W7;

L17:
    ISW1 = ISW2;
    ISW2 = false;
    IMIN = 4;

    // DO 40 I = 5, n
    for (int I = 5; I <= n; I++) {
        IAUS = I - IMIN;
        W4 = C16(BIG, 0);
        W5 = X[I-1];
        W7 = X[I] - X[I-1];
        if (ZERO_F(W7)) goto L_24;
        W4 = C16(1.0, 0) / W7;
        if (std::real(W1) == BIG) goto L25;
        W6 = W4 - W1;
        if (APXABS(W6) > 1.0e-12 * APXABS(W4)) goto L22;
        ISW2 = true;
        if (NONZER(W6)) goto L22;
        W5 = C16(BIG, 0);
        W6 = W1;
        if (std::real(W2) != BIG) goto L28;
        goto L26;
    L22:
        W5 = X[I-1] + C16(1.0, 0) / W6;
        if (APXABS(W5) >= 1.0e-10 * APXABS(X[I-1])) goto L_24;
        if (NONZER(W5)) goto L36_OUT;
    L_24: // label 24
        W7 = W5 - W2;
        if (NONZER(W7)) goto L27;
    L25:
        W6 = C16(BIG, 0);
    L26:
        ISW2 = false;
        X[IAUS] = W2;
        goto L39;
    L27:
        W6 = W1 + C16(1.0, 0) / W7;
    L28:
        if (!ISW1) goto L33;
        if (std::real(W2) == BIG) goto L32;
        W7 = W5/(W2-W5) + T/(W2-T) + X[I-2]/(X[I-2]-W2);
        if (ZERO_F(W7 + C16(1.0, 0))) goto L38;
        X[IAUS] = W7 * W2 / (C16(1.0, 0) + W7);
        goto L39;
    L32:
        X[IAUS] = W5 + T - X[I-2];
        goto L39;
    L33:
        W7 = W6 - W3;
        if (ZERO_F(W7)) goto L38;
        X[IAUS] = W2 + C16(1.0, 0) / W7;
        if (APXABS(X[IAUS]) >= 1.0e-10 * APXABS(W2)) goto L37;
        if (NONZER(X[IAUS])) goto L36_OUT;
    L37:
        if (std::real(W2) != BIG) goto L39;
    L38:
        X[IAUS] = C16(BIG, 0);
        IMIN = I;
    L39:
        W1 = W4;
        T  = W2;
        W2 = W5;
        W3 = W6;
        ISW1 = ISW2;
        ISW2 = false;
    } // end DO 40

    n -= IMIN;
    goto L_OUTER;

L36_OUT: // label 36
    n = IAUS - 1;
    goto L_OUTER;
}

// GRDSET: moved to grdset_translated.cpp

// INIT8 + entries: initialize allocator objects
void INIT8(int I, double X8)
{
    int LOC = LOCPTRS.Z[I];
    int LEN = LENGTH.LENG[I];
    double* ptr = &ALLOC(LOC);
    for (int k = 0; k < LEN; k++) ptr[k] = X8;
}

void INIT4(int I, float x4)
{
    int LOC = LOCPTRS.Z[I];
    int LEN = LENGTH.LENG[I];
    float* ptr = &ALLOC4(ALLOCS.FACFR4 * (LOC - 1) + 1);
    int n4 = LEN * (8 / sizeof(float));
    for (int k = 0; k < n4; k++) ptr[k] = x4;
}

void INIT4i(int I, int i4)
{
    int LOC = LOCPTRS.Z[I];
    int LEN = LENGTH.LENG[I];
    int* ptr = &ILLOC(ALLOCS.FACFR4 * (LOC - 1) + 1);
    int n4 = LEN * (8 / sizeof(int));
    for (int k = 0; k < n4; k++) ptr[k] = i4;
}

void INIT2(int I, int16_t x2)
{
    int LOC = LOCPTRS.Z[I];
    int LEN = LENGTH.LENG[I];
    int16_t* ptr = &ILLOC2(ALLOCS.FACFR4 * (LOC - 1) + 1);
    int n2 = LEN * (8 / sizeof(int16_t));
    for (int k = 0; k < n2; k++) ptr[k] = x2;
}

void INIT1(int I, int8_t x1)
{
    int LOC = LOCPTRS.Z[I];
    int LEN = LENGTH.LENG[I];
    char* ptr = &CHLLOC(ALLOCS.FACFR4 * (LOC - 1) + 1);
    int n1 = LEN * 8;
    for (int k = 0; k < n1; k++) ptr[k] = x1;
}

void INTRCF(int NCHN, int JP, int JT, int LX, int LDEL,
            int NUMLIS, int* LIS, int LSTEP, int LXMAX,
            double* SR, double* SI, int LOFF1, int ISDIM1, int LISW,
            int LSTRT, int LMAX_arg, double* SOUT, float* SOUT4,
            int LOFF2, int ISDIM2, int LBASE2, int FOURSW,
            double& BIGBOY, int& LPEAK, double* CMPLXL, double* CMPLXS)
{
    // Uses continued fraction to interpolate S-matrices
    // Translated from source.f L22815-22955

    double RE, HIM, F, G, SIZE;
    complex16 CFLI, CTEMP;

    BIGBOY = 0;
    LPEAK = LSTRT;

    // Skip interpolation if LSTEP = 1
    if (LSTEP == 1) goto L800;

    {
        // Fit continued fraction to S-matrices
        int ICOUNT = 0;
        for (int ILI = 1; ILI <= NUMLIS; ILI++) {
            int LI = LIS[ILI];  // 1-based
            if (LI < LXMAX) continue;
            ICOUNT++;
            // CMPLXL and CMPLXS are REAL*8(2,NUMLIS) treated as complex
            CMPLXL[2*(ICOUNT-1)]   = (double)LI;
            CMPLXL[2*(ICOUNT-1)+1] = 0.0;
            int II = ILI;
            if (LISW) II = LI + 1;
            // SR(LOFF1, II) → SR[(LOFF1-1) + (II-1)*ISDIM1]
            CMPLXS[2*(ICOUNT-1)]   = SR[(LOFF1-1) + (II-1)*ISDIM1];
            CMPLXS[2*(ICOUNT-1)+1] = SI[(LOFF1-1) + (II-1)*ISDIM1];
        }

        // CCNFRC: complex continued fraction setup
        CCNFRC(ICOUNT, reinterpret_cast<complex16*>(CMPLXL),
                        reinterpret_cast<complex16*>(CMPLXS));

        // Compute S for LSTRT to LMAX using continued fraction
        int ILI = 1;
        for (int LI = LSTRT; LI <= LMAX_arg; LI++) {
            CFLI = complex16((double)LI, 0.0);
            CCONTF(ICOUNT, reinterpret_cast<complex16*>(CMPLXL),
                            reinterpret_cast<complex16*>(CMPLXS),
                            CFLI, CTEMP);
            RE = std::real(CTEMP);
            HIM = std::imag(CTEMP);
            SIZE = RE*RE + HIM*HIM;

            // Check against exact values when available
L500:       if (LIS[ILI] - LI < 0) goto L510;
            if (LIS[ILI] - LI == 0) goto L550;
            goto L780;
L510:       if (ILI == NUMLIS) goto L550;
            ILI++;
            goto L500;
L550:       {
                int II = ILI;
                if (LISW) II = LI + 1;
                F = SR[(LOFF1-1) + (II-1)*ISDIM1];
                G = SI[(LOFF1-1) + (II-1)*ISDIM1];
                SIZE = F*F + G*G;
                ILI++;
                if (LI < LXMAX) goto L775;
                if ((F-RE)*(F-RE) + (G-HIM)*(G-HIM) > 1.0e-6 * SIZE)
                    std::printf(" **** WARNING: FOR CHAN, 2*JP, 2*JT, LX, LO-LI, LI =%5d%5d%5d%5d%5d%5d\n"
                                "     CONT. FRAC. =%13.5G%13.5G     EXACT =%13.5G%13.5G\n",
                                NCHN, JP, JT, LX, LDEL, LI, RE, HIM, F, G);
            }
L775:       RE = F;
            HIM = G;
L780:       if (FOURSW) goto L785;
            // SOUT(1, LOFF2, LI-LBASE2) → SOUT[0 + 2*(LOFF2-1) + 2*ISDIM2*(LI-LBASE2-1)]
            SOUT[0 + 2*(LOFF2-1) + 2*ISDIM2*(LI-LBASE2-1)] = RE;
            SOUT[1 + 2*(LOFF2-1) + 2*ISDIM2*(LI-LBASE2-1)] = HIM;
            goto L790;
L785:       // SOUT4(1, LOFF2, LI-LBASE2)
            SOUT4[0 + 2*(LOFF2-1) + 2*ISDIM2*(LI-LBASE2-1)] = (float)RE;
            SOUT4[1 + 2*(LOFF2-1) + 2*ISDIM2*(LI-LBASE2-1)] = (float)HIM;
L790:       if (SIZE < BIGBOY) continue;
            BIGBOY = SIZE;
            LPEAK = LI;
        }
    }
    goto L850;

    // No interpolation (LSTEP = 1), just transfer and scan
L800:
    {
        int LBASE1 = LIS[1] - 1;
        if (LISW) LBASE1 = -1;
        for (int LI = LSTRT; LI <= LMAX_arg; LI++) {
            F = SR[(LOFF1-1) + (LI-LBASE1-1)*ISDIM1];
            G = SI[(LOFF1-1) + (LI-LBASE1-1)*ISDIM1];
            if (FOURSW) goto L835;
            SOUT[0 + 2*(LOFF2-1) + 2*ISDIM2*(LI-LBASE2-1)] = F;
            SOUT[1 + 2*(LOFF2-1) + 2*ISDIM2*(LI-LBASE2-1)] = G;
            goto L840;
L835:       SOUT4[0 + 2*(LOFF2-1) + 2*ISDIM2*(LI-LBASE2-1)] = (float)F;
            SOUT4[1 + 2*(LOFF2-1) + 2*ISDIM2*(LI-LBASE2-1)] = (float)G;
L840:       SIZE = F*F + G*G;
            if (SIZE < BIGBOY) continue;
            BIGBOY = SIZE;
            LPEAK = LI;
        }
    }

L850:
    return;
}

void IXSORT(int* IX, double* X, int N)
{
    // Simple insertion sort on X, permuting IX
    for (int I = 2; I <= N; I++) {
        double KEY = X[IX[I]];
        int JVAL = IX[I];
        int J = I - 1;
        while (J >= 1 && X[IX[J]] > KEY) {
            IX[J + 1] = IX[J];
            J--;
        }
        IX[J + 1] = JVAL;
    }
}

void LINLSQ(int IFUNC, int N, double* XVALS, double* SVALS, double& CVAL,
             double& AVAL, double& B, double& CHI, int PBUGSW)
{ CHI = 0.0; /* translate from source.f L23407-23509 */ }

void MEBDEF(int IORD, int I12, int IBN, double* XNS, int IBC, double* XCS,
            int MEORD, int LXX2, int* LXX2S, int IZ,
            double RNUC, double RC, double DELTAE, int J2, int J1, int K2, int K1,
            double& RMENUC, double& RMECOU, int& IRTN, int IPT)
{ RMENUC = 0.0; RMECOU = 0.0; IRTN = 0; /* translate from source.f L26066-26335 */ }

void MEBROT(int IORD, int I12, int IBN, double* XNS, int IBC, double* XCS,
            int MEORD, int LXX2, int* LXX2S, int IZ,
            double RNUC, double RC, double DELTAE, int J2, int J1, int K2, int K1,
            double& RMENUC, double& RMECOU, int& IRTN)
{ RMENUC = 0.0; RMECOU = 0.0; IRTN = 0; /* translate from source.f L26337-26567 */ }

void MEBVIB(int IORD, int I12, int IBN, double* XNS, int IBC, double* XCS,
            int MEORD, int LXX2, int* LXX2S, int IZ,
            double RNUC, double RC, double DELTAE, int J2, int J1, int K2, int K1,
            double& RMENUC, double& RMECOU, int& IRTN)
{ RMENUC = 0.0; RMECOU = 0.0; IRTN = 0; /* translate from source.f L26569-26889 */ }

// MUELCO → anapow_translated.cpp

// PARPRT → parprt_translated.cpp

// PHSPRT → phsprt_translated.cpp

// PRBPRT: moved to prbprt_translated.cpp

void PRTCHI(int& IRET)
{ IRET = 0; /* translate from source.f L29731-30040 */ }

void PRTDER(int& IRET)
{ IRET = 0; /* translate from source.f L30042-30229 */ }

void RCASYM(int L, double ETA, double RHO, int IPRN, double SIGL,
            double* ZETA, double* PHI, double* DZETA,
            double* F, double* FP, double* G, double* GP,
            double* Z, double* DZSQ, double* S, double* ZINV,
            double EPS, int& NMAX, int& NTZ, int& ISIG)
{ NMAX = 0; NTZ = 0; ISIG = 0; /* translate from source.f L30438-30607 */ }

void REQUIV(double& RC, double& RC0, int IZ, int IA)
{ /* translate from source.f L30859-30924 */ }

double RTXLNX(double X)
{ if (X <= 0.0) return 0.0; return std::sqrt(X * std::log(X)); }

void SAVEHS()
{ /* translate from source.f L30960-31176 */ }

void SAVECM()
{ /* ENTRY of SAVEHS */ }

void USEHS()
{ /* ENTRY of SAVEHS */ }

void SETINT(int N)
{
    // Set up integer workspace — minimal implementation
    INTWRK.MAXINT = N;
}

void SETLOG(int N)
{
    // Set up log factorial table
    if (N <= LOGFAC.MAXLF) return;
    int ISTART = LOGFAC.MAXLF + 1;
    LOGFAC.MAXLF = N;
    for (int I = ISTART; I <= N; I++) {
        LOGFAC.LF[I + 1] = LOGFAC.LF[I] + std::log((double)I);
    }
}

double TCOEF(int ITYP, int LAS, int L, int JP, int JSP)
{ /* translate from source.f L34242-34323 */ return 1.0; }

// Alias for translated files using TCOEF_5
double TCOEF_5(int ITYP, int LAS, int L, int JP, int JSP)
{ return TCOEF(ITYP, LAS, L, JP, JSP); }

void TEMP1()
{ /* translate from source.f L34341-34573 */ }

void TEMP2() { }

void TEMP3() { }

// File-scope arrays shared between SETVSQ and VCSQ12
double VCSQ12_R1[4] = {}, VCSQ12_R2[4] = {}, VCSQ12_VC0[4] = {};
double VCSQ12_A[4] = {}, VCSQ12_B[4] = {}, VCSQ12_C[4] = {};
double VCSQ12_D[4] = {}, VCSQ12_E[4] = {}, VCSQ12_F[4] = {};
double VCSQ12_X[4] = {}, VCSQ12_Y[4] = {};

void VCSQ12(double RVAL, double& X, int K)
{
    // Coulomb potential between two uniform spheres (Poling et al. Phys.Rev.C13,648,1976)
    // Uses coefficients precomputed by SETVSQ for sphere pair K (1-indexed).
    // Three regions: outside (point-charge), fully inside (quadratic), partial overlap (polynomial/R).
    double R1 = VCSQ12_R1[K];  // larger radius (enforced by SETVSQ)
    double R2 = VCSQ12_R2[K];  // smaller radius
    if (RVAL >= R1 + R2) {
        // Outside both spheres: point-charge Coulomb
        X = VCSQ12_VC0[K] / RVAL;
    } else if (RVAL <= R1 - R2) {
        // Fully inside larger sphere: quadratic formula
        X = VCSQ12_X[K] + VCSQ12_Y[K] * RVAL * RVAL;
    } else {
        // Partial overlap: polynomial / R
        double R = RVAL;
        X = VCSQ12_A[K] + R*(VCSQ12_B[K] + R*(VCSQ12_C[K] + R*(VCSQ12_D[K]
              + R*(VCSQ12_E[K] + R*R*VCSQ12_F[K]))));
        X = X / R;
    }
    if (K <= 2 && !std::isfinite(X))
        std::fprintf(stderr, "VCSQ12_NAN K=%d RVAL=%.5e R1=%.5e R2=%.5e VC0=%.5e X=%.5e\n",
            K, RVAL, R1, R2, VCSQ12_VC0[K], X);
}

void SETVSQ(double RR1, double RR2, int IZ1, int IZ2, int K)
{
    // ENTRY of VCSQ12 — translate from source.f L34823-34861
    // Computes coefficients for Coulomb potential between two uniform spheres
    VCSQ12_R1[K] = RR1;
    VCSQ12_R2[K] = RR2;
    VCSQ12_VC0[K] = IZ1 * IZ2 * (CNSTNT.HBARC / CNSTNT.AFINE);
    if (VCSQ12_R1[K] < VCSQ12_R2[K]) {
        double TEM = VCSQ12_R2[K]; VCSQ12_R2[K] = VCSQ12_R1[K]; VCSQ12_R1[K] = TEM;
    }
    double R1SQ = VCSQ12_R1[K] * VCSQ12_R1[K];
    double R2SQ = VCSQ12_R2[K] * VCSQ12_R2[K];
    double R13 = R1SQ * VCSQ12_R1[K];
    if (VCSQ12_R2[K] != 0.0) {
        double R14 = R13 * VCSQ12_R1[K];
        double R15 = R14 * VCSQ12_R1[K];
        double R16 = R15 * VCSQ12_R1[K];
        double R23 = R2SQ * VCSQ12_R2[K];
        double R24 = R23 * VCSQ12_R2[K];
        double R25 = R24 * VCSQ12_R2[K];
        double R26 = R25 * VCSQ12_R2[K];
        double R12C = R13 * R23;
        double V_val = VCSQ12_VC0[K] / R12C;
        VCSQ12_A[K] = (V_val * (R16 - 9*R14*R2SQ + 16*R12C - 9*R1SQ*R24 + R26)) / 32;
        VCSQ12_B[K] = -(3*V_val * (R15 - 5*R13*R2SQ - 5*R1SQ*R23 + R25)) / 20;
        VCSQ12_C[K] = (9*V_val * (R14 - 2*R1SQ*R2SQ + R24)) / 32;
        VCSQ12_D[K] = -(V_val * (R13 + R23)) / 4;
        VCSQ12_E[K] = (3*V_val * (R1SQ + R2SQ)) / 32;
        VCSQ12_F[K] = -V_val / 160;
    }
    if (VCSQ12_R1[K] == 0.0) return;
    VCSQ12_X[K] = 0.3 * (5*R1SQ - R2SQ) * (VCSQ12_VC0[K] / R13);
    VCSQ12_Y[K] = -0.5 * (VCSQ12_VC0[K] / R13);
}

double second() { return 0.0; } // timing stub
