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
{
    // Construction of kinematical factors W(1,2)
    // See the notes, Eqn 21
    // source.f L143-165
    const double& RT4PI = CNSTNT.RT4PI;

    if (L1 < 0 || L2 < 0) {
        AKREST = 0.0;
        return;
    }
    double RME = CLEBSH(L1, LXX2, 0, 0, L2, 0)
               * SIXJ(J1, L1, J, L2, J2, LXX2)
               * std::pow(-1.0, (J1 + L1 + J) / 2)
               / RT4PI
               * std::sqrt((LXX2 + 1.0) * (L1 + 1.0));
    AKREST = RME;
}

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
           double* RMENUC, double* RMECOU, double* RNUC, int* IZPT,
           double& RMEGN, double& RMEGC, int& IPT1)
{
    // Computes the <JA, JBIGA|| G(LX) ||JB, JBIGB>
    // See notes equations 47, 70, 46, 71
    // source.f L13083-13177
    const double& PI = CNSTNT.PI;
    double* TABLE = FACTRL.FACTBL;  // 1-based: TABLE(I+2) = sqrt((I+1)!)

    // Statement function F (local)
    auto F = [&](int L1, int L2, int L_arg) -> double {
        return std::sqrt((L1 + 1.0) * (L2 + 1.0) / (16.0 * PI * (L_arg + 1.0)))
             * CLEBSH(L1, L2, 0, 0, L_arg, 0);
    };

    int LX = LXX2 / 2;
    double XN = 1.0;
    double XC = 1.0;

    // Check for order of coupling
    IPT1 = IPT;
    if (IORDER == 2) goto L300;

    // First and second order at the same vertex
L200:
    {
        double SGNFAC = +1.0;
        if (IPT == 1) SGNFAC = std::pow(-1.0, LX);
        RMEGN = SGNFAC * RNUC[IPT] * XN * RMENUC[IORDER];
        RMEGC = SGNFAC * 4.0 * PI / ((LXX2 + 1) * 3.0 * IZPT[IPT])
              * XC * RMECOU[IORDER];
        return;
    }

    // Second-order excitation
L300:
    IPT1 = IPT + 2;
    if (IPT > 2) goto L500;

    // Non-mutual second-order excitation (projectile or target)
    // A few extra factors on top of first order for nuclear
    // Also double if LX1 != LX2
    if (LXX2S[1] != LXX2S[2]) XC = 2.0;
    XN = XC * RNUC[IPT] * F(LXX2S[1], LXX2S[2], LXX2);
    goto L200;

    // Second order mutual excitation of target and projectile
L500:
    IPT1 = 5;
    {
        double SGNFAC = std::pow(-1.0, LXX2S[1] / 2);
        RMEGN = SGNFAC * 2.0 * F(LXX2S[1], LXX2S[2], LXX2)
              * RNUC[1] * RNUC[2] * RMENUC[1] * RMENUC[2];

        RMEGC = 0.0;
        if (LXX2 != LXX2S[1] + LXX2S[2]) return;
        RMEGC = SGNFAC * 4.0 * PI * CNSTNT.RT4PI
              * TABLE[LXX2 + 2] / (TABLE[LXX2S[1] + 2] * TABLE[LXX2S[2] + 2])
              / (3.0 * (LXX2 + 1) * IZPT[1] * IZPT[2])
              * RMECOU[1] * RMECOU[2];
    }
}

void CUPSPN(int& IRTN, int LXX2, int IPT, int* SCHN,
            int* LXX2S, int* JNUC, double RMEGN, double RMEGC,
            double& XN, double& XC)
{
    // Calculates channel-spin coupling matrix elements
    // source.f L13179-13292
    // SCHN(2) is INTEGER array (1-based)
    // JNUC(2,2) is INTEGER array in column-major: JNUC(IPT,I)
    //   Fortran JNUC(2,2) stored column-major as flat[0..3]:
    //   JNUC(1,1)=flat[0], JNUC(2,1)=flat[1], JNUC(1,2)=flat[2], JNUC(2,2)=flat[3]
    //   Access: JNUC(i,j) = flat[(j-1)*2 + (i-1)]
    #define JNUC_F(i,j) JNUC[((j)-1)*2 + ((i)-1)]

    double RME;

    // Check for mutual excitation
    if (IPT == 3) goto L500;

    // First and second order at the same vertex
    // L200:
    {
        RME = std::sqrt((SCHN[1] + 1.0) * (SCHN[2] + 1))
            * SIXJ(SCHN[1], SCHN[2], LXX2,
                   JNUC_F(IPT, 2), JNUC_F(IPT, 1), JNUC_F(3 - IPT, 1))
            * std::pow(-1.0, (LXX2 + JNUC_F(1, IPT) + JNUC_F(2, IPT) + SCHN[3 - IPT]) / 2);
    }

    goto L300;

L500:
    // Second order mutual excitation of target and projectile
    {
        // int MUL = LXX2S[1] + LXX2S[2]; // unused
        // IPT1 = 5; // unused locally

        RME = std::sqrt((SCHN[1] + 1) * (SCHN[2] + 1.0) * (LXX2 + 1.0))
            * WIG9J(SCHN[1],    SCHN[2],    LXX2,
                    JNUC_F(1,1), JNUC_F(1,2), LXX2S[1],
                    JNUC_F(2,1), JNUC_F(2,2), LXX2S[2]);
    }

L300:
    XN = RME * RMEGN;
    XC = RME * RMEGC;
    IRTN = 1;

    #undef JNUC_F
}

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
        goto L37;
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
    // Quicksort on X, permuting IX — 1-based arrays
    // source.f L22957-23063
    int LV[17], UV[17]; // 1-based pushdown stack (max depth 16)
    int P, LP, UP, IY;
    double Y;

    // Initialize IX
    for (int I = 1; I <= N; I++)
        IX[I] = I;

    // Initialize the pushdown list
    LV[1] = 1;
    UV[1] = N;
    P = 1;

    // Partition next segment
L10:
    if (P < 1) return;
L20:
    if ((UV[P] - LV[P]) >= 1) goto L30;
    P = P - 1;
    goto L10;

    // Choose bound
L30:
    LP = LV[P] - 1;
    UP = UV[P];
    Y  = X[UP];
    IY = IX[UP];

    // Move lower pointer
L40:
    if ((UP - LP) < 2) goto L70;
    LP = LP + 1;
    if (X[LP] <= Y) goto L40;
    X[UP]  = X[LP];
    IX[UP] = IX[LP];

    // Move upper pointer
L50:
    if ((UP - LP) < 2) goto L60;
    UP = UP - 1;
    if (X[UP] >= Y) goto L50;
    X[LP]  = X[UP];
    IX[LP] = IX[UP];
    goto L40;

    // Finish up
L60:
    UP = UP - 1;
L70:
    X[UP]  = Y;
    IX[UP] = IY;
    if ((UP - LV[P]) >= (UV[P] - UP)) goto L80;
    LV[P+1] = LV[P];
    UV[P+1] = UP - 1;
    LV[P]   = UP + 1;
    goto L90;

L80:
    LV[P+1] = UP + 1;
    UV[P+1] = UV[P];
    UV[P]   = UP - 1;
L90:
    P = P + 1;
    goto L20;
}

void LINLSQ(int IFUNC, int N, double* XVALS, double* SVALS, double& CVAL,
             double& AVAL, double& B, double& CHI, int PBUGSW)
{
    double G1=0, G2=0, G11=0, G12=0, G22=0;
    double SCAL=1, X, S, F2_v, F, WT, DELTA, A, C;

    if (IFUNC==2 || IFUNC==4) SCAL = 0.5*(SVALS[1]+SVALS[2]);
    if (IFUNC==4) SCAL = std::pow(0.5*(XVALS[1]+XVALS[2]), B) / SCAL;

    for (int I=1; I<=N; I++) {
        X = XVALS[I];
        S = SVALS[I];
        switch (IFUNC) {
        case 1: F2_v = -X;             F = DLOG(DABS(S));       break;
        case 2: F2_v = DEXP(B*X);      F = 1.0/S;              break;
        case 3: F2_v = DLOG(DABS(X));   F = DLOG(DABS(S));      break;
        case 4: F2_v = X;               F = S/std::pow(DABS(X), B); break;
        case 5: F2_v = std::pow(DABS(X), B); F = S;             break;
        }
        F = SCAL * F;
        WT = 1.0 / (F*F);
        G1  = G1  + WT*F;
        G2  = G2  + WT*F*F2_v;
        G11 = G11 + WT;
        G12 = G12 + WT*F2_v;
        G22 = G22 + WT*F2_v*F2_v;
    }

    DELTA = G11*G22 - G12*G12;
    A = (G22*G1 - G12*G2) * (1.0/(SCAL*DELTA));
    C = (G11*G2 - G12*G1) * (1.0/(SCAL*DELTA));

    CHI = N - (A*G1 + C*G2) * SCAL;
    AVAL = A;
    CVAL = C;
    if (PBUGSW) std::printf(" LINLSQ:%16.9G%14.6G%14.6G%16.9G%11.2G\n", B, AVAL, CVAL, CHI, SCAL);
}

void MEBDEF(int IORD, int I12, int IBN, double* XNS, int IBC, double* XCS,
            int MEORD, int LXX2, int* LXX2S, int IZ,
            double RNUC, double RC, double DELTAE, int J2, int J1, int K2, int K1,
            double& RMENUC, double& RMECOU, int& IRTN, int IPT)
{
    // Deformed rotational reduced matrix element to and from beta
    // source.f L26066-26335
    const double& PI = CNSTNT.PI;
    auto& Z = LOCPTRS.Z;
    int NOTDEF = *reinterpret_cast<int*>(&INTRNL.NOTDEF);  // Fortran integer overlay

    static const char* BETNAM[3] = { "", "BETAP", "BETAT" }; // 1-based

    double CLEBS[21];  // 1-based: CLEBS[1]..CLEBS[20]
    double XC, XN, AME1, AME2, AME2X, AMEUSE, BETDEF;
    double X2, X3, A, BETA;

    IRTN = 0;

    // Cannot proceed if IORD != 1
    if (IORD != 1) goto L900;

    {
        int LX = LXX2 / 2;
        double TJ1P1 = J1 + 1;

        // Coefficients for RME to beta (rotational model)
        // AME1 converts < ||ALPHA(COULOMB)|| > to < ||E(LX)|| > in sharp cutoff model
        AME1 = 3.0 * IZ * std::pow(RC / 10.0, LX) / (4.0 * PI);

        // AME2 converts beta(effective) to < ||ALPHA(COULOMB)|| >
        AME2 = CLEBSH(J1, LXX2, K1, K2 - K1, J2, K2) * std::sqrt(TJ1P1);
        if (K1 == K2) goto L50;
        if (K2 != 0 && K1 != 0) goto L920;
        AME2 = std::sqrt(2.0) * AME2;
L50:
        if (AME2 == 0.0) goto L900;
        AMEUSE = AME1 * AME2;
        AME2X = AME2;

        if (IBC == -99) goto L190;

        // Get the power series coefficients that relate beta to beta(effective)
        // for the long-range part of the Coulomb interaction
        // BETA(EFFECTIVE) = BETA + X2*BETA**2 + X3*BETA**3
        X3 = 0.0;

        // In the following loop L = 1/2 LAMBDA
        for (int L = 0; L <= LX; L++) {
            CLEBS[L + 1] = CLEBSH(LXX2, LXX2, 0, 0, 4 * L, 0);
            X3 = X3 + std::pow(CLEBS[L + 1], 4) / (4 * L + 1);
        }

        X2 = 0.0;
        if ((LX & 1) == 0)  // .NOT. btest(LX, 0) => LX is even
            X2 = 0.50 * (LX + 2) * std::sqrt((LXX2 + 1) / (4.0 * PI))
               * std::pow(CLEBS[LX / 2 + 1], 2);
        X3 = (LX + 2) * (LX + 1) * (LXX2 + 1.0) * (LXX2 + 1.0) / (24.0 * PI) * X3;

        XC = XCS[I12];

        // Now get the beta Coulomb if necessary
        switch (IBC) {
        case 1: goto L200;
        case 2: goto L110;
        case 3: goto L150;
        case 4: goto L130;
        default: goto L930;
        }

        // <|| ALPHA ||>  get beta effective
L110:
        XC = XC / AME2;
        goto L160;

        // <LX,0||R.M.E.||0>
L130:
        AME2X = 1.0;

        // <J2,K2||R.M.E.||J1,K1>  get beta effective
L150:
        XC = XC / (AME1 * AME2X);

        // Now get beta from beta(effective) by three iterations of Newton-Raphson
        // XC has the desired beta(effective)
L160:
        BETA = XC;
        for (int I = 1; I <= 3; I++) {
            A = 1.0 + BETA * (2.0 * X2 + 3.0 * X3 * BETA);
            BETA = BETA + (XC - BETA * (1.0 + BETA * (X2 + BETA * X3))) / A;
        }
        XC = BETA;
        goto L200;

        // Nothing about Coulomb defined; use zero
L190:
        XC = 0.0;
        XCS[I12] = 0.0;
        RMECOU = 0.0;
        goto L210;

        // BETA is in XC; get everything from it
        // First the R.M.E.
L200:
        XCS[I12] = XC;
        RMECOU = AMEUSE * XC * (1.0 + XC * (X2 + XC * X3));

        // Advance to nuclear first order
L210:
        XN = XNS[I12];
        BETDEF = 0.0;
        if (IBN != -99 && XN == 0.0) goto L300;

        // Locate the definitions
        {
            int L1 = INTGER.IPARM4;
            if (L1 == NOTDEF) L1 = 2;
            int IBETA = NAMLOC(BETNAM[IPT]);
            if (IBETA == 0) goto L910;
            int L2 = LX - L1 + 1;
            if (L2 <= LENGTH.LENG[IBETA] && L2 >= 1)
                BETDEF = ALLOC(Z[IBETA] - 1 + L2);
        }

        if (IBN == -99) goto L260;
        switch (IBN) {
        case 1: goto L300;
        case 2: goto L250;
        default: goto L930;
        }

        // <||ALPHA||>
L250:
        XN = XN / AME2;
        goto L300;

        // Get nuclear from Coulomb if possible
L260:
        if (IBC == -99) goto L270;
        XN = RC / RNUC * XC;
        goto L300;

        // Nuclear beta is not defined; use deformed value if possible
L270:
        XNS[I12] = BETDEF;
        RMENUC = 1.0;
        goto L310;

        // XN is now BETAN or 0 if it could not be found
L300:
        XNS[I12] = XN;
        RMENUC = 1.0;
        if (BETDEF != 0.0) RMENUC = XN / BETDEF;
        if (XN == 0.0) RMENUC = 0.0;
L310:
        RMENUC = AME2 * RMENUC;
        return;
    }

L900:
    RMECOU = 0.0;
    RMENUC = 0.0;
    return;

L910:
    std::printf("\n*** \"%s\" MUST BE DEFINED FOR DEFORMED MODEL\n", BETNAM[IPT]);
    goto L925;

L920:
    std::printf("\n*** PRESENT ROTATIONAL MODEL DOES NOT WORK FOR K1, K2 = %4d/2%4d/2\n", K1, K2);
L925:
    IRTN = 1;
    goto L900;

L930:
    std::printf("\n*** INVALID IBC OR IBN IN MEBDEF*%12d%12d\n", IBC, IBN);
    goto L925;
}

void MEBROT(int IORD, int I12, int IBN, double* XNS, int IBC, double* XCS,
            int MEORD, int LXX2, int* LXX2S, int IZ,
            double RNUC, double RC, double DELTAE, int J2, int J1, int K2, int K1,
            double& RMENUC, double& RMECOU, int& IRTN)
{
    // Converts rotational reduced matrix element to and from beta
    // source.f L26337-26567
    const double& PI = CNSTNT.PI;
    double XC, XN, AME1, AME2, AME2X, AME3, AME4, AMEUSE, SGNFAC;

    IRTN = 0;

    // Cannot proceed if both deformations are not defined
    if (IBN == -99 && IBC == -99) goto L900;

    {
        int LX = LXX2 / 2;
        double TJ1P1 = J1 + 1;

        // Coefficients for RME to beta (rotational model)
        // AME1 converts < ||ALPHA|| > to < ||E(LX)|| > in sharp cutoff model
        AME1 = 3.0 * IZ * std::pow(RC / 10.0, LX) / (4.0 * PI);

        // AME2 converts beta or beta1*beta2 to < ||ALPHA|| >
        AME2 = CLEBSH(J1, LXX2, K1, K2 - K1, J2, K2) * std::sqrt(TJ1P1);
        if (K1 == K2) goto L50;
        if (K2 != 0 && K1 != 0) goto L920;
        AME2 = std::sqrt(2.0) * AME2;
L50:
        if (AME2 == 0.0) goto L900;
        AMEUSE = AME1 * AME2;
        AME2X = AME2;

        if (IORD == 2) goto L500;

        // AME3 is second-order contribution to R.M.E.
        AME3 = 0.0;
        if (MEORD == 1 && IORD == 1)
            AME3 = 0.25 * (LX + 2) * std::sqrt((LXX2 + 1) / PI)
                 * std::pow(CLEBSH(LXX2, LXX2, 0, 0, LXX2, 0), 2);

        // First order
        XC = XCS[I12];

        // Now get the beta Coulomb if necessary
        switch (IBC) {
        case 1: goto L200;
        case 2: goto L110;
        case 3: goto L150;
        case 4: goto L130;
        default: goto L930;
        }

        // <|| ALPHA ||>
L110:
        XC = XC / AME2;
        goto L200;

        // <LX,0||R.M.E.||0>
L130:
        AME2X = 1.0;
        // fall through

        // <J2,K2||R.M.E.||J1,K1>
        // We are solving M.E. = A1*A2*BETA*(1 + A3*BETA)
L150:
        XC = XC / (AME1 * AME2X);
        if (AME3 != 0.0)
            XC = (std::sqrt(1.0 + 4.0 * AME3 * XC) - 1.0) / (2.0 * AME3);

        // BETA is in XC; get everything from it.
        // First the R.M.E. (model dependent)
L200:
        XCS[I12] = XC;
        RMECOU = AMEUSE * XC;

        // Advance to nuclear first order
        XN = XNS[I12];
        if (IBN == -99) goto L260;
        switch (IBN) {
        case 1: goto L300;
        case 2: goto L250;
        default: goto L930;
        }

        // <||ALPHA||>
L250:
        XN = XN / AME2;
        goto L300;

        // Get nuclear from Coulomb
L260:
        XN = RC / RNUC * XC;

        // XN is now BETAN
L300:
        XNS[I12] = XN;
        RMENUC = AME2 * XN;
        return;

        // Second order
        // First Coulomb
L500:
        XC = XCS[2];

        // AME4 converts < ALPHA ALPHA > to < E(LX) >
        AME4 = AME1 * (LX + 2) * std::sqrt((LXX2S[1] + 1.0) * (LXX2S[2] + 1.0)
             / (16.0 * PI * (LXX2 + 1.0)))
             * CLEBSH(LXX2S[1], LXX2S[2], 0, 0, LXX2, 0);
        AME2 = AME2
             * CLEBSH(LXX2S[1], LXX2S[2], K2 - K1, 0, LXX2, K2 - K1);

        switch (IBC) {
        case 1: goto L520;
        case 2: goto L600;
        case 3: goto L540;
        default: goto L930;
        }

        // BETAC1 and BETAC2, get ALPHA X ALPHA
L520:
        XC = AME2 * XCS[1] * XCS[2];
        goto L600;

        // RME E(LX), get ALPHA X ALPHA
L540:
        XC = XC / AME4;

        // Now we have ALPHA X ALPHA, get RME E(LX)
L600:
        RMECOU = AME4 * XC;

        // Go on to nuclear
        if (IBN == -99) goto L760;
        switch (IBN) {
        case 1: goto L720;
        case 2: goto L740;
        default: goto L930;
        }

        // BETAN1, BETAN2
L720:
        RMENUC = AME2 * XNS[1] * XNS[2];
        goto L800;

        // ALPHA X ALPHA given
L740:
        RMENUC = XNS[2];
        goto L800;

        // Unknown, use Coulomb
L760:
        RMENUC = std::pow(RC / RNUC, 2) * XC;

L800:
        return;
    }

L900:
    RMECOU = 0.0;
    RMENUC = 0.0;
    return;

L920:
    std::printf("\n*** PRESENT ROTATIONAL MODEL DOES NOT WORK FOR K1, K2 = %4d/2%4d/2\n", K1, K2);
L925:
    IRTN = 1;
    goto L900;

L930:
    std::printf("\n*** INVALID IBC OR IBN IN MEBROT*%12d%12d\n", IBC, IBN);
    goto L925;
}

void MEBVIB(int IORD, int I12, int IBN, double* XNS, int IBC, double* XCS,
            int MEORD, int LXX2, int* LXX2S, int IZ,
            double RNUC, double RC, double DELTAE, int J2, int J1, int K2IN, int K1IN,
            double& RMENUC, double& RMECOU, int& IRTN)
{
    // Converts vibrational reduced matrix element to and from beta
    // source.f L26569-26889
    const double& PI = CNSTNT.PI;
    double XC, XN, AME1, AME2, AME4, SGNFAC;

    IRTN = 0;

    // Cannot proceed if both deformations are not defined
    if (IBN == -99 && IBC == -99) goto L900;

    {
        int LX = LXX2 / 2;

        // Arrange things so that K2 > K1
        // Fill in undefined numbers of phonons
        int K1 = K1IN;
        int K2 = K2IN;
        if (K1 == 0 && J1 > 0) K1 = 2;
        if (K2 == 0 && J2 > 0) K2 = 2;
        if (K1 >= 0) goto L20;
        if (K2 < 0) goto L940;
        K1 = 2 * IORD;
        if (K2 == 2) K1 = 6 - K1;
        goto L40;
L20:
        if (K2 >= 0) goto L40;
        K2 = 2 * IORD;
        if (K1 == 2) K2 = 6 - K2;

L40:
        int KK1, KK2, JJ1, JJ2;
        if (K2 < K1) goto L45;
        KK1 = K1 / 2;
        KK2 = K2 / 2;
        JJ1 = J1;
        JJ2 = J2;
        SGNFAC = +1.0;
        goto L50;
L45:
        KK1 = K2 / 2;
        KK2 = K1 / 2;
        JJ1 = J2;
        JJ2 = J1;
        SGNFAC = std::pow(-1.0, (J2 - J1) / 2);

L50:
        {
            double TJ2P1 = JJ2 + 1;

            // Coefficients for RME to beta (vibrational model)
            // AME1 converts < ||ALPHA|| > to < ||E(LX)|| > in sharp cutoff model
            AME1 = 3.0 * IZ * std::pow(RC / 10.0, LX) / (4.0 * PI);

            // AME2 converts beta or beta1*beta2 to < ||ALPHA|| >
            // If this is a case we don't know about, AME2 = 0
            AME2 = 0.0;

            if (IORD == 2) goto L400;

            // First order
            if (KK2 >= 3) goto L100;
            if (KK2 != KK1 + 1) goto L100;
            if (KK1 == 0 && JJ1 != 0) goto L100;
            switch (KK2) {
            case 1: goto L70;
            case 2: goto L80;
            default: goto L100;
            }

            // < 1 phonon || LX || 0 >   LX = J2 is already known
L70:
            AME2 = 1.0;
            goto L100;

            // < 2 phonons; J || LX || 1 phonon >
L80:
            AME2 = std::sqrt(TJ2P1 / (LXX2 + 1));
            if (LXX2 != JJ1) goto L100;
            if (JJ2 % 4 != 0) goto L99;
            AME2 = std::sqrt(2.0) * AME2;
            goto L100;

L99:
            AME2 = 0.0;
L100:
            AME2 = SGNFAC * AME2;
            XC = XCS[I12];

            // Now get the beta Coulomb if necessary and possible
            switch (IBC) {
            case 1: goto L170;
            case 2: goto L110;
            case 3: goto L150;
            case 4: goto L130;
            default: goto L930;
            }

            // <|| ALPHA ||>  convert to < ||E(LX)|| >
L110:
            XC = XC * AME1;
            goto L150;

            // <LX||R.M.E.||0>   get beta
L130:
            XC = XC / AME1;
            goto L170;

            // <J2,K2||R.M.E.||J1,K1>  get beta if possible
L150:
            RMECOU = XC;
            XC = 0.0;
            if (AME2 != 0.0) XC = RMECOU / (AME1 * AME2);
            goto L200;

            // BETA is in XC; get everything from it
L170:
            if (AME2 == 0.0) goto L920;
            RMECOU = AME1 * AME2 * XC;

L200:
            XCS[I12] = XC;

            // Advance to nuclear first order
            XN = XNS[I12];
            if (IBN == -99) goto L260;
            switch (IBN) {
            case 1: goto L270;
            case 2: goto L250;
            default: goto L930;
            }

            // <||ALPHA||>
L250:
            RMENUC = XN;
            XN = 0.0;
            if (AME2 != 0.0) XN = RMENUC / AME2;
            goto L300;

            // Get nuclear alpha from Coulomb RME
L260:
            XN = RC / RNUC * RMECOU / AME1;
            goto L250;

            // XN is now BETAN
L270:
            if (AME2 == 0.0) goto L920;
            RMENUC = AME2 * XN;

L300:
            XNS[I12] = XN;
            return;

            // Second order
            // Get conversion factors
L400:
            if (KK2 >= 3) goto L500;
            if (KK1 == 0 && KK2 == 2) goto L410;
            if (KK1 == 1 && KK2 == 1) goto L420;
            goto L500;

            // < 2 phonon || (LX1 LX2)LX || 0 >   LX = J2 is known
L410:
            AME2 = std::sqrt(TJ2P1 / ((LXX2S[1] + 1.0) * (LXX2S[2] + 1.0)));
            if (LXX2S[1] == LXX2S[2]) AME2 = std::sqrt(2.0) * AME2;
            goto L500;

            // < 1 phonon || (LX1 LX2)LX || 1 phonon >
L420:
            if (LXX2S[1] == JJ1 && LXX2S[2] == JJ2) goto L425;
            if (LXX2S[1] == JJ2 && LXX2S[2] == JJ1) goto L425;
            goto L500;

L425:
            AME2 = std::pow(-1.0, LXX2S[1] / 2) * std::sqrt((LXX2 + 1.0)
                 / ((LXX2S[1] + 1.0) * (LXX2S[2] + 1.0)));
            if (JJ2 == JJ1) AME2 = 2.0 * AME2;

L500:
            AME2 = SGNFAC * AME2;
            AME4 = AME1 * (LX + 2) * std::sqrt((LXX2S[1] + 1.0) * (LXX2S[2] + 1.0)
                 / (16.0 * PI * (LXX2 + 1.0)))
                 * CLEBSH(LXX2S[1], LXX2S[2], 0, 0, LXX2, 0);

            // Now do Coulomb
            XC = XCS[2];
            switch (IBC) {
            case 1: goto L520;
            case 2: goto L600;
            case 3: goto L540;
            default: goto L930;
            }

            // BETAC1 and BETAC2, get ALPHA X ALPHA
L520:
            if (AME2 == 0.0) goto L920;
            XC = AME2 * XCS[1] * XCS[2];
            goto L600;

            // RME E(LX), get ALPHA X ALPHA
L540:
            XC = XC / AME4;

            // Now we have ALPHA X ALPHA, get RME E(LX)
L600:
            RMECOU = AME4 * XC;

            // Go on to nuclear
            if (IBN == -99) goto L760;
            switch (IBN) {
            case 1: goto L720;
            case 2: goto L740;
            default: goto L930;
            }

            // BETAN1, BETAN2
L720:
            if (AME2 == 0.0) goto L920;
            RMENUC = AME2 * XNS[1] * XNS[2];
            goto L800;

            // ALPHA X ALPHA given
L740:
            RMENUC = XNS[2];
            goto L800;

            // Unknown, use Coulomb
L760:
            RMENUC = std::pow(RC / RNUC, 2) * XC;

L800:
            return;
        }
    }

L900:
    RMECOU = 0.0;
    RMENUC = 0.0;
    return;

L920:
    std::printf("\n*** PRESENT VIBRATIONAL MODEL DOES NOT WORK"
                " FOR J2, K2, LX1, LX2, LX, J1, K1, ORDER =\n"
                "      %4d/2%4d/2%4d/2%4d/2%4d/2%4d/2%4d/2%4d\n",
                J2, K2IN, LXX2S[1], LXX2S[2], LXX2, J1, K1IN, IORD);
L925:
    IRTN = 1;
    goto L900;

L930:
    std::printf("\n*** INVALID IBC OR IBN IN MEBVIB*%12d%12d\n", IBC, IBN);
    goto L925;

L940:
    std::printf("\n**** CANNOT HAVE INDEFINITE NUMBER OF PHONONS IN BOTH STATES.\n");
    goto L925;
}

// MUELCO → anapow_translated.cpp

// PARPRT → parprt_translated.cpp

// PHSPRT → phsprt_translated.cpp

// PRBPRT: moved to prbprt_translated.cpp

void PRTCHI(int& IRET)
{ IRET = 0; /* translate from source.f L29731-30040 */ }

void PRTDER(int& IRET)
{ IRET = 0; /* translate from source.f L30042-30229 */ }


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

void SETLOG(int MAXI)
{
    static char LOGNAM[9] = "LOGFACS ";
    int ILOGFC = NAMLOC(LOGNAM);
    extern void DUMMY1();
    DUMMY1();
    if (MAXI <= LOGFAC.MAXLF) {
        if (ILOGFC == 0) return;
        goto L300;
    }

    {
        LOGFAC.MAXLF = MAXI + 10;
        ILOGFC = NALLOC(LOGNAM, LOGFAC.MAXLF + 1);
        int LLOGFC = LOCPTRS.Z[ILOGFC];
        ALLOC(LLOGFC) = 0;
        double X = 0;
        double DI = 1;
        for (int I = 1; I <= LOGFAC.MAXLF; I++) {
            X = X + std::log(DI);
            DI = DI + 1;
            ALLOC(LLOGFC + I) = X;
        }
    }

L300:
    {
        int LLOGFC = LOCPTRS.Z[ILOGFC];
        double* alloc_base = &ALLOC(1);
        double* lf_base = &LOGFAC.LF[0];
        LOGFAC.LFBIAS = (int)(alloc_base - lf_base) + LLOGFC - 2;
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
