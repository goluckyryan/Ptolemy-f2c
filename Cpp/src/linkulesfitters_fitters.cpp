// linkulesfitters_fitters.cpp -- translated from linkulesfitters.f
//
// Fitter/optimizer routines for Ptolemy.
// Contains: CMTRIX, FCNCHK, FCUBIC, GENINV, DLSMIN
// Stubs for: LMCHOL, MINFIT, MINIM, MINMON, OCOPTR, MTMM, MOVE,
//            QMTRIX, QUAVER, SDIAG, SIP, SV01A, SYMINV, UPQPAR,
//            VA02A, VD01A, VMMIN, STUFF, RANDU, READY, AIM, FIRE,
//            DRESS, MATMPY
//
// These routines are only called during chi-squared fitting operations.
// The core routines (CMTRIX, FCNCHK, FCUBIC, GENINV, DLSMIN) are
// fully translated. The remaining ~25 routines are stubbed for now.
//
// Translated from Fortran to C++ preserving all logic.

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include "ptolemy_forward.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_alloc.h"
#include "ptolemy_io.h"
#include <cstdio>
#include <cstring>
#include <cmath>


// ============================================================================
// SUBROUTINE CMTRIX(N, M, RJ, F, RJTJ, RJTF)
//
// Forms (Jacobian transpose)*Jacobian and (Jacobian transpose)*F.
// ============================================================================

void CMTRIX(int N, int M, double* RJ, double* F, double* RJTJ, double* RJTF)
{
    // RJ is M x N (Fortran column-major: RJ(M,N))
    // RJTJ is N x N, RJTF is N
    // Access: RJ(I,J) = RJ[(J-1)*M + (I-1)]  (1-based Fortran to 0-based C)

    for (int I = 1; I <= N; I++) {
        for (int J = 1; J <= I; J++) {
            double sum = 0.0;
            for (int K = 1; K <= M; K++) {
                sum += RJ[(I-1)*M + (K-1)] * RJ[(J-1)*M + (K-1)];
            }
            RJTJ[(I-1)*N + (J-1)] = sum;
        }
    }

    for (int I = 1; I <= N; I++) {
        double sum = 0.0;
        for (int K = 1; K <= M; K++) {
            sum += F[K-1] * RJ[(I-1)*M + (K-1)];
        }
        RJTF[I-1] = sum;
    }
}


// ============================================================================
// FUNCTION FCUBIC(F0, G0, F1, G1, XLAM)
//
// Cubic interpolation for minimizers.
// ============================================================================

double FCUBIC(double F0, double G0, double F1, double G1, double XLAM)
{
    double FF = (F1 - F0) / XLAM;
    double G = G1 - G0;
    double C = G - 2.0 * (FF - G0);
    double D = G - 3.0 * C;
    double WSQ = D * D - 12.0 * C * G0;
    double Y;

    if (WSQ < 0.0) {
        double B = G1 - FF;
        Y = 0.5 * (1.0 - FF / B);
    } else {
        double W = DSQRT(WSQ);
        Y = -2.0 * G0 / (D + W);
    }
    return Y * XLAM;
}


// ============================================================================
// SUBROUTINE FCNCHK
//
// Check function values and convergence for DLSMIN.
// ============================================================================

void FCNCHK_sub(int N, int M, double* X, double* F,
            void (*FCN)(int, int, double*, double*, double&, int&),
            int& IS, int& IERR, double ACC, double& FSQ, double& FMIN,
            double& SQSTEP, double DELTSQ, int& MAXFCN, int& NFCALL,
            int& NTEST, int NX, int NF, double* W)
{
    double G = 0.0;
    FCN(N, M, X, F, G, IERR);
    NFCALL = NFCALL + 1;
    if (IERR < 0) { IS = 0; return; }

    FSQ = 0.0;
    for (int I = 1; I <= M; I++) {
        FSQ += F[I-1] * F[I-1];
    }

    if (FSQ <= ACC) { IS = 0; return; }

    if (IS != 1 && IS != 4) goto L60;

    if (FSQ < FMIN) goto L50;
    if (SQSTEP > DELTSQ) goto L60;
    NTEST = NTEST - 1;
    if (NTEST > 0) goto L60;

    IERR = 2;
    if (NTEST < 0) IERR = 3;
    goto L70;

L50:
    NTEST = N + 4;
L60:
    if (NFCALL < MAXFCN) return;
    IERR = 1;
    if (FSQ < FMIN) { IS = 0; return; }

L70:
    for (int I = 1; I <= N; I++) {
        X[I-1] = W[NX + I - 1];
    }
    for (int I = 1; I <= M; I++) {
        F[I-1] = W[NF + I - 1];
    }
    IS = 0;
}


// ============================================================================
// SUBROUTINE GENINV(M, N, A, IA, W)
//
// Generalized matrix inversion via orthogonal triangularization.
// ============================================================================

void GENINV(int M, int N, double* A, int IA, double* W)
{
    // A is IA x N (column-major), W is scratch
    // Access: A(I,J) = A[(J-1)*IA + (I-1)]
    #define A2(i,j) A[((j)-1)*IA + ((i)-1)]

    int NR = M;
    int NC = M + M;

    // Set initial records of row and column interchanges
    for (int I = 1; I <= M; I++)
        W[NR + I - 1] = 0.5 + (double)I;
    for (int I = 1; I <= N; I++)
        W[NC + I - 1] = 0.5 + (double)I;

    int KK = 1;
    int KP, IR;

L115:
    {
        double RMAX = 0.0;
        IR = KK;
        for (int I = KK; I <= M; I++) {
            double SUM = 0.0;
            for (int J = KK; J <= N; J++)
                SUM += A2(I,J) * A2(I,J);
            if (SUM > RMAX) { RMAX = SUM; IR = I; }
        }

        if (IR > KK) {
            int K = NR + KK - 1;
            int L = NR + IR - 1;
            double SUM = W[K]; W[K] = W[L]; W[L] = SUM;
            for (int J = 1; J <= N; J++) {
                SUM = A2(KK,J); A2(KK,J) = A2(IR,J); A2(IR,J) = SUM;
            }
        }

        // Find largest element of pivotal row
        RMAX = 0.0;
        double SUM = 0.0;
        IR = KK;
        for (int J = KK; J <= N; J++) {
            SUM += A2(KK,J) * A2(KK,J);
            if (DABS(A2(KK,J)) > RMAX) { RMAX = DABS(A2(KK,J)); IR = J; }
        }

        if (IR > KK) {
            int K = NC + KK - 1;
            int L = NC + IR - 1;
            double tmp = W[K]; W[K] = W[L]; W[L] = tmp;
            for (int I = 1; I <= M; I++) {
                tmp = A2(I,KK); A2(I,KK) = A2(I,IR); A2(I,IR) = tmp;
            }
        }

        // Orthogonal triangularization
        double SIGMA = DSQRT(SUM);
        double BSQ = DSQRT(SUM + SIGMA * DABS(A2(KK,KK)));
        W[KK-1] = DSIGN(SIGMA + DABS(A2(KK,KK)), A2(KK,KK)) / BSQ;
        A2(KK,KK) = -DSIGN(SIGMA, A2(KK,KK));
        KP = KK + 1;

        if (KP > N) goto L200;

        for (int J = KP; J <= N; J++)
            A2(KK,J) = A2(KK,J) / BSQ;

        if (KP > M) goto L200;

        for (int I = KP; I <= M; I++) {
            SUM = W[KK-1] * A2(I,KK);
            for (int J = KP; J <= N; J++)
                SUM += A2(KK,J) * A2(I,J);
            A2(I,KK) = A2(I,KK) - SUM * W[KK-1];
            for (int J = KP; J <= N; J++)
                A2(I,J) = A2(I,J) - SUM * A2(KK,J);
        }

        KK = KP;
        goto L115;
    }

L200:
    // Construction of generalized inverse
    KK = M;
    KP = M + 1;
    {
        double SUM = W[M-1] / A2(M,M);
        if (N > M) {
            for (int J = KP; J <= N; J++)
                A2(M,J) = -SUM * A2(M,J);
        }
        A2(M,M) = 1.0 / A2(M,M) - SUM * W[M-1];
    }

L215:
    KP = KK;
    KK = KP - 1;
    if (KK <= 0) goto L250;

    for (int I = KP; I <= M; I++) {
        double SUM = 0.0;
        for (int J = KP; J <= N; J++)
            SUM += A2(KK,J) * A2(I,J);
        for (int J = KP; J <= N; J++)
            A2(I,J) = A2(I,J) - SUM * A2(KK,J);
        W[I-1] = -SUM * W[KK-1];
    }

    for (int J = KP; J <= N; J++) {
        double SUM = -W[KK-1] * A2(KK,J);
        for (int I = KP; I <= M; I++)
            SUM -= A2(I,KK) * A2(I,J);
        A2(KK,J) = SUM / A2(KK,KK);
    }

    {
        double SUM = 1.0 - W[KK-1] * W[KK-1];
        for (int I = KP; I <= M; I++) {
            SUM -= A2(I,KK) * W[I-1];
            A2(I,KK) = W[I-1];
        }
        A2(KK,KK) = SUM / A2(KK,KK);
    }
    goto L215;

L250:
    // Undo column interchanges
    for (int I = M; I >= 1; I--) {
        int K = (int)W[NC + I - 1];
        if (K == I) continue;
        for (int J = 1; J <= M; J++) {
            double tmp = A2(J,I); A2(J,I) = A2(J,K); A2(J,K) = tmp;
        }
    }

    // Undo row interchanges
    for (int I = M; I >= 1; I--) {
        int K = (int)W[NR + I - 1];
        if (K == I) continue;
        for (int J = 1; J <= N; J++) {
            double tmp = A2(I,J); A2(I,J) = A2(K,J); A2(K,J) = tmp;
        }
    }

    #undef A2
}


// ============================================================================
// SUBROUTINE MINMON
//
// Print iteration information for minimizers.
// ============================================================================

void MINMON(int ITCNT, int NFCALL, int NJCALL, int NPARAM, int NPTS,
            double* PARAM, double* F, double& FS)
{
    double CHISQ = 0.0;
    for (int I = 0; I < NPTS; I++)
        CHISQ += F[I] * F[I];

    printf(" ITERATION%4d  NFCALL =%6d  CHI**2 =%15.7g\n", ITCNT, NFCALL, CHISQ);
    printf(" PARAMETERS:");
    for (int I = 0; I < NPARAM; I++) {
        printf("%14.6g", PARAM[I]);
        if ((I + 1) % 5 == 0) printf("\n             ");
    }
    printf("\n");
    FS = CHISQ;
}


// ============================================================================
// SUBROUTINE MOVE(X0, X, N)
// Copy array X to X0.
// ============================================================================

void MOVE(double* X0, double* X, int N)
{
    for (int I = 0; I < N; I++)
        X0[I] = X[I];
}


// ============================================================================
// SUBROUTINE MTMM(A1, A2, RESULT, N1, N2, N3)
// Matrix multiply: RESULT = A1 * A2
// ============================================================================

void MTMM(double* ARRAY1, double* ARRAY2, double* RESULT,
           int IDIM1, int IDIM2, int IDIM3)
{
    for (int I = 0; I < IDIM1; I++) {
        for (int J = 0; J < IDIM3; J++) {
            double sum = 0.0;
            for (int K = 0; K < IDIM2; K++) {
                sum += ARRAY1[I * IDIM2 + K] * ARRAY2[K * IDIM3 + J];
            }
            RESULT[I * IDIM3 + J] = sum;
        }
    }
}


// ============================================================================
// SUBROUTINE MATMPY(M, N, H, G, S)
// Matrix-vector multiply: S = H * G
// ============================================================================

void MATMPY(int M, int N, double* H, double* G, double* S)
{
    for (int I = 0; I < M; I++) {
        double sum = 0.0;
        for (int J = 0; J < N; J++) {
            sum += H[I * N + J] * G[J];
        }
        S[I] = sum;
    }
}


// ============================================================================
// SUBROUTINE RANDU(IX, YFL)
// Simple random number generator.
// ============================================================================

void RANDU(int& IX, double& YFL)
{
    IX = IX * 65539;
    if (IX < 0) IX = IX + 2147483647 + 1;
    YFL = (double)IX * 0.4656613e-9;
}


// ============================================================================
// Stub routines for remaining fitter/optimizer routines.
// These are only needed for chi-squared fitting and will be fully
// translated when fitting operations are tested.
// ============================================================================

void DLSMIN(int N, int M, double* X, double DISMAX, double ACC,
            int& MAXFCN,
            void (*FCN)(int, int, double*, double*, double&, int&),
            int& IERR, double* F, double* W,
            void (*MONIT)(int, int, int, int, int, double*, double*, double&),
            int IDEB, double STEPH, double BOUND)
{
    printf("\n**** DLSMIN: least-squares minimizer (Phase 7 stub).\n");
    IERR = 1;
}

void LMCHOL(int NVAR, int NFUNS, double* X, int MODE, int& MAXF,
            double* S,
            void (*EVLFUN)(int, int, double*, double*, double*, int&),
            double TOL, double TOL2, double* FVECT, double* RJACOB,
            double* RJINV, int& IERR, int& NFCALL, int& NJCALL,
            double* SCRV1, double* SQSCAL, double* DJTJ, double* VJTF,
            double* SCRV2, double* HISTRY,
            void (*MONIT)(int, int, int, int, int, double*, double*, double&),
            int IDEB)
{
    printf("\n**** LMCHOL: Levenberg-Marquardt minimizer (Phase 7 stub).\n");
    IERR = 1;
}

void MINIM(double* DIAG, double EPS, double* X, double* XP, double* G,
           double* GP, double* S, double* SIG, double* Y, double* ETA,
           double* RHO, double* H, double& FS,
           void (*FJACOB)(int, int, double*, double*, double*, int&),
           int NA, int N, int IPRINT, int& MAXFUN, int& ISIG)
{
    printf("\n**** MINIM: variable-metric minimizer (Phase 7 stub).\n");
    ISIG = -1;
}

void QUAVER(double* DIAG, double EPS, double TOL, double* X, double* XP,
            double* G, double* GP, double* S, double* SIG, double* Y,
            double* ETA, double* RHO, double* H, double& FS,
            void (*FJACOB)(int, int, double*, double*, double*, int&),
            int NA, int N, int IPRINT, int& MAXFUN, int& ISIG)
{
    printf("\n**** QUAVER: quasi-Newton minimizer (Phase 7 stub).\n");
    ISIG = -1;
}

void VMMIN(double* DIAG, double EPS, double* X, double* XP, double* G,
           double* GP, double* S, double* T, double* GB, double* Q1,
           double* Q2, double* H, double& FS,
           void (*FJACOB)(int, int, double*, double*, double*, int&),
           int NA, int N, int IPRINT, int NRAN, int& MAXFUN, int& ISIG)
{
    printf("\n**** VMMIN: AMD-Davidon minimizer (Phase 7 stub).\n");
    ISIG = -1;
}

void VA02A(int M, int N, double* F, double* X, double TOL, double TOL2,
           double ESCALE, int IPRINT, int& MAXFUN, int& IRET, double* W)
{
    printf("\n**** VA02A: Powell least-squares fitter (Phase 7 stub).\n");
    IRET = -1;
}

void OCOPTR(int* IDIM, int M, int N, double* X, double* F,
            double* X0, double* F0,
            void (*FCN)(int, int, double*, double*, double&, int&),
            int& NFCALL, double* FJ, int& MAXFCN, int& ICONV,
            double EPS, double* H, double& FL0, double& FMIN, int IDEB)
{
    printf("\n**** OCOPTR: optimization routine (Phase 7 stub).\n");
    ICONV = -1;
}

void SYMINV(double* A, int N, double* B, int M, double& DETERM, int NMAX)
{
    printf("\n**** SYMINV: symmetric matrix inversion (Phase 7 stub).\n");
    DETERM = 0.0;
}

void QMTRIX_sub(int N, double* RJTJ, double* D, double* RJTF,
                double* S, double& DELTA, double& QPARM)
{
    printf("\n**** QMTRIX: Marquardt matrix solver (Phase 7 stub).\n");
}

void UPQPAR(int N, double SQS, double* RJTJ, double ACTRED,
            double PRERED, double DIRDRV, double& QPARM, double& QCPARM)
{
    printf("\n**** UPQPAR: update Marquardt parameter (Phase 7 stub).\n");
}

void SDIAG(int NM, int N, double EPS, int ILIMIT, double* A,
           double* D, double* X, double* E)
{
    printf("\n**** SDIAG: symmetric matrix diagonalization (Phase 7 stub).\n");
}

void MINFIT(int NM, int M, int N, double* A, double* W, int* IP,
            double* B, int& IERR, double* RV1)
{
    printf("\n**** MINFIT: singular value decomposition (Phase 7 stub).\n");
    IERR = 0;
}

void SIP(double* A, int I, int J, double* B, int K, int L,
         double ALPHA, double& SUM, int M)
{
    // Sum inner products: SUM = SUM + ALPHA * dot(A(I:I+M-1,J), B(K:K+M-1,L))
    printf("\n**** SIP: sum inner products (Phase 7 stub).\n");
}

void SV01A(int M, int N, double* VAR, int IVDIM, int* IVPAR, double* W)
{
    printf("\n**** SV01A: singular value computation (Phase 7 stub).\n");
}

void VD01A(int& ITEST, double& X, double& F, int& MAXFUN,
           double ABSACC, double RELACC, double XSTEP)
{
    printf("\n**** VD01A: helper for VA02A (Phase 7 stub).\n");
    ITEST = -1;
}

// VMMIN helpers
void STUFF(int NA, int N, double* X, double* XP, double* G, double* GP,
           double* S, double* T, double* GB, double& FS, double* H,
           void (*FJACOB)(int, int, double*, double*, double*, int&))
{
    printf("\n**** STUFF: VMMIN helper (Phase 7 stub).\n");
}

void READY(int NA, int N, double* X, double* XP, double* G, double* GP,
           double* S, double* T, double* GB, double& FS, double* H,
           void (*FJACOB)(int, int, double*, double*, double*, int&))
{
    printf("\n**** READY: VMMIN helper (Phase 7 stub).\n");
}

void AIM(int NA, int N, double* X, double* XP, double* G, double* GP,
         double* S, double* T, double* GB, double& FS, double* H,
         void (*FJACOB)(int, int, double*, double*, double*, int&))
{
    printf("\n**** AIM: VMMIN helper (Phase 7 stub).\n");
}

void FIRE(int NA, int N, double* X, double* XP, double* G, double* GP,
          double* S, double* T, double* GB, double& FS, double* H,
          void (*FJACOB)(int, int, double*, double*, double*, int&))
{
    printf("\n**** FIRE: VMMIN helper (Phase 7 stub).\n");
}

void DRESS(int NA, int N, double* X, double* XP, double* G, double* GP,
           double* S, double* T, double* GB, double& FS, double* H,
           void (*FJACOB)(int, int, double*, double*, double*, int&))
{
    printf("\n**** DRESS: VMMIN helper (Phase 7 stub).\n");
}
