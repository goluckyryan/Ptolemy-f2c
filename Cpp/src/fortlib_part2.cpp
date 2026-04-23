// fortlib_part2.cpp — second half of fortlib.f translation
// INTRPC, LAGBC, LGRECR, LGROOT, LAGUER, LIN, LINBSV, LINSAV,
// LINTRI, LSQPOL, MATINV, PLMSUB, PLSUB, SPLNCB, SYSERR, YLMSUB

#include "ptolemy_commons.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_forward.h"
#include <cstdio>
#include <cmath>
#include <cstdlib>

// ==========================================================================
// INTRPC — cubic spline interpolation (lines 2312-2419)
// ==========================================================================
void INTRPC(int NCUBIC, double* XCUBES, double* AS, double* BS, double* CS,
            double* DS, int NPTS, double* XS, double* YS)
{
    double INF = 1.0e300;

    double XSIGN = DSIGN(1.0, XCUBES[2] - XCUBES[1]);

    // Fake starting conditions that will cause direct branch to label_200
    double XPREV = INF;
    double XNEXT = INF;
    double XBASE = 0.0;
    double A = 0.0, B = 0.0, C = 0.0, D = 0.0;
    int N = 0;

    for (int I = 1; I <= NPTS; I++) {
        double X = XS[I];
        double XC = XSIGN * X;

    label_100:
        if (XC < XNEXT) goto label_200;

        // Go to the next spline
        XPREV = XNEXT;
        XBASE = XNEXT * XSIGN;
        N = N + 1;
        XNEXT = XCUBES[N + 1] * XSIGN;
        // When we are at the last spline, allow it to go to infinity
        if (N == NCUBIC - 1) XNEXT = INF;

    label_150:
        A = AS[N];
        B = BS[N];
        C = CS[N];
        D = DS[N];
        // See if X is inside this new spline
        goto label_100;

    label_200:
        if (XC >= XPREV) goto label_300;

        // The XS array is not monotonic increasing -- go back to start
        // Also initialize for first time
        N = 1;
        XPREV = -INF;
        XNEXT = XCUBES[2] * XSIGN;
        XBASE = XCUBES[1];
        goto label_150;

    label_300:
        // The set of splines to use has been found, use it
        double DEL = X - XBASE;
        YS[I] = A + DEL * (B + DEL * (C + DEL * D));
    // label_399: continue
    }

    return;
}

// ==========================================================================
// LAGBC — Laguerre recurrence coefficients (lines 2564-2580)
// ==========================================================================
void LAGBC(int NN, double ALF, double* A, double* B, double* C)
{
    for (int N = 1; N <= NN; N++) {
        double EN = N;
        A[N] = 1.0 / EN;
        B[N] = ALF + EN + EN - 1.0;
        C[N] = (ALF + EN - 1.0) / EN;
    }
    return;
}

// ==========================================================================
// LGRECR — Laguerre polynomial recurrence (lines 2539-2563)
// ==========================================================================
void LGRECR(double& PN, double& DPN, double& PN1, double X,
            int NN, double ALF, double* A, double* B, double* C)
{
    double P1 = 1.0;
    double P = ALF + 1.0 - X;
    double DP1 = 0.0;
    double DP = -1.0;

    for (int J = 2; J <= NN; J++) {
        double Q = A[J] * (B[J] - X) * P - C[J] * P1;
        double DQ = A[J] * (B[J] - X) * DP - C[J] * DP1 - A[J] * P;
        P1 = P;
        P = Q;
        DP1 = DP;
        DP = DQ;
    }

    PN = P;
    DPN = DP;
    PN1 = P1;
    return;
}

// ==========================================================================
// LGROOT — Newton iteration for Laguerre root (lines 2515-2538)
// ==========================================================================
void LGROOT(double& X, int NN, double ALF, double& DPN, double& PN1,
            double* A, double* B, double* C, double EPS)
{
    double P, DP, D;

    for (int ITER = 1; ITER <= 10; ITER++) {
        LGRECR(P, DP, PN1, X, NN, ALF, A, B, C);
        D = P / DP;
        X = X - D;
        if (DABS(D / X) < EPS) goto label_3;
    }

label_3:
    DPN = DP;
    return;
}

// ==========================================================================
// LAGUER — Gauss-Laguerre quadrature (lines 2421-2514)
// ==========================================================================
void LAGUER(int NN, double* X, double* W, double ALF, double EPS,
            double& CSX, double& CSW, double& TSX, double& TSW,
            double* A, double* B, double* C)
{
    double FN = NN;
    double XT, R1, R2, RATIO, FI, DPN, PN1;

    LAGBC(NN, ALF, A, B, C);

    CSX = 0.0;
    CSW = 0.0;
    double CC = DGAMMA(ALF + 1.0);
    TSX = FN * (FN + ALF);
    TSW = CC;

    for (int J = 2; J <= NN; J++) {
        CC = CC * C[J];
    }

    for (int I = 1; I <= NN; I++) {
        // Arithmetic IF: IF (I - 2) 2, 4, 5
        if (I - 2 < 0) {
            // label_2: smallest zero
            XT = (1.0 + ALF) * (3.0 + 0.920 * ALF) /
                 (1.0 + 2.40 * FN + 1.80 * ALF);
        } else if (I - 2 == 0) {
            // label_4: second zero
            XT = XT + (15.0 + 6.250 * ALF) /
                 (1.0 + 0.90 * ALF + 2.50 * FN);
        } else {
            // label_5: all other zeros
            FI = I - 2;
            R1 = (1.0 + 2.550 * FI) / (1.90 * FI);
            R2 = 1.260 * FI * ALF / (1.0 + 3.50 * FI);
            RATIO = (R1 + R2) / (1.0 + 0.30 * ALF);
            XT = XT + RATIO * (XT - X[I - 2]);
        }

        // label_6:
        LGROOT(XT, NN, ALF, DPN, PN1, A, B, C, EPS);
        X[I] = XT;
        W[I] = -CC / DPN / PN1;
        CSX = CSX + XT;
        CSW = CSW + W[I];
    }

    return;
}

// ==========================================================================
// LIN — Linear equation solver, no saving (lines 2639-2690)
// ==========================================================================
void LIN(double* A, int NMAX, int NDIM, double* B, int MMAX,
         double POWER, double& DET)
{
    // Validate inputs
    if ((NMAX <= 0) || (MMAX <= 0) || (NDIM < NMAX)) goto label_900;

    if (NMAX > 1) goto label_100;

    // NMAX == 1 case
    for (int M = 1; M <= MMAX; M++) {
        // B(1,M) = B(1,M)/A(1,1)
        B[1 + (M - 1) * NDIM] = B[1 + (M - 1) * NDIM] / A[1 + (1 - 1) * NDIM];
    }
    DET = A[1 + (1 - 1) * NDIM]; // A(1,1)
    return;

label_100:
    {
        int NMAXM1 = NMAX - 1;

        // Forward elimination (triangularize)
        for (int K = 1; K <= NMAXM1; K++) {
            int KP1 = K + 1;
            for (int J = KP1; J <= NMAX; J++) {
                // PIV = -A(J,K)/A(K,K)
                double PIV = -A[J + (K - 1) * NDIM] / A[K + (K - 1) * NDIM];
                for (int M = 1; M <= MMAX; M++) {
                    // B(J,M) = B(J,M) + PIV*B(K,M)
                    B[J + (M - 1) * NDIM] = B[J + (M - 1) * NDIM]
                                           + PIV * B[K + (M - 1) * NDIM];
                }
                for (int I = J; I <= NMAX; I++) {
                    // A(I,J) = A(I,J) + PIV*A(I,K)
                    A[I + (J - 1) * NDIM] = A[I + (J - 1) * NDIM]
                                           + PIV * A[I + (K - 1) * NDIM];
                }
            }
        }

        // Back substitution
        for (int M = 1; M <= MMAX; M++) {
            // B(NMAX,M) = B(NMAX,M)/A(NMAX,NMAX)
            B[NMAX + (M - 1) * NDIM] = B[NMAX + (M - 1) * NDIM]
                                      / A[NMAX + (NMAX - 1) * NDIM];
            for (int II = 1; II <= NMAXM1; II++) {
                int I = NMAX - II;
                double SUM = 0.0;
                int IP1 = I + 1;
                for (int J = IP1; J <= NMAX; J++) {
                    // SUM = SUM + A(J,I)*B(J,M)
                    SUM = SUM + A[J + (I - 1) * NDIM] * B[J + (M - 1) * NDIM];
                }
                // B(I,M) = (B(I,M) - SUM)/A(I,I)
                B[I + (M - 1) * NDIM] = (B[I + (M - 1) * NDIM] - SUM)
                                        / A[I + (I - 1) * NDIM];
            }
        }

        // Determinant computation
        if ((POWER != 1.0) && (POWER != 0.0)) goto label_400;
        DET = A[1 + (1 - 1) * NDIM]; // A(1,1)
        for (int I = 2; I <= NMAX; I++) {
            // DET = DET*A(I,I)
            DET = DET * A[I + (I - 1) * NDIM];
        }
        goto label_500;

    label_400:
        DET = 1.0;
        for (int I = 1; I <= NMAX; I++) {
            // DET = DET*DSIGN(DABS(A(I,I))**POWER, A(I,I))
            double AII = A[I + (I - 1) * NDIM];
            DET = DET * DSIGN(std::pow(DABS(AII), POWER), AII);
        }

    label_500:
        return;
    }

label_900:
    std::printf("\n $$$ LIN/LINSAV -- INVALID INPUT...%10d%10d%10d\n",
                NMAX, NDIM, MMAX);
    SYSERR();
    return;
}

// ==========================================================================
// LINBSV — Linear solver preserving inhomogeneous term (ENTRY in Fortran)
// ==========================================================================
void LINBSV(double* A, int NMAX, int NDIM, double* BSV, int MMAX,
            double* B, double POWER, double& DET)
{
    if ((NMAX <= 0) || (MMAX <= 0) || (NDIM < NMAX)) goto label_900;

    // Copy BSV to B
    for (int M = 1; M <= MMAX; M++) {
        for (int I = 1; I <= NMAX; I++) {
            // B(I,M) = BSV(I,M)
            B[I + (M - 1) * NDIM] = BSV[I + (M - 1) * NDIM];
        }
    }

    LIN(A, NMAX, NDIM, B, MMAX, POWER, DET);
    return;

label_900:
    std::printf("\n $$$ LIN/LINSAV -- INVALID INPUT...%10d%10d%10d\n",
                NMAX, NDIM, MMAX);
    SYSERR();
    return;
}

// ==========================================================================
// LINSAV — Linear solver preserving both matrix and inhomogeneous term
// ==========================================================================
void LINSAV(double* ASV, int NMAX, int NDIM, double* BSV, int MMAX,
            double* A, double* B, double POWER, double& DET)
{
    if ((NMAX <= 0) || (MMAX <= 0) || (NDIM < NMAX)) goto label_900;

    // Copy ASV to A
    for (int I = 1; I <= NMAX; I++) {
        for (int J = 1; J <= NMAX; J++) {
            // A(J,I) = ASV(J,I)
            A[J + (I - 1) * NDIM] = ASV[J + (I - 1) * NDIM];
        }
    }

    // Fall through to LINBSV
    LINBSV(A, NMAX, NDIM, BSV, MMAX, B, POWER, DET);
    return;

label_900:
    std::printf("\n $$$ LIN/LINSAV -- INVALID INPUT...%10d%10d%10d\n",
                NMAX, NDIM, MMAX);
    SYSERR();
    return;
}

// ==========================================================================
// LINTRI — Linear solver with triangular rearrangement
// A(I,J) = A(J,I+1) for I >= J, then copies BSV→B and solves
// ==========================================================================
void LINTRI(double* A, int NMAX, int NDIM, double* BSV, int MMAX,
            double* B, double POWER, double& DET)
{
    if ((NMAX <= 0) || (MMAX <= 0) || (NDIM < NMAX)) goto label_900;

    // Rearrange: A(I,J) = A(J,I+1) for I >= J
    for (int I = 1; I <= NMAX; I++) {
        for (int J = 1; J <= I; J++) {
            // A(I,J) = A(J,I+1)
            A[I + (J - 1) * NDIM] = A[J + ((I + 1) - 1) * NDIM];
        }
    }

    // Go to label 3 in LINSAV, which is the BSV copy + LIN call
    LINBSV(A, NMAX, NDIM, BSV, MMAX, B, POWER, DET);
    return;

label_900:
    std::printf("\n $$$ LIN/LINSAV -- INVALID INPUT...%10d%10d%10d\n",
                NMAX, NDIM, MMAX);
    SYSERR();
    return;
}

// ==========================================================================
// LSQPOL — Least squares polynomial fit (lines 2712-2824)
// ==========================================================================
void LSQPOL(double* X, double* Y, double* W, double* RESID, int NSUB,
            double* SUM, int LSUB, double* A, double* B, int MSUB,
            int NMAX, int MMAX)
{
    // XPOWER overlays F402 common block (100 doubles of scratch space)
    // F402 is reused as raw double workspace here
    double* XPOWER = reinterpret_cast<double*>(&F402) - 1; // 1-based access

    // EQUIVALENCE (K1,DETERM),(K2,I2),(TERM,POLY) — separate variables in C++
    double DETERM;
    double TERM, POLY;
    int K1, K2, I2;

    int N = NSUB;
    int L = LSUB;
    int M = MSUB;
    int M1 = M + 1;
    int M3 = M + M + M;
    int M31 = M3 - 1;
    int M41 = M31 + M;

    // Scale X into (-1,1) to prevent overflow or underflow
    double XMAX = 0.0;
    for (int K = 1; K <= N; K++) {
        XMAX = DMAX1(XMAX, DABS(X[K]));
    }
    double XSCALE = 1.0 / XMAX;
    for (int K = 1; K <= N; K++) {
        X[K] = XSCALE * X[K];
    }

    // Formation and inversion of system of normal equations
    for (K2 = M1; K2 <= M41; K2++) {
        XPOWER[K2] = 0.0;
    }

    for (K1 = 1; K1 <= N; K1++) {
        TERM = W[K1];
        for (K2 = M1; K2 <= M31; K2++) {
            XPOWER[K2] = TERM + XPOWER[K2];
            TERM = X[K1] * TERM;
        }
    }

    for (int I = 1; I <= M; I++) {
        for (int J = 1; J <= M; J++) {
            K2 = I + J + M - 1;
            // A(I,J) = XPOWER(K2)
            A[I + (J - 1) * MMAX] = XPOWER[K2];
        }
    }

    for (int J = 1; J <= L; J++) {
        for (int K = 1; K <= N; K++) {
            // TERM = W(K)*Y(K,J)
            TERM = W[K] * Y[K + (J - 1) * NMAX];
            for (K2 = M3; K2 <= M41; K2++) {
                XPOWER[K2] = TERM + XPOWER[K2];
                TERM = X[K] * TERM;
            }
        }
        for (int I = 1; I <= M; I++) {
            K2 = I + M31;
            // B(I,J) = XPOWER(K2)
            B[I + (J - 1) * MMAX] = XPOWER[K2];
            if (J != L) XPOWER[K2] = 0.0;
        }
    }

    // Call MATINV
    MATINV(A, M, B, L, DETERM, MMAX);

    // Evaluation of residuals
    for (int J = 1; J <= L; J++) {
        SUM[J] = 0.0;
        for (int K = 1; K <= N; K++) {
            POLY = 0.0;
            for (I2 = 1; I2 <= M; I2++) {
                int I = M1 - I2;
                // POLY = X(K)*POLY + B(I,J)
                POLY = X[K] * POLY + B[I + (J - 1) * MMAX];
            }
            // RESID(K,J) = POLY - Y(K,J)
            RESID[K + (J - 1) * NMAX] = POLY - Y[K + (J - 1) * NMAX];
            // SUM(J) = SUM(J) + W(K)*RESID(K,J)**2
            double residval = RESID[K + (J - 1) * NMAX];
            SUM[J] = SUM[J] + W[K] * residval * residval;
        }
    }

    // Un-scale the coefficients and X
    for (int J = 1; J <= L; J++) {
        for (I2 = 2; I2 <= M; I2++) {
            for (int I = I2; I <= M; I++) {
                // B(I,J) = XSCALE*B(I,J)
                B[I + (J - 1) * MMAX] = XSCALE * B[I + (J - 1) * MMAX];
            }
        }
    }
    for (int K = 1; K <= N; K++) {
        X[K] = XMAX * X[K];
    }

    return;
}

// ==========================================================================
// MATINV — Matrix inversion with linear equation solution (lines 2826-2948)
// Gauss-Jordan elimination with full pivot searching
// ==========================================================================
void MATINV(double* A, int N, double* B, int M, double& DETERM, int NMAX)
{
    // PIVOT is REAL*4 in Fortran, INDEX is INTEGER
    // Use F402.PIVOT as double (stores float-like values; works for the algorithm)
    // Use F402.INDEX as int
    double* PIVOT = F402.PIVOT; // 1-based, [1..50]
    int*    INDEX = F402.INDEX; // 1-based, [1..50]

    double DETMAX = 1.0e300;
    double DETMIN = 1.0e-300;

    int IROW = 0, ICOLUM = 0;
    double AMAX, TEMP, SWAP, T;
    int IDET;

    // Initialize determinant and pivot element array
    DETERM = 1.0;
    IDET = 0;
    for (int I = 1; I <= N; I++) {
        PIVOT[I] = 0.0;
        INDEX[I] = 0;
    }

    // Perform successive pivot operations (grand loop)
    for (int I = 1; I <= N; I++) {

        // Search for pivot element
        AMAX = 0.0;
        for (int J = 1; J <= N; J++) {
            if (PIVOT[J] != 0.0) goto label_105;
            for (int K = 1; K <= N; K++) {
                if (PIVOT[K] != 0.0) goto label_100;
                // TEMP = DABS(A(J,K))
                TEMP = DABS(A[J + (K - 1) * NMAX]);
                if (TEMP < AMAX) goto label_100;
                IROW = J;
                ICOLUM = K;
                AMAX = TEMP;
            label_100:;
            }
        label_105:;
        }

        INDEX[I] = 4096 * IROW + ICOLUM;
        int J = IROW;
        // AMAX = A(J,ICOLUM)
        AMAX = A[J + (ICOLUM - 1) * NMAX];
        DETERM = AMAX * DETERM;

        if (DABS(DETERM) < DETMAX) goto label_130;
        DETERM = DETERM * DETMIN;
        IDET = IDET + 1;
        goto label_140;

    label_130:
        if (DABS(DETERM) > DETMIN) goto label_140;
        DETERM = DETERM * DETMAX;
        IDET = IDET - 1;

    label_140:
        // Return if matrix is singular (zero pivot)
        if (DETERM == 0.0) goto label_600;

        PIVOT[ICOLUM] = AMAX;

        // Interchange rows to put pivot element on diagonal
        if (IROW == ICOLUM) goto label_260;
        DETERM = -DETERM;
        for (int K = 1; K <= N; K++) {
            // SWAP = A(J,K); A(J,K) = A(ICOLUM,K); A(ICOLUM,K) = SWAP
            SWAP = A[J + (K - 1) * NMAX];
            A[J + (K - 1) * NMAX] = A[ICOLUM + (K - 1) * NMAX];
            A[ICOLUM + (K - 1) * NMAX] = SWAP;
        }
        if (M <= 0) goto label_260;
        for (int K = 1; K <= M; K++) {
            // SWAP = B(J,K); B(J,K) = B(ICOLUM,K); B(ICOLUM,K) = SWAP
            SWAP = B[J + (K - 1) * NMAX];
            B[J + (K - 1) * NMAX] = B[ICOLUM + (K - 1) * NMAX];
            B[ICOLUM + (K - 1) * NMAX] = SWAP;
        }

    label_260:
        // Divide pivot row by pivot element
        // A(ICOLUM,ICOLUM) = 1.0
        A[ICOLUM + (ICOLUM - 1) * NMAX] = 1.0;
        for (int K = 1; K <= N; K++) {
            // A(ICOLUM,K) = A(ICOLUM,K)/AMAX
            A[ICOLUM + (K - 1) * NMAX] = A[ICOLUM + (K - 1) * NMAX] / AMAX;
        }
        if (M <= 0) goto label_380;
        for (int K = 1; K <= M; K++) {
            // B(ICOLUM,K) = B(ICOLUM,K)/AMAX
            B[ICOLUM + (K - 1) * NMAX] = B[ICOLUM + (K - 1) * NMAX] / AMAX;
        }

    label_380:
        // Reduce non-pivot rows
        for (int J2 = 1; J2 <= N; J2++) {
            if (J2 == ICOLUM) continue; // GO TO 550
            // T = A(J2,ICOLUM)
            T = A[J2 + (ICOLUM - 1) * NMAX];
            // A(J2,ICOLUM) = 0.0
            A[J2 + (ICOLUM - 1) * NMAX] = 0.0;
            for (int K = 1; K <= N; K++) {
                // A(J2,K) = A(J2,K) - A(ICOLUM,K)*T
                A[J2 + (K - 1) * NMAX] = A[J2 + (K - 1) * NMAX]
                                        - A[ICOLUM + (K - 1) * NMAX] * T;
            }
            if (M <= 0) continue; // GO TO 550
            for (int K = 1; K <= M; K++) {
                // B(J2,K) = B(J2,K) - B(ICOLUM,K)*T
                B[J2 + (K - 1) * NMAX] = B[J2 + (K - 1) * NMAX]
                                        - B[ICOLUM + (K - 1) * NMAX] * T;
            }
        } // 550 CONTINUE
    } // end grand loop on I

label_600:
    // Interchange columns after all pivot operations
    for (int I = 1; I <= N; I++) {
        int I1 = N + 1 - I;
        int K = INDEX[I1] / 4096;
        ICOLUM = INDEX[I1] - 4096 * K;
        if (K == ICOLUM) goto label_710;
        for (int J = 1; J <= N; J++) {
            // SWAP = A(J,K); A(J,K) = A(J,ICOLUM); A(J,ICOLUM) = SWAP
            SWAP = A[J + (K - 1) * NMAX];
            A[J + (K - 1) * NMAX] = A[J + (ICOLUM - 1) * NMAX];
            A[J + (ICOLUM - 1) * NMAX] = SWAP;
        }
    label_710:;
    }

    // Store IDET in PIVOT(1) — Fortran stores int into REAL*4
    PIVOT[1] = static_cast<double>(IDET);
    return;
}

// ==========================================================================
// PLMSUB — Associated Legendre functions P(L,M)(X) (lines 2949-3072)
// ==========================================================================
void PLMSUB(int LMAX, int MMAXin, double X, double* PRAY)
{
    int MMAX = MMAXin;
    double ROOT;

    if (MMAX == 0) goto label_200;

    ROOT = std::sqrt(std::fabs(1.0 - X * X));
    if (std::fabs(X) > 1.0) ROOT = -ROOT;

label_200:
    PRAY[1] = 1.0;

    int I = 1;
    int J = 1;
    int IST, IEND;

    for (int M = 0; M <= MMAX; M++) {
        if (M >= LMAX) return;
        J = I;

        // Generate the P(L=M+1, M) for this M
        I = I + 1;
        PRAY[I] = (2 * M + 1) * X * PRAY[I - 1];
        int MP2 = M + 2;
        if (LMAX < MP2) goto label_450;

        // Now do L = M+2, M+3, ..., LMAX for this M
        {
            double DL = MP2;
            double TEMP1 = (2 * DL - 1) * X;
            double TEMP2 = M - 1 + DL;
            double TEMP3 = DL - M;
            IST = I + 1;
            IEND = I + 1 + LMAX - MP2;

            for (I = IST; I <= IEND; I++) {
                PRAY[I] = (TEMP1 * PRAY[I - 1] - TEMP2 * PRAY[I - 2]) / TEMP3;
                TEMP1 = TEMP1 + (2.0 * X);
                TEMP2 = TEMP2 + 1.0;
                TEMP3 = TEMP3 + 1.0;
            }
        }

        // label_400:
        I = IEND;

    label_450:
        if (M == MMAX) return;

        // Generate the next P(L=M+1, M+1)
        // J points to the last P(M, M)
        I = I + 1;
        PRAY[I] = -(2 * M + 1) * ROOT * PRAY[J];
    } // 499 CONTINUE

    return;
}

// ==========================================================================
// PLSUB — Legendre polynomials P(L,0)(X) (ENTRY in Fortran)
// ==========================================================================
void PLSUB(int LMAX, double X, double* PRAY)
{
    PLMSUB(LMAX, 0, X, PRAY);
}

// ==========================================================================
// SPLNCB — Natural cubic spline (lines 4946-5133)
// ==========================================================================
void SPLNCB(int N, double* X, double* Y, double* B, double* C, double* D)
{
    // Natural cubic by default
    // IOPT = 1

    if (N <= 1) return;

    int NM1 = N - 1;
    int IBASE;
    int I, I1, J;
    double H, F, G, E, GBY3, EPSIM1, RIM1B3;

    // Long strings of constant values can cause underflow problems;
    // remove such strings from the start.
    for (I = 1; I <= NM1; I++) {
        B[I] = 0.0;
        C[I] = 0.0;
        D[I] = 0.0;
        if (1.0e15 * DABS(Y[I + 1] - Y[I]) > DABS(Y[I])) goto label_80;
    }
    I = NM1;

label_80:
    IBASE = I - 1;
    for (I1 = I; I1 <= NM1; I1++) {
        J = N + I - I1;
        B[J] = 0.0;
        C[J] = 0.0;
        D[J] = 0.0;
        if (1.0e15 * DABS(Y[J - 1] - Y[J]) > DABS(Y[J])) goto label_90;
    }

label_90:
    {
        int NUSE = J - IBASE;
        if (NUSE <= 1) return;
        NM1 = NUSE - 1;

        H = X[IBASE + 2] - X[IBASE + 1];
        F = (Y[IBASE + 2] - Y[IBASE + 1]) / H;
        if (NUSE == 2) goto label_800;

        for (I = 2; I <= NM1; I++) {
            G = H;
            H = X[IBASE + I + 1] - X[IBASE + I];
            E = F;
            F = (Y[IBASE + I + 1] - Y[IBASE + I]) / H;
            GBY3 = G / 3.0;
            D[IBASE + I - 1] = GBY3 * B[IBASE + I - 1];
            EPSIM1 = G + H;
            RIM1B3 = F - E;

            // label_300:
            B[IBASE + I] = 1.0 / ((2.0 / 3.0) * EPSIM1 - GBY3 * D[IBASE + I - 1]);
            C[IBASE + I] = RIM1B3 - D[IBASE + I - 1] * C[IBASE + I - 1];
        }
        D[IBASE + NM1] = 0.0;

        // Back substitution
        for (I1 = 2; I1 <= NM1; I1++) {
            I = NM1 + 2 - I1;
            C[IBASE + I] = B[IBASE + I] * C[IBASE + I]
                         - D[IBASE + I] * C[IBASE + I + 1];
        }

        // Compute B, C, D from the tridiagonal solution
        // label_600:
        for (I = 1; I <= NM1; I++) {
            H = X[IBASE + I + 1] - X[IBASE + I];
            D[IBASE + I] = (C[IBASE + I + 1] - C[IBASE + I]) / (3.0 * H);
            B[IBASE + I] = (Y[IBASE + I + 1] - Y[IBASE + I]) / H
                         - (H * D[IBASE + I] + C[IBASE + I]) * H;
        }

        // Define the N'th cubic as the N-1'th cubic
        D[IBASE + NUSE] = D[IBASE + NM1];
        B[IBASE + NUSE] = B[IBASE + NM1]
                        + (2 * C[IBASE + NM1] + 3 * D[IBASE + NM1] * H) * H;

        return;

    label_800:
        // Special case for N = 2
        D[IBASE + 1] = 0.0;
        B[IBASE + 1] = F;
        return;
    }
}

// ==========================================================================
// SYSERR — Fatal error handler (lines 5135-5146)
// ==========================================================================
void SYSERR()
{
    std::printf("\n\n\n$*$*$*$*$*$  SYSERR CALLED     "
                "$*$*$*$*$*$*$*$*$*$*$*$*$*$*$*$*$*$*$*$*\n");
    std::exit(9876);
}

// ==========================================================================
// YLMSUB — Spherical harmonics Y(L,M)(cos theta, phi=0) (lines 5148-5362)
// ==========================================================================
void YLMSUB(int LMAX, int MMAX, double Z, double* YLMRAY)
{
    double ORT4PI = 0.28209479177387814;
    int IZERO = 0;

    double ZSQ, ROOT, DM, DL, A, ALAST;
    int I, J, MP2, LMMAX, MMAX1, LMAX1, L1;

    // Check the arguments for validity
    ZSQ = Z * Z;
    if (LMAX < 0 || MMAX < 0 || ZSQ > 1.00000000000200e+0)
        goto label_120;

    goto label_200;

label_120:
    std::printf("\n$*$*$* INVALID INPUT TO YLMSUB:%15d%15d%15d%15d%20.10g\n",
                LMAX, MMAX, I, ROOTIS.IMXDIM, Z);
    SYSERR();

label_200:
    if (MMAX != 0) ROOT = DSQRT(DABS(1.0 - ZSQ));

    YLMRAY[1] = ORT4PI;
    if (LMAX == 0) return;

    I = 1;

    // The major loop is on M
    LMMAX = MIN0(MMAX, LMAX - 1);

    MMAX1 = ROOTIS.IMXDIM / 2 - 3;
    LMAX1 = MIN0(LMAX, MMAX1);

    DM = 0.0;

    for (int M = IZERO; M <= LMMAX; M++) {

        J = I;
        I = I + 1;
        MP2 = M + 2;
        if (M > MMAX1) goto label_500;

        // Generate the Y(L=M+1, M) for this M
        // ROOTI is 1-based: ROOTI(k) = sqrt(k-1)
        YLMRAY[I] = ROOTIS.ROOTI[4 + 2 * M] * Z * YLMRAY[I - 1];
        if (LMAX < MP2) goto label_850;

        // Now do L = M+2, M+3, ..., LMAX for this M
        ALAST = 1.0 / ROOTIS.ROOTI[2 * MP2];

        for (int L = MP2; L <= LMAX1; L++) {
            I = I + 1;
            A = ROOTIS.ROOTI[L - M + 1] * ROOTIS.ROOTI[L + M + 1]
              / (ROOTIS.ROOTI[2 * L] * ROOTIS.ROOTI[2 * L + 2]);
            YLMRAY[I] = (Z * YLMRAY[I - 1] - ALAST * YLMRAY[I - 2]) / A;
            ALAST = A;
        }

        if (LMAX1 >= LMAX) goto label_850;
        L1 = LMAX1 + 1;
        goto label_600;

    label_500:
        // Large LMAX: compute required square roots every time
        // Generate the Y(L=M+1, M) for this M
        YLMRAY[I] = DSQRT(3.0 + DM + DM) * Z * YLMRAY[I - 1];
        if (LMAX < MP2) goto label_850;
        L1 = MP2;

    label_600:
        // Now do L = L1, L1+1, ..., LMAX for this M
        DL = L1;
        ALAST = DSQRT((DL - 1.0 - DM) * (DL - 1.0 + DM)
                     / ((DL + DL - 1.0) * (DL + DL - 3.0)));

        for (int L = L1; L <= LMAX; L++) {
            I = I + 1;
            A = DSQRT((DL - DM) * (DL + DM)
                     / ((DL + DL + 1.0) * (DL + DL - 1.0)));
            YLMRAY[I] = (Z * YLMRAY[I - 1] - ALAST * YLMRAY[I - 2]) / A;
            ALAST = A;
            DL = DL + 1.0;
        }

    label_850:
        if (M == MMAX) return;

        // Generate the next Y(L=M+1, M+1)
        // J points to the last Y(M, M)
        I = I + 1;
        if (M > MMAX1) goto label_880;
        YLMRAY[I] = -(ROOTIS.ROOTI[4 + 2 * M] / ROOTIS.ROOTI[3 + 2 * M])
                   * ROOT * YLMRAY[J];
        goto label_890;

    label_880:
        YLMRAY[I] = -(DSQRT(3.0 + DM + DM) / DSQRT(2.0 + DM + DM))
                   * ROOT * YLMRAY[J];

    label_890:
        DM = DM + 1.0;
    } // 899 CONTINUE

    return;
}
