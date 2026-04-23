// gfortran_stuff.cpp — translated from gfortran_stuff.f
//
// Contains:
//   CMPUTR  — returns machine name and hostname
//   GRAB    — allocates memory from the ALLOC common block pool
//   DGAMMA  — complete gamma function with ENTRY points DGAMAP, DLGAMA, DLNGAM, DIGAMA
//

#include "ptolemy_commons.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_alloc.h"
#include "ptolemy_io.h"
#include <cstdio>
#include <cstring>
#include <cstdint>
#include <unistd.h>

// ============================================================================
// SUBROUTINE CMPUTR ( MACH, hostname )
// ============================================================================
void CMPUTR(char* MACH, char* hostname) {
    //
    // Returns machine name and hostname
    //
    std::strcpy(MACH, "Intel Pentium, Redhat Linux");

    std::memset(hostname, ' ', 1);
    hostname[0] = '\0';
    // Fortran hostnm intrinsic — use POSIX gethostname
    char hbuf[256];
    int ierr = gethostname(hbuf, sizeof(hbuf));
    if (ierr == 0) {
        std::strcpy(hostname, hbuf);
    }
}

// ============================================================================
// SUBROUTINE GRAB ( S, I4OR8, IBASE, NWORDS )
// ============================================================================
//
// THE SPECIFIED AMOUNT OF MEMORY IS OBTAINED.
//
// FOR USE IN A ONE-DIMENSIONAL ARRAY  'S'.... S[IBASE] IS THE FIRST
// WORD OF THE ALLOCATED CORE.
//
//     DESCRIPTION OF PARAMETERS
//  S    -   THE ARRAY RELATIVE TO WHICH 'GRABBED' CORE IS ADDRESSED
//  I4OR8  -   AN I*4 VARIABLE WHOSE VALUE IS THE NUMBER OF BYTES PER
//             WORD IN THE ARRAY 'S'.  1,2,4,8,16 AND 32 ARE VALID
//  IBASE - OUTPUT:  THE INDEX OF THE FIRST WORD OF THE 'GRABBED' CORE
//  NWORDS -   THE NUMBER OF I4OR8-BYTE WORDS ALLOCATED
//
void GRAB(double* S, int I4OR8, int& IBASE, int& NWORDS) {
    //
    // COMMON /ALLOCS/ I1, I2, I3, intjunk, ALLOC( com_mem_size )
    // In C++: ALLOCS.FACFR4 (=I1), ALLOCS.FACFR2 (=I2), ALLOCS.FACFR1 (=I3),
    //         ALLOCS.intjunk, ALLOCS.ALLOC_pool[]
    //
    constexpr int com_mem_size = COM_MEM_SIZE;
    int MAXWRD = com_mem_size;

    // LOC(S) — address of S in bytes
    intptr_t LOCS = reinterpret_cast<intptr_t>(S);

    //
    // FIXED-SIZE COMMON MEMORY
    //
    NWORDS = (MAXWRD * 8 / I4OR8) - 10;
    printf("0FIXED COMMON BLOCK VERSION OF PTOLEMY, "
           " ALLOCATOR SIZE =%9d WORDS.\n\n", NWORDS);

    //
    // In the fixed-common version we know that ALLOC above is offset
    // by 2 double words from the real ALLOC which will be at the
    // beginning of /ALLOCS/.  We must space over the leading
    // two words.  We remove 10 words for safety.
    //
    intptr_t ITSAT = LOCS + 2 * 8;

    //
    // If the LOC function is working, one may use
    //
    intptr_t ITSAT2 = reinterpret_cast<intptr_t>(&ALLOCS.ALLOC_pool[1]);
    printf(" Check of LOC: %ld %ld\n", (long)ITSAT, (long)ITSAT2);

    //
    // HERE WE HAVE nwords words STARTING at byte LOCATION ITSAT
    //
    //  <------------ ITSAT BYTES ---------> <-------- IGOT BYTES ----->
    //
    // 0                S[1]                S[IBASE]        S[IBASE+NWORDS]
    //  <-loc(S) BYTES->
    //
    IBASE = (int)((ITSAT - LOCS) / I4OR8) + 1;

    printf(" CURRENT ALLOCATOR SIZE =%11d words"
           "          BIASES = %12ld%12ld%12d\n",
           NWORDS, (long)ITSAT, (long)LOCS, IBASE);
}

// ============================================================================
// FUNCTION DGAMMA with ENTRY points DGAMAP, DLGAMA, DLNGAM, DIGAMA
// ============================================================================
//
// Complete Gamma Function and related functions.
//
// Entry Points (all are double precision):
//
//   DGAMMA:  Complete Gamma Function.
//   DGAMAP:  First derivative of the Complete Gamma Function.
//   DLGAMA, DLNGAM:  Natural logarithm of DABS( GAMMA(X) ).
//   DIGAMA:  Digamma, or Psi, function.
//

enum DGammaEntry {
    ENTRY_DGAMMA = 1,
    ENTRY_DGAMAP = 2,
    ENTRY_DLGAMA = 3,
    ENTRY_DIGAMA = 4
};

static double dgamma_impl(double X, DGammaEntry entry) {
    //
    // Local variables (shared across all entry paths)
    //
    bool LNEG, LEXP;
    int NE, N, I;
    double Z, H, Y, G, D, T, ZI, ZIS;

    static const char NAME[4][7] = { "DGAMMA", "DGAMAP", "DLGAMA", "DIGAMA" };
    static const double PI    = 3.141592653589793240e0;
    static const double R2PILN = 0.918938533204672420e0;

    //
    // ELIM is the biggest number s.t. DEXP(ELIM) doesn't overflow.
    //
    static const double ELIM = 695.0e0;

    // Set entry point
    NE = static_cast<int>(entry);

    // - - - - - - - - - - DGAMMA, DGAMAP, DLGAMA enter at 100 - - - - - - - -
    // - - - - - - - - - - DIGAMA enters at 500 - - - - - - - - - - - - - - - -

    if (NE == 4) goto label_500;

    // Log Gamma code common to DGAMMA, DLGAMA, and DLNGAM.

    // Check for argument < 0.5 and transform.

label_100:
    Z = X;
    LEXP = false;
    LNEG = Z < 0.50e0;
    if (LNEG) Z = 1.0e0 - Z;

    // If Z < 12, we transform it down into the range 0.5 <= Z < 1.5
    // (A&S 6.1.15) and use a 19-term series expansion of 1/Gamma(1+Z)
    // (from A&S 6.1.34).

    if (Z >= 12.0e0) goto label_200;
    H = 1.0e0;
    N = (int)(Z - 0.50e0);
    if (N == 0) goto label_120;
    for (I = 1; I <= N; I++) {
        Z = Z - 1.0e0;
        H = H * Z;
    }

label_120:
    Y = Z - 1.0e0;
    G = -.000020134854780700e0;
    if (DABS(Y) > 0.10e0)
        G = G
            + Y * (-.00000125049348210e0
            + Y * ( 0.00000113302723200e0
            + Y * (-.00000020563384170e0
            + Y * ( 0.00000000611609500e0
            + Y * ( 0.00000000500200750e0
            + Y * (-.00000000118127460e0
            + Y * ( 0.00000000010434270e0
            )))))));
    G = 1.0e0
        + Y * ( 0.577215664901532900e0
        + Y * (-.655878071520253800e0
        + Y * (-.042002635034095200e0
        + Y * ( 0.166538611382291500e0
        + Y * (-.042197734555544300e0
        + Y * (-.009621971527877000e0
        + Y * ( 0.007218943246663000e0
        + Y * (-.001165167591859100e0
        + Y * (-.000215241674114900e0
        + Y * ( 0.000128050282388200e0
        + Y * G ))))))))));
    G = H / G;
    LEXP = true;
    goto label_250;

    // Calculate log gamma using 6 terms of the Bernoulli expansion.
    // For |X| > 164, we can get 17 decimal places with 2 terms.

label_200:
    ZI = 1.0e0 / Z;
    ZIS = ZI * ZI;
    G = (-1.0e0 / 360.0e0);
    if (Z < 164.0e0)
        G = G
            + ZIS * ((1.0e0 / 1260.0e0)
            + ZIS * ((-1.0e0 / 1680.0e0)
            + ZIS * ((1.0e0 / 1188.0e0)
            + ZIS * ((-691.0e0 / 312840.0e0)
            ))));
    G = R2PILN - Z + (Z - 0.50e0) * DLOG(Z)
        + ZI * ((1.0e0 / 12.0e0)
        + ZIS * G);

    // If X < 0.5, transform the result using A&S 6.1.17.

label_250:
    if (!LNEG) goto label_300;
        T = DSIN(PI * Z);
        if (T != 0.0e0) goto label_260;
            printf("\n **** ERROR IN %6s:  ARGUMENT = %25.16G ****\n", NAME[NE - 1], X);
            printf(" **** ARGUMENT IS A NEGATIVE INTEGER OR ZERO ****\n");
label_260:
        if (LEXP) goto label_270;
            G = DLOG(PI / DABS(T)) - G;
            LNEG = T < 0;
            goto label_300;
label_270:
        G = PI / (G * T);
        LNEG = false;

    // We now have log gamma (or gamma if LEXP = .TRUE.).
    // Branch depends on the entry point used.

label_300:
    if (NE == 3) goto label_30;
    if (LEXP) {
        if (NE == 1) goto label_10;
        if (NE == 2) goto label_500;
    }
    if (G <= ELIM) goto label_350;
        printf("\n **** ERROR IN %6s:  ARGUMENT = %25.16G ****\n", NAME[NE - 1], X);
        printf(" **** EXPONENT OVERFLOW ****\n");
label_350:
    G = DEXP(G);
    if (LNEG) G = -G;
    if (NE == 1) goto label_10;

    // Digamma code common to DGAMAP, DIGAMA.

    // Check for argument < 0.5 and transform.

label_500:
    Z = X;
    LNEG = Z < 0.50e0;
    if (LNEG) Z = 1.0e0 - Z;

    // Transform to Z >= 9 if necessary.

    H = 0.0e0;
    if (Z >= 9.0e0) goto label_600;
    N = (int)(10.0e0 - Z);
    for (I = 1; I <= N; I++) {
        H = H + 1.0e0 / Z;
        Z = Z + 1.0e0;
    }

    // Calculate Digamma using 8 terms of the Bernoulli expansion.

label_600:
    ZI = 1.0e0 / Z;
    ZIS = ZI * ZI;
    D = DLOG(Z) - 0.50e0 * ZI - H
        - ZIS * ((1.0e0 / 12.0e0)
        + ZIS * ((-1.0e0 / 120.0e0)
        + ZIS * ((1.0e0 / 252.0e0)
        + ZIS * ((-1.0e0 / 240.0e0)
        + ZIS * ((1.0e0 / 132.0e0)
        + ZIS * ((-691.0e0 / 32760.0e0)
        + ZIS * ((1.0e0 / 12.0e0)
        + ZIS * ((-3617.0e0 / 8160.0e0)
        ))))))));

    // Modify the result if the argument was negative.

    if (!LNEG) goto label_700;
        T = DSIN(PI * Z);
        if (T != 0.0e0) goto label_690;
            printf("\n **** ERROR IN %6s:  ARGUMENT = %25.16G ****\n", NAME[NE - 1], X);
            printf(" **** ARGUMENT IS A NEGATIVE INTEGER OR ZERO ****\n");
label_690:
        D = D + PI * DCOS(PI * Z) / T;
label_700:
    if (NE == 4) goto label_40;
    goto label_20;

    // Return points for each entry

label_10:
    return G;

label_20:
    return G * D;

label_30:
    if (LEXP) G = DLOG(G);
    return G;

label_40:
    return D;
}

// Inline wrapper functions for each entry point

double DGAMMA(double X) {
    return dgamma_impl(X, ENTRY_DGAMMA);
}

double DGAMAP(double X) {
    return dgamma_impl(X, ENTRY_DGAMAP);
}

double DLGAMA(double X) {
    return dgamma_impl(X, ENTRY_DLGAMA);
}

double DLNGAM(double X) {
    return dgamma_impl(X, ENTRY_DLGAMA);
}

double DIGAMA(double X) {
    return dgamma_impl(X, ENTRY_DIGAMA);
}
