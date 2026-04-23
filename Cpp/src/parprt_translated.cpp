// parprt_translated.cpp — PARPRT subroutine
// Translated from source.f lines 28139-28241
// Prints one potential term row; suppresses when VPAR=0 (except I=3,4 spin-orbit).

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include "ptolemy_alloc.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_forward.h"
#include <cstdio>
#include <cmath>

// Fortran Gw.d format: prints val in a field of width w with d significant digits.
// When 0.1 <= |val| < 10^d: uses F(w-4).(d-e) followed by 4 trailing spaces,
// where e = number of digits before the decimal point (0 if |val| < 1).
// Otherwise: uses E(w).(d-1) format.
void print_G(int w, int d, double val)
{
    double absval = std::fabs(val);
    double limit  = std::pow(10.0, (double)d);
    if (absval == 0.0 || (absval >= 0.1 && absval < limit)) {
        int e = 0;
        if (absval >= 1.0)
            e = (int)std::floor(std::log10(absval)) + 1;
        int dec = d - e;
        if (dec < 0) dec = 0;
        // F(w-4).dec + 4 trailing spaces = total w chars
        std::printf("%*.*f    ", w - 4, dec, val);
    } else {
        // Fortran E notation: 0.dddddE±ee (exponent is 1 higher than C's)
        int exp_f = (int)std::floor(std::log10(absval)) + 1;
        double mant = val / std::pow(10.0, (double)exp_f);
        std::printf("%*.*fE%+03d", w - 4, d, mant, exp_f);
    }
}

// PARPRT: print parameters for one potential term.
//
// I        term index in LNKADR (0 = no linkule possible).
// NAME     18-char label for the term (right-padded with spaces).
// VPAR     coupling constant.
// RPAR     radius.
// APAR     diffuseness.
// RPAR0    radius parameter (R/R0MASS).
// P5       optional 5th parameter (printed only if nonzero).
// P5NAME   8-char label for P5.
// IRET     error flag; set to -1 on linkule error.
//
// Suppression rules (no LINKULE case):
//   if VPAR == 0 and I != 3 and I != 4: return (suppress)
//   if VPAR == 0 and (I == 3 or 4) and P5 == 0: return (suppress)
//
// FORMAT 113: ' ', a18, G14.5, T37, F7.4, T49, F7.4, T67, F7.4 [, 5X, A8, ' =', G13.5]
// Column layout (1-based, Fortran carriage control ' ' at col 1):
//   Col  1    : ' '
//   Cols 2-19 : NAME (a18)
//   Cols 20-33: VPAR (G14.5 = F10.dec + 4 spaces)
//   T37 fill  : cols 34-36 (3 spaces)
//   Cols 37-43: RPAR (F7.4)
//   T49 fill  : cols 44-48 (5 spaces)
//   Cols 49-55: APAR (F7.4)
//   T67 fill  : cols 56-66 (11 spaces)
//   Cols 67-73: RPAR0 (F7.4)
//   5X        : cols 74-78 (5 spaces)
//   Cols 79-86: P5NAME (A8)
//   Cols 87-88: ' ='
//   Cols 89-101: P5 (G13.5)

void PARPRT(int I, const char* NAME, double VPAR, double RPAR, double APAR,
            double RPAR0, double P5, const char* P5NAME, int& IRET)
{
    // If a linkule exists for this term, call it to print its own info.
    if (I > 0 && I <= NUMLINKULES && LNKBLK.LNKADR[I][3] != 0) {
        // Blank line before (FORMAT '0'), then header
        std::printf("\n%.18s POTENTIAL IS BEING COMPUTED BY THE %.8s LINKULE:\n",
                    NAME, (char*)&LNKBLK.LNKADR[I][1]);
        // Full LINKUL call with mode=2 (print) would go here.
        // Linkule printing not yet fully implemented — fall through to standard print.
        std::printf(" \n");
        return;
    }

    // No linkule: print 4 or 5 parameters, unless VPAR is zero.
    if (VPAR == 0.0) {
        if (I != 3 && I != 4) return;   // suppress non-spin-orbit zero terms
        if (P5 == 0.0) return;           // suppress spin-orbit when both VSO and TAU are 0
    }

    // Leading space (Fortran carriage control ' ') + NAME (a18)
    std::printf(" %-18.18s", NAME);
    // VPAR in G14.5
    print_G(14, 5, VPAR);
    // 3 fill spaces to reach T37, then RPAR in F7.4
    std::printf("   %7.4f", RPAR);
    // 5 fill spaces to reach T49, then APAR in F7.4
    std::printf("     %7.4f", APAR);
    // 11 fill spaces to reach T67, then RPAR0 in F7.4
    std::printf("           %7.4f", RPAR0);

    if (P5 != 0.0) {
        // 5X + P5NAME (A8) + ' =' + G13.5
        std::printf("     %-8.8s =", P5NAME);
        print_G(13, 5, P5);
    }
    std::printf("\n");
}
