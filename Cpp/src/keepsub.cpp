// keepsub.cpp -- translated from keepsub.f
//
// KEEP  - writes a Speakeasy KEEP object to a file (Unix version)
// KEEPCK - same as KEEP on Unix (extra checking only on IBM)
//
// Original: 5/1/83 JEG; 3/20/04 s.p. - Unix version for Ptolemy
//
// Translated from Fortran to C++ preserving all logic and labels.

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include "ptolemy_forward.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_io.h"
#include <cstdio>
#include <cstring>

// ============================================================================
// Card buffer — 80 bytes matching Fortran EQUIVALENCE
// EQUIVALENCE (CARD1(1),CARD(1)), (CARD2(1),CARD(1)),
//             (CARD4(1),CARD(1)), (RCARD(1),CARD(1))
// ============================================================================
struct KeepCardBuf {
    union {
        char     card1[80];      // CHARACTER*1 CARD1(80)
        double   rcard[10];      // REAL*8 RCARD(10)
        int32_t  card4[20];      // INTEGER*4 CARD4(20)
        int16_t  card2[40];      // INTEGER*2 CARD2(40)
    };
    // Access as 8-byte CHARACTER*8 word (1-based Fortran index)
    char* word(int i) { return &card1[(i-1)*8]; }
};

// ============================================================================
// KEEP subroutine
//
// Writes an object to a file in Speakeasy KEEP format.
// OBJ is the data (treated as 8-byte words).
// IDIM is 0-based in C++ (IDIM[0] = first dimension).
// ============================================================================

void KEEP(char* OBJ, int KIND, int IFAM, int NDIM, int* IDIM,
          char* LIB, char* NAME, int& IER)
{
    // implicit real*8 (a-h, o-z)

    // The sizes (in bytes) and boundaries as a function of Speakeasy kind
    // 1-based: index by KIND (1..9)
    static const int ISIZE[]  = {0, 4, 8, 0, 16, 4, 8, 2, 1, 1};
    static const int IBOUND[] = {0, 4, 8, 0,  8, 4, 8, 2, 1, 1};

    KeepCardBuf card;
    char FILE_arr[258];   // path buffer

    // NAME handling: NAME2/NAME1 equivalenced
    char NAME2[8];
    char* NAME1 = NAME2;

    char TIME[8];  std::memset(TIME, ' ', 8);
    char DATE[8];  std::memset(DATE, ' ', 8);

    int32_t ZFFA3 = static_cast<int32_t>(0xFFA30000u);
    int32_t VER   = 4061000;

    int IDCB, ICODE;
    int NELS, NN, IFACT, LHEAD, NTWORD, NWORD, NLEFT, IOUT, LMAX;

    // --- Both KEEP and KEEPCK jump here (label 20) ---
    // On Unix, KEEPCK is equivalent to KEEP.

    IER = 0;
    std::memcpy(NAME2, NAME, 8);

    // Validate inputs
    if (KIND <= 0 || KIND >= 10) goto L900;
    if (KIND == 3) goto L900;
    if (NDIM < 0 || NDIM > 15) goto L900;

    NELS = 1;
    if (NDIM == 0) goto L30;
    for (int I = 0; I < NDIM; I++) {
        if (IDIM[I] < 0) goto L900;
        NELS = NELS * IDIM[I];
    }
L30:
    if (KIND == 4) NELS = 2 * NELS;

    // Construct the full path name of the file
    {
        int K, L;
        for (K = 1; K <= 247; K++) {
            L = K;
            if (LIB[K-1] == ' ') goto L48;
            FILE_arr[K-1] = LIB[K-1];
        }
        // Loop completed without finding blank
        L = 247;
L48:
        if (FILE_arr[L-2] != '/' && L != 1) {
            FILE_arr[L-1] = '/';
        } else {
            L = L - 1;
        }
        int J;
        K = 0;
        for (J = 1; J <= 8; J++) {
            K = J;
            if (NAME1[J-1] == ' ') goto L66;
            FILE_arr[L+J-1] = NAME1[J-1];
        }
        K = 9;
L66:
        int pos = L + K - 1;
        FILE_arr[pos]   = '.';
        FILE_arr[pos+1] = 'k';
        FILE_arr[pos+2] = 'e';
        FILE_arr[pos+3] = 'e';
        FILE_arr[pos+4] = 'p';
        FILE_arr[pos+5] = '\0';
    }

    // CALL SWOPEN(IDCB,ICODE,FFILE,1,0)
    {
        int oversw = 1, append_val = 0;
        swopen(&IDCB, &ICODE, FILE_arr, &oversw, &append_val);
    }
    if (ICODE != 0) goto L920;

    // Set up the first card (label 200)
    std::memset(&card, 0, 80);
    card.card4[0]  = ZFFA3;          // CARD4(1)
    card.card4[1]  = 1;              // CARD4(2)
    card.card4[2]  = 0;              // CARD4(3)
    card.card4[3]  = 1;              // CARD4(4)
    card.card4[4]  = VER;            // CARD4(5)
    card.card4[5]  = 0;              // CARD4(6)
    card.card4[5]  = 0;              // CARD4(6) (set twice in Fortran)
    card.card4[6]  = ZFFA3;          // CARD4(7)
    card.card4[7]  = 4;              // CARD4(8)
    card.card4[8]  = 0;              // CARD4(9)
    card.card4[9]  = 2;              // CARD4(10)
    std::memcpy(card.word(6), TIME, 8);  // CARD(6) = TIME
    std::memcpy(card.word(7), DATE, 8);  // CARD(7) = DATE
    card.card4[14] = ZFFA3;          // CARD4(15)
    card.card4[15] = 1;              // CARD4(16)
    card.card4[16] = 0;              // CARD4(17)

    NN = 0;
    IFACT = 8 / IBOUND[KIND];
    if (((NELS / IFACT) * IFACT) != NELS) NN = 1;
    LHEAD = 2 + NDIM / 2;
    NTWORD = (NELS * IBOUND[KIND]) / 8 + LHEAD + NN + 1;
    card.card4[17] = NTWORD;              // CARD4(18)
    std::memcpy(card.word(10), NAME2, 8); // CARD(10) = NAME

    // Write out the first card
    {
        int numchr = 80;
        swrite(&IDCB, &ICODE, card.card1, &numchr);
    }
    if (ICODE != 0) goto L920;

    // Set up the second card information
    std::memset(&card, 0, 80);
    card.card1[2] = static_cast<char>(LHEAD);            // CARD1(3)
    card.card1[3] = static_cast<char>(0);                // CARD1(4)
    card.card2[2] = static_cast<int16_t>(KIND);          // CARD2(3)
    card.card1[6] = static_cast<char>(IFAM);             // CARD1(7)
    card.card1[7] = static_cast<char>(NDIM);             // CARD1(8)
    card.card1[8] = static_cast<char>(ISIZE[KIND]);      // CARD1(9)
    card.card1[9] = static_cast<char>(IBOUND[KIND]);     // CARD1(10)

    if (NDIM != 0) {
        for (int K = 0; K < NDIM; K++) {
            card.card4[K+3] = IDIM[K];   // CARD4(K+4) in Fortran 1-based
        }
    }

    IOUT = 0;
    NWORD = NTWORD - LHEAD - 1;
    NLEFT = NWORD;

    // Loop to output the data, starting somewhere inside the second card (label 350)
L350:
    if (NLEFT == 0) goto L400;
    LMAX = MIN0(NLEFT, 10 - LHEAD);
    for (int K = 1; K <= LMAX; K++) {
        // CARD(K+LHEAD) = OBJ(K+IOUT)
        std::memcpy(card.word(K + LHEAD), &OBJ[(K + IOUT - 1) * 8], 8);
    }
    IOUT += LMAX;
    LHEAD += LMAX;
    NLEFT -= LMAX;

    if (LHEAD < 10) goto L400;

    // Write out the card
    {
        int numchr = 80;
        swrite(&IDCB, &ICODE, card.card1, &numchr);
    }
    if (ICODE != 0) goto L920;

    // Subsequent cards are fully data
    LHEAD = 0;
    goto L350;

    // Write out the last card, close and return (label 400)
L400:
    // CARD(LHEAD+1) = ENDDAT
    std::memcpy(card.word(LHEAD + 1), "ENDOFDAT", 8);
    LHEAD = LHEAD + 1;
    if (LHEAD > 10) goto L420;
    for (int I = LHEAD; I <= 10; I++) {
        card.rcard[I-1] = 0.0;   // RCARD(I) in Fortran 1-based
    }
L420:
    {
        int numchr = 80;
        swrite(&IDCB, &ICODE, card.card1, &numchr);
    }
    if (ICODE != 0) goto L920;
    swclse(&IDCB);
    return;

L920:
    IER = 1;
    return;

L900:
    IER = -1;
    return;
}

// ============================================================================
// KEEPCK — equivalent to KEEP on Unix
// (Extra PDS checking is IBM-only, not applicable here)
// ============================================================================

void KEEPCK(char* OBJ, int KIND, int IFAM, int NDIM, int* IDIM,
            char* LIB, char* NAME, int& IER)
{
    KEEP(OBJ, KIND, IFAM, NDIM, IDIM, LIB, NAME, IER);
}
