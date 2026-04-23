// keptsub.cpp -- translated from keptsub.f
//
// KEPT   - reads a Speakeasy KEEP object from a file (Unix version)
// KEPTCN - continues a KEPT read after MBYTE=0 call
// KEPTCL - cancels a KEPT read after MBYTE=0 call
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
// ============================================================================
struct KeptCardBuf {
    union {
        char     card1[80];      // CHARACTER*1 CARD1(80)
        double   rcard[10];      // REAL*8 RCARD(10)
        int32_t  card4[20];      // INTEGER*4 CARD4(20)
        int16_t  card2[40];      // INTEGER*2 CARD2(40)
    };
    char* word(int i) { return &card1[(i-1)*8]; }
};

// ============================================================================
// Dispatch enum for KEPT / KEPTCN / KEPTCL
// ============================================================================
enum class KeptEntry { KEPT, KEPTCN, KEPTCL };

// ============================================================================
// Internal dispatch — shared static state persists across KEPT/KEPTCN/KEPTCL
// ============================================================================
static void KEPT_dispatch(KeptEntry entry,
    char* OBJ,
    int* KIND_ptr, int* IFAM_ptr, int* NDIMS_ptr, int* IDIM,
    char* LIB, char* NAME,
    int* IER_ptr,
    int MBYTE_arg)
{
    // SAVE'd state (persists between KEPT/KEPTCN/KEPTCL calls)
    static int IOFF, NWDS;
    static int IDCB;
    static KeptCardBuf card;

    // Local state
    int ICALL;
    int IER = 0;
    int MBYTE;
    int ICODE, NUMCHR;

    int32_t ZFFA3 = static_cast<int32_t>(0xFFA30000u);

    // ---------------------------------------------------------------
    // Entry point dispatch
    // ---------------------------------------------------------------
    if (entry == KeptEntry::KEPTCN) {
        // ENTRY KEPTCN ( OBJ, IER2, MBYTE2 )
        ICALL = 2;
        MBYTE = MBYTE_arg;
        goto L350;
    }

    if (entry == KeptEntry::KEPTCL) {
        // ENTRY KEPTCL
        ICALL = 3;
        goto L990;
    }

    // --- KEPT entry ---
    ICALL = 1;
    IER = 0;

    {
        // NAME handling: NAME2/NAME1 equivalenced
        char NAME2[8];
        char* NAME1 = NAME2;
        std::memcpy(NAME2, NAME, 8);

        // Construct file path
        char FILE_arr[258];
        int K, L;
        for (K = 1; K <= 247; K++) {
            L = K;
            if (LIB[K-1] == ' ') goto L48;
            FILE_arr[K-1] = LIB[K-1];
        }
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
        {
            int pos = L + K - 1;
            FILE_arr[pos]   = '.';
            FILE_arr[pos+1] = 'k';
            FILE_arr[pos+2] = 'e';
            FILE_arr[pos+3] = 'e';
            FILE_arr[pos+4] = 'p';
            FILE_arr[pos+5] = '\0';
        }

        // CALL SROPEN(IDCB,ICODE,FFILE)
        sropen(&IDCB, &ICODE, FILE_arr);
        if (ICODE != 0) goto L915;

        // Read the first card
        {
            int maxchr = 80, rtype = 1;
            srread(&IDCB, &ICODE, card.card1, &NUMCHR, &maxchr, &rtype);
        }
        if (ICODE == 3) goto L910;
        if (ICODE != 0) goto L920;

        // Process 1st card of the KEEP file
        if (card.card4[0] != ZFFA3 || card.card4[14] != ZFFA3) goto L900;
        if (card.card4[1] != 1) goto L900;
        if (card.card4[2] != 0) goto L900;
        if (card.card4[5] != 0) goto L900;

        if (card.card4[15] != 1) goto L905;
        NWDS = card.card4[17];   // CARD4(18)

        // Read the second card
        {
            int maxchr = 80, rtype = 1;
            srread(&IDCB, &ICODE, card.card1, &NUMCHR, &maxchr, &rtype);
        }
        if (ICODE == 3) goto L910;
        if (ICODE != 0) goto L920;

        {
            int LHEAD_val = static_cast<unsigned char>(card.card1[2]);  // ICHAR(CARD1(3))
            NWDS = NWDS - LHEAD_val - 1;

            *KIND_ptr  = card.card2[2];                                  // KIND=CARD2(3)
            *IFAM_ptr  = static_cast<unsigned char>(card.card1[6]);      // IFAM=ICHAR(CARD1(7))
            *NDIMS_ptr = static_cast<unsigned char>(card.card1[7]);      // NDIMS=ICHAR(CARD1(8))

            // int ISIZE_val  = static_cast<unsigned char>(card.card1[8]);   // not used further
            // int IBOUND_val = static_cast<unsigned char>(card.card1[9]);   // not used further

            if (*NDIMS_ptr != 0) {
                for (K = 0; K < *NDIMS_ptr; K++) {
                    IDIM[K] = card.card4[K+3];   // CARD4(K+4) in 1-based
                }
            }

            // Skip to end of the header if necessary (label 140)
            int ISKIP = LHEAD_val / 10;
            if (ISKIP != 0) {
                for (int MM = 1; MM <= ISKIP; MM++) {
                    int maxchr = 80, rtype = 1;
                    srread(&IDCB, &ICODE, card.card1, &NUMCHR, &maxchr, &rtype);
                    if (ICODE == 3) goto L910;
                    if (ICODE != 0) goto L920;
                    LHEAD_val = LHEAD_val - 10;
                }
            }
            // Label 170
            IOFF = LHEAD_val;
        }
    }

    // Label 300 — ready to process the data
    IER = 0;
    if (ICALL == 1) *IER_ptr = 0;
    if (MBYTE_arg == 0) {
        if (ICALL == 1) *IER_ptr = IER;
        return;
    }
    MBYTE = MBYTE_arg;
    goto L350;

    // ---------------------------------------------------------------
    // Data reading loop (label 350)
    // ---------------------------------------------------------------
L350:
    IER = 0;
    if (NWDS == 0) goto L990;
    {
        int MAXWRD = MBYTE / 8;
        if (MAXWRD < NWDS) IER = -4;
        if (MAXWRD == 0) goto L990;
        NWDS = MIN0(MAXWRD, NWDS);
    }
    {
        int L_data = 0;

        // Process the data on this card (label 330)
L330:
        {
            int NDO = MIN0(NWDS, 10 - IOFF);
            for (int I = 1; I <= NDO; I++) {
                // OBJ(L+I) = CARD(IOFF+I)
                std::memcpy(&OBJ[(L_data + I - 1) * 8], card.word(IOFF + I), 8);
            }
            L_data += NDO;
            NWDS -= NDO;
            IOFF = 0;
        }
        if (NWDS == 0) goto L990;
        {
            int maxchr = 80, rtype = 1;
            srread(&IDCB, &ICODE, card.card1, &NUMCHR, &maxchr, &rtype);
        }
        if (ICODE == 3) goto L910;
        if (ICODE != 0) goto L920;
        goto L330;
    }

    // ---------------------------------------------------------------
    // Error messages
    // ---------------------------------------------------------------
L900:
    IER = -1;
    goto L990;
L905:
    IER = -2;
    goto L990;
L910:
    IER = 9;
    goto L990;
L915:
    IER = 1;
    goto L990;
L920:
    IER = 5;
    goto L990;
    // L930: IER = -3; goto L990;  (Speakeasy-III path, not reached)

    // ---------------------------------------------------------------
    // Close and return (label 990)
    // ---------------------------------------------------------------
L990:
    srclse(&IDCB);
    if (ICALL == 1 && IER_ptr) *IER_ptr = IER;
    if (ICALL == 2 && IER_ptr) *IER_ptr = IER;
    return;
}

// ============================================================================
// Public interface functions
// ============================================================================

void KEPT(char* OBJ, int& KIND, int& IFAM, int& NDIMS, int* IDIM,
          char* LIB, char* NAME, int& IER1, int MBYTE1)
{
    KEPT_dispatch(KeptEntry::KEPT, OBJ, &KIND, &IFAM, &NDIMS, IDIM,
                  LIB, NAME, &IER1, MBYTE1);
}

void KEPTCN(char* OBJ, int& IER2, int MBYTE2)
{
    KEPT_dispatch(KeptEntry::KEPTCN, OBJ, nullptr, nullptr, nullptr, nullptr,
                  nullptr, nullptr, &IER2, MBYTE2);
}

void KEPTCL()
{
    KEPT_dispatch(KeptEntry::KEPTCL, nullptr, nullptr, nullptr, nullptr, nullptr,
                  nullptr, nullptr, nullptr, 0);
}
