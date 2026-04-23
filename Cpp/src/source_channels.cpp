// source_channels.cpp — Remaining stubs
// BASLBL → baslbl_translated.cpp, SETFG → setfg_translated.cpp
// CHANEL → chanel_translated.cpp, SETCHN → setchn_translated.cpp (when ready)

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include "ptolemy_alloc.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_forward.h"
#include "ptolemy_io.h"
#include <cstdio>
#include <cstring>
#include <cmath>

void BASCPL(int& IRET)
{ IRET = 0; std::printf(" BASCPL: stub\n"); }

void CCHAN(int& IRET)
{ IRET = 0; std::printf(" CCHAN: stub\n"); }

void CCOUP(int& IRET)
{ IRET = 0; std::printf(" CCOUP: stub\n"); }

// CHNSCN — Scans an input line for channel specification
// Translated from source.f lines 8267-8523
void CHNSCN(char8* GUY, double* ESTARS, int* NODVLS, int* LVALS, int* JVALS,
            int* IPARYS, int& IRTN)
{
    // Character constants
    const char SLASH = '/';
    const char PLUS  = '+';
    const char MINUS = '-';
    const char EQUAL = '=';
    const char BLANK = ' ';
    const char LEFTP = '(';

    int I, II;
    int INCHST;
    bool PREDSW, EREDSW;
    int ITOKEN, IVALU, IMESS;
    char STPCHR;
    char8 CVALU;
    double VALU, DUMVAL;
    char dummy[4] = {};   // unused MESS parameter for MSCAN

    IRTN    = 0;
    INCHST  = INPBUF.INCH;
    PREDSW  = false;
    EREDSW  = false;

    // Initialise all 3 slots (1-based)
    char8 BLANK8;
    std::memset(BLANK8.data, ' ', 8);
    for (I = 1; I <= 3; I++) {
        GUY[I]    = BLANK8;
        NODVLS[I] = NOTDEF_INT;
        LVALS[I]  = NOTDEF_INT;
        IPARYS[I] = 0;
        JVALS[I]  = NOTDEF_INT;
        ESTARS[I] = INTRNL.UNDEF;
    }

    I = 1;

// -----------------------------------------------------------------------
//  LOOP ON REACTION PARTICIPANTS
// -----------------------------------------------------------------------
L100:
    MSCAN(-20, ITOKEN, GUY[I], DUMVAL, IVALU, dummy, IMESS, STPCHR);
    if (IMESS > 5)  goto L9004;
    if (IMESS == 0) goto L9005;
    INPBUF.INCH = INPBUF.INCH + 1;   // advance past stop char

    // If stopped by '(' → process excited-state / J specification
    if (STPCHR != LEFTP) goto L400;

// -----------------------------------------------------------------------
//  PROCESS EXCITED STATE INFO  (E*, J, parity, level descriptor)
// -----------------------------------------------------------------------
L200:
    MSCAN(18, ITOKEN, CVALU, VALU, IVALU, dummy, IMESS, STPCHR);
    if (ITOKEN == -9) goto L9900;
    if (ITOKEN == -8) goto L9007;
    if (ITOKEN != -2) goto L230;

    // It is E*
    if (ESTARS[I] != INTRNL.UNDEF) {
        std::printf("0**** ATTEMPTING TO ENTER 2 VALUES FOR E*\n");
        goto L9900;
    }
    ESTARS[I] = VALU;
    goto L200;

L230:
    if (ITOKEN != -5) goto L270;

    // It is J (half-integer token)
    if (JVALS[I] != NOTDEF_INT) {
        std::printf("0**** ATTEMPTING TO ENTER 2 VALUES FOR J\n");
        goto L9900;
    }
    JVALS[I] = 2 * IVALU;
    if (STPCHR != SLASH) goto L250;

    // Bypass /2: advance past '/', then verify the '2'
    INPBUF.INCH = INPBUF.INCH + 1;
    JVALS[I] = IVALU;
    MSCAN(18, ITOKEN, CVALU, VALU, IVALU, dummy, IMESS, STPCHR);
    if (ITOKEN == -9) goto L9900;
    if (ITOKEN == -8) goto L9007;
    if (!(ITOKEN == -5 && IVALU == 2)) {
        std::printf("0**** INVALID J FIELD\n");
        goto L9900;
    }

L250:
    // Parity may be given by terminal + or -
    if (STPCHR == PLUS)  { IPARYS[I] = +1; INPBUF.INCH++; goto L200; }
    if (STPCHR == MINUS) { IPARYS[I] = -1; INPBUF.INCH++; goto L200; }
    goto L200;

L270:
    if (ITOKEN != -3) goto L280;

    // It is a nodes-and-L level descriptor (e.g. 2P3/2-)
    if (LVALS[I] != NOTDEF_INT) {
        std::printf("0**** ATTEMPTING TO ENTER TWO LEVEL DESCRIPTORS\n");
        goto L9900;
    }
    LVALS[I]  = IVALU;
    NODVLS[I] = (int)VALU;
    goto L200;

L280:
    // Is it the closing ')'?
    if (ITOKEN != 11) {
        std::printf("0**** UNRECOGNIZABLE EXCITED STATE DESCRIPTOR (TOKEN =%4d)\n", ITOKEN);
        goto L9900;
    }
    // Done: bypass the ')' already consumed by MSCAN; set STPCHR to char after it
    STPCHR         = INPBUF.IBUF[INPBUF.INCH];
    INPBUF.INCH    = INPBUF.INCH + 1;

// -----------------------------------------------------------------------
//  READY TO PROCESS THE NEXT FIELD (or end of channel spec)
// -----------------------------------------------------------------------
L400:
    INPBUF.INCH = INPBUF.INCH - 1;
    if (PREDSW && EREDSW) goto L800;

    // Skip blanks
L420:
    INPBUF.INCH = INPBUF.INCH + 1;
    if (STPCHR != BLANK) goto L450;
    if (INPBUF.INCH >= INPBUF.NOCH) goto L460;
    STPCHR = INPBUF.IBUF[INPBUF.INCH];
    goto L420;

    // Only '+' and '=' are valid delimiters inside the specification
L450:
    if (STPCHR == EQUAL) goto L500;
    if (STPCHR == PLUS)  goto L550;

    // End of channel spec: STPCHR will be reconsidered by CONTRL
    INPBUF.INCH = INPBUF.INCH - 1;
L460:
    if (PREDSW) goto L800;
    std::printf("0**** CHANNEL SPECIFICATION REQUIRES TWO"
                " NUCLEAR SYMBOLS SEPARATED BY A \"+\".\n");
    goto L9900;

// -----------------------------------------------------------------------
//  '=' SIGN — bound-state result
// -----------------------------------------------------------------------
L500:
    if (EREDSW) goto L9020;
    EREDSW = true;
    // Always store result of '=' in 3rd slot
    if (I == 2) {
        I = 3;
    } else {
        GUY[3]    = GUY[1];
        JVALS[3]  = JVALS[1];
        ESTARS[3] = ESTARS[1];
        // I stays the same (1), re-scan GUY[1]
    }
    goto L100;

// -----------------------------------------------------------------------
//  '+' SIGN
// -----------------------------------------------------------------------
L550:
    if (PREDSW) goto L9022;
    PREDSW = true;
    I = I + 1;
    goto L100;

// -----------------------------------------------------------------------
//  SUCCESS
// -----------------------------------------------------------------------
L800:
    IRTN = 2;
    if (EREDSW) IRTN = 3;
    return;

// -----------------------------------------------------------------------
//  ERROR MESSAGES
// -----------------------------------------------------------------------
L9004:
    std::printf("0**** TOO MANY CHARACTERS IN SYMBOL AND ATOMIC MASS.\n");
    goto L9900;
L9005:
    std::printf("0**** UNEXPECTED DELIMITOR OR END OF LINE ENCOUNTERED"
                " WHEN LOOKING FOR ELEMENT SYMBOL.\n");
    goto L9900;
L9007:
    std::printf("0**** COMPLETE SPECIFICATION MUST BE ON SAME LINE AS"
                " CHANNEL KEYWORD.\n");
    goto L9900;
L9020:
    std::printf("0**** TWO \"=\" SIGNS ARE NOT ALLOWED\n");
    goto L9900;
L9022:
    std::printf("0**** TWO \"+\" SIGNS ARE NOT ALLOWED.\n");
    goto L9900;

L9900:
    {
        int INCHM1 = INPBUF.INCH - 1;
        std::printf("      WAS PROCESSING FIELD%2d OF THE SPECIFICATION\n"
                    "      SPECIFICATION UP TO THE ERROR IS: ", I);
        for (II = INCHST; II <= INCHM1 && II <= IBUFSIZE; II++)
            std::printf("%c", INPBUF.IBUF[II]);
        std::printf("\n");
    }
    IRTN = -1;
    return;
}

// CLRCHN → clrchn_translated.cpp

void CTRAN(int& IRET)
{ IRET = 0; std::printf(" CTRAN: stub\n"); }

// SETCHN → setchn_translated.cpp
