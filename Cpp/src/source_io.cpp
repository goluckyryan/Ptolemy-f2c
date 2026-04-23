// source_io.cpp — I/O, parsing, input buffer handling
// Translated from source.f: NXINT+entries, MSCAN, DATAIN, GETIT, GETNUM, LSTKEY, QVSCAN, DEFINE

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include "ptolemy_alloc.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_forward.h"
#include "ptolemy_io.h"

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <cctype>

// Forward declarations for NXINT-family entry point functions
// and allocator entry points used here
extern int GIVEAL(int IWRDS);
extern int CHOPIT(int IWRDS);

// ============================================================================
// MSCAN — Field scanner for input
// ============================================================================
void MSCAN(int KEYEX, int& ITOKEN, char8& CVALU, double& VALU, int& IVALU,
           char* MESS, int& IMESS, char& STOP)
{
    // Local character constants
    const char DOLLAR = '$';
    const char BLANK = ' ';
    const char DOT = '.';
    const char STAR = '*';
    const char NUM0 = '0';
    const char NUM9 = '9';
    const char LETTRA = 'A';
    const char LETTRZ = 'Z';
    const char LOWERA = 'a';
    const char LOWERZ = 'z';
    const char LETTRE = 'E';
    const char LETTRD = 'D';

    static const char LCHARS[10] = { ' ', 'S', 'P', 'D', 'F', 'G', 'H', 'I', 'J', 'K' };
    int NUMLS = 9;

    char8 WORD, BLANK8(" ");
    char WORD1[9]; // 1-based
    int IWORD;

    static const char CHARS[22] = { ' ', '%', '#', '&', ':', ';',
        '-', '+', '/', '*', '(',
        ')', '=', '@', '<', '>',
        '|', '~', '?', '"', '.',
        ',' };

    static const char DELIMS[5] = { ' ', '+', '-', '\'', '$' };

    char ICH;
    double VAL, POW, SIGN, P2, VALU1;
    int ICNT, IDEC, INPHAS, ICOM, ISG, IPOW;
    int INCHST, I;

    // Presets
    INCHST = INPBUF.INCH;
    SIGN = +1.0;
    VAL = 0.0;
    ICNT = 0;
    INPBUF.INCH = INPBUF.INCH - 1;
    IMESS = 0;
    ICOM = 0;
    IDEC = 0;

    // Phases: 1=nothing, 2=alpha, 3=numeric
    INPHAS = 1;

    // Get next char
L100:
    INPBUF.INCH = INPBUF.INCH + 1;
    if (INPBUF.INCH > INPBUF.NOCH) goto L710;
    ICH = INPBUF.IBUF[INPBUF.INCH];

    // Is a literal coming in
    if (IMESS > 0) goto L149;

    // Is a comment coming in
    if (ICOM == 0) goto L101;
    if (ICH != DOLLAR) goto L100;
    ICOM = 0;
    goto L100;

    // Blank?
L101:
    if (ICH == BLANK) {
        switch (INPHAS) {
            case 1: goto L100;
            case 2: goto L210;
            case 3: goto L340;
        }
    }

    // Number?
    if (ICH >= NUM0 && ICH <= NUM9) goto L300;

    // Alphabetic?
    if ((ICH >= LETTRA && ICH <= LETTRZ) || ICH == STAR || ICH == '_') goto L200;

    // Convert lowercase to uppercase
    if (ICH < LOWERA || ICH > LOWERZ) goto L102;
    ICH = (char)(ICH + 'A' - 'a');
    goto L200;

    // Is it a period
L102:
    if (ICH == DOT) goto L320;

    // All other special characters stop a field in progress
    switch (INPHAS) {
        case 1: goto L105;
        case 2: goto L210;
        case 3: goto L340;
    }

    // Check for + - ' $
L105:
    for (I = 1; I <= 4; I++) {
        if (ICH == DELIMS[I]) goto L140;
    }

    // Look for tokens
    for (ITOKEN = 1; ITOKEN <= 21; ITOKEN++) {
        if (ICH == CHARS[ITOKEN]) goto L590;
    }
    goto L9000;

L140:
    switch (I) {
        case 1: goto L141;
        case 2: goto L142;
        case 3: goto L143;
        case 4: goto L144;
    }
L141:
    SIGN = +1.0;
    ITOKEN = 7;
    if (KEYEX < 0) goto L590;
    goto L302;
L142:
    SIGN = -1.0;
    ITOKEN = 6;
    if (KEYEX < 0) goto L590;
    goto L302;
L143:
    IMESS = 1;
    goto L100;
L144:
    ICOM = 1;
    goto L100;

    // Collecting a literal
L149:
    MESS[IMESS] = ICH;
    IMESS = IMESS + 1;
    if (IMESS > 80) goto L9006;
    if (ICH != DELIMS[3]) goto L100;
    IMESS = IMESS - 2;
    MESS[IMESS + 1] = BLANK;
    ITOKEN = -4;
    INPBUF.INCH = INPBUF.INCH + 1;
    goto L602;

    // Alphabetic character check phase
L200:
    switch (INPHAS) {
        case 1: goto L201;
        case 2: goto L202;
        case 3: goto L330;
    }
    // Start
L201:
    WORD = BLANK8;
    std::memcpy(WORD1 + 1, WORD.data, 8);
    IWORD = 0;
L202:
    if (IWORD < 8) { WORD1[IWORD + 1] = ICH; IWORD = IWORD + 1; }
    INPHAS = 2;
    goto L100;
    // Done
L210:
    IMESS = IWORD;
    std::memset(WORD.data, ' ', 8);
    for (int k = 1; k <= IWORD && k <= 8; k++) WORD.data[k-1] = WORD1[k];
    CVALU = WORD;
    ITOKEN = -1;
    goto L600;

    // Numeric
L300:
    if (KEYEX < 0) goto L200;
    switch (INPHAS) {
        case 1: goto L301;
        case 2: goto L202;
        case 3: goto L303;
    }
    // Start
L301:
    IDEC = 0;
    VAL = (double)(ICH - NUM0);
    ICNT = 1;
L302:
    POW = 1.0;
    INPHAS = 3;
    goto L100;
    // More
L303:
    VAL = 10.0 * VAL + (double)(ICH - NUM0);
    if (IDEC != 0) POW = POW * 10.0;
    ICNT = ICNT + 1;
    goto L100;

    // Decimal pt.
L320:
    ITOKEN = 20;
    if (KEYEX < 0) goto L590;
    switch (INPHAS) {
        case 1: goto L321;
        case 2: goto L9002;
        case 3: goto L310;
    }

    // Decimal pt encountered in numeric field
L310:
    if (IDEC != 0) goto L9002;
    IDEC = 1;
    goto L100;
L321:
    IDEC = 1;
    VAL = 0.0;
    goto L302;

    // Alphabetic in numeric field - check E for exponents
L330:
    if (ICH == LETTRE) goto L332;
    if (IDEC == 0) goto L400;
    if (ICH != LETTRD) goto L9012;

    // E or D format
L332:
    ISG = 1;
    IDEC = 1;
    IPOW = 0;
    INPBUF.INCH = INPBUF.INCH + 1;
    if (INPBUF.INCH > INPBUF.NOCH) goto L9008;
    ICH = INPBUF.IBUF[INPBUF.INCH];
    if (ICH == DELIMS[1]) goto L334;
    if (ICH != DELIMS[2]) goto L335;
    ISG = -1;
L334:
    INPBUF.INCH = INPBUF.INCH + 1;
    if (INPBUF.INCH > INPBUF.NOCH) goto L3360;
    ICH = INPBUF.IBUF[INPBUF.INCH];
L335:
    if (ICH < NUM0 || ICH > NUM9) goto L336;
    IPOW = IPOW * 10 + (ICH - NUM0);
    if (IPOW > 300) goto L9004;
    goto L334;

L3360:
    ;
L336:
    if (IPOW == 0) goto L340;
    P2 = std::pow(10.0, (double)IPOW);
    if (ISG < 0) goto L337;
    POW = POW / P2;
    goto L340;
L337:
    POW = POW * P2;
    goto L340;
L340:
    ITOKEN = -2;
    VALU1 = VAL / POW;
    VALU = VALU1 * SIGN;
    IVALU = (int)VALU;
    if (IDEC == 0) ITOKEN = -5;
    if (ICNT == 0) goto L9010;
    goto L600;

    // Alphabetic found in integer numeric field - level specification
L400:
    VALU = SIGN * VAL - 1.0;
    for (I = 1; I <= NUMLS; I++) {
        if (ICH == LCHARS[I]) goto L420;
    }
    goto L9012;
L420:
    IVALU = I - 1;
    ITOKEN = -3;
    INPBUF.INCH = INPBUF.INCH + 1;
    goto L600;

    // Special character - should we skip over it?
L590:
    if (ITOKEN > IABS(KEYEX)) goto L100;
    INPBUF.INCH = INPBUF.INCH + 1;

    // Field has been found and processed
L600:
    ;
L602:
    if (INPBUF.INCH >= 1 && INPBUF.INCH <= IBUFSIZE)
        STOP = INPBUF.IBUF[INPBUF.INCH];
    else
        STOP = ' ';
    return;

    // End of line encountered
L710:
    switch (INPHAS) {
        case 1: goto L701;
        case 2: goto L210;
        case 3: goto L340;
    }
L701:
    ITOKEN = -8;
    return;

    // Error messages
L9000:
    std::printf("\n**** UNRECOGNIZABLE CHARACTER '%c' (HEX IS %02X).\n", ICH, (unsigned char)ICH);
    goto L9019;
L9002:
    std::printf("\n**** INVALID DECIMAL POINT.\n");
    goto L9019;
L9004:
    std::printf("\n**** INVALID OR TOO LARGE EXPONENT FIELD.\n");
    goto L9019;
L9006:
    std::printf("\n**** FIELD IS TOO LONG.\n");
    goto L9019;
L9008:
    std::printf("\n**** E OR D MUST BE FOLLOWED BY NUMBER\n");
    goto L9019;
L9010:
    std::printf("\n**** MISPLACED DECIMAL POINT OR PLUS OR MINUS SIGN.\n");
    goto L9019;
L9012:
    std::printf("\n**** INVALID ALPHABETIC IN NUMERIC FIELD\n");

L9019:
    std::printf(" FIELD UP TO THE ERROR IS: ");
    for (int ii = INCHST; ii <= INPBUF.INCH && ii <= IBUFSIZE; ii++)
        std::printf("%c", INPBUF.IBUF[ii]);
    std::printf("\n");
    ITOKEN = -9;
    return;
}


// ============================================================================
// NXINT + entries (NXVAL, NXVALF, NXHINT, NXWORD, FITKEY, NEWCD)
// Uses enum dispatch pattern
// ============================================================================

enum class NxintEntry {
    NXINT_main,   // KEY=2
    NXVAL,        // KEY=3
    NXVALF,       // KEY=6
    NXHINT,       // KEY=5
    NXWORD,       // KEY=8
    FITKEY,       // KEY=1
    NEWCD         // special
};

// Internal implementation
// Returns: 0=normal, 1/2/3/4=alternate returns
// For NXINT: INTARG is output int
// For NXVAL/NXVALF: VALARG is output double
// For NXWORD/FITKEY: CVARG is output char8
static int NXINT_impl(NxintEntry entry, int& INTARG, double& VALARG, char8& CVARG)
{
    // SAVE variables
    static int KEYEX_save = 11;

    int KEY, ICD, INT_val, IMESS, INCHST, IINT;
    double VAL;
    char8 CVALU;
    char STOP;
    char IBUFN[2]; // IBUFN(1)

    // NEWCD special case
    if (entry == NxintEntry::NEWCD) {
        INPBUF.NOCH = 0;
        INPBUF.INCH = INPBUF.NOCH + 1;
        return 0;
    }

    // Set KEY based on entry point
    int KEYEX = KEYEX_save;
    switch (entry) {
        case NxintEntry::NXINT_main: KEY = 2; break;
        case NxintEntry::NXVAL:      KEY = 3; break;
        case NxintEntry::NXVALF:     KEY = 6; KEYEX = 12; break;
        case NxintEntry::NXHINT:     KEY = 5; break;
        case NxintEntry::NXWORD:     KEY = 8; break;
        case NxintEntry::FITKEY:     KEY = 1; KEYEX = 12; break;
        default: KEY = 2; break;
    }

    // Main scanner loop
L50:
    INCHST = INPBUF.INCH;
    MSCAN(KEYEX, ICD, CVALU, VAL, INT_val, IBUFN, IMESS, STOP);
    if (ICD != -8) goto L100;

    // Time to read a new card
    {
        char line[201];
        if (std::fgets(line, sizeof(line), stdin) == nullptr) {
            // EOF
            return 1;
        }
        int len = (int)std::strlen(line);
        if (len > 0 && line[len-1] == '\n') { line[--len] = '\0'; }
        if (len > IBUFSIZE) len = IBUFSIZE;
        for (int i = 1; i <= len; i++) INPBUF.IBUF[i] = line[i-1];
        for (int i = len + 1; i <= IBUFSIZE; i++) INPBUF.IBUF[i] = ' ';
        // Find NOCH (last non-blank)
        for (INPBUF.NOCH = IBUFSIZE; INPBUF.NOCH >= 1; INPBUF.NOCH--) {
            if (INPBUF.IBUF[INPBUF.NOCH] != ' ') break;
        }
        if (INPBUF.NOCH < 1) INPBUF.NOCH = 1;
    }
    INPBUF.INCH = 1;
    if (SWITCH.IECHO != 0) {
        std::printf("0INPUT... ");
        for (int i = 1; i <= INPBUF.NOCH; i++) std::printf("%c", INPBUF.IBUF[i]);
        std::printf("\n");
    }
    goto L50;

L100:
    if (ICD == 0) goto L50;
    // Computed GOTO based on KEY: 1=810, 2=200, 3=300, 4=550, 5=500, 6=310, 7=2000, 8=800, 9=2000, 10=860
    switch (KEY) {
        case 1: goto L810;
        case 2: goto L200;
        case 3: goto L300;
        case 4: goto L550;
        case 5: goto L500;
        case 6: goto L310;
        case 7: goto L2000;
        case 8: goto L800;
        case 9: goto L2000;
        case 10: goto L860;
    }

    // NXINT: integer requested
L200:
    if (ICD != -2 && ICD != -5) goto L2000;
    if (VAL != (double)INT_val) goto L2000;
    if (STOP == '/') goto L2000;
    INTARG = INT_val;
    KEYEX_save = 11;
    return 0;

    // NXVAL
L300:
    if (ICD != -2 && ICD != -5) goto L2000;
    if (STOP == '/') goto L2000;
    VALARG = VAL;
    KEYEX_save = 11;
    return 0;

    // NXVALF (KEY=6)
L310:
    KEYEX = 11;
    KEYEX_save = 11;
    if (ICD == 11) return 4;
    goto L300;

    // NXHINT: J or JP values
L500:
    if (ICD != -5) goto L2000;
    INTARG = INT_val * 2;
    if (STOP != '/') goto L570;
    IINT = INTARG / 2;
    INTARG = IINT;
    INPBUF.INCH = INPBUF.INCH + 1;
    KEY = 4;
    goto L50;
L550:
    if (ICD != -5 || INT_val != 2) goto L2050;
    INTARG = IINT;
    // Skip trailing sign
L570:
    if (STOP == '+' || STOP == '-') INPBUF.INCH = INPBUF.INCH + 1;
    KEYEX_save = 11;
    return 0;

    // NXWORD
L800:
    if (ICD != -1) goto L2000;
    CVARG = CVALU;
    KEYEX_save = 11;
    return 0;

    // FITKEY
L810:
    KEYEX = 11;
    KEYEX_save = 11;
    CVARG = CVALU;
    if (ICD == -1) return 0;
    if (ICD == 10) return 1;
    if (ICD == 11) return 3;
    if (ICD == 12) return 2;
    INPBUF.INCH = INCHST;
    return 4;

L860:
    INTARG = ICD;
    KEYEX_save = 11;
    return 0;

    // Not found; rescan same field next time
L2000:
    if (ICD == 5) return 2;
    INPBUF.INCH = INCHST;
    if (ICD < 0) return 1;
    return 3;

L2050:
    std::printf("\n**** NUMBER FOLLOWING / MUST BE 2%10d\n", INT_val);
    return 1;
}

// Wrapper functions
void NXINT(int& INTARG, int& iret) {
    double dummy_v = 0.0; char8 dummy_c;
    iret = NXINT_impl(NxintEntry::NXINT_main, INTARG, dummy_v, dummy_c);
}

int NXVAL(double& VALARG) {
    int dummy_i = 0; char8 dummy_c;
    return NXINT_impl(NxintEntry::NXVAL, dummy_i, VALARG, dummy_c);
}

int NXVALF(double& VALARG) {
    int dummy_i = 0; char8 dummy_c;
    return NXINT_impl(NxintEntry::NXVALF, dummy_i, VALARG, dummy_c);
}

int NXHINT(int& INTARG) {
    double dummy_v = 0.0; char8 dummy_c;
    return NXINT_impl(NxintEntry::NXHINT, INTARG, dummy_v, dummy_c);
}

int NXWORD(char* CVARG) {
    int dummy_i = 0; double dummy_v = 0.0;
    char8 cv;
    int ret = NXINT_impl(NxintEntry::NXWORD, dummy_i, dummy_v, cv);
    std::memcpy(CVARG, cv.data, 8);
    return ret;
}

int FITKEY(char* CVARG) {
    int dummy_i = 0; double dummy_v = 0.0;
    char8 cv;
    int ret = NXINT_impl(NxintEntry::FITKEY, dummy_i, dummy_v, cv);
    std::memcpy(CVARG, cv.data, 8);
    return ret;
}

void NEWCD() {
    int dummy_i = 0; double dummy_v = 0.0; char8 dummy_c;
    NXINT_impl(NxintEntry::NEWCD, dummy_i, dummy_v, dummy_c);
}


// ============================================================================
// DATAIN — data input parsing
// ============================================================================
void DATAIN(int& IRET)
{
    static const char8 KEYS[19] = {
        char8(""), // 0 unused
        char8("ELAB"), char8("ECM"), char8("ANGLE"), char8("MBERROR"), char8("PERCENTE"),
        char8("POLARIZA"), char8("SIGMA"), char8("SIGMATOR"), char8("WEIGHT"), char8("RENORMAL"),
        char8("ANGLESHI"), char8("LABANGLE"), char8("CMANGLE"), char8("CMSIGMA"), char8("LABSIGMA"),
        char8("ERROR"), char8("CHANNEL"), char8("RENORM")
    };
    int NUMKEY = 18;
    char8 NAME8("DATA");
    char8 AKEY;

    int LENPRE = 5, LENKEY = 5, LENHED = 30;

    int KEY[5] = { 0, 0, 0, 0, 0 }; // 1-based

    IRET = 0;
    int IANGLE = 0, ISIGMA = 0, IERROR = 0, IPOL = 0, ICNT = 0, IAT = 0;
    int NPTSE = 0, NPTS = 0;
    int IELAB = 0;
    double ELAB_val, WEIGH, RNORM, ASHIFT, TEMP, ANGLAS;
    int NTHISE, LENDAT;
    int IORDER;
    bool MORESW;

    int I = NAMLOC(NAME8.data);
    if (I != 0) LOCPTRS.Z[I] = -LOCPTRS.Z[I];
    int IWRDS;
    int IDATA = GIVEAL(IWRDS);
    NAMCOM.NAMES[IDATA] = NAME8;
    int LDATA = LOCPTRS.Z[IDATA];
    int LDATA4 = ALLOCS.FACFR4 * (LDATA - 1) + 1;
    ALLOC(LDATA) = 0.0;
    ILLOC(LDATA4 + 2) = LENPRE;
    ILLOC(LDATA4 + 3) = LENKEY;
    ILLOC(LDATA4 + 4) = LENHED;
    int LEKEY = LDATA4 + LENPRE - 1;

    { char8 ak; int r = FITKEY(ak.data);
      if (r == 1) goto L90;
      goto L900;
    }

L90:
    IORDER = 0;
    { int LEHED = LEKEY + LENKEY + 1;
      int LADATA = LEHED + LENHED;
      ALLOC4(LEHED + 1) = 1.0f;
      ALLOC4(LEHED + 2) = 1.0f;
      ALLOC4(LEHED + 3) = 0.0f;
      IELAB = 0;

L100:
      { int r = NXWORD(AKEY.data);
        if (r == 1) goto L800;
        if (r == 2 || r == 3) goto L900;
      }
      ANGLAS = 0.0;

      // Match keyword
      for (I = 1; I <= NUMKEY; I++) {
          if (KEYS[I] == AKEY) {
              switch (I) {
                  case 1: goto L200; case 2: goto L210; case 3: goto L250;
                  case 4: goto L300; case 5: goto L310; case 6: goto L400;
                  case 7: goto L450; case 8: goto L470; case 9: goto L500;
                  case 10: goto L510; case 11: goto L520; case 12: goto L260;
                  case 13: goto L265; case 14: goto L460; case 15: goto L465;
                  case 16: goto L320; case 17: goto L550; case 18: goto L510;
              }
          }
      }
      goto L900;

      // ELAB
L200: IELAB = 0; goto L220;
L210: IELAB = 1;
L220: { int r = NXVAL(ELAB_val); if (r != 0) goto L900; }
      goto L100;

      // ANGLE
L250: KEY[4] = 0;
L255: IORDER++; IANGLE = IORDER; goto L100;
L260: KEY[4] = 1; goto L255;
L265: KEY[4] = 0; goto L255;

      // MBERROR / % ERROR / ERROR
L300: KEY[1] = 0;
L305: IORDER++; IERROR = IORDER; goto L100;
L310: KEY[1] = 1; goto L305;
L320: KEY[1] = 2; goto L305;

      // POLARIZATION
L400: IORDER++; IPOL = IORDER; KEY[3] = 1; goto L100;

      // SIGMA / CMSIGMA / LABSIGMA / SIGMATOR
L450: KEY[2] = -1;
L455: IORDER++; ISIGMA = IORDER; goto L100;
L460: KEY[2] = 1; goto L455;
L465: KEY[2] = 2; goto L455;
L470: KEY[2] = 0; goto L455;

      // WEIGHT
L500: { int r = NXVAL(WEIGH); if (r != 0) goto L900; }
      ALLOC4(LEHED + 1) = (float)WEIGH; goto L100;

      // RENORMALIZATION
L510: { int r = NXVAL(RNORM); if (r != 0) goto L900; }
      ALLOC4(LEHED + 2) = (float)RNORM; goto L100;

      // ANGLESHIFT
L520: { int r = NXVAL(ASHIFT); if (r != 0) goto L900; }
      ALLOC4(LEHED + 3) = (float)ASHIFT; goto L100;

      // CHANNEL
L550: CHANEL(IRET);
      if (IRET < 0) goto L100;
      if (IRET == 0) goto L950;
      goto L100;

      // Number encountered
L800: IAT = IAT + 4;
      if (IERROR == 0 || IANGLE == 0 || ISIGMA == 0) goto L930;
      if (KEY[2] == -1) KEY[2] = KEY[4] + 1;
      NTHISE = 0;
      LENDAT = 3;
      if (IPOL != 0) LENDAT = 4;

L850: for (I = 1; I <= LENDAT; I++) {
          { int r = NXVALF(TEMP);
            if (r == 1) goto L880;
            if (r == 2 || r == 3) goto L910;
            if (r == 4) goto L885;
          }
          if (I == IANGLE) goto L862;
          if (I == ISIGMA) goto L864;
          if (I == IERROR) goto L866;
          if (I == IPOL) goto L868;
          goto L945;

L862:     if (TEMP < 0.0 || TEMP > 180.0) goto L940;
          if (TEMP < 1.0)
              std::printf(" **** WARNING:  ANGLE (%13.5g DEGREES) IS UNUSUALLY SMALL.\n", TEMP);
          if (TEMP < ANGLAS)
              std::printf(" **** WARNING:  ANGLE (%6.2f DEGREES) IS LESS THAN PREVIOUS ANGLE.\n", TEMP);
          ANGLAS = TEMP;
          ALLOC4(LADATA) = (float)TEMP;
          NTHISE++;
          NPTS++;
          goto L860;
L864:     ALLOC4(LADATA + 2) = (float)TEMP; goto L860;
L866:     ALLOC4(LADATA + 1) = (float)TEMP; goto L860;
L868:     ALLOC4(LADATA + 3) = (float)TEMP; goto L860;
L860:     ;
      }
      LADATA = LADATA + LENDAT;
      goto L850;

      // Finished this group -- more
L880: MORESW = true; goto L890;
      // ) encountered
L885: MORESW = false;

L890: if (I != 1) goto L920;
      for (I = 1; I <= 4; I++) ILLOC(LEKEY + I) = KEY[I];
      ILLOC(LEKEY + 5) = LENDAT;
      ILLOC(LEHED) = NTHISE;
      ALLOC4(LEHED + 4) = (float)ELAB_val;
      ILLOC(LEHED + 5) = INTGER.IZP;
      ILLOC(LEHED + 6) = INTGER.IZT;
      ALLOC4(LEHED + 7) = (float)FLOAT_common.AMP;
      ALLOC4(LEHED + 8) = (float)FLOAT_common.AMT;
      {
          int* REACT4 = reinterpret_cast<int*>(HEDCOM.REACT);
          for (I = 9; I <= 19; I++) ILLOC(LEHED + I) = REACT4[I - 8];
      }
      ILLOC(LEHED + 20) = (int)JBLOCK.JSP;
      ILLOC(LEHED + 21) = (int)JBLOCK.JST;
      ILLOC(LEHED + 22) = 0;
      ILLOC(LEHED + 23) = 0;
      if (FLOAT_common.EXSPT[1] != INTRNL.UNDEF) {
          ILLOC(LEHED + 22) = 1;
          ILLOC(LEHED + 23) = 3;
          if (FLOAT_common.AMP != FLOAT_common.AMT || INTGER.IZP != INTGER.IZT
              || JBLOCK.JSP != JBLOCK.JST || FLOAT_common.EXSPT[1] != FLOAT_common.EXSPT[2])
              goto L897;
          if ((int)JBLOCK.JSP == 0) ILLOC(LEHED + 22) = 2;
          ILLOC(LEHED + 23) = ((int)JBLOCK.JSP % 2) + 1;
      }
L897: ALLOC4(LEHED + 25) = (float)FLOAT_common.AMXGPT[1];
      ALLOC4(LEHED + 26) = (float)FLOAT_common.AMXGPT[2];
      ALLOC4(LEHED + 27) = (float)FLOAT_common.EXSPT[1];
      ALLOC4(LEHED + 28) = (float)FLOAT_common.EXSPT[2];
      ILLOC(LEHED + 29) = IELAB;
      NPTSE++;
      LEKEY = LEKEY + LENKEY + LENHED + LENDAT * NTHISE;
      if (MORESW) goto L90;

      ILLOC(LDATA4) = NPTS;
      ILLOC(LDATA4 + 1) = NPTSE;
      IWRDS = (LEKEY - LDATA4) / ALLOCS.FACFR4 + 1;
      { int il = CHOPIT(IWRDS); (void)il; }
      IRET = 1;
      return;
    } // end scope

L900:
    std::printf(" DATA NOT ENCLOSED IN PARENTHESIS OR KEYWORD MISSPELLED\n");
    std::printf(" DATA KEYWORDS ARE:\n");
    for (int i = 1; i <= NUMKEY; i++)
        std::printf(" %.8s", KEYS[i].data);
    std::printf("\n");
    return;
L910:
    std::printf(" ERROR IN READING NUMERICAL DATA VALUES\n");
    return;
L920:
    std::printf(" INPUT ERROR - NUMBER OF INPUT NUMBERS MUST BE %5d*N\n", LENDAT);
    return;
L930:
    std::printf(" ORDER OF ANGLE,ERROR AND SIGMA MUST BE SPECIFIED\n");
    return;
L940:
    std::printf("\n**** INVALID ANGLE (%13.5g DEGREES); CHECK ORDER OF INPUT.\n", TEMP);
    return;
L945:
    std::printf("\n*** A DATA-ENTERING KEYWORD HAS BEEN SPECIFIED TWICE OR AN INTERNAL ERROR OCCURRED: %5d%5d%5d%5d%5d\n",
        I, IANGLE, ISIGMA, IERROR, IPOL);
L950:
    return;
}


// ============================================================================
// GETIT — get/read named Speakeasy object
// ============================================================================
void GETIT(int& IRET)
{
    static const int BYTES[10] = { 0, 4, 8, 8, 16, 4, 8, 2, 1, 1 }; // 1-based
    static const char8 KINDS[10] = {
        char8(""), char8("INTEGER"), char8("REAL*8"), char8("REAL*8"), char8("COMPLEX"),
        char8("REAL*4"), char8("NAME LIT"), char8("INTGER*2"), char8("INTGER*1"), char8("CHAR LIT")
    };
    char8 MEMBER;
    char8 MYKEEP("MYKEEP");
    char8 OBNAME;

    { int r = NXWORD(OBNAME.data); if (r != 0) goto L900; }
    { int r = NXWORD(MEMBER.data); if (r != 0) goto L900; }

    std::printf(" GET/KEPT needs to be coded for unix!!!!\n");
    goto L950;

    // ... (rest of GETIT is dead code in the unix version)

L900:
    std::printf("\n**** GET REQUIRES TWO NAMES AS OPERANDS.\n");
L950:
    IRET = 0;
    return;
}


// ============================================================================
// GETNUM — get number of an object from input
// Returns: 0=normal, 1=alternate return (error)
// ============================================================================
int GETNUM(int& NUM)
{
    char8 NAME;

    // Try to get an integer
    { int r; double dv; char8 dc;
      r = NXINT_impl(NxintEntry::NXINT_main, NUM, dv, dc);
      if (r == 0) goto L_found;
      if (r == 2) goto L900;
      // r==1 or r==3: try by name
      goto L200;
    }

L_found:
    // Number found, check if defined
    if (NUM > USAGE.IMAX || NUM <= 0) goto L150;
    if (LOCPTRS.Z[NUM] > 0 && LENGTH.LENG[NUM] > 0) return 0;
L150:
    std::printf("\n***** NUMBERED OBJECT%8d IS NOT DEFINED *****\n\n", NUM);
    NUM = 0;
    return 1;

    // Try to find by name
L200:
    { int r = NXWORD(NAME.data); if (r != 0) goto L900; }
    NUM = NAMLOC(NAME.data);
    if (NUM != 0) return 0;
    std::printf("\n******  %.8s IS NOT A DEFINED NAME. ******\n\n", NAME.data);
    NUM = 0;
    return 1;

L900:
    std::printf("\n***** INPUT ERROR - EXPECTING NAME OR NUMBER *****\n");
    NUM = 0;
    return 1;
}


// ============================================================================
// DEFINE — define a named object from input values
// ============================================================================
void DEFINE(char8 ANAME, int& ILEN, int& ISIZE)
{
    // Get the object number if it presently exists
    int IDOLD = NAMLOC(ANAME.data);

    // Get maximum amount possible from allocator
    int IDATA = GIVEAL(ISIZE);
    NAMCOM.NAMES[IDATA] = ANAME;
    int LDATA = LOCPTRS.Z[IDATA];
    ILEN = 0;
    double X;

L50:
    { int r = NXVAL(X);
      if (r == 1) goto L310;
      if (r == 2) goto L300;
      if (r == 3) goto L310;
    }
    ILEN++;
    if (ILEN > ISIZE + 1) goto L50;
    if (ILEN <= ISIZE) goto L100;
    std::printf("\n**** EXCEEDED ALLOCATOR WHEN DEFINING %.8s ROOM FOR ONLY%6d VALUES.\n", ANAME.data, ILEN);
    goto L50;

L100:
    ALLOC(LDATA - 1 + ILEN) = X;
    goto L50;

    // End of list -- backup for semicolon
L300:
    INPBUF.INCH = INPBUF.INCH - 1;
L310:
    { int ixxxx = CHOPIT(MIN0(ILEN, ISIZE)); (void)ixxxx; }

    // If object previously existed, assign it the same number
    if (IDOLD == 0) return;
    int LOLD = LOCPTRS.Z[IDOLD];
    LOCPTRS.Z[IDOLD] = LDATA;
    LOCPTRS.Z[IDATA] = -LOLD;
    X = ALLOC(LDATA - 1);
    ALLOC(LDATA - 1) = ALLOC(LOLD - 1);
    ALLOC(LOLD - 1) = X;
    LENGTH.LENG[IDATA] = LENGTH.LENG[IDOLD];
    LENGTH.LENG[IDOLD] = MIN0(ILEN, ISIZE);
    return;
}


// ============================================================================
// QVSCAN — Q-value / reaction scan
// ============================================================================
void QVSCAN(char8* GUY, double* ESTARS, int* NODVLS, int* LVALS,
            int* JVALS, int* IPARYS, int& IRTN)
{
    static const char STOPS[5] = { ' ', '(', ',', ')', ' ' }; // 1-based
    const char SLASH = '/';
    const char PLUS = '+';
    const char MINUS = '-';
    const char BLANK = ' ';
    const char COMMA = ',';
    char8 BLANK8(" ");
    char8 CVALU;
    double VALU, DUMMY;
    int ITOKEN, IVALU, IMESS;
    char STPCHR;
    int STPTYP = 0;
    int IP;

    IRTN = 0;
    int INCHST = INPBUF.INCH;

    // Loop on the reaction participants
    for (int I = 1; I <= 4; I++) {
        GUY[I] = BLANK8;
        NODVLS[I] = NOTDEF_INT;
        LVALS[I] = NOTDEF_INT;
        JVALS[I] = NOTDEF_INT;
        ESTARS[I] = INTRNL.UNDEF;
        IPARYS[I] = 0;

        MSCAN(-20, ITOKEN, GUY[I], VALU, IVALU, INPBUF.IBUF, IMESS, STPCHR);
        if (IMESS > 5) goto L9004;
        if (IMESS == 0) goto L9005;
        INPBUF.INCH = INPBUF.INCH + 1;

        // If stopped by ( and not first field, process excited state
        if (STPCHR != STOPS[1] || I == 1) goto L400;

        // Get E* and J
L200:   MSCAN(18, ITOKEN, CVALU, VALU, IVALU, INPBUF.IBUF, IMESS, STPCHR);
        if (ITOKEN == -9) goto L9900;
        if (ITOKEN == -8) goto L9007;
        if (ITOKEN != -2) goto L230;

        // It is E*
        if (ESTARS[I] != INTRNL.UNDEF) {
            std::printf("\n**** ATTEMPTING TO ENTER 2 VALUES FOR E*\n");
            goto L9900;
        }
        ESTARS[I] = VALU;
        goto L200;

L230:   if (ITOKEN != -5) goto L270;
        // It is J
        if (JVALS[I] != NOTDEF_INT) {
            std::printf("\n**** ATTEMPTING TO ENTER 2 VALUES FOR J\n");
            goto L9900;
        }
        JVALS[I] = 2 * IVALU;
        if (STPCHR != SLASH) goto L250;
        // Bypass /2
        INPBUF.INCH = INPBUF.INCH + 1;
        JVALS[I] = IVALU;
        MSCAN(18, ITOKEN, CVALU, VALU, IVALU, INPBUF.IBUF, IMESS, STPCHR);
        if (ITOKEN == -9) goto L9900;
        if (ITOKEN == -8) goto L9007;
        if (ITOKEN == -5 && IVALU == 2) goto L250;
        std::printf("\n**** INVALID J FIELD\n");
        goto L9900;

        // Bypass terminal + or -
L250:   IP = +1;
        if (STPCHR == PLUS) goto L255;
        IP = -1;
        if (STPCHR != MINUS) goto L200;
L255:   IPARYS[I] = IP;
        INPBUF.INCH = INPBUF.INCH + 1;
        goto L200;

L270:   if (ITOKEN != -3) goto L280;
        // Level descriptor
        if (LVALS[I] != NOTDEF_INT) {
            std::printf("\n**** ATTEMPTING TO ENTER TWO LEVEL DESCRIPTORS\n");
            goto L9900;
        }
        LVALS[I] = IVALU;
        NODVLS[I] = (int)VALU;
        goto L200;

L280:   if (ITOKEN == 11) goto L290;
        std::printf("\n**** UNRECOGNIZABLE EXCITED STATE DESCRIPTOR (TOKEN =%4d)\n", ITOKEN);
        goto L9900;

L290:   STPCHR = INPBUF.IBUF[INPBUF.INCH];
        INPBUF.INCH = INPBUF.INCH + 1;

        // Bypass blanks
L400:   if (STPCHR != BLANK || I == 2 || I == 4) goto L500;
        MSCAN(-20, ITOKEN, CVALU, DUMMY, IVALU, INPBUF.IBUF, IMESS, STPCHR);
        STPCHR = INPBUF.IBUF[INPBUF.INCH - 1];
        if (IMESS != 0) goto L9006;

L500:   if ((I == 1 || I == 3) && STPCHR != STOPS[I]) goto L9001;
        STPTYP = 2;
        if (STPCHR == BLANK || STPCHR == COMMA) STPTYP = 1;
        if (I == 2 && STPTYP != 1) goto L9002;
    }

    if (STPTYP != 1 && INPBUF.INCH <= INPBUF.NOCH) goto L9003;
    return;

L9001: std::printf("\n**** MISSING PARENTHESIS.\n"); goto L9900;
L9002: std::printf("\n**** INVALID DELIMITOR AFTER INCOMING PROJECTILE.\n"); goto L9900;
L9003: std::printf("\n**** ONLY 4 PARTICIPANTS IN REACTION ALLOWED.\n"); goto L9900;
L9004: std::printf("\n**** TOO MANY CHARACTERS IN SYMBOL AND ATOMIC MASS.\n"); goto L9900;
L9005: std::printf("\n**** UNEXPECTED DELIMITOR OR END OF LINE ENCOUNTERED WHEN LOOKING FOR ELEMENT SYMBOL.\n"); goto L9900;
L9006: std::printf("\n**** UNEXPECTED FIELD ENCOUNTERED WHEN LOOKING FOR DELIMITOR\n"); goto L9900;
L9007: std::printf("\n**** COMPLETE REACTION MUST BE ON SAME LINE AS REACTION KEYWORD.\n"); goto L9900;

L9900:
    {
        int INCHM1 = INPBUF.INCH - 1;
        std::printf("      WAS PROCESSING FIELD%2d OF THE REACTION\n", 0);
        std::printf("      REACTION UP TO THE ERROR IS: ");
        for (int II = INCHST; II <= INCHM1 && II <= IBUFSIZE; II++)
            std::printf("%c", INPBUF.IBUF[II]);
        std::printf("\n");
    }
    IRTN = -1;
    return;
}


// ============================================================================
// LSTKEY — list keywords and their values
// ============================================================================
void LSTKEY(int VALSW, int NUMFLT, char8* FLTWRD, int NUMINT, char8* INTWRD,
            int NUMJ, char8* JWRDS, int NUMSWW, char8* SWTWRD, int* SWTNUM,
            int* SWTSET, int NUMKEY_arg, char8* KEYWRD, int NUMALI, char8* ALIAI)
{
    char8 NOMORE("*NOMORE*");
    char8 BLANK8(" ");
    char8 NAMES[251]; // 1-based
    char8 SYNS[251];  // 1-based
    int INDEX[251];   // 1-based
    char VALSTR[251][17]; // values as strings, 1-based

    int K, N, I;

    for (I = 1; I <= 250; I++) {
        NAMES[I] = BLANK8;
        std::memset(VALSTR[I], ' ', 16); VALSTR[I][16] = '\0';
        INDEX[I] = I;
    }

    // Collect float keywords
    K = 0;
    for (I = 1; I <= NUMFLT; I++) {
        if (FLTWRD[I] == NOMORE) continue;
        K++;
        NAMES[K] = FLTWRD[I];
        if (!VALSW) continue;
        double fval = FLOAT_arr(I);
        if (fval == INTRNL.UNDEF) {
            std::snprintf(VALSTR[K], 17, "  NOT DEFINED   ");
            continue;
        }
        double X = DABS(fval);
        if (X == 0.0)
            std::snprintf(VALSTR[K], 17, "%10.0f", fval);
        else if (X < 0.01 || X > 1.0e6)
            std::snprintf(VALSTR[K], 17, "%14.6e", fval);
        else if (X + 1.0e-6 < 1.0)
            std::snprintf(VALSTR[K], 17, "%10.6f", fval);
        else if (std::fmod(X + 1.0e-8, 1.0) < 1.0e-7)
            std::snprintf(VALSTR[K], 17, "%10.0f", fval);
        else if (X < 100.0)
            std::snprintf(VALSTR[K], 17, "%10.5f", fval);
        else if (X < 10000.0)
            std::snprintf(VALSTR[K], 17, "%10.3f", fval);
        else
            std::snprintf(VALSTR[K], 17, "%10.1f", fval);
    }

    // Collect integer keywords
    int* INTGER_arr_ptr = &INTGER.L;
    for (I = 1; I <= NUMINT; I++) {
        if (INTWRD[I] == NOMORE) continue;
        K++;
        NAMES[K] = INTWRD[I];
        if (!VALSW) continue;
        int ival = INTGER_arr_f(I);
        if (ival == NOTDEF_INT) {
            std::snprintf(VALSTR[K], 17, "  NOT DEFINED   ");
            continue;
        }
        std::snprintf(VALSTR[K], 17, "%10d", ival);
    }

    // Collect J-block keywords
    for (I = 1; I <= NUMJ; I++) {
        if (JWRDS[I] == NOMORE) continue;
        K++;
        NAMES[K] = JWRDS[I];
        if (!VALSW) continue;
        double jval_d = JBLOCK_arr_f(I);
        if (jval_d == INTRNL.NOTDEF) {
            std::snprintf(VALSTR[K], 17, "  NOT DEFINED   ");
            continue;
        }
        int JVAL = (int)jval_d;
        if (JVAL % 2 != 0)
            std::snprintf(VALSTR[K], 17, "%8d/2", JVAL);
        else
            std::snprintf(VALSTR[K], 17, "%10d", JVAL / 2);
    }

    // Sort and print data-entering keywords
    N = K;
    ALSORT(NAMES, INDEX, N);
    if (!VALSW) {
        std::printf("\n DATA ENTERING KEYWORDS ARE:\n\n");
        for (I = 1; I <= N; I++) {
            std::printf(" %.8s  ", NAMES[INDEX[I]].data);
            if (I % 6 == 0) std::printf("\n");
        }
        if (N % 6 != 0) std::printf("\n");
    } else {
        std::printf("\n DATA ENTERING KEYWORDS AND THEIR VALUES ARE:\n\n");
        for (I = 1; I <= N; I++) {
            std::printf(" %.8s %.16s  ", NAMES[INDEX[I]].data, VALSTR[INDEX[I]]);
            if (I % 3 == 0) std::printf("\n");
        }
        if (N % 3 != 0) std::printf("\n");
    }

    // Switch settings
    K = 0;
    int* SWITCH_arr_ptr = &SWITCH.IASYMP;
    for (I = 1; I <= NUMSWW; I++) {
        if (SWTWRD[I] == NOMORE) continue;
        if (VALSW) {
            if (SWITCH_arr_ptr[SWTNUM[I] - 1] != SWTSET[I]) continue;
        }
        K++;
        NAMES[K] = SWTWRD[I];
    }
    for (I = 1; I <= K; I++) INDEX[I] = I;
    ALSORT(NAMES, INDEX, K);
    if (!VALSW) {
        std::printf("\n THE FOLLOWING KEYWORDS MAY BE USED TO SELECT OPTIONS:\n\n");
    } else {
        std::printf("\n THE FOLLOWING OPTIONS HAVE BEEN CHOSEN:\n\n");
    }
    for (I = 1; I <= K; I++) {
        std::printf(" %.8s  ", NAMES[INDEX[I]].data);
        if (I % 6 == 0) std::printf("\n");
    }
    if (K > 0 && K % 6 != 0) std::printf("\n");

    if (VALSW) return;

    // Control keywords
    K = 0;
    for (I = 1; I <= NUMKEY_arg; I++) {
        if (KEYWRD[I] == NOMORE) continue;
        K++;
        NAMES[K] = KEYWRD[I];
    }
    for (I = 1; I <= K; I++) INDEX[I] = I;
    ALSORT(NAMES, INDEX, K);
    std::printf("\n THE FOLLOWING ARE COMMAND KEYWORDS:\n\n");
    for (I = 1; I <= K; I++) {
        std::printf(" %.8s  ", NAMES[INDEX[I]].data);
        if (I % 6 == 0) std::printf("\n");
    }
    if (K > 0 && K % 6 != 0) std::printf("\n");

    // Aliases
    K = 0;
    for (I = 1; I <= NUMALI; I += 2) {
        K++;
        NAMES[K] = ALIAI[I];
        SYNS[K] = ALIAI[I + 1];
    }
    for (I = 1; I <= K; I++) INDEX[I] = I;
    ALSORT(NAMES, INDEX, K);
    std::printf("\n THE FOLLOWING ARE SYNONYMS FOR STANDARD PTOLEMY KEYWORDS\n\n");
    std::printf(" SYNONYM   KEYWORD       SYNONYM   KEYWORD       SYNONYM   KEYWORD\n\n");
    for (I = 1; I <= K; I++) {
        std::printf(" %.8s  %.8s     ", NAMES[INDEX[I]].data, SYNS[INDEX[I]].data);
        if (I % 3 == 0) std::printf("\n");
    }
    if (K > 0 && K % 3 != 0) std::printf("\n");

    return;
}
