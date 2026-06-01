// contrl_part1.cpp — First half of CONTRL subroutine
// Translated from source.f lines 9349–9997
// Declarations, initialization, keyword search, and LISTKEYS/KEYWORDS handlers
//
// The function brace is NOT closed here — continues in contrl_part2.cpp.

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include "ptolemy_alloc.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_forward.h"
#include "ptolemy_io.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <algorithm>

extern double second();

void CONTRL(int& IGOTO, int& JGOTO, int& IPRM1, int& IPRM2, int& IPRM3, double& PRM1)
{
    // ======================================================================
    // Local variables (CHARACTER*8, etc.)
    // ======================================================================
    char8 KYWORD, WORD;
    char BLANK1 = ' ';

    // Part 2 variables
    auto& z = LOCPTRS.Z;
    int IALL = 0;
    int nxret = 0, VALSW_int = 0, JUMP_target = 0;

    int VALSW;   // LOGICAL

    // Aliases for COMMON block members used frequently
    double& UNDEF   = INTRNL.UNDEF;
    double& NOTDEF  = INTRNL.NOTDEF;
    int& ICHANB     = INTRNL.ICHANB;
    int& ICHANW     = INTRNL.ICHANW;
    int& ICNTT      = INTRNL.ICNTT;
    int& IDONE      = INTRNL.IDONE;
    int& ISTRIP     = INTRNL.ISTRIP;
    int& IHSAVE     = INTRNL.IHSAVE;
    int& IPOT       = INTRNL.IPOT;
    int& IWFN       = INTRNL.IWFN;
    int& IEXCIT     = INTRNL.IEXCIT;
    double& E       = FLOAT_common.E;
    double& ECM     = FLOAT_common.ECM;
    double& ELAB    = FLOAT_common.ELAB;
    double& V       = FLOAT_common.V;
    double& STEPSZ  = FLOAT_common.STEPSZ;
    double& WRITES_val = FLOAT_common.WRITES;
    double& ASYMPT  = FLOAT_common.ASYMPT;
    double& BNDASY  = FLOAT_common.BNDASY;
    int& NEXT       = SWITCH.NEXT;
    int& PROBLM     = SWITCH.PROBLM;
    int& IBSPAS     = SWITCH.IBSPAS;
    int& ITSO       = SWITCH.ITSO;
    int& IPRINT     = INTGER.IPRINT;
    int& ISIZE      = INTGER.ISIZE;
    int& LX         = INTGER.LX;
    int& MAXLF      = LOGFAC.MAXLF;
    char* REACT     = HEDCOM.REACT;
    char* HEADER    = HEDCOM.HEADER;
    float& TBEGIN   = FTIME.TBEGIN;
    int& NSPTSW     = USAGE.NSPTSW;

    // ASSIGNED GOTO variable
    static int JUMP = 0;

    // Local scalars
    int I, II, III, IDUM, NUM1, NUM2, IVAL, IRET, INEED, ITWO, LLLL;
    int ISTP1, ISTP2, IARE, ISEND, LLL, KEYRET;
    int NWRITE, NSTEPS;
    int LOCC;
    int keylen;
    double VAL, TEMP, AME, RRR, TIME;

    // ======================================================================
    // DATA keyword tables — FLOAT keywords (153 entries)
    // ======================================================================
    static const int NUMFLT = 153;
    static char8 FLTWRD[154]; // 1-based: FLTWRD[1]..FLTWRD[153]
    static char8 FLT1[101];   // 1-based
    static char8 FLT2[54];    // 1-based
    static bool data_init = false;
    if (!data_init) {
        data_init = true;

        // FLT1: indices 1..100
        const char* flt1_data[100] = {
            "A",        "ACCURACY","AI",       "ANGLESTE","ANGLEMIN",
            "ANGLEMAX", "ASO",     "ASOI",     "ASYMPTOP","DELTAVK",
            "DWCUTOFF", "E",       "ECM",      "ELAB",    "E*A",
            "E*B",      "E*BIGA",  "E*BIGB",   "E*X",     "GAMMADIF",
            "GAMMASUM", "TAU",     "TAUI",     "M",       "MA",
            "MB",       "MBIGA",   "MBIGB",    "MX",      "MXCA",
            "MXCB",     "MXCBIGA", "MXCBIGB",  "MXCX",    "MXCGA",
            "MXCGB",    "MXCGBIGA","MXCGBIGB", "MXCGX",   "MP",
            "MT",       "Q",       "R",        "R0",      "RC",
            "RC0",      "RI",      "RI0",      "RSO",     "RSO0",
            "RSOI",     "RSOI0",   "SPAM",     "SPAMP",   "SPAMT",
            "STEPSIZE", "STEP1R",  "STEP1I",   "SUMMAX",  "SUMMID",
            "SUMMIN",   "V",       "VI",       "VSO",     "VSOI",
            "FITMULTI", "WRITESTE","ASI",      "RSI",     "RSI0",
            "VSI",      "MIDMULT", "STEPSPER", "R0E",     "RI0E",
            "AE",       "AIE",     "VE",       "VIE",     "FITACCUR",
            "LMINMULT", "SPFACP",  "SPFACT",   "R0ESQ",   "RI0ESQ",
            "AESQ",     "AIESQ",   "VESQ",     "VIESQ",   "DERIVSTE",
            "FITRATIO", "RCP",     "RCT",      "RC0P",    "RC0T",
            "SUMPOINT", "INELASAC","COULOMBM", "LMAXMULT","E*P"
        };
        for (int k = 0; k < 100; k++) FLT1[k+1] = char8(flt1_data[k]);

        // FLT2: indices 1..53
        const char* flt2_data[53] = {
            "E*T",      "REALPOWE","IMAGPOWE", "BOUNDASY","SCATASYM",
            "VTR",      "VTRI",    "VTL",      "VTLI",    "VTPR",     "VTPRI",
            "PARAM1",   "PARAM2",  "PARAM3",   "PARAM4",  "PARAM5",
            "PARAM6",   "PARAM7",  "PARAM8",   "PARAM9",  "PARAM10",
            "PARAM11",  "PARAM12", "PARAM13",  "PARAM14", "PARAM15",
            "PARAM16",  "PARAM17", "PARAM18",  "PARAM19", "PARAM20",
            "RTR",      "RTRI",    "RTL",      "RTLI",    "RTPR",     "RTPRI",
            "RTR0",     "RTRI0",   "RTL0",     "RTLI0",   "RTPR0",    "RTPRI0",
            "ATR",      "ATRI",    "ATL",      "ATLI",    "ATPR",     "ATPRI",
            "MXCGP",    "MXCGT",   "PHIMID",   "GAMMAPHI"
        };
        for (int k = 0; k < 53; k++) FLT2[k+1] = char8(flt2_data[k]);

        // Build combined FLTWRD from FLT1 and FLT2
        for (int k = 1; k <= 100; k++) FLTWRD[k] = FLT1[k];
        for (int k = 1; k <= 53; k++)  FLTWRD[100+k] = FLT2[k];
    }

    // ======================================================================
    // INTEGER keywords (54 entries)
    // ======================================================================
    static const int NUMINT = 54;
    static char8 INTWRD[55]; // 1-based
    static bool int_init = false;
    if (!int_init) {
        int_init = true;
        const char* intw[54] = {
            "L",        "LMAX",     "LMIN",     "LOOKSTEP","LX",
            "MAPSUM",   "MAPDIF",   "MAXITER",  "NAITKEN", "NBACK",
            "NCOSINE",  "NODES",    "NPDIF",    "NPHIADD", "NPPHI",
            "NPSUM",    "PRINT",    "SIZE",     "ZA",      "ZB",
            "ZBIGA",    "ZBIGB",    "ZX",       "ZP",      "ZT",
            "LSTEP",    "LBACK",    "MAXLEXTR", "LMAXADD", "MAXFUNCT",
            "REINITIA", "FITMODE",  "NUMRANDO", "VPOLYORD","NPCOULOM",
            "MAXCOULI", "LMINSUB",
            "PARA",     "PARB",     "PARBIGA",  "PARBIGB", "PARX",
            "PARP",     "PART",     "PARITY",
            "IPARAM1",  "IPARAM2",  "IPARAM3",  "IPARAM4", "IPARAM5",
            "BLOCKSIZ", "MAPPHI",   "NFIRSTOF", "LEPSBACK"
        };
        for (int k = 0; k < 54; k++) INTWRD[k+1] = char8(intw[k]);
    }

    // ======================================================================
    // J-value keywords (12 entries)
    // ======================================================================
    static const int NUMJ = 12;
    static char8 JWRDS[13]; // 1-based
    static bool j_init = false;
    if (!j_init) {
        j_init = true;
        const char* jw[12] = {
            "J",  "JA", "JB", "JBIGA", "JBIGB",
            "JX", "JP", "SP", "ST",
            "JPARAM1", "JPARAM2", "JPARAM3"
        };
        for (int k = 0; k < 12; k++) JWRDS[k+1] = char8(jw[k]);
    }

    // ======================================================================
    // Switch keywords (112 entries mapping to 34 switches)
    // ======================================================================
    static const int NUMSW = 34;
    static const int NUMSWW = 112;
    static char8 SWTWRD[113]; // 1-based
    static bool sw_init = false;
    if (!sw_init) {
        sw_init = true;
        const char* sww[112] = {
            "SKIPASYM","DOASYMPT",
            "CHECKASY","NOCHECKA","ECHO",     "NOECHO",
            "ELASTIC", "NOELASTI","FITE",     "FITV",    "USEPROJE","USETARGE",
            "REALWAVE","COMPLEXW","INTERACT", "BATCH",
            "LINEAR",  "QUADRATI","NUCONLY",  "USECOULO","USESIMPC","USECORE",
            "PROJECTI","TARGET",  "INCOMING", "OUTGOING",
            "GRIDSETU","RADIALIN","LINTERPO", "CROSSSEC",
            "DOTEMP1", "DOTEMP2", "DOTEMP3",
            "POWELL65","ROCORD",  "OCOPTR",   "LMGENV",  "DERIVCHE",
            "CHECKDER","LMCHOL",  "MINIM",    "QUAVER",  "DAVIDON", "FITTER10",
            "FITTER11","FITTER12","FITTER13",
            "EPOWERS", "EPOWER",  "EINVERSE",
            "SAVEMXCR","NOSAVEMX","SAVEF",    "NOSAVEF",
            "LABANGLE","CMANGLES","BOUNDSTA", "SCATTERI",
            "NZRDWBA", "DWBA",    "FITELAST", "INELASTI","COUPLE2",
            "CCINELAS",
            "TRUEMASS","INTEGERM","FLOATNUM",
            "R0DEFAUL","R0TARGET","R0SUM",    "R0MTOT",
            "ONEPASS", "MULTIPAS","NEXTPASS",
            "SW1ON",   "SW1OFF",  "SW2ON",    "SW2OFF",  "SW3ON",   "SW3OFF",
            "SCRFILE", "NOSCRFIL","BFACTORS", "NOBFACTO",
            "BLOCKREL","RELATIVE","ABSOLUTE",
            "WOODSSAX","LPOWER",  "LPOWER2",  "WKBEXTRA",
            "LEXTRAP5","LEXTRAP6",
            "SUBBORN", "NOSUBBOR","ME1STORD", "ME2NDORD",
            "FORMDERI","FORMASIS",
            "NOWRITEW","WRITEWAV","WRITEITE", "WRITEELA","WRTELAIT",
            "GALILEO", "FITZGERA","POTNONRE", "POTRELAT",
            "NOSAVESI","SAVESITE","SAVEPITE", "SAVESPIT"
        };
        for (int k = 0; k < 112; k++) SWTWRD[k+1] = char8(sww[k]);
    }

    // SWTNUM: which switch (1..34) each keyword sets
    static const int SWTNUM[113] = { 0, // index 0 unused
        1, 1,
        2, 2, 3, 3,
        4, 4, 5, 5, 6, 6,
        7, 7, 8, 8,
        9, 9, 10, 10, 10, 10,
        11, 11, 11, 11,
        11, 11, 11, 11,
        11, 11, 11,
        12, 12, 12, 12, 12,
        12, 12, 12, 12, 12, 12,
        12, 12, 12,
        13, 13, 13,
        14, 14, 15, 15,
        16, 16,
        17, 17, 17, 17, 17, 17, 17, 17,
        18, 18, 18,
        19, 19, 19, 19,
        20, 20, 20,
        21, 21, 22, 22, 23, 23,
        24, 24, 25, 25,
        26, 26, 26,
        27, 27, 27, 27, 27, 27,
        28, 28, 29, 29,
        30, 30,
        31, 31, 31, 31, 31,
        32, 32, 33, 33,
        34, 34, 34, 34
    };

    // SWTSET: what value to set the switch to
    static const int SWTSET[113] = { 0, // index 0 unused
        1, 0,
        1, 0, 1, 0,
        1, 0, 0, 1, 1, 2,
        1, 0, 1, 0,
        1, 0, 1, 2, 0, 3,
        1, 2, 3, 4,
        10, 11, 12, 13,
        30, 31, 32,
        1, 2, 3, 4, 5,
        5, 6, 7, 8, 9, 10,
        11, 12, 13,
        0, 0, 1,
        1, 0, 1, 0,
        1, 0, 5, 6,
        20, 20, 21, 22, 23,
        24,
        0, 1, 2,
        0, 1, 2, 3,
        0, 1, 2,
        1, 0, 1, 0, 1, 0,
        1, 0, 0, 1,
        0, 2, 1,
        1, 2, 3, 4, 5, 6,
        0, 1, 0, 1,
        0, 1,
        0, 1, 2, 11, 12,
        0, 1, 0, 1,
        0, 1, 2, 3
    };

    // ======================================================================
    // Control keywords (62 entries)
    // ======================================================================
    static const int NUMKEY = 62;
    static char8 KEYWRD[63]; // 1-based
    static bool key_init = false;
    if (!key_init) {
        key_init = true;
        const char* kw[62] = {
            "LISTKEYS","RETURN",  "RENAME",  "REACTION","HEADER",
            "PRINTFIT","KEEP",    "CLEAR",   "RESET",   "PRINTWAV",
            "DUMPR8",  "DUMPR4",  "NSCATALO","FREE",    "NUMRNUM",
            "CALLTSO", "BETA",    "UNDEFINE","PARAMETE","KEYWORDS",
            "DUMPI2",  "NSSTATUS","ALLOCATE","BIMULT",  "NRNLIMS",
            "WRITENS", "BETACOUL","SAVEHS",  "USEHS",   "DUMPALL",
            "COPY",    "FIT",     "DATA",    "CHANNEL", "NEWPAGE",
            "BELX",    "REALPOTE","IMAGPOTE","REALSOPO","IMAGSOPO",
            "COULOMBP","WAVEFUNC","GET",     "REALTRPO","IMAGTRPO",
            "REALTLPO","IMAGTLPO","REALTPRP","IMAGTPRP","DUMPI4",
            "DUMPHEX", "CHANNELS","COUPLING","TRANSIST","TIMELIMI",
            "DUMP",    "CCDUMP",  "BETAP",   "BETAT",   "SIPOTENT",
            "WRITENS8","WRITENS4"
        };
        for (int k = 0; k < 62; k++) KEYWRD[k+1] = char8(kw[k]);
    }

    // ======================================================================
    // Alias table (148 entries, paired: input word / internal word)
    // ======================================================================
    static const int NUMALI = 148;
    static char8 ALIAI[149]; // 1-based
    static bool ali_init = false;
    if (!ali_init) {
        ali_init = true;
        const char* ali[148] = {
            "W",       "VI",
            "VR",      "V",       "RR0",      "R0",      "AR",       "A",       "RR",       "R",
            "END",     "RETURN",  "QUIT",     "RETURN",  "STOP",     "RETURN",
            "NAMES",   "NSCATALO","BE2",      "BELX",    "BE",       "BELX",
            "BETANUCL","BETA",    "LMINSUBT", "LMINSUB", "TITLE",    "HEADER",
            "GENAUIGK","ACCURACY","WINKELMA", "ANGLEMAX",
            "WINKELMI","ANGLEMIN","WINKELSC", "ANGLESTE",
            "DWENDPUN","DWCUTOFF","PASSENGE", "FITACCUR",
            "PASSENVE","FITRATIO","LMAXHINZ", "LMAXADD",
            "LMAXVERM","LMAXMULT","LMINVERM", "LMINMULT",
            "LMINABZI","LMINSUB", "LSCHRITT", "LSTEP",   "KNOTEN",   "NODES",
            "SCHRITTL","STEPSIZE","SCHRITTE", "STEPSPER",
            "SCHREIBS","WRITESTE","ELASTISC", "ELASTIC",
            "EINGANG", "INCOMING","LABWINKE", "LABANGLE",
            "KERNEINZ","NUCONLY", "AUSGANG",  "OUTGOING",
            "WURFGESC","PROJECTI","REELLEWE", "REALWAVE",
            "STREUUNG","SCATTERI","ZIELSCHE", "TARGET",
            "MITCOULO","USECOULO","MITWURFG", "USEPROJE",
            "MITZIELS","USETARGE","MITKERNG", "USECORE",
            "ZUTEILEN","ALLOCATE","KANAL",    "CHANNEL",  "BESTIMMU","DATA",
            "ANPASSEN","FIT",     "UEBERSCH", "HEADER",
            "REAKTION","REACTION","VERLASSE", "QUIT",
            "NAMEN",   "NSCATALO","SCHLUESS", "KEYWORDS","ENDE",     "QUIT",
            "PROJEKTI","PROJECTI","MITPROJE", "USEPROJE",
            "DRUCKEN", "PRINT",   "IMAGHOCH", "IMAGPOWE",
            "REELLEHO","REALPOWE","GEBUNDAS", "BOUNDASY",
            "STEUASYM","SCATASYM","SUMMEPUN", "SUMPOINT",
            "EINLAGER","SAVEHS",  "GEBRAUCH", "USEHS",
            "DUMP8",   "DUMPR8",  "NSDUMP",   "DUMPR8",   "DUMPINTE","DUMPI4",
            "DUMPOCTA","DUMPHEX", "DUMPHEXA", "DUMPHEX",
            "DUMP4",   "DUMPR4",  "BETAPROJ", "BETAP",
            "BETATARG","BETAT",   "MAXITERA", "MAXITER",
            "TRANSITI","TRANSIST","TSO",       "INTERACT"
        };
        for (int k = 0; k < 148; k++) ALIAI[k+1] = char8(ali[k]);
    }

    // ======================================================================
    // JGOTO dispatch — resume after returning to main to get something done
    //   GO TO ( 10, 990, 5, 750, 550, 740 ), JGOTO
    // ======================================================================
    switch (JGOTO) {
        case 1: goto L10;
        case 2: goto L990;
        case 3: goto L5;
        case 4: goto L750;
        case 5: goto L550;
        case 6: goto L740;
        default: break;
    }

    // ======================================================================
    // Initialize on first call (label 5)
    // ======================================================================
L5:
    DEFALT(NUMFLT, NUMINT, NUMJ, NUMSW, FALSE_F);
    // debug removed
    NEWCD();

    // ======================================================================
    // Search for next keyword (label 10)
    //   CALL NXWORD( KYWORD, *925, *600, *925 )
    //   return 0=normal, 1=alt1(→925), 2=alt2(→600), 3=alt3(→925)
    // ======================================================================
L10:
    {
        int r = NXWORD(KYWORD.data);
        // debug removed
        if (r == 1 || r == 3) goto L925;
        if (r == 2) goto L600;
    }

    // ======================================================================
    // Look up keyword in FLOAT table (label 20)
    // ======================================================================
L20:
    // debug removed
    for (I = 1; I <= NUMFLT; I++) {
        if (FLTWRD[I] != KYWORD) continue;
        if (NXVAL(VAL) != 0) goto L910;
        FLOAT_arr(I) = VAL;
        goto L10;
    }

    // Look up in INTEGER table
    for (I = 1; I <= NUMINT; I++) {
        if (INTWRD[I] != KYWORD) continue;
        { int iret_nx; NXINT(IVAL, iret_nx); if (iret_nx != 0) goto L915; }
        INTGER_arr_f(I) = IVAL;
        // keyword stored
        goto L10;
    }

    // Look up in J-value table
    for (I = 1; I <= NUMJ; I++) {
        if (JWRDS[I] != KYWORD) continue;
        if (NXHINT(IVAL) != 0) goto L920;
        JBLOCK_arr_f(I) = IVAL;
        goto L10;
    }

    // Look up in switch table
    for (I = 1; I <= NUMSWW; I++) {
        if (SWTWRD[I] != KYWORD) continue;
        SWITCH_arr()[SWTNUM[I]-1] = SWTSET[I];
        goto L10;
    }

    // ======================================================================
    // If it is a control keyword, do it now
    // ======================================================================
    USAGE.NSPTSW = ((INTGER.IPRINT / 100000) % 10) > 1 ? TRUE_F : FALSE_F;

    for (I = 1; I <= NUMKEY; I++) {
        if (KEYWRD[I] != KYWORD) continue;
        //   GO TO (110, 140, 270, 500, 510,    I=1..5
        //          150, 290, 280, 130, 160,    I=6..10
        //          236, 231, 200, 260, 420,    I=11..15
        //          350, 360, 480, 380, 120,    I=16..20
        //          233, 210, 250, 400, 430,    I=21..25
        //          300, 360, 530, 540, 220,    I=26..30
        //          320, 492, 496, 505, 470,    I=31..35
        //          360, 560, 560, 560, 560,    I=36..40
        //          560, 560, 295, 560, 560,    I=41..45
        //          560, 560, 560, 560, 232,    I=46..50
        //          234, 497, 498, 499, 570,    I=51..55
        //          230, 340, 360, 360, 560,    I=56..60
        //          300, 305                    I=61..62
        //                                  ), I
        switch (I) {
            case  1: goto L110;
            case  2: goto L140;
            case  3: goto L270;
            case  4: goto L500;
            case  5: goto L510;
            case  6: goto L150;
            case  7: goto L290;
            case  8: goto L280;
            case  9: goto L130;
            case 10: goto L160;
            case 11: goto L236;
            case 12: goto L231;
            case 13: goto L200;
            case 14: goto L260;
            case 15: goto L420;
            case 16: goto L350;
            case 17: goto L360;
            case 18: goto L480;
            case 19: goto L380;
            case 20: goto L120;
            case 21: goto L233;
            case 22: goto L210;
            case 23: goto L250;
            case 24: goto L400;
            case 25: goto L430;
            case 26: goto L300;
            case 27: goto L360;
            case 28: goto L530;
            case 29: goto L540;
            case 30: goto L220;
            case 31: goto L320;
            case 32: goto L492;
            case 33: goto L496;
            case 34: goto L505;
            case 35: goto L470;
            case 36: goto L360;
            case 37: goto L560;
            case 38: goto L560;
            case 39: goto L560;
            case 40: goto L560;
            case 41: goto L560;
            case 42: goto L560;
            case 43: goto L295;
            case 44: goto L560;
            case 45: goto L560;
            case 46: goto L560;
            case 47: goto L560;
            case 48: goto L560;
            case 49: goto L560;
            case 50: goto L232;
            case 51: goto L234;
            case 52: goto L497;
            case 53: goto L498;
            case 54: goto L499;
            case 55: goto L570;
            case 56: goto L230;
            case 57: goto L340;
            case 58: goto L360;
            case 59: goto L360;
            case 60: goto L560;
            case 61: goto L300;
            case 62: goto L305;
            default: break;
        }
    }

    // ======================================================================
    // Is the word a synonym?  If so start the whole search over again.
    // ======================================================================
    for (I = 1; I <= NUMALI; I += 2) {
        if (ALIAI[I] != KYWORD) continue;
        KYWORD = ALIAI[I+1];
        goto L20;
    }

    // ======================================================================
    // Misspelled keyword
    // ======================================================================
    std::printf("\n**** MISSPELLED KEYWORD: %.8s\n", KYWORD.data);
    if (SWITCH.ITSO == 1) goto L980;
    JUMP = 955;
    VALSW = FALSE_F;
    goto L125;

    // ======================================================================
    // LISTKEYS  (keyword list without values)
    // ======================================================================
L110:
    VALSW = FALSE_F;
    JUMP = 10;
    goto L125;

    // ======================================================================
    // KEYWORDS  (list words and values)
    // ======================================================================
L120:
    VALSW = TRUE_F;
    JUMP = 10;
L125:
    LSTKEY(VALSW, NUMFLT, FLTWRD, NUMINT, INTWRD,
           NUMJ, JWRDS, NUMSWW, SWTWRD,
           const_cast<int*>(SWTNUM), const_cast<int*>(SWTSET),
           NUMKEY, KEYWRD, NUMALI, ALIAI);
    // GO TO JUMP, ( 10, 955 )
    if (JUMP == 10)  goto L10;
    if (JUMP == 955) goto L955;
    goto L10;

// -----------------------------------------------------------------------
//     RESET - RESET KEYWORDS AND ALLOCATOR
// -----------------------------------------------------------------------
L130:
    DEFALT(NUMFLT, NUMINT, NUMJ, NUMSW, TRUE_F);
    // debug removed
    goto L280;

// -----------------------------------------------------------------------
//     RETURN
// -----------------------------------------------------------------------
L140:
    KEYRET = 1;
    IGOTO = 0;
    if (IHSAVE > 0) {
        // REWIND 1 — no-op in C++ (file-based I/O not used)
    }
    if (ITSO == 1) return;
    goto L215;

// -----------------------------------------------------------------------
//     PRINTFIT
// -----------------------------------------------------------------------
L150:
    std::printf(" ENERGY =%16.7G     V =%16.7G\n", E, V);
    goto L10;

// -----------------------------------------------------------------------
//     PRINTWAVE
// -----------------------------------------------------------------------
L160:
    std::printf("    STEP       R          POT          WAVEFUNCTION\n");
    if (WRITES_val == UNDEF) goto L175;
    if (IPOT == 0) goto L170;
    NWRITE = (int)(WRITES_val / STEPSZ + 0.5);
    if (NWRITE < 1) NWRITE = 1;
    NSTEPS = (int)(BNDASY / STEPSZ + 1.5);
    if (ASYMPT != UNDEF) NSTEPS = (int)(ASYMPT / STEPSZ + 1.5);
    for (III = 1; III <= NSTEPS; III += NWRITE) {
        RRR = (III - 1) * STEPSZ;
        std::printf(" %6d%12.5f%15.5G%15.5G\n",
            III, RRR, ALLOC(z[IPOT] + III - 1), ALLOC(z[IWFN] + III - 1));
    }
    goto L10;

L170:
    std::printf("\n**** PRINTWAVE CAN BE USED ONLY AFTER BOUNDSTATE\n");
    goto L10;

L175:
    std::printf("\n**** CANNOT---- WRITESTEP MUST BE SET NON-ZERO \n");
    goto L10;

// -----------------------------------------------------------------------
//     NSCATALOG
// -----------------------------------------------------------------------
L200:
    if (ICNTT != 0) NSCAT(IDUM);
    goto L10;

// -----------------------------------------------------------------------
//     NSSTATUS
// -----------------------------------------------------------------------
L210:
    KEYRET = 0;
L215:
    if (ICNTT != 0) NSSTAT(IDUM);
    if (KEYRET != 1) goto L10;
    TIME = second() - TBEGIN;
    std::printf("-TOTAL PTOLEMY CPU TIME =%7.2f SECONDS.\n", (float)TIME);
    return;

// -----------------------------------------------------------------------
//     DUMPALL
// -----------------------------------------------------------------------
L220:
    NSDUMP(0, 1);
    goto L10;

// -----------------------------------------------------------------------
//     DUMPXX <NAME OR NUMBER>
//       WHERE XX IS ONE OF  8, R8, 4, R4, I4, I2, HEX, HEXADECIMAL,
//          OCTAL, INTEGER OR UNSPECIFIED.
//          NSDUMP IS ALSO THE SAME AS DUMP.
// -----------------------------------------------------------------------
L230:
    I = 1;
    goto L240;
L231:
    I = 2;
    goto L240;
L232:
    I = 3;
    goto L240;
L233:
    I = 4;
    goto L240;
L234:
    I = 5;
    goto L240;
L236:
    I = 6;

L240:
    if (GETNUM(IDUM) != 0) goto L980;
    NSDUMP(IDUM, I);
    goto L10;

// -----------------------------------------------------------------------
//     ALLOCATE <NAME> <COMPLETE LIST OF VALUES>
// -----------------------------------------------------------------------
L250:
    if (NXWORD(WORD.data) != 0) goto L930;
    if (ICNTT == 0) ISTART(ISIZE);
    ICNTT = 1;
    DEFINE(WORD, NUM1, ISIZE);  // must pass ISIZE (allocator size), not NUM2 — NUM2 is uninitialized here
    std::printf("%8d NUMBERS READ TO DEFINE OBJECT %.8s\n", NUM1, WORD.data);
    if (NUM1 > ISIZE) goto L950;
    goto L10;

// -----------------------------------------------------------------------
//     FREE <NAME OR NUMBER>
// -----------------------------------------------------------------------
L260:
    if (GETNUM(NUM1) != 0) goto L980;
    z[NUM1] = -z[NUM1];
    goto L10;

// -----------------------------------------------------------------------
//     RENAME  <NAME OR NUMBER>  <NEW NAME>
// -----------------------------------------------------------------------
L270:
    if (GETNUM(NUM1) != 0) goto L950;
    if (NXWORD(WORD.data) != 0) goto L930;
    NAMCOM.NAMES[NUM1] = WORD;
    goto L10;

// -----------------------------------------------------------------------
//     CLEAR
// -----------------------------------------------------------------------
L280:
    if (ICNTT != 0) ICLEAR(IDUM);
    MAXLF = 0;
    goto L10;

// -----------------------------------------------------------------------
//     KEEP  <NAME>  <MEMBER NAME>
// -----------------------------------------------------------------------
L290:
    KEEPIT();
    // Cap: KEEP NOT AVAILABLE message removed (KEEPIT handles it)
    goto L10;

// -----------------------------------------------------------------------
//     GET  <PTOLEMY NAME>   <MEMBER NAME>
// -----------------------------------------------------------------------
L295:
    if (ICNTT == 0) ISTART(ISIZE);
    ICNTT = 1;
    GETIT(IRET);
    if (IRET < 0) goto L950;
    if (IRET == 0) goto L950;
    goto L10;

// -----------------------------------------------------------------------
//     WRITENS, WRITENS8, WRITENS4 <NAME OR NUM>
//     WRITES TO FILE 15
// -----------------------------------------------------------------------
L300:
    keylen = 1;
    goto L310;
L305:
    keylen = 2;

L310:
    if (GETNUM(ISEND) != 0) goto L980;
    LLL = keylen * LENGTH.LENG[ISEND];
    std::printf(" NUMBERED STORAGE OBJECT %3d with%7d elements"
                " WRITTEN TO FT15 WITH 5G14.5 FORMAT\n", ISEND, LLL);
    LOCC = keylen * (z[ISEND] - 1) + 1;
    if (LLL > 0) {
        if (keylen == 1) {
            for (III = 1; III <= LLL; III++) {
                std::printf(" %14.5G", ALLOC(LOCC + III - 1));
                if (III % 5 == 0 || III == LLL) std::printf("\n");
            }
        }
        if (keylen == 2) {
            for (III = 1; III <= LLL; III++) {
                std::printf(" %14.5G", (double)ALLOC4(LOCC + III - 1));
                if (III % 5 == 0 || III == LLL) std::printf("\n");
            }
        }
    } else {
        std::printf("  Nothing to write!!\n");
    }
    goto L10;

// -----------------------------------------------------------------------
//     COPY  <NAME1 OR NUM1>  <NAME2>
// -----------------------------------------------------------------------
L320:
    if (GETNUM(NUM1) != 0) goto L900;
    INEED = LENGTH.LENG[NUM1];
    if (NXWORD(WORD.data) != 0) goto L900;
    // IS THE "TO" AREA ALREADY DEFINED
    NUM2 = NAMLOC(WORD.data);
    if (NUM2 == 0) goto L325;
    // YES, IS IT THE CORRECT SIZE
    LLL = LENGTH.LENG[NUM2];
    if (INEED == LLL) goto L326;
    // NO, FREE AND MAKE A NEW ONE
    std::printf(" WARNING:  LENGTH OF %.8s IS PRESENTLY%6d"
                " BUT COPY IS CHANGING IT TO%6d\n", WORD.data, LLL, INEED);
    z[NUM2] = -z[NUM2];
L325:
    NUM2 = NALLOC(WORD.data, INEED);
L326:
    for (I = 1; I <= INEED; I++) {
        ALLOC(z[NUM2] - 1 + I) = ALLOC(z[NUM1] - 1 + I);
    }
    goto L10;

// -----------------------------------------------------------------------
//     CCDUMP
// -----------------------------------------------------------------------
L340:
    CCDUMP();
    goto L10;

// -----------------------------------------------------------------------
//     CALLTSO  KEYWORD ENTERED
// -----------------------------------------------------------------------
L350:
    IGOTO = 3;
    JGOTO = 1;
    return;

// -----------------------------------------------------------------------
//     DEFINE A KEYWORD THAT HAS A VARIABLE NUMBER OF VALUES.
//     THE VALUES ARE STORED IN THE ALLOCATOR WITH THE KEYWORD NAME.
//
//     PRESENT CASES ARE BETA, BETAP, BETAT, BETACOULOMB & BELX
// -----------------------------------------------------------------------
L360:
    WORD = KYWORD;
    if (ICNTT == 0) ISTART(ISIZE);
    ICNTT = 1;
    DEFINE(WORD, NUM1, ISIZE);  // must pass ISIZE (allocator size), not NUM2
    if (NUM1 > NUM2) goto L950;
    goto L10;

// -----------------------------------------------------------------------
//     PARAMETERSET KEYWORD IS USED TO INITIALIZE PARAMETER PACKAGES
// -----------------------------------------------------------------------
L380:
    PARAM(IGOTO);
    // debug removed
    if (IGOTO < 0) goto L10;
    if (IGOTO == 0) goto L900;
    goto L10;

// -----------------------------------------------------------------------
//     BIMULT <NAME OR NUM>  <NAME OR NUM>
//     OVERLAP OF TWO WAVEFUNCTIONS
// -----------------------------------------------------------------------
L400:
    if (GETNUM(NUM1) != 0) goto L980;
    if (GETNUM(NUM2) != 0) goto L980;
    INEED = LENGTH.LENG[NUM1];
    ITWO = LENGTH.LENG[NUM2];
    if (INEED != ITWO)
        std::printf("\n*** WARNING: LENGTHS NOT EQUAL, THEY ARE"
                    "%8d%8d.  THE SHORTER WILL BE USED.\n", INEED, ITWO);
    INEED = std::min(INEED, ITWO);
    IALL = IALLOC(INEED);
    for (II = 1; II <= INEED; II++) {
        ALLOC(z[IALL] + II - 1) = ALLOC(z[NUM1] + II - 1) * ALLOC(z[NUM2] + II - 1);
    }
    std::printf(" HAVE MULTIPLIED TWO %5d-ELEMENT ARRAYS %3d AND"
                "%3d AND PUT THE RESULT IN %4d\n", INEED, NUM1, NUM2, IALL);
    goto L10;

// -----------------------------------------------------------------------
//     NUMRNUM <NAME OR NUM>  <POWER>  <NAME OR NUM>
//     MATRIX ELEMENT OF A POWER OF R
// -----------------------------------------------------------------------
L420:
    if (GETNUM(NUM1) != 0) goto L980;
    NXINT(IARE, nxret);
    if (nxret != 0) goto L900;
    if (GETNUM(NUM2) != 0) goto L980;
    ISTP1 = 1;
    INEED = LENGTH.LENG[NUM1];
    ITWO = LENGTH.LENG[NUM2];
    if (INEED != ITWO)
        std::printf("\n*** WARNING: LENGTHS NOT EQUAL, THEY ARE"
                    "%8d%8d.  THE SHORTER WILL BE USED.\n", INEED, ITWO);
    ISTP2 = std::min(INEED, ITWO);
    goto L435;

// -----------------------------------------------------------------------
//      NRNLIMS <NAME OR NUM> <POWER> <NAME OR NUM> <START INDEX>
//                <END INDEX>
//     PARTIAL MATRIX ELEMENT OF R**N
// -----------------------------------------------------------------------
L430:
    if (GETNUM(NUM1) != 0) goto L980;
    NXINT(IARE, nxret);
    if (nxret != 0) goto L900;
    if (GETNUM(NUM2) != 0) goto L980;
    if (NXVAL(TEMP) != 0) goto L900;
    ISTP1 = (int)(TEMP / STEPSZ + 1.5);
    if (NXVAL(TEMP) != 0) goto L900;
    ISTP2 = (int)(TEMP / STEPSZ + 1.5);

L435:
    std::printf(" AM CALCULATING ME BETWEEN OBJECTS%3d AND%3d OF R**"
                "%1d FROM STEPS%5d TO %5d\n", NUM1, NUM2, IARE, ISTP1, ISTP2);
    AME = 0.0;
    LLLL = std::min(LENGTH.LENG[NUM1], LENGTH.LENG[NUM2]);
    if (LLLL >= ISTP2) goto L445;
    TEMP = STEPSZ * (LLLL - 1);
    std::printf("\n**** WARNING:  OBJECTS EXTEND ONLY TO R =%10.2f"
                " FMS; WE WILL STOP AT END OF OBJECTS.\n", TEMP);
    ISTP2 = LLLL;
L445:
    for (II = ISTP1; II <= ISTP2; II++) {
        RRR = (II - 1) * STEPSZ;
        TEMP = 1.0;
        if (IARE != 0) TEMP = std::pow(RRR, IARE);
        AME = AME + TEMP * STEPSZ * ALLOC(z[NUM1] + II - 1) * ALLOC(z[NUM2] + II - 1);
    }
    std::printf("\n        MATRIX ELEMENT =%14.5E\n\n", AME);
    goto L10;

// -----------------------------------------------------------------------
//     NEWPAGE
// -----------------------------------------------------------------------
L470:
    std::printf("\f");  // form feed
    goto L10;

// -----------------------------------------------------------------------
//     UNDEFINE  < DATA KEYWORD >
// -----------------------------------------------------------------------
L480:
    if (NXWORD(WORD.data) != 0) goto L900;
    for (I = 1; I <= NUMFLT; I++) {
        if (FLTWRD[I] != WORD) continue;
        FLOAT_arr(I) = UNDEF;
        goto L10;
    }
    for (I = 1; I <= NUMINT; I++) {
        if (INTWRD[I] != WORD) continue;
        INTGER_arr_f(I) = (int)NOTDEF;
        goto L10;
    }
    for (I = 1; I <= NUMJ; I++) {
        if (JWRDS[I] != WORD) continue;
        JBLOCK_arr_f(I) = NOTDEF;
        goto L10;
    }
    std::printf("\n**** %.8s IS NOT A DATA-DEFINING KEYWORD.\n", WORD.data);
    goto L950;

// -----------------------------------------------------------------------
//   FIT BEING SPECIFIED--- READ INPUT
// -----------------------------------------------------------------------
L492:
    if (ICNTT == 0) ISTART(ISIZE);
    ICNTT = 1;
    FITINP(IRET, NUMFLT, FLTWRD, NUMALI, ALIAI);
    PROBLM = 21;
    goto L535;

// -----------------------------------------------------------------------
//     DATA ( ... )  READS IN DATA FOR FITTER
// -----------------------------------------------------------------------
L496:
    if (ICNTT == 0) ISTART(ISIZE);
    ICNTT = 1;
    DATAIN(IRET);
    goto L535;

// -----------------------------------------------------------------------
//     CHANNELS (...)  DEFINES CHANNS FOR COUPLED CHAN
// -----------------------------------------------------------------------
L497:
    if (ICNTT == 0) ISTART(ISIZE);
    ICNTT = 1;
    CCHAN(IRET);
    goto L535;

// -----------------------------------------------------------------------
//     COUPLING (...)  DEFINES THE COUPLING SCHEME AND SETS
//                     THE CCINELAS SWITCH
// -----------------------------------------------------------------------
L498:
    if (ICNTT == 0) ISTART(ISIZE);
    ICNTT = 1;
    CCOUP(IRET);
    PROBLM = 24;
    goto L535;

// -----------------------------------------------------------------------
//     TRANSITIONS (...)  DEFINES THE TRANSITIONS TO USE
// -----------------------------------------------------------------------
L499:
    if (ICNTT == 0) ISTART(ISIZE);
    ICNTT = 1;
    CTRAN(IRET);
    goto L535;

// -----------------------------------------------------------------------
//     REACTION  208PB(16O, N15)BI209(1/2- 3.23)  TYPE INPUT
// -----------------------------------------------------------------------
L500:
    REACTN(IRET);
    goto L535;

// -----------------------------------------------------------------------
//     CHANNEL   P + 208PB = 209BI(7/2- .9)
// -----------------------------------------------------------------------
L505:
    CHANEL(IRET);
    goto L535;

// -----------------------------------------------------------------------
//     HEADER < REST OF LINE IS HEADER >
// -----------------------------------------------------------------------
L510:
    for (I = 1; I <= 65; I++) {
        HEDCOM.HEADER[I] = ' ';
    }
    for (I = INPBUF.INCH; I <= INPBUF.NOCH; I++) {
        if (INPBUF.IBUF[I] != ' ') goto L520;
    }
    goto L528;
L520:
    NUM2 = std::min(INPBUF.NOCH, I + 64);
    II = 0;
    for (int INCH_loc = I; INCH_loc <= NUM2; INCH_loc++) {
        if (INPBUF.IBUF[INCH_loc] != ' ') II = INCH_loc;
    }

    // RIGHT JUSTIFY HEADER FOR PRINTING ON RIGHT SIDE OF PAGE
    for (int INCH_loc = I; INCH_loc <= II; INCH_loc++) {
        HEDCOM.HEADER[65 - II + INCH_loc] = INPBUF.IBUF[INCH_loc];
    }
L528:
    INPBUF.INCH = INPBUF.NOCH + 1;
    goto L10;

// -----------------------------------------------------------------------
//     SAVEHS
// -----------------------------------------------------------------------
L530:
    SAVEHS();
    // Fortran: CALL SAVEHS( IRET ) — IRET set by SAVEHS
    // C++ stub SAVEHS() has no args; assume IRET=1 (success)
    IRET = 1;
L535:
    if (IRET < 0) goto L10;
    if (IRET == 0) goto L950;
    goto L10;

// -----------------------------------------------------------------------
//     USEHS
// -----------------------------------------------------------------------
L540:
    IRET = 1;
L545:
    // Fortran: CALL USEHS(NUMFLT,FLOAT,NUMINT,INTGER,NUMJ,JBLOCK,NUMSW,SWITCH,IRET)
    // C++ stub USEHS() has no args; IRET set manually
    USEHS();
    // IRET not modified by stub; keep existing value
    if (IRET < 0) goto L10;
    if (IRET == 0) goto L950;
    goto L547;

    // ARE WE STILL PROCESSING BOUND STATE
L547:
    JGOTO = 5;
    if (NEXT <= 2) goto L720;
    NEXT = 0;
    // INDICATE BOTH BOUND STATES DONE
    IDONE = 3;
    if (IBSPAS > 1) IDONE = 15;
    goto L10;

// -----------------------------------------------------------------------
//     DONE WITH A READIN BOUND STATE; GET NEXT SET OF COMMON BLOCKS
// -----------------------------------------------------------------------
L550:
    IRET = 2;
    goto L545;

// -----------------------------------------------------------------------
//     REALPOTEN <NAME>  ,  IMAGPOTEN <NAME>   , ETC.
// -----------------------------------------------------------------------
L560:
    GETLNK(KYWORD, IRET);
    if (IRET < 0) goto L950;
    if (IRET == 0) goto L950;
    goto L10;

// -----------------------------------------------------------------------
//     TIMELIMIT (VAX/VMS ONLY)
// -----------------------------------------------------------------------
L570:
    if (NXVAL(TEMP) != 0) goto L900;
    std::printf("\n**** TIMELIMIT IS NOT IN THIS VERSION.  IGNORED.\n");
    goto L10;

// -----------------------------------------------------------------------
//     ; ENCOUNTERED.  DO A STAGE OF THE CALCULATION
// -----------------------------------------------------------------------
L600:
    NSPTSW = (((IPRINT / 100000) % 10) > 0) ? TRUE_F : FALSE_F;
    if (ICNTT == 0) ISTART(ISIZE);
    ICNTT = 1;

    // NUMBERED STORAGE INITIALIZED
    // NOW ON OUR WAY

    // NOW DETERMINE WHAT TO DO FOR THE ;
    if (NEXT == 0) NEXT = PROBLM;
    if (NEXT != 0) goto L700;

    // HE DIDN'T SPECIFY; EITHER FINISH A DWBA OR DO A STAND ALONE CHANNEL.
    if (IDONE > 0) goto L650;

    // IT IS STAND ALONE
    if (E != UNDEF) goto L630;
    if (ELAB != UNDEF || ECM != UNDEF) goto L635;
    std::printf("\n**** E OR ELAB OR ECM MUST BE SPECIFIED.\n");
    goto L950;
L630:
    NEXT = 5;
    if (E < 0.0) goto L700;
L635:
    NEXT = 6;
    goto L700;

    // ASSUME (FOR THE MOMENT) THAT ALL 2-BODY CHANNELS ARE
    // DEFINED AND START A TRANSFER OR INELASTIC SCATTERING CALCULATION
L650:
    NEXT = 20;
    if (ISTRIP == 0) NEXT = 22;
    PROBLM = NEXT;

    // NEXT SAYS WHAT TO DO NEXT
L700:
    if (NEXT >= 10) goto L770;

    // INITIALIZE A CHANNEL AND ITS POTENTIAL
    SETCHN(SWITCH.NEXT, IRET);
    if (IRET == 0) goto L995;
    SETPOT(IRET);
    if (IRET == 0) goto L995;

    // NOW DO BOUND STATE OR SCATTERING.  UPON RETURN CLEAR MANY DEFINITIONS
    JGOTO = 4;
L720:
    IGOTO = 1;
    switch (NEXT) {
        case 1: goto L722;
        case 2: goto L724;
        case 3: goto L726;
        case 4: goto L728;
        case 5: goto L730;
        case 6: goto L732;
        case 7: goto L900;
        case 8: goto L736;
        default: goto L900;
    }
    // PROJECTILE
L722:
    ICHANB = 1;
    return;
    // TARGET
L724:
    ICHANB = 2;
    return;
    // INCOMING
L726:
    ICHANW = 1;
    IGOTO = 5;
    return;
    // OUTGOING
L728:
    ICHANW = 2;
    IGOTO = 5;
    return;
L730:
    ICHANB = ICHANB + 1;
    return;
L732:
    ICHANW = ICHANW + 1;
    IGOTO = 13;
    return;

    // ********   THIS METHOD OF GETTING TO FITTER NOT USED
L736:
    ICHANW = 1;
    IGOTO = 2;
    return;

    // WE HAVE RETURNED FROM A PASS OF A MULTIPASS CALCULATION.
    // PREPARE FOR NEW BOUND STATE PARAMETERS.
    // (Label 740 is entered via JGOTO=6 from the top of CONTRL)
    // Label 740 is not in contrl_label switch since it's only
    // reached from JGOTO dispatch in part1. We include it here
    // for completeness since it uses labels in this part.
    // In the merged function it would be reachable.

    // NOTE: The following labels (740, 750, 760, 765, 770, 780, 800,
    // 810, 850, 860, 870, 900..998) are reached via fall-through or
    // goto from code within this same function body.

    {
        // Unreachable entry — these are only reached via goto within CONTRL
        goto L_skip_internal_labels;

L740:
        NEXT = 0;
        IDONE = 31;
        IBSPAS = 2;
        LX = (int)NOTDEF;
        IRET = 3;
        if (IHSAVE == 2) goto L545;
        goto L760;

        // WE HAVE RETURNED FROM A BOUND STATE OR SCATTERING STATE
L750:
        IRET = 1;
        if (IPRM1 == 0) goto L760;
        if (NEXT <= 4) IDONE = IDONE | (1 << (NEXT - 1));
        if (E > 0.0) goto L760;
        IPOT = IPRM2;
        IWFN = IPRM3;
        if (IHSAVE == 1) {
            // Fortran: CALL SAVECM(NUMFLT,FLOAT,NUMINT,INTGER,NUMJ,JBLOCK,NUMSW,SWITCH,IRET)
            SAVECM();
            // IRET not modified by stub; keep IRET=1
        }
L760:
        CLRCHN(NEXT);
        if (IPRM1 == 0) goto L995;
        if (IRET == 0) goto L950;
L765:
        NEXT = 0;
        goto L10;

        // WE ARE NOW DOING THINGS SUBSEQUENT TO THE CHANNEL INPUT
L770:
        JGOTO = 2;
        if (NEXT == 21) goto L780;
        if (NEXT < 30) goto L800;

        // THIS IS A DOTEMP CALL
        IGOTO = NEXT - 20;
        return;

        // START ELASTIC SCATTERING FIT
L780:
        IGOTO = 2;
        JGOTO = 4;
        return;

        // CARRY ON WITH THE DWBA OR C.C. IF ALL CHANNELS WERE DEFINED.
L800:
        if (ISTRIP == 0) goto L805;
        if (PROBLM == 0) PROBLM = 20;
        if (IBSPAS > 0) JGOTO = 6;
        if ((IDONE % 16) == 15) goto L810;
        std::printf("\n**** PROJECTILE, TARGET, INCOMING AND OUTGOING"
                    " CHANNELS MUST ALL BE DEFINED\n");
        goto L950;
L805:
        if (PROBLM == 0) PROBLM = 22;
        if (((IDONE / 4) % 4) == 3) goto L810;
        std::printf("\n**** INCOMING AND OUTGOING CHANNELS MUST BE"
                    " DEFINED FOR INELASTIC EXCITATION.\n");
        goto L950;

L810:
        if (NEXT >= 20) goto L850;

        // HE IS DOING IT IN PIECES AND WE DON'T CHECK MUCH
        IDONE = IDONE | 16;
        IGOTO = NEXT - 4;
        if (NEXT == 10) goto L870;

        // IF CROSSSECTIONS IS SPECIFIED, IT IMPLIES LINTERPOLATION ALSO
        // UNLESS THAT HAS JUST BEEN DONE.
        if (NEXT == 13 && ((IDONE / 32) % 2) == 0) IGOTO = 15;
        if (NEXT == 12) IDONE = IDONE | 32;
        return;

        // DO THE WHOLE THING; DON'T ALLOW REPEATS (EXCEPT MULTIPASS)
L850:
        if (((IDONE / 16) % 2) == 0) goto L860;
        if (IBSPAS > 0) goto L855;
        std::printf("\n**** YOU ARE REQUESTING THINGS OUT OF ORDER\n");
        goto L950;
L855:
        IBSPAS = 2;
        IDONE = IDONE % 16;
L860:
        IGOTO = 14;

        // FOR SAVEHS SAVE COMMONS AT ENTRY TO GRDSET
L870:
        IRET = 1;
        if (IHSAVE == 1) {
            // Fortran: CALL SAVECM(NUMFLT,FLOAT,NUMINT,INTGER,NUMJ,JBLOCK,NUMSW,SWITCH,IRET)
            SAVECM();
            // IRET not modified by stub; keep IRET=1
        }
        if (IRET == 0) goto L950;
        return;

L_skip_internal_labels:
        (void)0; // placeholder
    }

// -----------------------------------------------------------------------
//     ERRORS DETECTED, FOR TSO KEEP GOING, FOR BATCH QUIT
// -----------------------------------------------------------------------
L900:
    std::printf("\n****         ------  INPUT ERROR  -------\n"
                " KEYWORD BEING PROCESSED WAS  %.8s\n\n", KYWORD.data);
    goto L950;

    // (label 910, 915, 920, 925, 930 are in part 1 for the keyword
    //  search section.  They are included here for completeness.)

L910:
    std::printf("\n**** A NUMERIC DATA VALUE MUST FOLLOW THE %.8s"
                " KEYWORD.\n", KYWORD.data);
    goto L950;

L915:
    std::printf("\n**** AN INTEGER DATA VALUE MUST FOLLOW THE %.8s"
                " KEYWORD.\n", KYWORD.data);
    goto L950;

L920:
    std::printf("\n**** AN INTEGER OR HALF-INTEGER DATA VALUE MUST "
                "FOLLOW THE %.8s KEYWORD.\n", KYWORD.data);
    goto L950;

L925:
    III = INPBUF.INCH - 1;
    std::printf("\n**** INVALID FIELD ENCOUNTERED WHEN LOOKING FOR A"
                " KEYWORD.\n");
    std::printf("      INPUT UP TO ERROR IS: ");
    for (I = 1; I <= III; I++) {
        std::printf("%c", INPBUF.IBUF[I]);
    }
    std::printf("\n");
    goto L950;

L930:
    std::printf("\n**** A NAME MUST FOLLOW THE %.8s KEYWORD.\n", KYWORD.data);
    goto L950;

L950:
    if (ITSO == 1) goto L980;

    // IN BATCH, GIVE A LIST OF EVERYTHING AND TERMINATE
    VALSW_int = TRUE_F;
    JUMP_target = 955;
    goto L125;
L955:
    if (ICNTT != 0) NSCAT(IDUM);
    if (IHSAVE > 0) {
        // REWIND 1 — no-op
    }
    std::exit(1230);

    // IN TSO GIVE HIM ANOTHER TRY
L980:
    std::printf(" ***** CONTINUING WITH NEXT INPUT LINE **** \n\n");
    NEWCD();
    goto L10;

// -----------------------------------------------------------------------
//     CHECK RETURN CODE AFTER A STEP OF THE CALCULATION
// -----------------------------------------------------------------------
    // Label 990 is entered from JGOTO=2 in the top of CONTRL (part1)
L990:
    if (IPRM1 != 0) goto L765;
L995:
    std::printf("\n**** ERROR DETECTED IN PREVIOUS STEP.\n\n");
    goto L950;

}   // end CONTRL

