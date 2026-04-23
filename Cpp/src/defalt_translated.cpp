// defalt_translated.cpp — DEFALT
// Translated from source.f lines 13737-14240
// Verbatim transliteration from Fortran to C++

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include "ptolemy_alloc.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_forward.h"
#include "ptolemy_io.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <ctime>

// External timing function
extern double second();

// ============================================================================
// SUBROUTINE DEFALT ( NUMFLT, NUMINT, NUMJ, NUMSW, RESTSW )
//
//     INITIALIZATION OF THE DWBA CODE
//
//     1/18/75 - SEPARATEDED FROM CONTRL - S.P.
//     1/18/75 - REARRANGED TO BE OVERLAYED - S. P.
//     2/10/75 - MOVE  NUCONLY AND JX
//     2/16/75 - LOAD IEWSZOVR, PRINT DATE
//     2/20/75 - MOST DEFAULTS SET TO 0 - S. P.
//     3/5/75 - NEW GRID DEFAULTS
//     4/15/75 - MAJOR REVISION FOR NEW COMMON BLOCKS
//     5/5/75 - STUFF FOR SETCHN
//     5/18/75 - RESTSW ADDED.
//     5/22/75 - VSI ADDED
//     5/30/75 - FIX FOR RESET
//     6/11/75 - SAVEHS/USEHS
//     8/7/75  - ADD R0E,VE,ETC TO COMMON
//     8/28/75 - LMINMULT, LMAXMULT, SPFACP SPFACT
//     9/3/75 - MAXLEXTRAP ADDED; USECOULOMB IS DEFAULT
//     9/25/75 - CHANGE DWCUT; REMOVE VCUT
//     9/30/75 - LMAXMULT BECOMES LMAXADD
//     10/1/75 - FITTING KEYWORDS
//     10/19/75 - ESQ PARAMETERS, TBEGIN
//     11/6/75 - DERIVSTEP, MAXFUNCTIONS;  LMCHOL IS DEFAULT FITTER
//     1/29/76 - FITMULTI, FITRATIO, FITMODE, NUMRANDO
//     2/26/76 - LMAXADD INCREASED FOR ELASTIC STUFF
//     6/1/76 - NVPOLY, SUMPTS ADDED
//     6/17/76 - LARGE FAKE AMDMLT
//     9/16/76 - SEPT. 76 VERSION; FITMUL FIX
//     12/26/76 - INELASTIC SCATTERING
//     4/15/77 - LMINSB, ALMXAD ADDED
//     5/8/77 - PRINT=10001; MAY 77 VERSION
//     9/24/77 - USECORE DEFAULT; SEPT 77 VERSION
//     10/16/77 - POWERS OF WOODS-SAXONS
//     11/20/77 - SEPARATE BOUND & SCAT ASYMP; DEC 77 VERSION
//     11/25/77 - LINKULES
//     3/19/78 - MARCH, 1978 VERSION
//     4/24/78 - APRIL VERSION FOR LI, LO REVERSALS - S.P.
//     4/28/78 - BIG AND SMALL STUFF IN /CNSTNT/ - S.P.
//     5/23/78 - JUNE 1978 VERSION - S.P.
//     12/7/78 - ADD TENSOR POTENTIALS - RPG
//     3/21/79 - APRIL 1979 VERSION - S.P.
//     5/6/79 - MAY 1979 VERSION; PARITIES AND LEVEL SPEC ARRAYS - S.P.
//     6/10/79 - CMS VERSION - S.P.
//     6/8/79 - JUNE 1979 VERSION; STANDARD PREFIXES, CHARACTER TYPE,
//        UNIVAC STUFF - RPG
//     6/18/79 - PROPERLY RESTORE IECHO FOR RESET - S.P.
//     9/10/79 - SEPTEMBER 1979 VERSION - RPG
//     11/3/79 - OCTOBER 1979 VERSION; TIME AND CPU VERSION STAMP - S.P.
//     11/15/79 - NOVEMBER 1979 VERSION; CMPUTR SUBROUTINE - S.P.
//     12/31/79 - DECEMBER 1979 VERSION; VAX STUFF - RPG
//     1/14/80 - FIX CVA IUNDEF/AUNDEF DATA; CHARACTER*1 BLANK1 - RPG
//     1/18/80 - CALL EDTDAT FOR CND; CHAACTER*8 DATE, THISDT - RPG
//     1/25/80 - CND CHARACTER*1 THEDAT, WRITE(6,-) TO PRINT - RPG
//     2/5/80 - INITIALIZE NBLKSZ, CVA CALL TLIMIT - RPG
//     2/19/80 - INITIALIZE MAPPHI, PHIMID, GAMPHI; FEB. 1980 VERSION - R
//     3/15/80 - MARCH 1980 VERSION FOR RPG MERGED BACK TO ANL - S.P.
//     5/2/80 - APRIL 1980 VERSION - C.C. WORKING FAIRLY WELL! S.P.
//     6/13/80 - CHAR*8 FOR CDC THEDAT, SHIFT FOR NOS - S.P.
//     7/10/80 - CWD, CHR PREFIXES, JULY 1980 VERSION - S.P.
//     7/17/80 - FIX CRAY NOTDEF, UNDEF VALUES - S.P.
//     8/25/80 - JULY 1980 VERSION, SOMEHOW LOST ON CDC - S.P.
//     9/16/80 - SEPT. 1980 VERSION; ZERO CCBLK; L-EXTRAP DEFAULT - S.P.
//     10/8/80 - OCTOBER 1980 VERSION FOR NEW BASCPL, ETC. - S.P.
//     12/17/80 - DECEMBER 1980 VERSION FOR COMPLETION OF ABOVE - S.P.
//     1/8/81 - JANUARY 1981 VERSION FOR NEW COUPLING SCHEME  - S.P.
//     6/24/81 - JUNE 1981 VERSION FOR REWORKS OF COUPLN - S.P.
//     8/28/81 - AUGUST 1981 VERSION FOR DEFORMED POTENS - S.P.
//     9/17/81 - SEPT 1981 VERSION; LEBACK DEFAULT - S.P.
//     10/12/81 - OCT 1981 VERSION FOR GIESSEN CHANGES - S.P.
//     2/11/82 - JANUARY 1982 VERSION
//     3/12/82 - MARCH 1982 VERSION
//     4/10/82 - APRIL 1982 VERSION
//     5/3/83 - MAY 1983 VERSION - VARIOUS SMALL FIXES - S.P.
//     7/20/83 - USE ISIZE=-500 FOR NOS - S.P
//     8/17/83 - AUGUST 1983 VERSION = FINAL AT MUNICH - S.P.
//     10/16/84 - OCT 1984 VERSION; SCTASY < 0 - S.P.
//     12/19/84 - DEC 1984 VERSION; ROOM FOR SIPOTENT - S.P.
//     5/10/85 - INITIALIZE PLM WORK ARRAY FOR CRA - S.P.
//     6/10/85 - JUNE, 1985 VERSION - S.P.
//     11/16/01 - Nov, 2001 version; new date routine
//
// ============================================================================

void DEFALT(int NUMFLT, int NUMINT, int NUMJ, int NUMSW, int RESTSW)
{
    //
    // EQUIVALENCE ( FLOAT(1), A ),  ( INTGER(1), L ),
    //   ( JBLOCK(1), J ),  ( SWITCH(1), IASYMP )
    //
    // Array aliases via inline functions from ptolemy_commons.h:
    //   FLOAT_arr()  => &FLOAT_common.A
    //   INTGER_arr() => &INTGER.L
    //   JBLOCK_arr() => &JBLOCK.J
    //   SWITCH_arr() => &SWITCH.IASYMP
    //
    // EQUIVALENCE ( VTEN(1), VTR )  => VTEN() inline
    //

    char BLANK1 = ' ';

    //     DATA  IUNDEF / Z"F0F0F0F0"/,  AUNDEF / Z"F0F0F0F0F0F0F0F0"/
    int IUNDEF;
    {
        unsigned int tmp = 0xF0F0F0F0u;
        std::memcpy(&IUNDEF, &tmp, sizeof(int));
    }
    double AUNDEF;
    {
        uint64_t tmp = 0xF0F0F0F0F0F0F0F0ull;
        std::memcpy(&AUNDEF, &tmp, sizeof(double));
    }

    char THEDAT[10];    // CHARACTER*9 + NUL
    std::memset(THEDAT, ' ', sizeof(THEDAT));
    THEDAT[9] = '\0';

    char THETIM[11];    // CHARACTER*10 + NUL
    std::memset(THETIM, ' ', sizeof(THETIM));
    THETIM[10] = '\0';

    char CPUID[41];
    std::memset(CPUID, ' ', sizeof(CPUID));
    CPUID[40] = '\0';

    char hostname[41];
    std::memset(hostname, ' ', sizeof(hostname));
    hostname[40] = '\0';

    //
    //     IS THIS THE FIRST CALL TO DEFALT
    //
    if (RESTSW) goto L10;
    //
    //     START TIMER
    //
    FTIME.TBEGIN = (float)second();
    //
    //  INITIALIZE INPUT
    //
    //     FIRST SET EVERYTHING TO UNDEFINED STATUS
    //
L10:
    {
        // NOTDEF = IUNDEF
        // In Fortran, NOTDEF is an integer-sized bit pattern stored in a double-word.
        // We replicate: write the bit pattern of IUNDEF into the first 4 bytes of NOTDEF.
        double nd_tmp = 0.0;
        std::memcpy(&nd_tmp, &IUNDEF, sizeof(int));
        INTRNL.NOTDEF = nd_tmp;
    }
    INTRNL.UNDEF = AUNDEF;
    //
    for (int I = 1; I <= NUMFLT; I++) {
        // FLOAT(I) = UNDEF
        FLOAT_arr(I) = INTRNL.UNDEF;
    }
    for (int I = 1; I <= NUMINT; I++) {
        // INTGER(I) = NOTDEF
        // Fortran stores integer bit pattern; replicate via memcpy
        int nd_int;
        std::memcpy(&nd_int, &INTRNL.NOTDEF, sizeof(int));
        INTGER_arr_f(I) = nd_int;
    }
    for (int I = 1; I <= NUMJ; I++) {
        // JBLOCK(I) = NOTDEF
        // JBLOCK is double-word array; store NOTDEF (a double containing IUNDEF pattern)
        int nd_int;
        std::memcpy(&nd_int, &INTRNL.NOTDEF, sizeof(int));
        // JBLOCK is double, but Fortran stores integer NOTDEF into it
        double jtmp = 0.0;
        std::memcpy(&jtmp, &nd_int, sizeof(int));
        JBLOCK_arr_f(I) = jtmp;
    }
    //
    //     SET ALL OPTIONS OFF
    //
    int ITSOSV = SWITCH.ITSO;
    int IECHSV = SWITCH.IECHO;
    for (int I = 1; I <= NUMSW; I++) {
        SWITCH_arr()[I - 1] = 0;
    }
    SWITCH.ITSO = ITSOSV;
    SWITCH.IECHO = IECHSV;
    //
    for (int I = 1; I <= 45; I++) {
        HEDCOM.REACT[I] = BLANK1;
    }
    for (int I = 1; I <= 65; I++) {
        HEDCOM.HEADER[I] = ' ';
    }
    //
    //
    //     FIND OUT ABOUT THE MACHINE AND OPERATING SYSTEM
    //
    CMPUTR(CPUID, hostname);
    //
    if (RESTSW) goto L110;
    //
    SWITCH.IECHO = 1;
    //     INDICATE NO ALLOCATOR YET.
    INTRNL.ICNTT = 0;
    //  IECHO MUST BE INITIALIZED BEFORE NEWCD IF FIRST INPUT... ISTOBEVISIBL
L110:
    if (SWITCH.ITSO == 0) std::printf("%49s%s\n", "", "P T O L E M Y");
    if (SWITCH.ITSO == 1) std::printf("1           P T O L E M Y\n");
    //
    //     GET CURRENT DATE
    //
    // L150:
    get_date(THEDAT);
    //
    // call date_and_time ( date8, thetim, zone, time_values )
    {
        time_t now = time(nullptr);
        struct tm* tm_info = localtime(&now);
        char date8[9];
        strftime(date8, sizeof(date8), "%Y%m%d", tm_info);
        strftime(THETIM, sizeof(THETIM), "%H%M%S.000", tm_info);
    }
    //
    //     **********************************************************
    //     *                                                        *
    //     *                STOP     THINK                          *
    //     *                                                        *
    //     *    IS IT TIME TO CHANGE THE DATE OF THE VERSION ??     *
    //     *                                                        *
    //     **********************************************************
    //
    // Trim CPUID and hostname for printing (remove trailing spaces)
    {
        // trim CPUID
        int len = 40;
        while (len > 0 && CPUID[len - 1] == ' ') len--;
        CPUID[len] = '\0';
        // trim hostname
        len = 40;
        while (len > 0 && hostname[len - 1] == ' ') len--;
        hostname[len] = '\0';
    }

    if (SWITCH.ITSO == 0) {
        std::printf(" April 2007 version     Computer: %s  %s          %s  %s\n\n",
                    CPUID, hostname, THEDAT, THETIM);
    }
    if (SWITCH.ITSO == 1) {
        std::printf(" February 2007 version     %s\n\n", THEDAT);
    }
    //
    //
    //     SETUP CONSTANTS
    //
    CNSTNT.PI = 3.14159265358979320;
    CNSTNT.RT4PI = std::sqrt(4 * CNSTNT.PI);
    CNSTNT.PIINV = 1.0 / CNSTNT.PI;
    CNSTNT.RADIAN = 180.0 / CNSTNT.PI;
    CNSTNT.DEGREE = 1.0 / CNSTNT.RADIAN;
    //
    //     FOLLOWING FROM  E. R. COHEN AND B. N. TAYLOR, J. PHYS. AND CHEM.
    //     REF. DATA, 2, 663 (1973).
    //
    CNSTNT.HBARC = 197.328580;
    CNSTNT.AMUMEV = 931.50160;
    CNSTNT.AFINE = 137.036040;
    //
    //     UNDERFLOW AND OVERFLOW LIMITS
    //
    CNSTNT.BIGEST = 1.E+300;
    //  bignum must fit in real*4
    CNSTNT.BIGNUM = 1.E+30;
    //
    CNSTNT.SMLEST = 1.0 / CNSTNT.BIGEST;
    CNSTNT.SMLNUM = 1.0 / CNSTNT.BIGNUM;
    CNSTNT.BIGLOG = std::log(CNSTNT.BIGNUM);
    //
    //     TO TRAP A TROLLOP TRAVELLING ON A TROLLY
    //     (I.E. RESET POINTER TO BOUND STATE POLENTIAL.)
    //
    INTRNL.IPOT = 0;
    //
    //     SETUP DEFAULT VALUES
    //
    FLOAT_common.DWCUT = 1.0E-3;
    INTGER.NCOSIN = 256;
    INTGER.NBACK = 4;
    FLOAT_common.DELTVK = 0.050;
    INTGER.MAXITR = 10;
    SWITCH.IFIT = 1;
    SWITCH.IVRTEX = 1;
    FLOAT_common.ACCURA = 1.0E-12;
    FLOAT_common.BNDASY = 20.0;
    FLOAT_common.SCTASY = -20.0;
    FLOAT_common.STEP1R = 1.0;
    FLOAT_common.STEP1I = 1.0;
    INTGER.ISIZE = 200000;
    FLOAT_common.STEPSZ = 0.10;
    //  DEFAULT VALUES FOR ANGLES
    FLOAT_common.ANGMIN = 0;
    FLOAT_common.ANGMAX = 90;
    FLOAT_common.ANGSTP = 1;
    INTGER.LSTEP = 1;
    INTGER.NPPHI = 10;
    INTGER.NPSUM = 15;
    INTGER.NPDIF = 10;
    INTGER.NVPOLY = 3;
    FLOAT_common.SUMPTS = 6;
    FLOAT_common.SPAMP = 1.0;
    FLOAT_common.SPAMT = 1.0;
    FLOAT_common.SPFACP = 1;
    FLOAT_common.SPFACT = 1;
    INTGER.NAITKN = 4;
    INTGER.IPRINT = 10001;
    if (SWITCH.ITSO == 1) INTGER.IPRINT = 01;
    INTGER.LOOKST = 250;
    INTRNL.ICHANB = 0;
    INTRNL.ICHANW = 0;
    INTRNL.IDONE = 0;
    INTRNL.ISTRIP = 0;
    INTGER.MAPSUM = 2;
    INTGER.MAPDIF = 1;
    INTGER.MAPPHI = 2;
    FLOAT_common.GAMSUM = 1;
    FLOAT_common.GAMDIF = 5;
    FLOAT_common.GAMPHI = 1.0E-6;
    FLOAT_common.AMDMLT = 2;
    FLOAT_common.PHIMID = 0.5;
    INTGER.NPHIAD = 4;
    SWITCH.NUCONL = 3;
    INTRNL.IHSAVE = 0;
    FLOAT_common.ALMXMT = 1.6;
    INTGER.LMAXAD = 30;
    FLOAT_common.ALMNMT = 0.6;
    INTGER.LMINSB = 20;
    INTGER.MAXLEX = 100;
    INTGER.LEBACK = 15;
    //
    WAVCOM.NWAVEF = 0;
    //     EXCITATION ENERGIES ARE 0 BY DEFAULT
    //     PARITY = 0 MEANS UNDEFINED
    //
    for (int I = 1; I <= 5; I++) {
        INTGER.PARITS[I] = 0;
        FLOAT_common.EXS[I] = 0;
    }
    INTGER.PARITY = 0;
    //
    INTGER.NBLKSZ = 5632;
    //
    //     VAX/VMS:  THE DEFAULT CPU TIME LIMIT IS ONE HOUR.
    //
    for (int I = 1; I <= 4; I++) {
        INTRNL.LSPECS[I] = *reinterpret_cast<int*>(&INTRNL.NOTDEF);
        INTRNL.NODESP[I] = *reinterpret_cast<int*>(&INTRNL.NOTDEF);
    }
    for (int I = 1; I <= 2; I++) {
        INTRNL.LSPCPT[I] = *reinterpret_cast<int*>(&INTRNL.NOTDEF);
        INTRNL.NODEPT[I] = *reinterpret_cast<int*>(&INTRNL.NOTDEF);
        INTGER.PARIPT[I] = 0;
        INTRNL.EBNDS[I] = INTRNL.UNDEF;
    }
    //
    //     POTENTIALS ARE 0 BY DEFAULT
    //
    FLOAT_common.V = 0;
    FLOAT_common.VI = 0;
    FLOAT_common.VSO = 0;
    FLOAT_common.VSOI = 0;
    FLOAT_common.TAU = 0;
    FLOAT_common.TAUI = 0;
    FLOAT_common.VSI = 0;
    for (int I = 1; I <= 6; I++) {
        VTEN()[I - 1] = 0;
    }
    //
    //  ENERGY DEPENDENCE ZERO BY DEFAULT
    //
    FLOAT_common.R0E = 0;
    FLOAT_common.RI0E = 0;
    FLOAT_common.AE = 0;
    FLOAT_common.AIE = 0;
    FLOAT_common.VE = 0;
    FLOAT_common.VIE = 0;
    FLOAT_common.R0ESQ = 0;
    FLOAT_common.RI0ESQ = 0;
    FLOAT_common.AESQ = 0;
    FLOAT_common.AIESQ = 0;
    FLOAT_common.VESQ = 0;
    FLOAT_common.VIESQ = 0;
    //
    //     POTENTIALS ARE LINEAR WOODS SAXONS
    //
    FLOAT_common.POWRL = 1;
    FLOAT_common.POWIM = 1;
    //
    //     INDICATE NO LINKULES IN USE.
    //
    LNKBLK.NUMLNK = NUMLINKULES;
    LNKBLK.IUNIQU = 100;
    for (int I = 1; I <= NUMLINKULES; I++) {
        LNKBLK.LNKADR[I][3] = 0;
    }
    //
    //     FITTING KEYWORDS
    //
    FLOAT_common.FITACC = 1.0E-3;
    FLOAT_common.FITMUL = 500;
    FLOAT_common.FITRAT = 1.0E-4;
    INTGER.MODEFT = 1;
    INTGER.NUMRAN = 0;
    INTGER.MAXFUN = 50;
    INTGER.IREINI = 0;
    SWITCH.IFITER = 6;
    FLOAT_common.DERIVS = 1.0E-6;
    //
    //     INELASTIC SCATTERING
    //
    FLOAT_common.ACCINE = 1.0E-5;
    FLOAT_common.COULML = 1;
    INTGER.NPCOUL = 8;
    INTGER.MXCOUL = 80;
    //
    //     COUPLED CHANNELS
    //
    SWITCH.IEXTYP = 1;
    CCBLK.ICHNDF = 0;
    CCBLK.ICHNCP = 0;
    CCBLK.ITRANS = 0;
    CCBLK.IBASDF = 0;
    INTGER.NFIROF = 1;
    //
    return;
}
