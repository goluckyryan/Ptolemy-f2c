// prbprt_part1.cpp — PRBPRT subroutine, Part 1 of 2
// Translated from source.f lines 28809–29422 (through statement 480)
// Makes consistency checks and prints problem summary.
//
// The function brace is NOT closed here — continues in prbprt_part2.cpp.

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

void PRBPRT(int& IRET)
{
    // =========================================================================
    // Local variable declarations
    // =========================================================================

    // Static DATA arrays (initialized once)
    //   INTEGER  IND(5) / 3, 1, 2, 4, 5 /
    static const int IND[6] = { 0, 3, 1, 2, 4, 5 };  // 1-based

    //   CHARACTER*8  PARTIC(5) / 'A', 'B', 'BIGA', 'BIGB', 'X' /
    static const char PARTIC[6][9] = {
        "",
        "A       ", "B       ", "BIGA    ", "BIGB    ", "X       "
    };

    //   CHARACTER*8  XWORDS(2) / ' ', ';   X' /
    static const char XWORDS[3][9] = {
        "",
        "        ", ";   X   "
    };

    //   CHARACTER*8  WORDS(4) / 'PROJECTI', 'LE', 'TARGET', ' ' /
    static const char WORDS[5][9] = {
        "",
        "PROJECTI", "LE      ", "TARGET  ", "        "
    };

    //   CHARACTER*8  BETNAM(4) / 'BETA', 'BETACOUL', 'BETARATS', 'BELX ' /
    static const char BETNAM[5][9] = {
        "",
        "BETA    ", "BETACOUL", "BETARATS", "BELX    "
    };

    //   CHARACTER*8 PWORD(3) / '     -1', 'UNKNOWN', '     +1' /
    static const char PWORD[4][9] = {
        "",
        "     -1 ", "UNKNOWN ", "     +1 "
    };

    // LOGICAL variables (Fortran LOGICAL)
    int PRNTSW;   // LOGICAL
    int BETSWS[5] = {0,0,0,0,0}; // LOGICAL, 1-based
    int TRANSW;   // LOGICAL
    int COUPSW;   // LOGICAL

    // Local scalars
    int I, II, IB, IPRNT;
    int JBPMIN, JBPMAX, JBTMIN, JBTMAX;
    int IPARIT;
    int JIN, JOUT;
    int LMININ;
    int ILVALU, LLVALU, L1, LEN;
    int LIBASE = 0;
    int NUMLXI;
    int JXMIN, JXMAX;
    double BETA, CHARGE, ATERM, TEMP;
    double EMUL, R2MASS, VC;
    double R2S[5] = {};  // 1-based local (Fortran EQUIVALENCE to RBNDS not needed in C++)
    int K;
    int IBETNR, LATERM, LBETA, LBETAC, LBETAR, LBELX_loc, LBETNR;
    int ITEST, NUMPAR;
    (void)LIBASE; // may be unused in transfer path

    // =========================================================================
    // EQUIVALENCE / COMMON aliases used in this routine
    // =========================================================================

    // AMS(I): 1-based pointer to AMA (AMA, AMB, AMBIGA, AMBIGB, AMX)
    // AMS(I) uses FLOAT_arr() to handle struct padding correctly
    // AMS_ptr is a 1-based pointer: AMS_ptr[1]=AMA, AMS_ptr[2]=AMB, etc.
    double* AMS_ptr = AMS();  // from ptolemy_commons.h

    // JS(I): In JBLOCK, the Fortran /JBLOCK/ has J then JS(5), JP, JSP, JST
    // JS(1)=JA, JS(2)=JB, JS(3)=JBIGA, JS(4)=JBIGB, JS(5)=JX
    // In the C++ header, JBLOCK.JS is double[6] (1-based).
    // Fortran EQUIVALENCE: JS(1)<=>JA, JS(2)<=>JB, JS(3)<=>JBIGA,
    //   JS(4)<=>JBIGB, JS(5)<=>JX
    // Access as JBLOCK.JS[I] for 1-based I.
    // For integer comparisons (JS(I).NE.NOTDEF), cast to int:
    //   (int)JBLOCK.JS[I] != NOTDEF_INT
    // But Fortran JS is DOUBLE — the EQUIVALENCE in Fortran source says
    // JS(1)=JA etc. and the NOTDEF check is on an integer-stored value.
    // In this code we follow the established pattern: compare as int.
    // Helper lambda for JS integer comparison:
    auto JS_NE_NOTDEF = [&](int i) -> bool {
        // JS values stored as double but compared as integer (Fortran quirk)
        return (int)JBLOCK.JS[i] != NOTDEF_INT;
    };

    // IZS: INTGER.IZS[I], 1-based (IZS[1]..IZS[5])
    int* IZS = INTGER.IZS;

    // PARITS: INTGER.PARITS[I], 1-based
    int* PARITS = INTGER.PARITS;

    // LBNDS equivalenced to LBP: LBNDS(1)=LBP, LBNDS(2)=LBT
    int* LBNDS = &FORMF.LBP;

    // NBNDS equivalenced to NBP: NBNDS(1)=NBP, NBNDS(2)=NBT
    int* NBNDS = &FORMF.NBP;

    // R2S(4) equivalenced to RBNDS(1): R2S(1)=RBNDS(1)..R2S(4) is RCBNDS overrun
    // From commons: inline double* R2S() { return &FORMF.RBNDS[1]; }
    // R2S(1)=RBNDS[1], R2S(2)=RBNDS[2], R2S(3)=RBNDS[3], R2S(4)=RCBNDS[1]
    // We use R2S() inline accessor.

    // LCRITS: KANDM.LCRITS[1], KANDM.LCRITS[2]
    // LNKAD2: LOCFIT.linkule.LNKAD2[I][J][K], 1-based

    // Convenient references to frequently used COMMON members
    double& UNDEF     = INTRNL.UNDEF;
    int&    ISTRIP    = INTRNL.ISTRIP;
    int&    IEXCIT    = INTRNL.IEXCIT;
    double& R0MASS    = INTRNL.R0MASS;
    double& ELAB      = FLOAT_arr(14);     // position 14 (before padding)
    double& AMA       = FLOAT_arr(25);     // position 25
    double& AMB       = FLOAT_arr(26);
    double& AMBIGA    = FLOAT_arr(27);
    double& AMBIGB    = FLOAT_arr(28);
    double& AMX       = FLOAT_arr(29);
    double& E         = FLOAT_common.E;
    double& R         = FLOAT_common.R;
    double& RI        = FLOAT_common.RI;
    double& R0        = FLOAT_common.R0;
    double& RI0       = FLOAT_common.RI0;
    double& RC        = FLOAT_common.RC;
    double& RC0       = FLOAT_common.RC0;
    double& RSI       = FLOAT_common.RSI;
    double& RSI0      = FLOAT_common.RSI0;
    double& ALMNMT    = FLOAT_common.ALMNMT;
    double& ALMXMT    = FLOAT_common.ALMXMT;
    double& V         = FLOAT_common.V;
    double& VI        = FLOAT_common.VI;
    double& VE        = FLOAT_common.VE;
    double& VIE       = FLOAT_common.VIE;
    double& VESQ      = FLOAT_common.VESQ;
    double& VIESQ     = FLOAT_common.VIESQ;
    double& AE        = FLOAT_common.AE;
    double& AIE       = FLOAT_common.AIE;
    double& AESQ      = FLOAT_common.AESQ;
    double& AIESQ     = FLOAT_common.AIESQ;
    double& R0E       = FLOAT_common.R0E;
    double& RI0E      = FLOAT_common.RI0E;
    double& R0ESQ     = FLOAT_common.R0ESQ;
    double& RI0ESQ    = FLOAT_common.RI0ESQ;
    double& RCP       = FLOAT_common.RCP;
    double& RCT       = FLOAT_common.RCT;
    // A and AI for deformation lengths
    double& A_flt     = FLOAT_common.A;
    double& AI_flt    = FLOAT_common.AI;
    // TVR, TVI, TAR, TAI from TEMPVS
    double& TVR       = TEMPVS.TVR;
    double& TVI       = TEMPVS.TVI;
    double& TAR       = TEMPVS.TAR;
    double& TAI       = TEMPVS.TAI;
    int&    PROBLM    = SWITCH.PROBLM;
    int&    IEPOW     = SWITCH.IEPOW;
    int&    NUCONL    = SWITCH.NUCONL;
    int&    IPRINT    = INTGER.IPRINT;
    int&    LX        = INTGER.LX;
    int&    LMIN      = INTGER.LMIN;
    int&    LMAX      = INTGER.LMAX;
    int&    LMINSB    = INTGER.LMINSB;
    int&    LMAXAD    = INTGER.LMAXAD;
    int&    LSTEP     = INTGER.LSTEP;
    int&    L         = INTGER.L;
    int&    LBP       = FORMF.LBP;
    int&    LBT       = FORMF.LBT;
    int&    NBP       = FORMF.NBP;
    int&    NBT       = FORMF.NBT;
    int&    JBP       = FORMF.JBP;
    int&    JBT       = FORMF.JBT;
    int&    LCRIT     = KANDM.LCRIT;
    int&    LXMIN     = INELCM.LXMIN;
    int&    LXMAX     = INELCM.LXMAX;
    int&    LXSTEP    = INELCM.LXSTEP;
    int&    NUMLX     = INELCM.NUMLX;
    int&    NUMLIS    = INELCM.NUMLIS;
    int&    ILIS      = INELCM.ILIS;
    int&    JPMIN     = INELCM.JPMIN;
    int&    JPMAX     = INELCM.JPMAX;
    int&    JTMIN     = INELCM.JTMIN;
    int&    JTMAX     = INELCM.JTMAX;
    int&    LSKIP     = INELCM.LSKIP;
    int&    ONELSW    = INELCM.ONELSW;
    int&    IDENSW    = INELCM.IDENSW;
    double& HBARC     = CNSTNT.HBARC;
    double& AFINE     = CNSTNT.AFINE;

    // JS references (Fortran JS(1)=JA, JS(2)=JB, JS(3)=JBIGA, JS(4)=JBIGB, JS(5)=JX)
    // In C++ JBLOCK.JS[1..5]
    // Individual named references matching Fortran EQUIVALENCE
    double& JA_ref    = JBLOCK.JS[1];  // JS(1) = JA
    double& JB_ref    = JBLOCK.JS[2];  // JS(2) = JB
    double& JBIGA_ref = JBLOCK.JS[3];  // JS(3) = JBIGA
    double& JBIGB_ref = JBLOCK.JS[4];  // JS(4) = JBIGB
    double& JX_ref    = JBLOCK.JS[5];  // JS(5) = JX

    // Cast macros for half-integer spin values stored as doubles
    // Fortran prints them as integer (e.g., FORMAT I5) — they hold 2*J
    // Use (int) casts when printing or comparing as integers
    #define JA_INT    ((int)JA_ref)
    #define JB_INT    ((int)JB_ref)
    #define JBIGA_INT ((int)JBIGA_ref)
    #define JBIGB_INT ((int)JBIGB_ref)
    #define JBP_INT   JBP
    #define JBT_INT   JBT
    #define JX_INT    ((int)JX_ref)

    // EXS: FLOAT_common.EXS[1..5], 1-based
    // LSKIPS, ISTATS from WAVCOM: WAVCOM.LSKIPS[1..2], WAVCOM.ISTATS[1..2]

    // LNKAD2(I,J,K): LOCFIT.linkule.LNKAD2[I][J][K], 1-based
    // Access pattern: LOCFIT.linkule.LNKAD2[i][j][k]

    // =========================================================================
    //      IRET = 1
    // =========================================================================
    IRET = 1;

    //      IPRNT = MOD( IPRINT, 10 )
    //      PRNTSW = IPRNT .GE. 1
    IPRNT = MOD(IPRINT, 10);
    PRNTSW = (IPRNT >= 1);

    //      IF ( PRNTSW )  write (6, 3) REACT, ELAB, HEADER
    //  3   FORMAT ( '1', T49, 'P T O L E M Y' /
    //    1  '   CONSTRUCTION OF THE INTEGRATION GRIDS', T70,
    //    2  'MERELY TO CONCEIVE OF SUCH THINGS MAKES THEM APPEAR RIDICULOUS'
    //    3   /
    //    4  '0', 45A1, 'ELAB =', F7.2, ' MEV', 5X, 65A1 / )
    //
    // FORMAT '1' = page break (just print content); '0' = blank line before
    if (PRNTSW) {
        // FORMAT '1' = literal '1', T49 = 48 spaces after '1', then title
        std::printf("1%47sP T O L E M Y\n", "");
        // '   CONSTRUCTION...', T70 = 29 fill spaces (text is 40 chars, T70 from col 41)
        std::printf("   CONSTRUCTION OF THE INTEGRATION GRIDS%29sMERELY TO CONCEIVE OF SUCH THINGS MAKES THEM APPEAR RIDICULOUS\n", "");
        // FORMAT '0' prefix, 45A1, ELAB, 65A1
        std::printf("0");
        for (int _i = 1; _i <= 45; _i++) std::printf("%c", HEDCOM.REACT[_i]);
        std::printf("ELAB =%7.2f MEV     ", ELAB);
        for (int _i = 1; _i <= 65; _i++) std::printf("%c", HEDCOM.HEADER[_i]);
        std::printf("\n\n");
    }

    // =========================================================================
    //      TRANSW = ISTRIP .NE. 0
    //      COUPSW = PROBLM .EQ. 24
    // =========================================================================
    TRANSW = (ISTRIP != 0);
    COUPSW = (PROBLM == 24);

    //      LXSTEP = 1
    //      IF ( TRANSW )  GO TO 40
    LXSTEP = 1;
    if (TRANSW) goto stmnt_40;

    //      IF ( .NOT. COUPSW )  LXSTEP = 2
    //      LBP = IABS(JB-JA)/2
    //      LBT = IABS(JBIGB-JBIGA)/2
    if (!COUPSW) LXSTEP = 2;
    LBP = std::abs(JB_INT - JA_INT) / 2;
    LBT = std::abs(JBIGB_INT - JBIGA_INT) / 2;

    // =========================================================================
    //  40  DO 49  I = 1, 4
    // =========================================================================
stmnt_40:
    for (I = 1; I <= 4; I++) {
        //         IF ( AMS(I) .NE. UNDEF )  GO TO 42
        if (AMS_ptr[I] != UNDEF) goto stmnt_42;
        //         write (6, 41)  PARTIC(I)
        // 41      FORMAT ( '0***** ERROR:  M', A4, ' IS NOT DEFINED.' )
        std::printf("\n***** ERROR:  M%.4s IS NOT DEFINED.\n", PARTIC[I]);
        IRET = 0;
stmnt_42:
        //         IF ( JS(I) .NE. NOTDEF )  GO TO 44
        if (JS_NE_NOTDEF(I)) goto stmnt_44;
        //         write (6, 43)  PARTIC(I)
        // 43      FORMAT ( '0***** ERROR:  J', A4, ' IS NOT DEFINED.' )
        std::printf("\n***** ERROR:  J%.4s IS NOT DEFINED.\n", PARTIC[I]);
        IRET = 0;
stmnt_44:
        //         IF ( IZS(I) .NE. NOTDEF )  GO TO 49
        if (IZS[I] != NOTDEF_INT) goto stmnt_49;
        //         write (6, 45)  PARTIC(I)
        // 45      FORMAT ( '0***** ERROR:  Z', A4, ' IS NOT DEFINED.' )
        std::printf("\n***** ERROR:  Z%.4s IS NOT DEFINED.\n", PARTIC[I]);
        IRET = 0;
stmnt_49:;
    }

    //      IF ( AMX .EQ. UNDEF )  AMX = DABS( AMA-AMB )
    if (AMX == UNDEF) AMX = std::abs(AMA - AMB);

    //      IF ( DABS(AMA+AMBIGA-AMB-AMBIGB) .LT. .3 )  GO TO 55
    if (std::abs(AMA + AMBIGA - AMB - AMBIGB) < 0.3) goto stmnt_55;
    //      write (6, 53) AMBIGA, AMA, AMB, AMBIGB
    // 53   FORMAT ( '0**** THE REACTION DOES NOT CONSERVE MASS:' /
    //    1  ' MBIGA(MA, MB)MBIGB =', 4G15.5 )
    std::printf("\n**** THE REACTION DOES NOT CONSERVE MASS:\n");
    std::printf(" MBIGA(MA, MB)MBIGB =%15.5G%15.5G%15.5G%15.5G\n",
                AMBIGA, AMA, AMB, AMBIGB);
    IRET = 0;

stmnt_55:
    //  55  IF ( IZS(1)+IZS(3) .EQ. IZS(2)+IZS(4) )  GO TO 60
    if (IZS[1] + IZS[3] == IZS[2] + IZS[4]) goto stmnt_60;
    //      write (6, 57) IZS(3), IZS(1), IZS(2), IZS(4)
    // 57   FORMAT ( '0**** THE REACTION DOES NOT CONSERVE CHARGE:' /
    //    1  ' ZBIGA(ZA, ZB)ZBIGB =', 4I8 )
    std::printf("\n**** THE REACTION DOES NOT CONSERVE CHARGE:\n");
    std::printf(" ZBIGA(ZA, ZB)ZBIGB =%8d%8d%8d%8d\n",
                IZS[3], IZS[1], IZS[2], IZS[4]);
    IRET = 0;

stmnt_60:
    //  60  IF ( DABS( AMX - DABS(AMA-AMB) )  .LT. .3 )  GO TO 65
    if (std::abs(AMX - std::abs(AMA - AMB)) < 0.3) goto stmnt_65;
    //      write (6, 63)
    // 63   FORMAT ( '0**** THE EXCHANGED MASS IS INCONSISTANT WITH THE ',
    //    1  'OTHER FOUR MASSES.' )
    std::printf("\n**** THE EXCHANGED MASS IS INCONSISTANT WITH THE OTHER FOUR MASSES.\n");
    IRET = 0;

stmnt_65:
    //  65  IF ( AMX .GT. .7 )  GO TO 90
    if (AMX > 0.7) goto stmnt_90;
    //      IF ( .NOT. TRANSW )  GO TO 70
    if (!TRANSW) goto stmnt_70;
    //      write (6, 68) AMX
    // 68   FORMAT ( '0**** THE EXCHANGED MASS MUST NOT BE 0 OR NEGATIVE:',
    //    1  G15.5 )
    std::printf("\n**** THE EXCHANGED MASS MUST NOT BE 0 OR NEGATIVE:%15.5G\n", AMX);
    IRET = 0;
    goto stmnt_90;

stmnt_70:
    //  70  IF ( IZS(5) .EQ. 0 )  GO TO 90
    if (IZS[5] == 0) goto stmnt_90;
    //      write (6, 73) IZS(5)
    // 73   FORMAT ( '0**** CANNOT TRANSFER CHARGE WITHOUT MASS TRANSFER:',
    //    1  ' AMX = 0, ZX =', I4 )
    std::printf("\n**** CANNOT TRANSFER CHARGE WITHOUT MASS TRANSFER: AMX = 0, ZX =%4d\n",
                IZS[5]);
    IRET = 0;

stmnt_90:
    //  90  CONTINUE
    //      IF ( COUPSW )  GO TO 400
    if (COUPSW) goto stmnt_400;
    //      IF ( .NOT. TRANSW )  GO TO 200
    if (!TRANSW) goto stmnt_200;

    // =========================================================================
    //      FOR TRANSFER CANNOT HAVE NUCLEAR CORE CORRECTIONS FOR LINKULES
    // =========================================================================
    //      IF ( LNKAD2(3,1,1) .EQ. 0  .AND.  LNKAD2(3,2,1) .EQ. 0 )
    //    1    GO TO 200
    if (LOCFIT.linkule.LNKAD2[3][1][1] == 0  &&
        LOCFIT.linkule.LNKAD2[3][2][1] == 0)
        goto stmnt_200;
    //      IF ( NUCONL .LT. 3 )  GO TO 200
    if (NUCONL < 3) goto stmnt_200;
    //      write (6, 103)
    // 103  FORMAT ( '0**** WARNING - NUCLEAR CORE-CORE CORRECTIONS CANNOT',
    //    1  ' BE COMPUTED FOR LINKULE-COMPUTED POTENTIALS.' /
    //    2  6X, 'ONLY COULOMB CORE-CORE CORRECTIONS WILL BE USED.' )
    std::printf("\n**** WARNING - NUCLEAR CORE-CORE CORRECTIONS CANNOT"
                " BE COMPUTED FOR LINKULE-COMPUTED POTENTIALS.\n");
    std::printf("      ONLY COULOMB CORE-CORE CORRECTIONS WILL BE USED.\n");
    //      NUCONL = 2
    NUCONL = 2;

    // =========================================================================
    //      CHECK THAT JBP AND JBT ARE INSIDE ALLOWED BOUNDS
    // =========================================================================
stmnt_200:
    // 200  JBPMIN=IABS(JA-JB)
    //      JBPMAX=JA+JB
    //      JBTMIN=IABS(JBIGA-JBIGB)
    //      JBTMAX=JBIGA+JBIGB
    JBPMIN = std::abs(JA_INT - JB_INT);
    JBPMAX = JA_INT + JB_INT;
    JBTMIN = std::abs(JBIGA_INT - JBIGB_INT);
    JBTMAX = JBIGA_INT + JBIGB_INT;

    //      IF ( .NOT. TRANSW )  GO TO 270
    if (!TRANSW) goto stmnt_270;

    // =========================================================================
    //      FOR TRANSFER, CHECK (JBP, JBT) PAIR
    // =========================================================================
    //      IF ( JBP .GE. JBPMIN  .AND.  JBP .LE. JBPMAX
    //    1    .AND.  JBT .GE. JBTMIN  .AND.  JBT .LE. JBTMAX )  GO TO 260
    if (JBP >= JBPMIN  &&  JBP <= JBPMAX  &&
        JBT >= JBTMIN  &&  JBT <= JBTMAX)
        goto stmnt_260;

    //      write (6, 253) JA, JB, JBIGA, JBIGB, JBPMIN, JBPMAX,
    //    &   JBTMIN, JBTMAX
    // 253  FORMAT ('0***** ERROR IN INPUT:' /
    //    1  ' JA, JB, JBIGA, JBIGB =', 4(I5, '/2') /
    //    2  ' RESULTS IN' /
    //    3  ' JBPMIN, JBPMAX, JBTMIN, JBTMAX =', 4(I5, '/2') )
    std::printf("\n***** ERROR IN INPUT:\n");
    std::printf(" JA, JB, JBIGA, JBIGB =%5d/2%5d/2%5d/2%5d/2\n",
                JA_INT, JB_INT, JBIGA_INT, JBIGB_INT);
    std::printf(" RESULTS IN\n");
    std::printf(" JBPMIN, JBPMAX, JBTMIN, JBTMAX =%5d/2%5d/2%5d/2%5d/2\n",
                JBPMIN, JBPMAX, JBTMIN, JBTMAX);
    IRET = 0;

    //      IF ( JBP .EQ. NOTDEF  .OR.  JBT .EQ. NOTDEF )  GO TO 255
    if (JBP == NOTDEF_INT  ||  JBT == NOTDEF_INT) goto stmnt_255;
    //      write (6, 254) JBP, JBT
    // 254  FORMAT ( ' BUT  JBP, JBT =', 2(I5, '/2'), '  AND ARE OUT OF ',
    //    1  'BOUNDS.' )
    std::printf(" BUT  JBP, JBT =%5d/2%5d/2  AND ARE OUT OF BOUNDS.\n",
                JBP, JBT);
    goto stmnt_260;

stmnt_255:
    //  255 write (6, 257) JBP, JBT
    // 257  FORMAT ( ' THUS THERE IS MORE THAN ONE POSSIBILITY FOR',
    //    1  ' JP FOR THE PROJECTILE OR TARGET BOUND STATE.' /
    //    2  ' IN SUCH AMBIGUOUS CASES JP MUST BE EXPLICITLY INDICATED.' /
    //    3  ' WE PRESENTLY HAVE  JP(PROJ) = ', I3, '/2,', 4X,
    //    4    'JP(TARG) = ', I3, '/2' /
    //    5  ' WHERE THE *** INDICATES THE AMBIGUOUS CASE.' )
    std::printf(" THUS THERE IS MORE THAN ONE POSSIBILITY FOR"
                " JP FOR THE PROJECTILE OR TARGET BOUND STATE.\n");
    std::printf(" IN SUCH AMBIGUOUS CASES JP MUST BE EXPLICITLY INDICATED.\n");
    std::printf(" WE PRESENTLY HAVE  JP(PROJ) = %3d/2,    JP(TARG) = %3d/2\n",
                JBP, JBT);
    std::printf(" WHERE THE *** INDICATES THE AMBIGUOUS CASE.\n");

    // =========================================================================
    //      FIND LX RANGE FOR TRANSFER AND MAKE NECESSARY CHECKS
    // =========================================================================
stmnt_260:
    // 260  LXMIN=IABS(JBP-JBT)
    //      LXMAX=JBP+JBT
    //      IF(LXMIN.LT.IABS(2*(LBP-LBT)))LXMIN=IABS(2*(LBP-LBT))
    //      IF(LXMAX.GT.(2*(LBP+LBT)))LXMAX=2*(LBP+LBT)
    LXMIN = std::abs(JBP - JBT);
    LXMAX = JBP + JBT;
    if (LXMIN < std::abs(2*(LBP-LBT))) LXMIN = std::abs(2*(LBP-LBT));
    if (LXMAX > 2*(LBP+LBT)) LXMAX = 2*(LBP+LBT);

    //      LXMAX=LXMAX/2
    //      LXMIN=LXMIN/2
    //      IPARIT = 1 - 2*MOD(LBP+LBT, 2)
    //      JPMIN = JBP
    //      JPMAX = JBP
    //      JTMIN = JBT
    //      JTMAX = JBT
    LXMAX = LXMAX / 2;
    LXMIN = LXMIN / 2;
    IPARIT = 1 - 2*MOD(LBP+LBT, 2);
    JPMIN = JBP;
    JPMAX = JBP;
    JTMIN = JBT;
    JTMAX = JBT;

    //      IF ( LXMAX .GE. LXMIN )  GO TO 300
    if (LXMAX >= LXMIN) goto stmnt_300;
    //      write (6, 263) JBP, JBT, LBP, LBT, LXMIN, LXMAX
    // 263  FORMAT ( '0**** INVALID BOUND STATE ANGULAR MOMENTA:' /
    //    1  ' JP, JT, LP, LT =', 2(I5, '/2'), 2I5 /
    //    2  ' RESULTS IN  LXMIN, LXMAX =', 2I8 )
    std::printf("\n**** INVALID BOUND STATE ANGULAR MOMENTA:\n");
    std::printf(" JP, JT, LP, LT =%5d/2%5d/2%5d%5d\n",
                JBP, JBT, LBP, LBT);
    std::printf(" RESULTS IN  LXMIN, LXMAX =%8d%8d\n", LXMIN, LXMAX);
    IRET = 0;
    goto stmnt_300;

    // =========================================================================
    //      SETUP LX RANGE FOR INELASTIC SCATTERING
    // =========================================================================
stmnt_270:
    // 270  IF ( EXS(4) .NE. EXS(3) )  GO TO 275
    if (FLOAT_common.EXS[4] != FLOAT_common.EXS[3]) goto stmnt_275;

    // PROJECTILE EXCITATION
    //      IEXCIT = 1
    //      LXMAX = JBPMAX
    //      LXMIN = JBPMIN
    //      JPMIN = JBPMIN
    //      JPMAX = JBPMAX
    //      JTMIN = 0
    //      JTMAX = 0
    //      JIN = JA
    //      JOUT = JB
    //      IPARIT = PARITS(1)*PARITS(2)
    IEXCIT = 1;
    LXMAX = JBPMAX;
    LXMIN = JBPMIN;
    JPMIN = JBPMIN;
    JPMAX = JBPMAX;
    JTMIN = 0;
    JTMAX = 0;
    JIN  = JA_INT;
    JOUT = JB_INT;
    IPARIT = PARITS[1] * PARITS[2];

    //      IF ( EXS(2) .NE. EXS(1) )  GO TO 280
    if (FLOAT_common.EXS[2] != FLOAT_common.EXS[1]) goto stmnt_280;
    //      write (6, 273)
    // 273  FORMAT ( '0**** E*B OR E*BIGB MUST BE DEFINED' )
    std::printf("\n**** E*B OR E*BIGB MUST BE DEFINED\n");
    IRET = 0;
    goto stmnt_280;

stmnt_275:
    // TARGET EXCITATION
    // 275  IEXCIT = 2
    //      LXMAX = JBTMAX
    //      LXMIN = JBTMIN
    //      JPMIN = 0
    //      JPMAX = 0
    //      JTMIN = JBTMIN
    //      JTMAX = JBTMAX
    //      JIN = JBIGA
    //      JOUT = JBIGB
    //      IPARIT = PARITS(3)*PARITS(4)
    IEXCIT = 2;
    LXMAX = JBTMAX;
    LXMIN = JBTMIN;
    JPMIN = 0;
    JPMAX = 0;
    JTMIN = JBTMIN;
    JTMAX = JBTMAX;
    JIN  = JBIGA_INT;
    JOUT = JBIGB_INT;
    IPARIT = PARITS[3] * PARITS[4];

    //      IF ( EXS(2) .EQ. EXS(1) )  GO TO 280
    if (FLOAT_common.EXS[2] == FLOAT_common.EXS[1]) goto stmnt_280;
    //      write (6, 278)
    // 278  FORMAT ( '0**** CANNOT SIMULTANEOUSLY EXCITE BOTH PROJECTILE',
    //    1  ' AND TARGET.' )
    std::printf("\n**** CANNOT SIMULTANEOUSLY EXCITE BOTH PROJECTILE AND TARGET.\n");
    IRET = 0;

stmnt_280:
    // 280  IF ( MOD(LXMIN, 2) .EQ. 0 )  GO TO 285
    if (MOD(LXMIN, 2) == 0) goto stmnt_285;
    //      write (6, 283) LXMIN, LXMAX, IEXCIT
    // 283  FORMAT ( '0**** INELASTIC EXCITATION MUST INVOLVE AN',
    //    1  ' INTEGER CHANGE IN NUCLEAR SPINS', 3I8 )
    std::printf("\n**** INELASTIC EXCITATION MUST INVOLVE AN"
                " INTEGER CHANGE IN NUCLEAR SPINS%8d%8d%8d\n",
                LXMIN, LXMAX, IEXCIT);
    IRET = 0;

stmnt_285:
    // 285  LXMIN = LXMIN/2
    //      LXMAX = LXMAX/2
    LXMIN = LXMIN / 2;
    LXMAX = LXMAX / 2;

    //      IF ( IPARIT .EQ. 0 )  GO TO 290
    if (IPARIT == 0) goto stmnt_290;
    {
        //      I = (IPARIT+3)/2
        //      LXMIN = LXMIN + MOD(I+LXMIN, 2)
        //      LXMAX = LXMAX - MOD(I+LXMAX, 2)
        //      LBNDS(IEXCIT) = LXMIN
        int Itmp = (IPARIT + 3) / 2;
        LXMIN = LXMIN + MOD(Itmp + LXMIN, 2);
        LXMAX = LXMAX - MOD(Itmp + LXMAX, 2);
        LBNDS[IEXCIT] = LXMIN;
    }

stmnt_290:
    // 290  IF ( LXMAX .EQ. LXMIN  .OR.  LX .NE. NOTDEF )  GO TO 300
    if (LXMAX == LXMIN  ||  LX != NOTDEF_INT) goto stmnt_300;
    //      IF ( IPARIT .NE. 0 )  GO TO 295
    if (IPARIT != 0) goto stmnt_295;
    //      write (6, 293)
    // 293  FORMAT ( '0', 120('*') /
    //    1  ' **** WARNING, THE PARITIES OF THE STATES MUST BE KNOWN',
    //    2    ' TO UNIQUELY DETERMIN THE MULTIPOLARITY OF THE EXCITATION' /
    //    3  ' **** PTOLEMY ASSUMES THAT THE CHANGE IN THE INTRINSIC PARITY',
    //    4    ' IS JUST  (-1)**|J(IN)-J(OUT)|.' /
    //    5  ' ', 120('*') )
    std::printf("\n");
    for (int _i = 0; _i < 120; _i++) std::printf("*");
    std::printf("\n");
    std::printf(" **** WARNING, THE PARITIES OF THE STATES MUST BE KNOWN"
                " TO UNIQUELY DETERMIN THE MULTIPOLARITY OF THE EXCITATION\n");
    std::printf(" **** PTOLEMY ASSUMES THAT THE CHANGE IN THE INTRINSIC PARITY"
                " IS JUST  (-1)**|J(IN)-J(OUT)|.\n");
    std::printf(" ");
    for (int _i = 0; _i < 120; _i++) std::printf("*");
    std::printf("\n");
    //      LXMAX = LXMAX - MOD( LXMIN+LXMAX, 2 )
    LXMAX = LXMAX - MOD(LXMIN + LXMAX, 2);
    goto stmnt_300;

stmnt_295:
    // 295  IF ( LXMAX .GT. LXMIN )  GO TO 300
    if (LXMAX > LXMIN) goto stmnt_300;
    //      write (6, 298) IPARIT, JIN, JOUT, LXMIN, LXMAX
    // 298  FORMAT ( '0**** ERROR:  PARITIES AND J-VALUES RESULT IN',
    //    1  ' IMPOSSIBLE EXCITATION:', I4, 2(I4,'/2'), 2I4 )
    std::printf("\n**** ERROR:  PARITIES AND J-VALUES RESULT IN"
                " IMPOSSIBLE EXCITATION:%4d%4d/2%4d/2%4d%4d\n",
                IPARIT, JIN, JOUT, LXMIN, LXMAX);
    IRET = 0;

    // =========================================================================
    //      DOES HE WANT ONLY ONE LX OR ALL OF THEM?
    // =========================================================================
stmnt_300:
    //  300 ONELSW = LX .NE. NOTDEF
    ONELSW = (LX != NOTDEF_INT);
    //      IF ( .NOT. ONELSW )  GO TO  340
    if (!ONELSW) goto stmnt_340;

    //      IF ( LXMIN .LE. LX  .AND.  LX .LE. LXMAX )  GO TO  330
    if (LXMIN <= LX  &&  LX <= LXMAX) goto stmnt_330;

    //      IF ( TRANSW )
    //    1  write (6,  303) JBP, JBT, LBP, LBT, LXMIN, LXMAX, LX
    // 303  FORMAT ( '0**** ERROR IN INPUT:' /
    //    1  ' JBP, JBT, LBP, LBT =', 2(I5, '/2'), 2I5 /
    //    2  ' AND RESULTS IN' /
    //    3  ' LXMIN, LXMAX =', 2I5 /
    //    4  ' BUT LX HAS BEEN INPUT AS', I5 / )
    if (TRANSW) {
        std::printf("\n**** ERROR IN INPUT:\n");
        std::printf(" JBP, JBT, LBP, LBT =%5d/2%5d/2%5d%5d\n",
                    JBP, JBT, LBP, LBT);
        std::printf(" AND RESULTS IN\n");
        std::printf(" LXMIN, LXMAX =%5d%5d\n", LXMIN, LXMAX);
        std::printf(" BUT LX HAS BEEN INPUT AS%5d\n\n", LX);
    }

    //      IF ( .NOT. TRANSW )  write (6, 307) JS(3), JS(1), JS(2), JS(4),
    //    1  IPARIT, LXMIN, LXMAX, LX
    // 307  FORMAT ( '0**** ERROR:  PARTICLE SPINS ARE', I5, '/2(',
    //    1    I3, '/2,', I3, '/2)', I3, '/2',
    //    1    4X, 'PARITY CHANGE =', I2 /
    //    1    ' **** WHICH RESULT IN',
    //    2  '  LXMIN, LXMAX =', 2I4, 5X, 'BUT LX =', I4,
    //    3    ' WAS SPECIFIED.' )
    if (!TRANSW) {
        std::printf("\n**** ERROR:  PARTICLE SPINS ARE%5d/2(%3d/2,%3d/2)%3d/2"
                    "    PARITY CHANGE =%2d\n",
                    JBIGA_INT, JA_INT, JB_INT, JBIGB_INT, IPARIT);
        std::printf(" **** WHICH RESULT IN  LXMIN, LXMAX =%4d%4d     BUT LX =%4d"
                    " WAS SPECIFIED.\n",
                    LXMIN, LXMAX, LX);
    }
    IRET = 0;

stmnt_330:
    //  330 LXMIN = LX
    //      LXMAX = LXMIN
    LXMIN = LX;
    LXMAX = LXMIN;

stmnt_340:
    // 340  NUMLX = LXMAX-LXMIN+1
    //      NUMLXI = (LXMAX-LXMIN)/2 + 1
    //      IF ( .NOT. TRANSW )  GO TO 400
    NUMLX = LXMAX - LXMIN + 1;
    NUMLXI = (LXMAX - LXMIN) / 2 + 1;
    if (!TRANSW) goto stmnt_400;

    // =========================================================================
    //      SUM OVER JX
    // =========================================================================
    //      JXMAX=2*LBP+JBP
    //      IF(JXMAX.GT.(2*LBT+JBT))JXMAX=2*LBT+JBT
    //      JXMIN=IABS(2*LBP-JBP)
    //      IF(JXMIN.LT.IABS(2*LBT-JBT))JXMIN=IABS(2*LBT-JBT)
    JXMAX = 2*LBP + JBP;
    if (JXMAX > 2*LBT + JBT) JXMAX = 2*LBT + JBT;
    JXMIN = std::abs(2*LBP - JBP);
    if (JXMIN < std::abs(2*LBT - JBT)) JXMIN = std::abs(2*LBT - JBT);

    //      IF ( JXMAX .GE. JXMIN )  GO TO  342
    if (JXMAX >= JXMIN) goto stmnt_342;
    //      write (6,  341) JBP, JBT, LBP, LBT, JXMIN, JXMAX
    //  341 FORMAT ( '0**** INVALID BOUND STATE ANGULAR MOMENTA:' /
    //    1  ' JP, JT, LP, LT =', 2(I5, '/2'), 2I5 /
    //    2  ' RESULTS IN  JXMIN, JXMAX =', 2(I7, '/2') )
    std::printf("\n**** INVALID BOUND STATE ANGULAR MOMENTA:\n");
    std::printf(" JP, JT, LP, LT =%5d/2%5d/2%5d%5d\n",
                JBP, JBT, LBP, LBT);
    std::printf(" RESULTS IN  JXMIN, JXMAX =%7d/2%7d/2\n", JXMIN, JXMAX);
    IRET = 0;

stmnt_342:
    //  342 IF ( JX .NE. NOTDEF )  GO TO  345
    if (JX_INT != NOTDEF_INT) goto stmnt_345;
    //      IF ( JXMIN .EQ. JXMAX )  GO TO  350
    if (JXMIN == JXMAX) goto stmnt_350;
    //      write (6,  343) JXMIN, JXMAX
    //  343 FORMAT ( '0***** ERROR: MORE THAN ONE JX IS POSSIBLE BUT JX WAS',
    //    1  ' NOT SET.  JXMIN, JXMAX =', 2(I5, '/2') / )
    std::printf("\n***** ERROR: MORE THAN ONE JX IS POSSIBLE BUT JX WAS"
                " NOT SET.  JXMIN, JXMAX =%5d/2%5d/2\n\n",
                JXMIN, JXMAX);
    IRET = 0;
    goto stmnt_400;

stmnt_345:
    //  345 IF ( JXMIN .LE. JX  .AND.  JX .LE. JXMAX )  GO TO  400
    if (JXMIN <= JX_INT  &&  JX_INT <= JXMAX) goto stmnt_400;
    //      write (6,  348) JXMIN, JXMAX, JX
    //  348 FORMAT ( '0**** ERROR, JXMIN, JXMAX, JX =', 3(I5, '/2') / )
    std::printf("\n**** ERROR, JXMIN, JXMAX, JX =%5d/2%5d/2%5d/2\n\n",
                JXMIN, JXMAX, JX_INT);
    IRET = 0;
    goto stmnt_400;

stmnt_350:
    //  350 JX = JXMIN
    JX_ref = (double)JXMIN;

    // =========================================================================
    //      COMPUTE L CRITICAL
    // =========================================================================
stmnt_400:
    // 400  LCRIT = (LCRITS(1)+LCRITS(2))/2
    LCRIT = (KANDM.LCRITS[1] + KANDM.LCRITS[2]) / 2;

    // =========================================================================
    //      IF AN ARRAY NAMED LVALUES IS DEFINED, WE USE IT
    // =========================================================================
    //      ILVALU = NAMLOC( 'LVALUES  ' )
    //      IF ( ILVALU .EQ. 0 )  GO TO 430
    ILVALU = NAMLOC("LVALUES  ");
    if (ILVALU == 0) goto stmnt_430;

    //      NUMLIS = LENG(ILVALU)
    //      ILIS = NALLOC( (NUMLIS+1)/FACFR4, 'LIS     ' )
    NUMLIS = LENGTH.LENG[ILVALU];
    ILIS   = NALLOC("LIS     ", (NUMLIS + 1) / ALLOCS.FACFR4);

    //      LLVALU = Z(ILVALU)-1
    //      LMIN = ALLOC(LLVALU+1)
    LLVALU = LOCPTRS.Z[ILVALU] - 1;
    LMIN = (int)ALLOC(LLVALU + 1);

    //      L1 = -1
    //      LSTEP = 1
    L1 = -1;
    LSTEP = 1;

    //      DO 419  I = 1, NUMLIS
    for (I = 1; I <= NUMLIS; I++) {
        //         L = ALLOC(LLVALU+1)
        L = (int)ALLOC(LLVALU + 1);
        //         IF ( L .GT. L1 )  GO TO 415
        if (L > L1) goto stmnt_415;
        //         write (6, 413) I, L, L1
        // 413     FORMAT ( '0*** LVALUES ARRAY MUST BE INCREASING:  ',
        //    1     'ELEMENT', I4, ' IS ', I5, ', PREVIOUS IS', I4 )
        std::printf("\n*** LVALUES ARRAY MUST BE INCREASING:  "
                    "ELEMENT%4d IS %5d, PREVIOUS IS%4d\n",
                    I, L, L1);
        return;
stmnt_415:
        //  415     ILLOC(FACFR4*Z(ILIS)-FACFR4+I) = L
        //          LSTEP = MAX0( LSTEP, L-L1 )
        //          L1 = L
        ILLOC(ALLOCS.FACFR4 * LOCPTRS.Z[ILIS] - ALLOCS.FACFR4 + I) = L;
        LSTEP = std::max(LSTEP, L - L1);
        L1 = L;
    }

    //      LMAX = L
    LMAX = L;
    goto stmnt_450;

    // =========================================================================
    //      NOW ESTIMATE LMIN, LMAX FROM L CRITICAL IF THEY WERE NOT GIVEN
    // =========================================================================
stmnt_430:
    // 430  LMININ = LMIN
    LMININ = LMIN;
    //      IF ( LMIN .NE. NOTDEF )  GO TO 440
    if (LMIN != NOTDEF_INT) goto stmnt_440;
    //      IF ( LMINSB .EQ. NOTDEF )  LMINSB = 10
    if (LMINSB == NOTDEF_INT) LMINSB = 10;
    //      LMIN = LCRIT*ALMNMT
    LMIN = (int)(LCRIT * ALMNMT);
    //      LMIN = MAX0( 0, MIN0( LMIN, LCRIT-LMINSB ) )
    LMIN = std::max(0, std::min(LMIN, LCRIT - LMINSB));

stmnt_440:
    // 440  IF ( LMAX .NE. NOTDEF )  GO TO 450
    if (LMAX != NOTDEF_INT) goto stmnt_450;
    //      LMAX = ALMXMT*LCRIT
    LMAX = (int)(ALMXMT * LCRIT);
    //      LMAX = MAX0( LMAX, LCRIT + MAX0(4*LSTEP, LMAXAD) )
    LMAX = std::max(LMAX, LCRIT + std::max(4*LSTEP, LMAXAD));

    // =========================================================================
    //      IF THERE ARE IDENTICAL PARTICLES IN EXIT OR ENTRANCE CHANNEL,
    //      MUST FIX UP CHOICE OF L'S NOW.
    // =========================================================================
stmnt_450:
    // 450  IDENSW = .FALSE.
    //      LSKIP = MAX0( LSKIPS(1), LSKIPS(2) )
    //      IF ( ISTATS(1) .EQ. 3 )  GO TO 460
    //      IF ( JA .EQ. 0 )  GO TO 455
    IDENSW = false;
    LSKIP = std::max(WAVCOM.LSKIPS[1], WAVCOM.LSKIPS[2]);
    if (WAVCOM.ISTATS[1] == 3) goto stmnt_460;
    if (JA_INT == 0) goto stmnt_455;

    // 452  write (6, 453)
    // 453  FORMAT ( '0**** ONLY SPIN-0 IDENTICAL PARTICLES MAY BE USED',
    //    1    ' FOR NOW.' )
    //      IRET = 0
    //      GO TO 500
stmnt_452:
    std::printf("\n**** ONLY SPIN-0 IDENTICAL PARTICLES MAY BE USED FOR NOW.\n");
    IRET = 0;
    goto stmnt_500;

stmnt_455:
    // 455  LIBASE = 0
    //      GO TO 470
    LIBASE = 0;
    goto stmnt_470;

stmnt_460:
    // 460  IF ( ISTATS(2) .EQ. 3 )  GO TO 480
    //      IF ( JB .GT. 0 )  GO TO 452
    if (WAVCOM.ISTATS[2] == 3) goto stmnt_480;
    if (JB_INT > 0) goto stmnt_452;
    //      LIBASE = MOD( LBP+LBT, 2 )
    LIBASE = MOD(LBP + LBT, 2);

stmnt_470:
    // 470  IDENSW = .TRUE.
    //      LSTEP = LSTEP + MOD(LSTEP, 2)
    IDENSW = true;
    LSTEP = LSTEP + MOD(LSTEP, 2);

    // =========================================================================
    //      MAKE LMIN AND LMAX LSTEP MULTIPLES
    // =========================================================================
stmnt_480:
    // 480  IF ( LMININ .EQ. NOTDEF )  LMIN = LSTEP*(LMIN/LSTEP)
    //      IF ( LSKIP .NE. 1 )  LMIN = IABS( LMIN-LIBASE )
    //      LMAX = LSTEP*((LMAX-LMIN+LSTEP-1)/LSTEP) + LMIN
    if (LMININ == NOTDEF_INT) LMIN = LSTEP * (LMIN / LSTEP);
    if (LSKIP != 1) LMIN = std::abs(LMIN - LIBASE);
    LMAX = LSTEP * ((LMAX - LMIN + LSTEP - 1) / LSTEP) + LMIN;

stmnt_500:;
    // =========================================================================
    // Part 2: Process spectroscopic or deformation factors (statement 500+)
    // =========================================================================
    if (TRANSW) goto stmnt_580;

    // Setup excitation potential parameters
    {
        EMUL = ELAB;
        if (IEPOW == 1) EMUL = 1.0 / ELAB;
        if (R0 != UNDEF) R = R0MASS * (R0 + EMUL*(R0E + EMUL*R0ESQ));
        if (RI0 != UNDEF) RI = R0MASS * (RI0 + EMUL*(RI0E + EMUL*RI0ESQ));
        if (RSI0 != UNDEF) RSI = R0MASS * RSI0;
        TAR = A_flt + EMUL*(AE + EMUL*AESQ);
        TAI = AI_flt + EMUL*(AIE + EMUL*AIESQ);
        TVR = V + EMUL*(VE + EMUL*VESQ);
        TVI = VI + EMUL*(VIE + EMUL*VIESQ);
        if (RC0 != UNDEF) RC = R0MASS * RC0;

        // Always use point and sphere Coulomb for effective excitation
        RCP = 0.0;
        RCT = RC;

        if (RI == UNDEF) RI = R;

        if (R != UNDEF && RC != UNDEF) goto stmnt_510;
        std::printf("\n**** BOTH R AND RC (OR R0 AND RC0) MUST BE"
                    " DEFINED TO COMPUTE DEFORMATION LENGTHS:\n"
                    "      R, RC = %8.3f %8.3f\n", R, RC);
        IRET = 0;
        return;

stmnt_510:
        R0 = R / R0MASS;
        RI0 = RI / R0MASS;
        RC0 = RC / R0MASS;
        VC = -3.0 * IZS[1] * IZS[3] * HBARC / AFINE;
        R2S[4] = VC;

        if (COUPSW) goto stmnt_600;

        // Deformation lengths use only target or projectile radius
        R2MASS = std::pow(AMS_ptr[2*IEXCIT - 1], 1.0/3.0);
        R2S[1] = R0 * R2MASS;
        R2S[2] = RI0 * R2MASS;
        R2S[3] = RC0 * R2MASS;

        // Setup IBETAS array
        BETA = 1.0;
        for (II = 1; II <= 4; II++) {
            IB = NAMLOC(BETNAM[II]);
            BETSWS[II] = (IB != 0);
            LEN = 0;
            if (IB == 0) goto do_stmnt_540;

            BETA = ALLOC(LOCPTRS.Z[IB]);
            LEN = LENGTH.LENG[IB];
            if (LEN == NUMLXI) goto do_stmnt_548;
            if (LEN <= 1) goto do_stmnt_540;

            std::printf("\n**** %.8s WAS DEFINED TO HAVE%4d"
                        " VALUES BUT THERE ARE%4d VALUES OF LX\n",
                        BETNAM[II], LEN, NUMLXI);
            IRET = 0;
            goto do_stmnt_548;

do_stmnt_540:
            IB = NALLOC(BETNAM[II], NUMLXI);
            for (I = 1; I <= NUMLXI; I++)
                ALLOC(LOCPTRS.Z[IB] - 1 + I) = BETA;

do_stmnt_548:
            INELCM.IBETAS[II] = IB;
        }

        CHARGE = (double)IZS[2*IEXCIT - 1];
        IBETNR = NALLOC("BETANRAT", NUMLXI);
        LBETA   = LOCPTRS.Z[INELCM.IBETAS[1]] - LXMIN/2;
        LBETAC  = LOCPTRS.Z[INELCM.IBETAS[2]] - LXMIN/2;
        LBETAR  = LOCPTRS.Z[INELCM.IBETAS[3]] - LXMIN/2;
        LBELX_loc = LOCPTRS.Z[INELCM.IBETAS[4]] - LXMIN/2;
        LBETNR  = LOCPTRS.Z[IBETNR] - LXMIN/2;

        // Initialize for Clebsch-Gordan
        DUMMY1();
        K = JIN % 2;

        // Make sure BETA(Coulomb) are defined
        if (BETSWS[2]) goto stmnt_560;
        if (BETSWS[4]) goto stmnt_555;

        std::printf("\n**** WARNING:  BETA COULOMB IS BEING CHOOSEN SUCH THAT:\n"
                    "                RC*BETACOULOMB = R*BETA"
                    "    WHERE RC AND R REFER ONLY TO THE EXCITED NUCLEUS\n");
        for (int LXv = LXMIN; LXv <= LXMAX; LXv += 2)
            ALLOC(LBETAC + LXv/2) = (R2S[1] / R2S[3]) * ALLOC(LBETA + LXv/2);
        goto stmnt_560;

stmnt_555:
        for (int LXv = LXMIN; LXv <= LXMAX; LXv += 2) {
            ALLOC(LBETAC + LXv/2) =
                (4.0*CNSTNT.PI / (3.0*CHARGE))
                * std::sqrt(ALLOC(LBELX_loc + LXv/2))
                / (std::pow(R2S[3]/10.0, (double)LXv)
                   * CLEBSH(JIN, 2*LXv, K, 0, JOUT, K));
        }

stmnt_560:
        if (BETSWS[1]) goto stmnt_570;
        if (BETSWS[2] || BETSWS[4]) goto stmnt_565;

        std::printf("\n");
        for (int _i = 0; _i < 120; _i++) std::printf("*");
        std::printf("\n **");
        for (int _i = 0; _i < 115; _i++) std::printf(" ");
        std::printf("**\n");
        std::printf(" **    WARNING:  BETA ASSUMED TO BE 1 --"
                    " MUST MULTIPLY CROSS SECTIONS BY SQUARE OF ACTUAL"
                    " BETA TO GET TRUE VALUES.");
        for (int _i = 0; _i < 19; _i++) std::printf(" ");
        std::printf("**\n ** ");
        for (int _i = 0; _i < 115; _i++) std::printf(" ");
        std::printf("**\n ");
        for (int _i = 0; _i < 120; _i++) std::printf("*");
        std::printf("\n");
        goto stmnt_570;

stmnt_565:
        for (int LXv = LXMIN; LXv <= LXMAX; LXv += 2)
            ALLOC(LBETA + LXv/2) = (R2S[3] / R2S[1]) * ALLOC(LBETAC + LXv/2);

stmnt_570:
        // Define auxiliary quantities
        for (int LXv = LXMIN; LXv <= LXMAX; LXv += 2) {
            TEMP = CLEBSH(JIN, 2*LXv, K, 0, JOUT, K);
            ALLOC(LBETAR + LXv/2) = ALLOC(LBETAC + LXv/2) * TEMP
                                  * std::pow(R2S[3], (double)LXv) / (2*LXv + 1.0);
            ALLOC(LBETNR + LXv/2) = ALLOC(LBETA + LXv/2) * TEMP;
            ALLOC(LBELX_loc + LXv/2) = std::pow(
                3.0 * CHARGE * ALLOC(LBETAC + LXv/2) * TEMP
                * std::pow(R2S[3]/10.0, (double)LXv) / (4.0*CNSTNT.PI),
                2.0);
        }

        // Final IBETAS definition: replace BETACOULOMB with BETANRAT
        INELCM.IBETAS[2] = IBETNR;
        goto stmnt_600;
    } // end inelastic block

    // =========================================================================
    // Statement 580: Process spectroscopic info for transfer
    // =========================================================================
stmnt_580:
    if (FLOAT_common.SPAMP == 1.0) FLOAT_common.SPAMP = std::sqrt(FLOAT_common.SPFACP);
    FLOAT_common.SPFACP = FLOAT_common.SPAMP * FLOAT_common.SPAMP;
    if (FLOAT_common.SPAMT == 1.0) FLOAT_common.SPAMT = std::sqrt(FLOAT_common.SPFACT);
    FLOAT_common.SPFACT = FLOAT_common.SPAMT * FLOAT_common.SPAMT;

    // Allocate and fill A-term array
    INELCM.IBETAS[1] = NALLOC("ATERM   ", LXMAX + 1);

    // Initialize for Racahs
    DUMMY1();
    LATERM = LOCPTRS.Z[INELCM.IBETAS[1]];
    TEMP = (JBIGB_ref + 1.0) / (JBIGA_ref + 1.0);
    if (ISTRIP == -1) TEMP = (JB_ref + 1.0) / (JA_ref + 1.0);
    TEMP = std::sqrt(TEMP);

    for (int LXv = LXMIN; LXv <= LXMAX; LXv++) {
        double racah_val = RACAH(2*LBT, JBT, 2*LBP, JBP, JX_INT, 2*LXv);
        ATERM = TEMP * std::sqrt(2*LXv + 1.0) * FLOAT_common.SPAMP * FLOAT_common.SPAMT
              * racah_val;
        ITEST = JX_INT - JBP + 2*(LBP + LBT);
        if (ISTRIP == -1)
            ITEST = JX_INT - JBT + JA_INT + JBIGA_INT - JB_INT - JBIGB_INT;
        ITEST = ITEST/2 + 1;
        if (ITEST % 2 != 0) ATERM = -ATERM;
        ALLOC(LATERM + LXv) = ATERM;
        if (IPRNT >= 4)
            std::printf(" LX =%3d     SIGN, ATERM =%4d%15.8G\n",
                        LXv, ITEST, ATERM);
    }

    // =========================================================================
    // Statement 600: Produce a nice summary of the reaction
    // =========================================================================
stmnt_600:
    if (!PRNTSW) goto stmnt_800;

    {
        NUMPAR = 4 + std::abs(ISTRIP);

        std::printf("0%20sSUMMARY OF THE REACTION\n", "");
        std::printf("0%23s", "");
        for (int ch = 1; ch <= 45; ch++) std::printf("%c", HEDCOM.REACT[ch]);
        std::printf("\n");
        std::printf("0%17sBIGA   (   A   ,    B   )  BIGB   %s\n",
                    "", XWORDS[NUMPAR - 3]);
        std::printf("0M (AMU)      ");
        for (int J = 1; J <= NUMPAR; J++)
            std::printf("%9.2f", AMS_ptr[IND[J]]);
        std::printf("\n");

        std::printf(" Z         ");
        for (int J = 1; J <= NUMPAR; J++)
            std::printf("%9d", IZS[IND[J]]);
        std::printf("\n");

        std::printf(" E* (MEV)       ");
        for (int J = 1; J <= NUMPAR; J++)
            std::printf("%9.4f", FLOAT_common.EXS[IND[J]]);
        std::printf("\n");

        std::printf(" J           ");
        for (int J = 1; J <= NUMPAR; J++)
            std::printf("%7d/2", (int)JBLOCK.JS[IND[J]]);
        std::printf("\n");

        std::printf(" PARITY        ");
        for (int J = 1; J <= NUMPAR; J++) {
            std::printf("%7.7s", PWORD[PARITS[IND[J]] + 2]);
            if (J < NUMPAR) std::printf("  ");
        }
        std::printf("\n");
    }

    if (!TRANSW) goto stmnt_680;

    // Bound state properties (transfer only)
    {
        std::printf("\n0     BOUND STATE PROPERTIES\n");
        std::printf("0           PROJECTILE   TARGET\n");
        std::printf(" E         %10.4f%10.4f          Q =%9.4f\n",
                    INTRNL.EBNDS[1], INTRNL.EBNDS[2], FLOAT_common.Q);
        std::printf(" JP      %8d/2%8d/2\n", JBP, JBT);
        std::printf(" L      %10d%10d\n", LBP, LBT);
        std::printf(" NODES  %10d%10d\n", NBNDS[0], NBNDS[1]);
        std::printf(" SPEC. AMP.%10.4f%10.4f\n",
                    FLOAT_common.SPAMP, FLOAT_common.SPAMT);
        std::printf(" SPEC. FACTOR%8.4f%10.4f\n",
                    FLOAT_common.SPFACP, FLOAT_common.SPFACT);
    }

    if (ISTRIP ==  1) std::printf("0THIS IS A STRIPPING REACTION\n");
    if (ISTRIP == -1) std::printf("0THIS IS A PICKUP REACTION\n");

    std::printf("0THE EXCHANGE POTENTIAL CONSISTS OF\n");
    if (SWITCH.IVRTEX == 1)
        std::printf("    THE NUCLEAR POTENTIAL AT THE PROJECTILE VERTEX\n");
    if (SWITCH.IVRTEX == 2)
        std::printf("    THE NUCLEAR POTENTIAL AT THE TARGET VERTEX\n");
    if (NUCONL == 0)
        std::printf("    AND THE COULOMB POTENTIAL AT THE SAME VERTEX.\n");
    if (NUCONL >= 2)
        std::printf("    AND THE COULOMB POTENTIAL INCLUDING CORE-CORE"
                    " CORRECTION TERMS\n");
    if (NUCONL == 3)
        std::printf("    AND THE CORE-CORE CORRECTIONS FROM THE REAL"
                    " PART OF THE NUCLEAR OPTICAL POTENTIAL\n");
    goto stmnt_750;

    // =========================================================================
    // Statement 680: Printout for inelastic scattering
    // =========================================================================
stmnt_680:
    if (COUPSW) goto stmnt_750;

    std::printf("\nTHE REACTION IS INELASTIC EXCITATION OF THE %s%s\n",
                WORDS[2*IEXCIT-1], WORDS[2*IEXCIT]);
    std::printf("\nLX         DEFORMATION PARAMETERS%21sB(E(LX))\n", "");
    std::printf("            NUCLEAR         COULOMB        (E**2 BARN**LX)\n");
    for (int LXv = LXMIN; LXv <= LXMAX; LXv += 2) {
        std::printf("%3d     %15.5G%15.5G%15.5G\n",
                    LXv,
                    ALLOC(LBETA + LXv/2),
                    ALLOC(LBETAC + LXv/2),
                    ALLOC(LBELX_loc + LXv/2));
    }
    std::printf("\nDEFORMATION RADII:\n");
    std::printf("\n  REAL%14.4f\n",     R2S[1]);
    std::printf("   IMAGINARY%9.4f\n", R2S[2]);
    std::printf("   COULOMB%11.4f\n",  R2S[3]);

    // =========================================================================
    // Statement 750: Critical-L summary and L-range printout
    // =========================================================================
stmnt_750:
    std::printf("\n0ESTIMATED CRITICAL L'S:    INCOMING =%4d"
                "     OUTGOING =%4d     AVERAGE =%4d\n",
                KANDM.LCRITS[1], KANDM.LCRITS[2], LCRIT);

    if (COUPSW) goto stmnt_760;

    if (ONELSW)
        std::printf("0     LX =%3d\n", LXMIN);
    if (!ONELSW)
        std::printf("0%5d  =< LX =<%3d\n", LXMIN, LXMAX);
    std::printf("+%23s%5d  =<  LIN, LOUT  =<%4d\n", "", LMIN, LMAX);
    goto stmnt_765;

stmnt_760:
    std::printf("0FOR COMPUTED S MATRICES:%6d =< J OR LIN =<%5d\n",
                LMIN, LMAX);

stmnt_765:
    if (ILVALU == 0)
        std::printf("+%58sLSTEP =%3d\n", "", LSTEP);
    if (ILVALU != 0)
        std::printf("+%58sUSING SPECIFIED VALUES WITH A MAXIMUM LSTEP OF%4d\n",
                    "", LSTEP);

    // =========================================================================
    // Statement 800: Q-value consistency warnings
    // =========================================================================
stmnt_800:
    if (!TRANSW) goto stmnt_810;
    if (std::abs(INTRNL.EBNDS[1] - INTRNL.EBNDS[2] - (double)ISTRIP * FLOAT_common.Q) > 0.001)
        std::printf("\n**** WARNING:  Q IS NOT EQUAL TO THE DIFFERENCE"
                    " OF THE BOUND STATE ENERGIES.\n");

stmnt_810:
    if (std::abs(FLOAT_common.ECM + FLOAT_common.Q - WAVCOM.ES[2]) > 0.001)
        std::printf("\n**** WARNING:  THE DIFFERENCE OF THE C.M. SCATTERING"
                    " ENERGIES IS NOT EQUAL TO Q.\n");

    KANDM.LOMOST = LMAX + LXMAX;
    INTRNL.LIMOST = LMAX;

    if (LCRIT == 0) LCRIT = LMAX / 2;

    return;

    #undef JA_INT
    #undef JB_INT
    #undef JBIGA_INT
    #undef JBIGB_INT
    #undef JBP_INT
    #undef JBT_INT
    #undef JX_INT
} // end PRBPRT
