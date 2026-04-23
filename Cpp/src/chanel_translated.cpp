// chanel_translated.cpp — CHANEL
// Translated from source.f lines 7953-8230
// Verbatim transliteration from Fortran to C++

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include "ptolemy_alloc.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_forward.h"
#include <cstdio>
#include <cstring>
#include <cmath>

// ============================================================================
// SUBROUTINE CHANEL ( IRET )
//
//     SCANS A CHANNEL SPECIFICATION
//
//     FORMS ARE
//        N15 + 209BI(9/2- 1.3)     ( FOR SCATTERING )
//        P + 208PB = BI209(7/2 .9)    ( FOR BOUND STATES )
//        O16 = N15 + P                ( ALTERNATIVE B.S. )
//     NOTE THAT PROJECTILE ALWAYS COMES FIRST.
//
//
//     AND STORES  Z'S,  M'S,  J'S,  E*'S   OF
//     THE 2 OR 3 PARTILCES INTO THE COMMON BLOCKS.
//     FOR BOUNDSTATES THE TOTAL J AND BINDING ENERGY ARE ALSO STORED.
//
//     9/22/75 - BASED ON REACTN - S. PIEPER
//     4/24/77 - STORE E* IN EXSPT
//     9/24/77 - RECOGNIZE E*=0; ALLOW LATTER DEF OF MASS EXCESS
//     12/10/77 - FIX ERROR IN ABOVE TO DEFINE J FOR B.S.
//     5/6/79 - G.S. MASS EXCESS, PARITIES, L & NODES FROM LEVEL - S.P.
//     12/31/79 - STANDARD MORTRAN - RPG
//     1/14/80 - CHARDATAVALUE IS NUMERIC, " IN FORMAT 313 - RPG
//     1/25/80 - CHARACTER SYMBOL, ALPHA - RPG
//
// ============================================================================

void CHANEL(int& IRET)
{
    // EQUIVALENCE arrays — pointers into COMMON blocks
    // AMPTS(2) <=> AMP, AMT  (contiguous doubles in FLOAT_common)
    double* AMPTS = &FLOAT_common.AMP;   // AMPTS[1]=AMP, AMPTS[2]=AMT (1-based)

    // IZPTS(2) <=> IZP, IZT  (contiguous ints in INTGER)
    int* IZPTS = &INTGER.IZP;            // IZPTS[1]=IZP, IZPTS[2]=IZT (1-based)

    // JSPTS(2) <=> JSP, JST  (contiguous doubles in JBLOCK)
    // Fortran declares INTEGER JSPTS(2) equivalenced to JSP (REAL*8).
    // This is a type-pun; in C++ we access JSP/JST as doubles directly.
    double* JSPTS = &JBLOCK.JSP;         // JSPTS[1]=JSP, JSPTS[2]=JST (1-based)

    // Local variables
    char8  GUY[4];         // 1-based: GUY[1]..GUY[3]
    char   SYMBOL[3];      // CHARACTER*2
    char   ALPHA[5];       // CHARACTER*4
    double ESTARS[4];      // 1-based: ESTARS[1]..ESTARS[3]
    double XCESES[4];      // 1-based: XCESES[1]..XCESES[3]
    int    NODVLS[4];      // 1-based
    int    LVALS[4];       // 1-based
    int    JVALS[4];       // 1-based
    int    IPARYS[4];      // 1-based
    char   BLANK1 = ' ';
    char   mychar;

    char   NAMES[4][17];   // 1-based: NAMES[1]..NAMES[3], each CHARACTER*16
    std::strcpy(NAMES[1], "PROJECTILE");
    std::strcpy(NAMES[2], "TARGET");
    std::strcpy(NAMES[3], "BOUND STATE");

    int    IZ = 0, IA = 0, IRTN, NUMSYM, I, N;
    int    INCHST;
    int    IAB = 0, IZB = 0;
    double AMXCES;
    int    JJ;
    double ESTAR;
    int    JVAL, IPARIT;

    //
    IRET = 0;
    //
    //     GET PAST ANY INITIAL DELIMITORS
    //
L100:
    mychar = INPBUF.IBUF[INPBUF.INCH];
    if ((mychar >= 'A' && mychar <= 'Z') ||
        (mychar >= 'a' && mychar <= 'z') ||
        (mychar >= '0' && mychar <= '9'))  goto L150;
    INPBUF.INCH = INPBUF.INCH + 1;
    if (INPBUF.INCH <= INPBUF.NOCH)  goto L100;
    std::printf("\n **** CHANNEL KEYWORD AND COMPLETE SPECIFICATION"
                " MUST BE ON ONE LINE.\n");
    return;
    //
    //     START OF SPECIFICATION FOUND.  BREAK IT INTO 2 OR 3 PIECES.
    //
L150:
    INCHST = INPBUF.INCH;
    CHNSCN(GUY, ESTARS, NODVLS, LVALS, JVALS, IPARYS, IRTN);
    if (IRTN < 0)  return;
    //
    //     NOW GET THE A AND Z VALUES
    //
    NUMSYM = IRTN;
    for (I = 1; I <= NUMSYM; I++) {
        AZCODE(GUY[I].data, IZ, IA, SYMBOL, IRTN);
        if (IRTN == 0)  goto L230;
        if (IRTN != -2)  std::printf("\n **** A SYMBOL HAS INCORRECT SYNTAX:  "
                                     "SYMBOL = %.8s\n", GUY[I].data);
        if (IRTN == -2)  std::printf("\n **** THE SYMBOL %.6s IS NOT A KNOWN ELEMENT.\n",
                                     GUY[I].data);
        return;
L230:
        if (I == 3)  goto L250;
        IZPTS[I-1] = IZ;
        AMPTS[I-1] = IA;
        INTRNL.NODEPT[I] = NODVLS[I];
        INTRNL.LSPCPT[I] = LVALS[I];
        JSPTS[I-1] = JVALS[I];
        FLOAT_common.EXSPT[I] = ESTARS[I];
        if (FLOAT_common.EXSPT[I] == INTRNL.UNDEF)  FLOAT_common.EXSPT[I] = 0;
        INTGER.PARIPT[I] = IPARYS[I];
    } // 239
    goto L300;
    //
L250:
    JBLOCK.J = JVALS[3];
    INTGER.PARITY = IPARYS[3];
    INTGER.NODES = NODVLS[3];
    INTGER.L = LVALS[3];
    IAB = IA;
    IZB = IZ;
    if (IAB == (int)(AMPTS[0] + AMPTS[1]) && IZ == IZPTS[0] + IZPTS[1])  goto L300;
    std::printf("\n **** BOUNDSTATE DOES NOT CONSERVE NUCLEON NUMBER"
                " OR CHARGE:\n"
                " ZP, ZT, Z(BOUND) =%8d%8d%8d\n"
                " MP, MT, M(BOUND) =%8.0f%8.0f%8d\n",
                IZPTS[0], IZPTS[1], IZ,
                AMPTS[0], AMPTS[1], IA);
    return;
    //
    //     NOW GET GROUND STATE MASS EXCESS AND DETERMIN IF THE NUCLEII
    //     EXIST.
    //
L300:
    for (I = 1; I <= NUMSYM; I++) {
        if (I == 3) {
            IA = IAB;
            IZ = IZB;
        } else {
            IA = (int)AMPTS[I-1];
            IZ = IZPTS[I-1];
        }
        N = IA - IZ;
        AMXCES = EXCESS(IZ, IA, IRTN);
        if (IRTN == 0)  goto L330;
        std::printf("\n **** THE NUCLEUS WITH A =%4d,  Z =%4d"
                    ",  N =%4d  IS NOT BOUND "
                    "ACCORDING TO PTOLEMY\"S MASS TABLE.\n", IA, IZ, N);
        goto L359;
L330:
        XCESES[I] = AMXCES;
        if (I <= 2)  FLOAT_common.AMXGPT[I] = AMXCES;
        if (I <= 2)  XCESES[I] = XCESES[I] + FLOAT_common.EXSPT[I];
        if (I == 3 && ESTARS[3] != INTRNL.UNDEF)
            XCESES[3] = XCESES[3] + ESTARS[3];
L359:
        ;
    } // 359
    //
    //     NOW GET J FOR EACH OF THE GROUND STATES
    //
    for (I = 1; I <= NUMSYM; I++) {
        //
        if (I == 3) {
            IA = IAB;
            IZ = IZB;
            JJ = (int)JBLOCK.J;
        } else {
            IA = (int)AMPTS[I-1];
            IZ = IZPTS[I-1];
            JJ = (int)JSPTS[I-1];
        }
        ESTAR = ESTARS[I];
        //
        //     IF E* HAS BEEN ENTERED THEN G.S. J IS MEANINGLESS
        //
        if (ESTAR > 1.0e-5 && ESTAR != INTRNL.UNDEF)  goto L459;
        N = IA - IZ;
        GSINFO(IZ, IA, JVAL, IPARIT, ALPHA, IRTN);
        if (IRTN == 0 && JVAL != NOTDEF_INT)  goto L430;
        //     NO ERROR IF J IS ALREADY DEFINED OR IS EXPLICITLY GIVEN
        //     IN THE REACTION.
        if (JJ != NOTDEF_INT)  goto L459;
        //
        //     DINEUTRON OR DIPROTON IS J = 0
        //
        if (IA != 2 || IZ == 1)  goto L420;
        JVAL = 0;
        goto L450;
L420:
        std::printf("\n **** WARNING:  THE GROUND STATE SPIN OF THE"
                    " NUCLEUS WITH A =%4d,  Z =%4d,  N =%4d"
                    "  IS NOT KNOWN TO PTOLEMY.\n", IA, IZ, N);
        goto L459;
        //
        //     IF J HAS BEEN DEFINED TO BE DIFFERENT FROM G.S. THEN NEED E*
        //
L430:
        if (JJ == NOTDEF_INT)  goto L450;
        if (JJ == JVAL || ESTAR != INTRNL.UNDEF)  goto L459;
        std::printf("\n **** WARNING:  GROUND STATE SPIN FOR %-10s"
                    " IS%3d/2.  YOU HAVE SPECIFIED J"
                    " =%3d/2, BUT HAVE NOT SPECIFIED E*\n",
                    NAMES[I], JVAL, JJ);
        goto L459;
        //
L450:
        if (I != 3)  goto L455;
        JBLOCK.J = JVAL;
        if (IPARIT != 0)  INTGER.PARITY = IPARIT;
        goto L459;
L455:
        JSPTS[I-1] = JVAL;
        if (IPARIT != 0)  INTGER.PARIPT[I] = IPARIT;
L459:
        ;
    } // 459
    //
    //     DEFINE BINDING ENERGY
    //
    if (NUMSYM == 3)
        FLOAT_common.E = XCESES[3] - XCESES[1] - XCESES[2];
    //
    //     IT WAS SUCCESFUL, COPY REACTION FOR HEADER
    //
    for (I = 1; I <= 45; I++) {
        HEDCOM.REACT[I] = BLANK1;
    } // 649
    N = std::min(45, INPBUF.INCH - INCHST);
    for (I = 1; I <= N; I++) {
        HEDCOM.REACT[I] = INPBUF.IBUF[INCHST - 1 + I];
    } // 659
    IRET = 1;
    return;
}
