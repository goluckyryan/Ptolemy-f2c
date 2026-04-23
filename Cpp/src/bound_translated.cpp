// bound_part1.cpp — BOUND subroutine, first half (source.f L4082-4580)
//
// Declarations, input validation, potential setup (Coulomb, central,
// spin-orbit), and linkule initialization.
//
// Translated verbatim from Fortran to C++.

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include "ptolemy_alloc.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_forward.h"
#include "ptolemy_io.h"
#include <cstdio>
#include <cstring>
#include <cmath>

void BOUND(int& IRET, int IVOUT, int IWF2)
{
//
//     BOUND STATE PROGRAM
//
//  CODED 2/22/73  D.GLOECKNER
//     1/20/75 - FIX NSTEP CALCULATION - S. PIEPER
//     ... (extensive history, see source.f)
//
//     NOTE...  THE CHANNEL NUMBER (1 OR 2) IS INPUT VIA  ICHANB BELOW
//
//     IFIT = 0  MEANS FIT E TO V.
//            1  MEANS FIT V TO E.
//
//     ON OUTPUT:
//     IRET IS 1 IF BOUND WAS SUCCESSFUL  OR
//          IS 0 IF BOUND WAS NOT SUCCESSFUL
//     IVOUT POINTS TO THE COMPLETE POTENTIAL
//     IWF2 POINTS TO THE WAVEFUNCTION
//

    // implicit real*8 (a-h, o-z)

    // --- Part 2 COMMON aliases (moved here to avoid goto-crosses-init) ---
    auto& ACCURA = FLOAT_common.ACCURA;
    auto& WRITES = FLOAT_common.WRITES;
    auto& Q      = FLOAT_common.Q;
    auto& IASYMP = SWITCH.IASYMP;
    auto& LINEAR = SWITCH.LINEAR;
    auto& MAXITR = INTGER.MAXITR;
    auto& R0MASS = INTRNL.R0MASS;
    auto* EBNDS  = INTRNL.EBNDS;
    auto& ISTRIP = INTRNL.ISTRIP;
    auto& J_val  = JBLOCK.J;
    auto& JST    = JBLOCK.JST;
    double RGET = 0, TEMP = 0, X = 0, XX = 0;
    double PPHI2 = 0, DP = 0, DPP = 0, XL = 0, XQ_loc = 0, XP_loc = 0, PHIP = 0;
    double ALL = 0;
    int NGO = 0, K = 0, IVK = 0, NMATCH = 0, NWRITE = 0;

    // --- Part 1 COMMON block aliases ---
    auto& PI      = CNSTNT.PI;
    auto& HBARC   = CNSTNT.HBARC;
    auto& AFINE   = CNSTNT.AFINE;
    auto& BIGNUM  = CNSTNT.BIGNUM;

    auto& A       = FLOAT_common.A;
    auto& ASO     = FLOAT_common.ASO;
    auto& ASYMPT  = FLOAT_common.ASYMPT;
    auto& DELTVK  = FLOAT_common.DELTVK;
    auto& E       = FLOAT_common.E;
    auto& TAU     = FLOAT_common.TAU;
    auto& AM      = FLOAT_common.AM;
    auto& AMP     = FLOAT_common.AMP;
    auto& AMT     = FLOAT_common.AMT;
    auto& R       = FLOAT_common.R;
    auto& R0      = FLOAT_common.R0;
    auto& RC      = FLOAT_common.RC;
    auto& RSO     = FLOAT_common.RSO;
    auto& STEPSZ  = FLOAT_common.STEPSZ;
    auto& V       = FLOAT_common.V;
    auto& VSO     = FLOAT_common.VSO;
    auto& STPSPR  = FLOAT_common.STPSPR;
    auto& POWRL   = FLOAT_common.POWRL;

    auto& L       = INTGER.L;
    auto& NODES   = INTGER.NODES;
    auto& IPRINT  = INTGER.IPRINT;
    auto& IZP     = INTGER.IZP;
    auto& IZT     = INTGER.IZT;
    auto& PARITY  = INTGER.PARITY;
    auto* PARIPT  = INTGER.PARIPT;  // 1-based

    auto& UNDEF   = INTRNL.UNDEF;
    int   NOTDEF  = NOTDEF_INT;
    auto& ICHANB  = INTRNL.ICHANB;
    auto& IHSAVE  = INTRNL.IHSAVE;
    auto& RATMAS  = INTRNL.RATMAS;

    auto& IFIT    = SWITCH.IFIT;
    auto& NUCONL  = SWITCH.NUCONL;

    auto& JP      = JBLOCK.JP;
    auto& JSP     = JBLOCK.JSP;

    auto& IBDS    = FORMF.IBDS;    // 1-based

    auto& IUNIQU  = LNKBLK.IUNIQU;
    // LNKADR accessed as LNKBLK.LNKADR[col][row] (1-based)

    auto* Z       = LOCPTRS.Z;     // 1-based

    // --- Local variables ---
    double T1, UK, STPSZ2, AFAC, DL2;
    double CI1, CI2, C0, VC, RVAL, VMUL;
    double ALS, BLS, VSOSV;
    double RCSAVE;
    double DELTA, DELTAV, DELTAK;
    double VST, UKST;
    double UKM, UKP, ETA;
    int NSTEPS, NMTOP, NMBOT;
    int IV1L, IV2L, IWF1, IWF3;
    int LV1, LV2, LVOUT, LWF1, LWF2, LWF3_loc;
    int ISEAR, ITER, ITERN;
    int NODEH;
    int LL, II, I;
    int N1, N2;
    int IRETUR;
    double DUMMY;
    int JUMP;  // for ASSIGN/GOTO dispatch

    // EQUIVALENCES from /FORMF/:
    //   DIMENSION LBNDS(2), NBNDS(2), JBNDS(2)
    //   EQUIVALENCE (LBNDS(1), LBP), (NBNDS(1), NBP), (JBNDS(1), JBP)
    int* LBNDS = &FORMF.LBP;  // LBNDS(1)=LBP, LBNDS(2)=LBT
    int* NBNDS = &FORMF.NBP;  // NBNDS(1)=NBP, NBNDS(2)=NBT
    int* JBNDS = &FORMF.JBP;  // JBNDS(1)=JBP, JBNDS(2)=JBT

    // CHARACTER*8  NAMES(2)  / 'PHIX', 'VBOUNDX' /
    char NAMES_arr[2][9];
    std::strncpy(NAMES_arr[0], "PHIX    ", 9);
    std::strncpy(NAMES_arr[1], "VBOUNDX ", 9);
    // NAMES1 overlays NAMES for character-level access
    // NAMES1(4,1) and NAMES1(7,2) will be set below

    // CHARACTER*8 CHANMS(2,3) / 'PROJECTI','LE','    TARG','ET',' ',' ' /
    // (used later, declare for completeness)
    char CHANMS[3][2][9];
    std::strncpy(CHANMS[0][0], "PROJECTI", 9);
    std::strncpy(CHANMS[0][1], "LE      ", 9);
    std::strncpy(CHANMS[1][0], "    TARG", 9);
    std::strncpy(CHANMS[1][1], "ET      ", 9);
    std::strncpy(CHANMS[2][0], "        ", 9);
    std::strncpy(CHANMS[2][1], "        ", 9);

    // DIMENSION UKS(3), ETAS(3), LWFS(3), UK2S(3), VS(3), DRIN(3),
    //   DROUT(3), PHIS(3), VALOUT(3), SUMOUT(3), SUMIN(3), NODCNT(3)
    double UKS[4], ETAS_loc[4], UK2S[4], VS_loc[4], DRIN[4];
    double DROUT[4], PHIS[4], VALOUT[4], SUMOUT[4], SUMIN[4];
    int LWFS[4], NODCNT[4];
    // EQUIVALENCE for LWFS and PHIS
    // LWF1 = LWFS(1), LWF2 = LWFS(2), LWF3 = LWFS(3)
    // PHI1 = PHIS(1), PHI2 = PHIS(2), PHI3 = PHIS(3)
    // (accessed via arrays below)

    int IPTYPS[5]; // 1-based
    IPTYPS[1] = 1; IPTYPS[2] = 3; IPTYPS[3] = 5; IPTYPS[4] = 6;

    // CHARACTER*8  PWORD(3) / '-1', 'UNKNOWN', '+1' /
    char PWORD[3][9];
    std::strncpy(PWORD[0], "-1      ", 9);
    std::strncpy(PWORD[1], "UNKNOWN ", 9);
    std::strncpy(PWORD[2], "+1      ", 9);

    // CHARACTER*1 CHNUM(9)
    char CHNUM[10]; // 1-based
    CHNUM[1] = '1'; CHNUM[2] = '2'; CHNUM[3] = '3'; CHNUM[4] = '4';
    CHNUM[5] = '5'; CHNUM[6] = '6'; CHNUM[7] = '7'; CHNUM[8] = '8';
    CHNUM[9] = '9';

    // CHARACTER*1 LINKID(4) / '*', '1', '0', '0' /
    char LINKID[5];
    LINKID[0] = '*'; LINKID[1] = '1'; LINKID[2] = '0'; LINKID[3] = '0';
    LINKID[4] = '\0';
    // character*8  linkfmtout — equivalenced to linkid(1)
    char linkfmtout[9];

    // LOGICAL  PRNTSW, PCONSW, SOSW
    bool PRNTSW, PCONSW, SOSW;

    char mychar;

//
//  IFIT=1  MEANS FIT V0
//  IFIT=0  MEANS FIT E
//
    ISEAR = 0;
//
L10:
    T1 = (float)dtime_();
    PCONSW = (IPRINT / 100) % 10 >= 2;
    PRNTSW = IPRINT % 10 >= 1;
//
//   INITIALIZE PARAMETERS
//
    if (ICHANB > 9) ICHANB = 1;
    mychar = CHNUM[ICHANB];
    NAMES_arr[0][3] = mychar;   // NAMES1(4,1) = mychar
    NAMES_arr[1][6] = mychar;   // NAMES1(7,2) = mychar
//
    SOSW = (VSO != 0.0 || TAU != 0.0);
//
//     60 IS A REASONABLE GUESS FOR V UNDER MANY CIRCUMSTANCES
//
    if (V == 0.0) V = 60.0;
//
//     TEST THE INPUT AND RETURN  0  COMPLETION CODE IF BAD
//
    IRET = 1;
    if (R > 0.0 && A > 0.0) goto L30;
    if (LNKBLK.LNKADR[1][3] > 0 || LNKBLK.LNKADR[6][3] > 0) goto L30;
    std::printf("\n **** R OR A HAS INVALID VALUE: %15.5G%15.5G\n", R, A);
    IRET = 0;
L30:
    if (!SOSW) goto L40;
    if (RSO > 0.0 && ASO > 0.0) goto L32;
    if (LNKBLK.LNKADR[3][3] > 0) goto L32;
    std::printf("\n **** RSO OR ASO HAS INVALID VALUE: %15.5G%15.5G\n", RSO, ASO);
    IRET = 0;
L32:
    if (JP != INTRNL.NOTDEF) goto L34;
    std::printf("\n **** JP MUST BE DEFINED FOR SPIN ORBIT FORCE.\n");
    IRET = 0;
L34:
    if (JSP != INTRNL.NOTDEF) goto L36;
    std::printf("\n **** WARNING:  SP WAS NOT DEFINED; "
                "IT IS ASSUMED TO BE 1/2 FOR THE SPIN-ORBIT FORCE.\n");
    JSP = 1;
L36:
    if ((int)JP <= (int)JSP + 2*L && (int)JP >= IABS((int)JSP - 2*L)
        && ((int)JP + (int)JSP) % 2 == 0) goto L40;
    std::printf("\n **** ERROR: L, SP AND JP DO NOT FORM A TRIANGLE OR "
                "THEY ARE MIXED INTEGER AND HALF-INTEGER: %5d%5d/2%5d/2\n",
                L, (int)JSP, (int)JP);
    IRET = 0;
L40:
    if (IFIT == 1 && V == 0.0) V = 50.0;
    if (IFIT == 0 && E == UNDEF) E = -15.0;
    if (E < 0.0 && V > 0.0 && AM > 0.0 && E != UNDEF) goto L50;
    std::printf("\n **** E OR V OR M (REDUCED MASS) HAS INVALID"
                " VALUE: %15.5G%15.5G%15.5G\n", E, V, AM);
    IRET = 0;
L50:
    if (R < ASYMPT || LNKBLK.LNKADR[1][3] > 0) goto L60;
    std::printf("\n **** R MUST BE LESS THAN ASYMTOPIA: %15.5G%15.5G\n", R, ASYMPT);
    IRET = 0;
L60:
    RCSAVE = RC;
    if (RC == UNDEF && IZP * IZT == 0) RC = 1.0;
    if (RC > 0.0 || LNKBLK.LNKADR[5][3] > 0) goto L70;
    std::printf("\n **** RC IS INVALID: %15.5G\n", RC);
    IRET = 0;
L70:
    if (LNKBLK.LNKADR[6][3] > 0) goto L75;
    if (L != NOTDEF && NODES != NOTDEF) goto L75;
    std::printf("\n **** BOTH L AND NODES MUST BE DEFINED.\n");
    IRET = 0;
L75:
    I = PARITY * PARIPT[1] * PARIPT[2];
    if (I == 0) goto L80;
    I = (I + 3) / 2 + L;
    if (I % 2 == 0) goto L80;
    std::printf("\n **** L AND PARITIES ARE INCOMPATABLE: %5d%5d%5d%5d\n",
                L, PARITY, PARIPT[1], PARIPT[2]);
    IRET = 0;
//
L80:
    ITER = 0;
    ITERN = 0;
//
    if (IRET == 0) goto L8875;
//
//  K IS
    UK = DSQRT(-2.0 * AM * E) / HBARC;
    if (STPSPR != UNDEF)
        STEPSZ = DMIN1(1.0 / UK, A) / STPSPR;
    STPSZ2 = STEPSZ * STEPSZ;
    AFAC = 2.0 * AM / (HBARC * HBARC);
    DL2 = (double)L * (L + 1);
//
    NSTEPS = (int)(ASYMPT / STEPSZ + 1.5);
//
//     WE DO THE MATCHING BEYOND BOTH THE REAL RADIUS
//     AND THE LAST NODE.
//
    NMTOP = (int)(DMAX1(R, A) / STEPSZ + 0.5);
    NMTOP = MAX0(NMTOP, 20);
    NMBOT = MAX0(NMTOP - 20, 3);
//
    IWF2 = 0;
    IVOUT = 0;
//
    NODEH = NODES;
    if (IHSAVE == 2) goto L900;
//
//  ASK FOR SPACE FOR POTENTIAL
    IV1L = NALLOC("BOUNDV1 ", NSTEPS);
    IV2L = NALLOC("BOUNDV2 ", NSTEPS);
    IVOUT = NALLOC(NAMES_arr[1], NSTEPS);
//     AND SPACE FOR THE WAVE FUNCTIONS
    IWF1  = NALLOC("BNDPHI1 ", NSTEPS);
    IWF2  = NALLOC(NAMES_arr[0], NSTEPS);
    IWF3  = NALLOC("BNDPHI3 ", NSTEPS);
//
//     INITIALIZE THE LINKULES FOR CENTRAL, COULOMB, SPIN-ORBIT,
//     AND WAVE FUNCTION.
//
    for (II = 1; II <= 4; II++) {
        I = IPTYPS[II];
        if (LNKBLK.LNKADR[I][3] == 0) goto L99;
        IUNIQU = IUNIQU + 1;
        std::snprintf(linkfmtout, sizeof(linkfmtout), "*%03d", IUNIQU);
        std::memcpy(LINKID, linkfmtout, 4);
        LINKUL(LNKBLK.LNKADR[I][3],
               *(char8*)&LNKBLK.LNKADR[I][1],
               &LNKBLK.LNKADR[I][5],
               I, 1,
               IRETUR, L, JP, 0.0, STEPSZ, NSTEPS, &DUMMY, &DUMMY, LINKID);
        if (IRETUR < 0) goto L295;
        LNKBLK.LNKADR[I][4] = IRET;
    L99:;
    }
//
//     NO MORE ALLOCATOR CALLS SO GET LOCATIONS FOR REST OF SUBROUTINE
//
    LV1      = Z[IV1L];
    LV2      = Z[IV2L];
    LVOUT    = Z[IVOUT];
    LWF1     = Z[IWF1];
    LWF2     = Z[IWF2];
    LWF3_loc = Z[IWF3];
//
//
//
//  AFINE IS FINE STRUCTURE CONSTANT
// NEED CI1,CI2 FOR COULOMB POTENTIAL INSIDE NUCLEAR RADIUS
//
    CI1 = IZP * IZT * HBARC * 1.5 / (AFINE * RC);
    CI2 = IZP * IZT * HBARC * 0.5 / (AFINE * RC * RC * RC);
//
//  C0 IS COULOMB OUTSIDE OF POTENTIAL
//
    C0 = IZP * IZT * HBARC / AFINE;
//
//
//  CALCULATE POTENTIAL AT ALL STEPS
//
//     THE ARRAYS ARE
//
//     ALLOC(LV1) = -(FIXED PART OF V) -- FIXED PART INCLUDES COULOMB
//          D      AND VSO PART OF SPIN ORBIT.
//     ALLOC(LV2) = -(VARIABLE PART OF V)/DEPTH -- INCLUDES REAL W.S.
//                  WELL AND TAU PART OF SPIN ORBIT.
//
//     TOTAL POTENTIAL IS
//        -( ALLOC(LV1) + V*ALLOC(LV2) )
//     ALLOC(LVOUT) = -(FIXED PART OF EFFECTIVE TRANSITION OPERATOR)
//
//     IN EACH CASE
//        ALLOC(L+I) = V( STEPSIZE*(I-1) )
//
//
//     GET THE COULOMB POTENTIAL.
//
    VMUL = -1.0;
    I = 5;
    LL = LV1;
    JUMP = 150;  // ASSIGN 150 TO JUMP
    if (LNKBLK.LNKADR[5][3] > 0) goto L290;
    RVAL = -STEPSZ;
    for (int II_c = 1; II_c <= NSTEPS; II_c++) {
        RVAL = RVAL + STEPSZ;
        if (RVAL >= RC) goto L120;
        VC = CI1 - CI2 * (RVAL * RVAL);
        goto L130;
    L120:
        VC = C0 / RVAL;
    L130:
        ALLOC(LV1 - 1 + II_c) = -VC;
    }
//
//     INITIALIZE THE FIXED PART OF THE EFFECTIVE TRANSITION OPERATOR.
//
L150:
    for (int II_t = 1; II_t <= NSTEPS; II_t++) {
        ALLOC(LVOUT - 1 + II_t) = 0.0;
        if (NUCONL != 1) ALLOC(LVOUT - 1 + II_t) = ALLOC(LV1 - 1 + II_t);
    }
//
//     NOW THE CENTRAL POTENTIAL
//
    VMUL = -1.0 / V;
    I = 1;
    LL = LV2;
    JUMP = 180;  // ASSIGN 180 TO JUMP
    if (LNKBLK.LNKADR[1][3] > 0) goto L290;
    WOODSX(NSTEPS, 0.0, STEPSZ, ALLOC_base(LV2), N1, N2,
           1, -1.0, R, A, POWRL);
//
//     NOW THE SPIN-ORBIT POTENTIAL (IF ANY)
//
L180:
    if (!SOSW) goto L300;
//     THIS IS SIGMA.L.  NOTE THAT JP, JSPS ARE DOUBLED.
    ALS = JP * (JP + 2) - JSP * (JSP + 2);
    ALS = (0.25 * ALS - DL2) / JSP;
    BLS = TAU * ALS;
    ALS = VSO * ALS;
    VSOSV = VSO;
    VSO = 1.0;
    VMUL = -1.0;
    I = 3;
    LL = LWF1;
    JUMP = 200;  // ASSIGN 200 TO JUMP
    if (LNKBLK.LNKADR[3][3] > 0) goto L290;
    WOODSX(NSTEPS, 0.0, STEPSZ, ALLOC_base(LWF1), N1, N2,
           2, -1.0, RSO, ASO, 1.0);
//     SPLIT SPIN-ORBIT INTO FIXED (VSO) AND VARIABLE (TAU) PARTS.
L200:
    for (int II_s = 1; II_s <= NSTEPS; II_s++) {
        ALLOC(LV1 - 1 + II_s) = ALLOC(LV1 - 1 + II_s) + ALS * ALLOC(LWF1 - 1 + II_s);
        ALLOC(LVOUT - 1 + II_s) = ALLOC(LVOUT - 1 + II_s) + ALS * ALLOC(LWF1 - 1 + II_s);
        ALLOC(LV2 - 1 + II_s) = ALLOC(LV2 - 1 + II_s) + BLS * ALLOC(LWF1 - 1 + II_s);
    }
    VSO = VSOSV;
    goto L300;
//
//     GET THE POTENTIAL FROM A LINKULE.
//
L290:
    {
        double neg1 = -1.0;
        LINKUL(LNKBLK.LNKADR[I][3],
               *(char8*)&LNKBLK.LNKADR[I][1],
               &LNKBLK.LNKADR[I][5],
               I, 3, IRETUR, L, JP,
               0.0, STEPSZ, NSTEPS, ALLOC_base(LL), &neg1, (char*)nullptr);
    }
    if (IRETUR < 0) goto L295;
    // GO TO JUMP, (150, 180, 200) — assigned GOTO dispatch
    if (JUMP == 150) goto L150;
    if (JUMP == 180) goto L180;
    if (JUMP == 200) goto L200;
L295:
    IRET = 0;
    return;
//
//
//     NOW CALCULATE THE WAVE FUNCTION.
//
//     IF A LINKULE DOES IT, SKIP MOST OF THE FOLLOWING.
//
L300:
    if (LNKBLK.LNKADR[6][3] > 0) goto L800;
//
//     WE FIND SOLUTIONS FOR THREE VALUES OF V OR K.  DELTA IS THE
//     RELATIVE SPACING TO USE.  AS APPROPRIATE, DELTAV OR DELTAK ARE
//     THE ABSOLUTE SPACINGS TO USE.  SOLUTIONS ARE FOUND FOR
//     V-DELTAV, V, V+DELTAV   OR   K-DELTAK, K, K+DELTAK
//
L3020:
    DELTA = DELTVK;
    if (PCONSW) std::printf(" %5d%16.8G%16.8G%16.8G\n", ITER, V, UK, DELTA);
    VST = V;
    UKST = UK;
//
//     SET K+, K, K- AND CORRESPONDING ETAS EACH TIME E CHANGES
//     (THIS IS FOR WHEN E AND NOT V IS BEING VARIED)
//
L3030:
    if (IFIT == 1) goto L301;
    DELTAK = UK * DELTA;
    DELTAV = 0.0;
    goto L302;
L301:
    DELTAV = V * DELTA;
    DELTAK = 0.0;
//
L302:
    UKM = UK - DELTAK;
    UKP = UK + DELTAK;
    UKS[1] = UKM;
    UKS[2] = UK;
    UKS[3] = UKP;
    for (int II_k = 1; II_k <= 3; II_k++) {
        UK2S[II_k] = UKS[II_k] * UKS[II_k];
        ETAS_loc[II_k] = IZP * IZT * AM / (HBARC * AFINE * UKS[II_k]);
    }
    ETA = ETAS_loc[2];
//
//     FIND ASYMPTOTIC FORM AS A WHITTAKER FUNCTION.
//


    // --- BOUND PART 2 ---

    // (Part 2 aliases and locals moved to top of function)

    // Sync LWFS with LWF1/LWF2/LWF3_loc (Fortran EQUIVALENCE)
    LWFS[1] = LWF1;
    LWFS[2] = LWF2;
    LWFS[3] = LWF3_loc;

//
//     FIND ASYMPTOTIC FORM AS A WHITTAKER FUNCTION.
//     IASYMP = 1 ( SKIPASYMP )  MEANS TO JUST USE 2 AND 1 AS THE
//            LAST TWO FUNCTION VALUES.  IF THIS IS ELECTED THEN
//            ASYMPTOPIA MUST BE EXTENDED BY ABOUT 5 FM TO GIVE
//            THE WAVEFUNCTION TIME TO HEAL TO THE CORRECT FORM.
//
    RGET = ASYMPT - STEPSZ;
//
//     WE NEED A NORM TERM THAT IS SMALL FOR SMALL RHO AND
//     ETA+L+1 FOR LARGE RHO.  WHEN THE WHITTAKER IS DIVIDED BY
//     GAMMA(THIS NORM) IT WILL GENERALLY NOT BE TOO LARGE OR SMALL
//
    TEMP = RGET * UK;
    TEMP = (L + ETA + 1) * std::pow(TEMP / (2.0 + TEMP), 2);
//
    for (II = 1; II <= 2; II++) {
        for (I = 1; I <= 3; I++) {
            if (I > 1 && IFIT == 1) goto L345;
            if (IASYMP == 1) goto L340;
            X = COULNG(L, ETAS_loc[I], UKS[I] * RGET, TEMP);
//
//     COULNG RETURNS THE WHITTAKER TIMES A FACTOR.  HERE WE
//     PATCH UP THE RATIO BETWEEN W(ASYMPT-STEPSZ) AND W(ASYMPT)
//
            if (II == 1) X = X *
                std::exp(UKS[I] * STEPSZ + ETAS_loc[I] * std::log(ASYMPT / RGET));
            goto L345;
    L340:
            X = (3 - II);
    L345:
            ALLOC(LWFS[I] + NSTEPS - 3 + II) = X;
        }
        RGET = RGET + STEPSZ;
    }
    if (PCONSW) std::printf(" ASYMPTOTIC FORM AT LAST TWO STEPS: %18.8G%18.8G\n",
        ALLOC(LWF2 + NSTEPS - 2), ALLOC(LWF2 + NSTEPS - 1));
//
//  NOW HAVE ASYMPTOTIC WAVEFUNCTIONS
//
//
//
//
//
//  STEP INTO MATCHING RADIUS
//
//     FIND THE MAXIMUM VALUE OF THE MINIMUM OVER V-DELTA, V, V+DELTA,
//     OF THE WAVEFUNCTION IN THE INNER 20 STEPS AND USE AS THE MATCHING
//     RADIUS.
//
//
//     WE COME HERE FOR EACH NEW VALUE OF V
//
L500:
    VS_loc[1] = V - DELTAV;
    VS_loc[2] = V;
    VS_loc[3] = V + DELTAV;
//
    XX = 0.0;
    NGO = NSTEPS - NMBOT + 1;
    for (I = 3; I <= NGO; I++) {
        II = NSTEPS - I + 1;
        ALL = DL2 / std::pow(STEPSZ * II, 2);
        TEMP = BIGNUM;
        for (IVK = 1; IVK <= 3; IVK++) {
            K = LWFS[IVK] + II;
            X = VS_loc[IVK] * ALLOC(LV2 + II) + ALLOC(LV1 + II);
            X = 2.0 + STPSZ2 * (UK2S[IVK] - AFAC * X + ALL);
            ALLOC(K - 1) = X * ALLOC(K) - ALLOC(K + 1);
            TEMP = DMIN1(TEMP, std::fabs(ALLOC(K)));
        }
        if (II > NMTOP) goto L529;
        if (TEMP <= XX) goto L529;
        NMATCH = II;
        XX = TEMP;
    L529:;
    }
//
//     COMPUTE OUTER MATCHING PARAMETERS -- DERIVATIVES, NORM
//     OF OUTER PART, VALUE AT MATCHING POINT AND NUMBER OF NODES.
//
    for (IVK = 1; IVK <= 3; IVK++) {
        K = LWFS[IVK] + NMATCH;
        DROUT[IVK] = (0.5 / STEPSZ) * (ALLOC(K + 1) - ALLOC(K - 1));
        VALOUT[IVK] = ALLOC(K);
        NODCNT[IVK] = 0;
        SUMIN[IVK] = 0.0;
        SUMOUT[IVK] = 0.5 * (std::pow(ALLOC(K), 2) + std::pow(ALLOC(LWFS[IVK] + NSTEPS - 1), 2));
        NGO = NMATCH + 3;
        for (II = NGO; II <= NSTEPS; II++) {
            K = LWFS[IVK] - 2 + II;
            SUMOUT[IVK] = SUMOUT[IVK] + std::pow(ALLOC(K), 2);
//     HERE WE COUNT THE NODE IF IT IS EXACTLY AT R(MATCH)
            if (ALLOC(K) * ALLOC(K - 1) <= 0.0) NODCNT[IVK] = NODCNT[IVK] + 1;
        }
    }
//
//  FIND POWER SERIES SOLUTION (KEEP FIRST TWO TERMS) IN ORDER
//  TO FIX FUNCTION AT FIRST STEP.  THEN STEP TO MATCHING RADIUS
//
//     FOLLOWING IS AN APPROXIMATE EXPRESION AND NEED NOT BE TOO
//     ACCURATE.  THIS EXPRESSION IS SUCH THAT FOR L = 0
//        PHI(R) = U(R)/R ---> 1   FOR R = 0.  (BEFORE NORMALIZATION)
//     WE USE THIS FACT ONLY FOR L = 0 SO WE CAN LATER ON EXTRACT
//     PHI(0).  REALLY THIS IS JUST FINESSE SINCE PHI(0) IS NEVER
//     SIGNIFICANT.  SINCE THIS VALUE IS SIGNIFICANT ONLY FOR L = 0,
//     WE ASSUME L = 0 AND THUS IGNORE SPIN-ORBIT FORCES.
//
//     INCLUDE IT.
//
    for (IVK = 1; IVK <= 3; IVK++) {
        K = LWFS[IVK];
        ALLOC(K) = 0.0;
//
//     POTENTIAL + K**2  AT R = STEPSZ.
//
        TEMP = -AFAC * (ALLOC(LV1 + 1) + VS_loc[IVK] * ALLOC(LV2 + 1))
            + UK2S[IVK];
        ALLOC(K + 1) = STEPSZ * (1.0 + STPSZ2 * TEMP / 6.0);
    }
//
//  STEP OUT TO MATCHING RADIUS
//
    for (II = 1; II <= NMATCH; II++) {
        ALL = DL2 / std::pow(STEPSZ * II, 2);
        for (IVK = 1; IVK <= 3; IVK++) {
            K = LWFS[IVK] + II;
            X = VS_loc[IVK] * ALLOC(LV2 + II) + ALLOC(LV1 + II);
            X = 2.0 + STPSZ2 * (UK2S[IVK] - AFAC * X + ALL);
            ALLOC(K + 1) = X * ALLOC(K) - ALLOC(K - 1);
            SUMIN[IVK] = SUMIN[IVK] + std::pow(ALLOC(K), 2);
            if (ALLOC(K) * ALLOC(K - 1) < 0.0) NODCNT[IVK] = NODCNT[IVK] + 1;
        }
    }
//
//     MATCH SOLUTIONS;  DERIVATIVES FOR INNER SIDE.
//
    for (IVK = 1; IVK <= 3; IVK++) {
        K = LWFS[IVK] + NMATCH;
        SUMIN[IVK] = SUMIN[IVK] - 0.5 * std::pow(ALLOC(K), 2);
//
//     SUMOUT MULTIPLIES THE OUTER SOLUTION TO NORMALIZE IT
//     SUMIN  MULTIPLIES THE INNER SOLUTION TO NORMALIZE IT
//
        XX = ALLOC(K) / VALOUT[IVK];
        X = 1.0 / std::sqrt(STEPSZ * (SUMIN[IVK] + XX * XX * SUMOUT[IVK]));
        SUMIN[IVK] = X;
        SUMOUT[IVK] = XX * X;
        VALOUT[IVK] = (XX * X) * VALOUT[IVK];
        DROUT[IVK] = (XX * X) * DROUT[IVK];
        DRIN[IVK] = X * (0.5 / STEPSZ) * (ALLOC(K + 1) - ALLOC(K - 1));
        PHIS[IVK] = DROUT[IVK] - DRIN[IVK];
    }
//
//     WE LOOK FOR A ZERO OF THE DIFFERENCE OF THE INNER AND OUTER
//     DERIVATIVES.
//
    NODEH = NODCNT[2];
    X = NMATCH * STEPSZ;
    if (PCONSW) std::printf(" DERIVATIVES AT R =%7.2f FM"
        "%57sV =%15.6G%15.6G%15.6G\n"
        "     INNER =%16.6G%16.6G%16.6G"
        "%35sK =%15.6G%15.6G%15.6G\n"
        "    OUTER =%16.6G%16.6G%16.6G"
        "%32sU(R) =%15.6G%15.6G%15.6G\n"
        "       DIF =%16.6G%16.6G%16.6G"
        "%31sNODES =%9d%6s%9d%6s%9d%6s\n",
        X, "", VS_loc[1], VS_loc[2], VS_loc[3],
        DRIN[1], DRIN[2], DRIN[3],
        "", UKS[1], UKS[2], UKS[3],
        DROUT[1], DROUT[2], DROUT[3],
        "", VALOUT[1], VALOUT[2], VALOUT[3],
        PHIS[1], PHIS[2], PHIS[3],
        "", NODCNT[1], "", NODCNT[2], "", NODCNT[3], "");
    if (std::fabs(PHIS[2]) < ACCURA) goto L600;
    PPHI2 = std::fabs(PHIS[2]);
    if (PCONSW) std::printf(" FOR ITER=%3d THE  FUNCTION IS%11.3G >%11.3G\n",
        ITER, PPHI2, ACCURA);
    if (ITER > MAXITR) goto L1001;
    DP = (PHIS[3] - PHIS[1]) / (2.0 * DELTA);
    DPP = (PHIS[3] - 2.0 * PHIS[2] + PHIS[1]) / (DELTA * DELTA);
    XL = -PHIS[2] / DP;
    XQ_loc = 0.0;
    X = XL;
//
//  IF DPP IS ZERO DO LINEAR APPROXIMATION
//
    if (std::fabs(DPP) < 1.0e-5 * std::fabs(DP)) goto L7398;
    { double CSQRT_val = (DP * DP) - 2.0 * PHIS[2] * DPP;
    if (CSQRT_val < 0.0) goto L7398;
//
    XQ_loc = (-DP + (DP / std::fabs(DP)) * std::sqrt(CSQRT_val)) / DPP;
    }
    if (LINEAR != 1) X = XQ_loc;
//
L7398:
    if (std::fabs(X) > 0.5) X = DSIGN(0.50, X);
    if (PCONSW) std::printf(" DP,DPP =%14.5G%14.5G\n"
        " LINEAR, QUADRATIC METHODS =%12.3G%12.3G"
        "     REL. STEP OF V OR K USED =%14.5G\n",
        DP, DPP, XL, XQ_loc, X);
//
//     FIND A NEW DELTAVK.  WE USE THE STEP LENGTH USED IN THIS
//     ITERATION UNLESS IT STAYED WITHIN THE THE PREVIOUS DELTA,
//     INWHICH CASE WE USE .5 OF IT.
//
    XX = DELTA;
    DELTA = std::fabs(X);
    if (DELTA < XX) DELTA = 0.5 * DELTA;
    if (DELTA < 1.0e-7) DELTA = 1.0e-7;
    DELTA = DMIN1(DELTA, 0.20);
    if (IFIT != 0) goto L3070;
    XX = UK * (1.0 + X);
    if (XX <= 0.0) UK = UK / 2.0;
    if (XX > 0.0) UK = XX;
    goto L3737;
L3070:
    XX = V * (1.0 + X);
    if (XX > 0.0) V = XX;
    if (XX <= 0.0) V = 0.5 * V;
L3737:
    ITER = ITER + 1;
    if (PCONSW) std::printf("\nITERATION%3d     V =%16.8G  K =%16.8G"
        "     DELTAVK =%12.4G\n", ITER, V, UK, DELTA);
    if (IFIT == 0) goto L3030;
    DELTAV = DELTA * V;
    goto L500;
//
//  MAXITRTER EXCEDED.  SEE IF WE HAVE RIGHT NUMBER OF NODES
//
L1001:
    if (NODES != NODEH) goto L6666;
    std::printf(" ******** NO SOLUTION FOUND *********\n");
//
//     RETURN AFTER A FAILURE
//
L1111:
    IRET = 0;
    goto L850;
//
//     HAVE CONVERGED ON A V OR UK - CHECK FOR PROPER NUMBER OF NODES
//
L600:
    PPHI2 = std::fabs(PHIS[2]);
    if (PCONSW) std::printf(" ON ITERATION %3d THE FUNCTION IS%11.3G <%11.3G\n",
        ITER, PPHI2, ACCURA);
    E = -UK2S[2] * HBARC * HBARC / (2.0 * AM);
    if (PCONSW) std::printf("\n   FIT LOCATED WITH V=%14.8f AND E=%14.8f\n\n",
        V, E);
//
//
//  RENORMALIZE WAVEFUNCTION AND CHECK NODES
//
    if (NODEH == NODES) goto L6041;
//
//     SOLUTION HAS WRONG NUMBER OF NODES - TRY ANOTHER V OR UK
//
L6666:
    std::printf(" FOR V =%8.2f,  E =%7.2f, WAVEFUNCTION HAS"
        "%3d NODES BUT%3d ARE DESIRED.\n", V, E, NODEH, NODES);
    if (ITERN <= MAXITR) goto L6668;
    std::printf("\n**** COULD NOT FIND A SOLUTION WITH%3d NODES.\n", NODES);
    goto L1111;
L6668:
    ITER = 0;
    ITERN = ITERN + 1;
//
//     WE USE
//
//       X = DESIRED NUMBER OF NODES
//       XQ = NUMBER OF NODES WE JUST GOT
//       XP = PREVIOUS NUMBER OF NODES
//       PHI1 = PRESENT V OR K
//       PHI2 = NEXT V OR K
//       PHIP = PREVIOUS V OR K
//
    X = NODES;
    XQ_loc = NODEH;
//
//     WE USE THE STARTING OR FINAL V OR K, WHICH EVER IS CLOSER
//     TO WHERE WE WANT TO GO.
//
    if (NODES < NODEH) goto L8130;
    V = DMAX1(V, VST);
    UK = DMIN1(UK, UKST);
    goto L8150;
L8130:
    V = DMIN1(V, VST);
    UK = DMAX1(UK, UKST);
L8150:
    PHIS[1] = V;
    if (IFIT == 0) PHIS[1] = UK * UK;
    if (ITERN <= 1) goto L8160;
    if (XP_loc != XQ_loc) goto L8200;
//
//     WE HAVE ONLY ONE POINT TO WORK WITH
//
L8160:
    if (IFIT == 0) goto L8170;
    PHIS[2] = (X + 1.0) * PHIS[1] / (XQ_loc + 1.0);
    goto L8300;
L8170:
    PHIS[2] = (1.0 + (XQ_loc - X) / (XQ_loc + 1.0)) * PHIS[1];
    goto L8300;
//
//     HAVE TWO POINTS - USE STRAIGHT LINE
//
L8200:
    PHIS[2] = PHIS[1] + (X - XQ_loc) * (PHIS[1] - PHIP) / (XQ_loc - XP_loc);
//
//     DO NOT ALLOW TO GO NEGATIVE
//
L8300:
    if (PHIS[2] <= 0.0) PHIS[2] = 0.5 * PHIS[1];
//
    PHIP = PHIS[1];
    XP_loc = XQ_loc;
    if (IFIT == 1) goto L8350;
    UK = std::sqrt(PHIS[2]);
    E = -(UK * HBARC) * (UK * HBARC) / (2.0 * AM);
    std::printf(" WILL TRY A NEW E: %14.5G\n", E);
    goto L3020;
L8350:
    V = PHIS[2];
    std::printf(" WILL TRY A NEW V:  %14.5G\n", V);
    goto L3020;
//
//     SOLUTION HAS PROPER NUMBER OF NODES - CLEAN UP
//
L6041:
    if (PCONSW) std::printf("  SOLUTION HAS %1d NODE(S). L=%2d JP=%2d/2\n",
        NODES, L, (int)JP);
    goto L850;
//
//
//     WAVE FUNCTION IS CALCULATED BY A LINKULE.
//     ARRAY1 IS THE ANSWER, ARRAY2 IS (-1/V)*THE VARIABLE POTENTIAL,
//     ID IS A POINTER TO -1*THE FIXED POTENTIAL.
//
L800:
    LINKUL(LNKBLK.LNKADR[6][3],
           *(char8*)&LNKBLK.LNKADR[6][1],
           &LNKBLK.LNKADR[6][5],
           6, 3,
           IRETUR, L, JP, 0.0, STEPSZ, NSTEPS,
           ALLOC_base(LWF2), ALLOC_base(LV2), (char*)&LV1);
    if (IRETUR < 0) goto L295;
    SUMIN[2] = 1.0;
    SUMOUT[2] = 1.0;
    NMATCH = NSTEPS / 2;
    NODEH = NODES;
//
//     NORMALIZE THE SOLUTION AND CONVERT TO   PHI = U/R
//
L850:
    RVAL = STEPSZ;
    NGO = NMATCH + 1;
    for (II = 1; II <= NGO; II++) {
        ALLOC(LWF2 + II) = SUMIN[2] * ALLOC(LWF2 + II) / RVAL;
        RVAL = RVAL + STEPSZ;
    }
    NGO = NGO + 2;
    for (II = NGO; II <= NSTEPS; II++) {
        ALLOC(LWF2 - 1 + II) = SUMOUT[2] * ALLOC(LWF2 - 1 + II) / RVAL;
        RVAL = RVAL + STEPSZ;
    }
//
//     FOR L = 0 GET WAVEFUNCTION AT ORIGIN BASED ON STARTING
//     VALUE DEFINED ABOVE
//
    if (L == 0) ALLOC(LWF2) = SUMIN[2];
//
//     FOR L=0 WAVEFUNCTIONS FROM LINKULES, WE EXTRAPOLATE, USING THE
//     FACT THAT THE DERIVATIVE = 0 AT R=0.
//
    if (L == 0 && LNKBLK.LNKADR[6][3] > 0)
        ALLOC(LWF2) = (4.0 * ALLOC(LWF2 + 1) - ALLOC(LWF2 + 2)) / 3.0;
//
//  FREE SECONDARY AREAS HERE
//
    Z[IWF1] = -Z[IWF1];
    Z[IWF3] = -Z[IWF3];
    Z[IV1L] = -Z[IV1L];
//
//     PRINT RESULTS NICELY
//
L900:
    if (!PRNTSW) goto L990;
    I = ICHANB;
    if (I > 2) I = 3;
    UK = std::sqrt(-2.0 * AM * E) / HBARC;
    ETA = IZP * IZT * AM / (HBARC * AFINE * UK);
    std::printf("\n0%8s%.8s%.2s BOUND STATE PARAMETERS\n",
        "", CHANMS[I-1][0], CHANMS[I-1][1]);
    std::printf("0E =%10.4f MEV     KAPPA =%8.5f\n", E, UK);
    {
    double Q1 = 0.0;
    double Q2 = 0.0;
    double Q3 = 0.0;
    }
    if (AMP == UNDEF || AMT == UNDEF) goto L920;
    std::printf(" PROJECTILE MASS =%7.2f AMU     "
        "TARGET MASS =%7.2f AMU     "
        "REDUCED MASS =%10.2f MEV/C**2\n", AMP, AMT, AM);
    goto L925;
L920:
    std::printf(" REDUCED MASS =%10.2f MEV/C**2\n", AM);
L925:
    std::printf(" L =%3d%9d NODES\n", L, NODEH);
    if ((int)JSP != NOTDEF_INT && (int)JST != NOTDEF_INT)
        std::printf(" PROJECTILE SPIN =%3d/2     "
            "TARGET SPIN =%3d/2\n", (int)JSP, (int)JST);
    if ((int)JP != NOTDEF_INT)
        std::printf(" J PROJECTILE =%3d/2\n", (int)JP);
    if ((int)J_val != NOTDEF_INT)
        std::printf(" TOTAL J =%3d/2\n", (int)J_val);
    std::printf(" PROJECTILE PARITY IS %.7s"
        "   TARGET PARITY IS %.7s"
        "   TOTAL PARITY IS %.7s\n",
        PWORD[PARIPT[1] + 1], PWORD[PARIPT[2] + 1], PWORD[PARITY + 1]);
    std::printf(" Z PROJECTILE =%4d     Z TARGET =%4d\n"
        "0POTENTIAL         COUPLING CONS.    RADIUS    DIFFUSENESS    RADIUS PARAM.\n",
        IZP, IZT);
    PARPRT(1, "REAL CENTRAL      ", V, R, A, R/R0MASS, POWRL, "   POWER", IRETUR);
    PARPRT(3, "REAL SPIN-ORBIT   ", VSO, RSO, ASO, RSO/R0MASS, TAU, "     TAU", IRETUR);
    PARPRT(5, "PT.&SPHERE COULOMB", ETA, RC, 0.0, RC/R0MASS, 0.0, "        ", IRETUR);
//
    TEMP = DMIN1(1.0 / UK, A) / STEPSZ;
    std::printf("0ASYMPTOPIA = %8.3f FM\n"
        " STEP SIZE =  %8.3f FM     %8.3f STEPS PER \"ASYMPTOTIC RANGE\"\n\n",
        ASYMPT, STEPSZ, TEMP);
//
//     IF THE WAVEFUNCTION IS CALCULATED BY A LINKULE, CALL IT TO
//     PRINT STUFF ABOUT ITSELF.
//
    if (LNKBLK.LNKADR[6][3] == 0) goto L990;
    std::printf("0 THE BOUND WAVEFUNCTION WAS COMPUTED BY THE %.8s LINKULE:\n",
        (char*)&LNKBLK.LNKADR[6][1]);
    LINKUL(LNKBLK.LNKADR[6][3],
           *(char8*)&LNKBLK.LNKADR[6][1],
           &LNKBLK.LNKADR[6][5],
           6, 2,
           IRETUR, L, JP, 0.0, STEPSZ, NSTEPS,
           ALLOC_base(LWF2), ALLOC_base(LV2), (char*)&LV1);
    if (IRETUR < 0) goto L295;
    std::printf("\n");
//
//
//     SAVE THE POTENTIAL AND OTHER STUFF FOR GRDSET
//
L990:
    if (ICHANB > 2 || ICHANB < 1) goto L6330;
    FORMF.IBDS[ICHANB] = IWF2;
    FORMF.NPOTS[ICHANB] = IVOUT;
    LBNDS[ICHANB - 1] = L;
    NBNDS[ICHANB - 1] = NODES;
    JBNDS[ICHANB - 1] = (int)JP;
    FORMF.BNDMAX[ICHANB] = ASYMPT;
    FORMF.BNDSTP[ICHANB] = STEPSZ;
    FORMF.NSTPBD[ICHANB] = NSTEPS;
    FORMF.RBNDS[ICHANB] = R;
    FORMF.ABNDS[ICHANB] = A;
    FORMF.RCBNDS[ICHANB] = RC;
    FORMF.BNDMAS[ICHANB] = AM;
    EBNDS[ICHANB] = E;
    FORMF.BRATMS[ICHANB] = RATMAS;
//
//     THE FOLLOWING WILL BE NEEDED FOR THE NEXT PASS (IF ANY).
//
    FORMF.RSOBDS[ICHANB] = RSO;
    FORMF.ASOBDS[ICHANB] = ASO;
    FORMF.VBNDS[ICHANB] = V;
    FORMF.VSOBDS[ICHANB] = VSO;
    FORMF.TAUBDS[ICHANB] = TAU;
    for (I = 1; I <= 12; I++) {
        for (II = 1; II <= 6; II++) {
            FORMF.LNKADB[I][ICHANB][II] = LNKBLK.LNKADR[I][II];
        }
    }
//
//     CAN WE COMPUTE THE Q OF THE REACTION
//
    if (EBNDS[3 - ICHANB] == UNDEF) goto L6330;
    Q = EBNDS[1] - EBNDS[2];
    if (ISTRIP == -1) Q = -Q;
//
L6330:
    if (IHSAVE == 2) goto L8875;
//
//     CONSTRUCT THE COMPLETE POTENTIAL FOR V PHI
//
    for (I = 1; I <= NSTEPS; I++) {
        ALLOC(LVOUT - 1 + I) = -(ALLOC(LVOUT - 1 + I)
            + V * ALLOC(LV2 - 1 + I));
    }
    Z[IV2L] = -Z[IV2L];
//
//     IF REQUESTED, PRINT THE WAVEFUNCTION
//
    RVAL = 0.0;
    if (WRITES == UNDEF) goto L8875;
    NWRITE = (int)DMAX1(WRITES / STEPSZ + 0.5, 1.0);
    std::printf(" PSI=F(R)/R\n"
        "    STEP      RADIUS        WAVEFUNCTION\n\n");
    for (II = 1; II <= NSTEPS; II += NWRITE) {
        RVAL = (II - 1) * STEPSZ;
        std::printf("%8d%12.5f%15.5G\n", II, RVAL, ALLOC(LWF2 + II - 1));
    }
//
L8875:
    RC = RCSAVE;
    T1 = (float)dtime_() - T1;
    if (PCONSW) std::printf(" BOUND STATE CPU TIME =%7.3f SECONDS.\n", T1);
    return;
//
//
//
}
