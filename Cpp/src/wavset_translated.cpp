// wavset_translated.cpp — WAVPOT, WAVSET
// Translated from source.f lines 36991-38020
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

// Forward declarations for routines called here but not in ptolemy_forward.h
// with full 5-argument SETVSQ signature
extern void SETVSQ(double RR1, double RR2, int IZ1, int IZ2, int K);

// ============================================================================
// SUBROUTINE WAVPOT ( IRET )
//
//     GETS L-DEPENDENT STUFF (COULOMB FCNS, S-ARRAYS) FOR WAVELJ
//
//     THIS IS THE SECOND PART OF THE OLD WAVEF.  IT SETS UP THE
//     ARRAYS FOR THE S-MATRICES AND FINDS THE COULOMB
//     BOUNDRY CONDITIONS.  POINTERS AND PARAMETERS ARE LEFT IN 2-ELEMENT
//     ARRAYS THAT ARE LABELED BY THE CHANNEL NUMBER ( ICHANW = 1, 2 ).
//
//     2/4/75 - MASSIVE RUPTURE OF WAVEF - S. P.
//     2/10/75 - COMPUTE COULOMB PHASE SHIFTS FOR INELDC
//     2/24/75 - INTEGRATE TO SEVERAL STEPS PAST LARGEST R NEEDED
//     2/25/75 - LOTS OF INPUT CHECKS; DON'T DIVIDE BY 0 DIFUSNESS
//     4/16/75 - NEW COMMON BLOCKS; REMOVE S.O. SINCE IT IS WRONG
//     4/20/75 - SPLIT INTO WAVSET AND WAVPOT.
//     5/21/75 - DO THINGS UPTO LIMOST
//     5/22/75 - POTENTIALS MOVED TO WAVELJ; STANSW
//     10/7/75 - POTENTIALS MOVED TO MAKPOT
//     10/15/75 - ELAB-DEPENDANCE OF STEPSIZE IN ELASTIC FITS.
//     3/23/76 - USE MANCHESTER COULOMB ROUTINE
//     4/15/77 DO NOT COMPUTE POTENTIAL WHEN NOT NEEDED
//     11/20/77 - AREA FOR MAKPOT TO CONSTRUCT POTENTIAL - S.P.
//                MUCH MOVED TO WAVSET.
//     4/24/78 - REVERSE LOMOST FOR LIMOST - S.P.
//     5/24/78 - INITIALIZE SMAT TO 0 - S.P.
//     8/1/79 - STANDARD PREFIXES, S-MATRIX INDICES - RPG
//     1/18/80 - SKIP COULOMB F & G FOR C.C. - S.P.
//     4/27/80 - COMPUTE FOR ONE VALUE OF L FOR C.C. - S.P.
//     5/14/82 - DO NOT DESTROY LX - HURTS DEFORMED FITS - S.P.
//
// WAVPOT is in source_potentials.cpp — not duplicated here
//
// ============================================================================
// WAVSET
//
// ============================================================================
void WAVSET(int& IRET, int ALONSW, int FITSW)
{
    // CHARACTER*8  CHANAM(2), STATWD(2), SOTNAM(2), NAMES(2,11)
    static const char CHANAM[2][9] = { "INCOMING", "OUTGOING" };
    static const char STATWD[2][9] = { "BOSONS  ", "FERMIONS" };
    static const char SOTNAM[2][9] = { "SPIN-ORB", "TENSOR  " };
    // CHARACTER*4  TENNAM(6)
    static const char TENNAM[6][5] = { "TR", "TRI", "TL", "TLI", "TP", "TPI" };
    // CHARACTER*8  NAMES(2,11) — for potential arrays
    static const char NAMES[2][11][9] = {
        { "REALVA  ", "IMAGVA  ", "REALVSOA", "IMAGVSOA", "CENTRIFA",
          "REALVTRA", "IMAGVTRA", "REALVTLA", "IMAGVTLA", "REALVTPA", "IMAGVTPA" },
        { "REALVB  ", "IMAGVB  ", "REALVSOB", "IMAGVSOB", "CENTRIFB",
          "REALVTRB", "IMAGVTRB", "REALVTLB", "IMAGVTLB", "REALVTPB", "IMAGVTPB" }
    };
    // CHARACTER*1  LINKID(4)
    char LINKID[4];
    // character*8  linkfmtout — equivalenced to linkid(1)
    char linkfmtout[9];

    LINKID[0] = '*';
    LINKID[1] = '1';
    LINKID[2] = '0';
    LINKID[3] = '0';

    int SOSW, SOISW, IMAGSW, TENSW;
    double EELAB, EEMUL, TVR_l, TVI_l, TAR_l, TAI_l, TVSOR_l, TVSOI_l;
    double UK, ALAMBD, ETA, H, RMAX, TEMP, TEMP2_l, TEMP3_l;
    int IPRNT, NSTEP, LSKIP, ISTAT, NSTPP1, LXSAVE;
    int IRTN, LC1, LC2;
    double DUMMY;
//
    IPRNT = 0;
    if (FITSW) goto L10;
    WAVCOM.NUMFIT = 0;
    WAVCOM.IOFFIT = 1;
    IPRNT = INTGER.IPRINT % 10;
L10:
    WAVCOM.PWAVSW = ((INTGER.IPRINT / 1000) % 10) >= 1;
    WAVCOM.PWBGSW = ((INTGER.IPRINT / 1000) % 10) >= 4;
//
    WAVCOM.STANSW = ALONSW;
//
//     SOME OF THE PARAMETERS CAN HAVE AN ENERGY DEPENDENCE WHICH
//     IS EXPRESSED IN TERMS OF THE LABORATORY ENERGY
//
    EELAB = (1 + INTRNL.RATMAS) * FLOAT_common.E;
    EEMUL = EELAB;
    if (SWITCH.IEPOW == 1) EEMUL = 1 / EELAB;
    TVR_l = FLOAT_common.V + FLOAT_common.VE * EEMUL + FLOAT_common.VESQ * EEMUL * EEMUL;
    TVI_l = FLOAT_common.VI + FLOAT_common.VIE * EEMUL + FLOAT_common.VIESQ * EEMUL * EEMUL;
    TAR_l = FLOAT_common.A + FLOAT_common.AE * EEMUL + FLOAT_common.AESQ * EEMUL * EEMUL;
    TAI_l = FLOAT_common.AI + FLOAT_common.AIE * EEMUL + FLOAT_common.AIESQ * EEMUL * EEMUL;
    // Store in TEMPVS common for use by MAKPOT
    TEMPVS.TVR = TVR_l;
    TEMPVS.TVI = TVI_l;
    TEMPVS.TAR = TAR_l;
    TEMPVS.TAI = TAI_l;
//
//     IS THERE AN IMAGINARY PART FOR THE POTENTIAL
//
    IMAGSW = TVI_l != 0;
//
//     IS THERE A SPIN-ORBIT POTENTIAL
//
    TVSOR_l = FLOAT_common.VSO + FLOAT_common.TAU * TVR_l;
    TVSOI_l = FLOAT_common.VSOI + FLOAT_common.TAUI * TVI_l;
    TEMPVS.TVSOR = TVSOR_l;
    TEMPVS.TVSOI = TVSOI_l;
    SOSW = TVSOR_l != 0;
    SOISW = TVSOI_l != 0;
    WAVCOM.SOSWT = SOSW || SOISW;
//
//     IS THERE A TENSOR POTENTIAL?
//
    WAVCOM.TCSWS[INTRNL.ICHANW] = FLOAT_common.VTR != 0 || FLOAT_common.VTRI != 0 ||
        FLOAT_common.VTP != 0 || FLOAT_common.VTPI != 0;
    TENSW = WAVCOM.TCSWS[INTRNL.ICHANW] || FLOAT_common.VTL != 0 || FLOAT_common.VTLI != 0;
    WAVCOM.SOSWS[INTRNL.ICHANW] = WAVCOM.SOSWT || TENSW;
//
//     TEST THE INPUT AND RETURN  0  COMPLETION CODE IF BAD
//
    IRET = 1;
    if (LNKBLK.LNKADR[1][3] != 0) goto L25;
    if (FLOAT_common.R > 0 && TAR_l > 0) goto L25;
    std::printf("0**** R OR A HAS INVALID VALUE:%15.5G%15.5G\n", FLOAT_common.R, TAR_l);
    IRET = 0;
L25:
    if (!IMAGSW) goto L30;
    if (LNKBLK.LNKADR[2][3] != 0) goto L30;
    if (FLOAT_common.RI > 0 && TAI_l > 0) goto L30;
    std::printf("0**** RI OR AI HAS INVALID VALUE:%15.5G%15.5G\n", FLOAT_common.RI, TAI_l);
    IRET = 0;
L30:
    if (!SOSW) goto L34;
    if (LNKBLK.LNKADR[3][3] != 0) goto L34;
    if (FLOAT_common.RSO > 0 && FLOAT_common.ASO > 0) goto L34;
    std::printf("0**** RSO OR ASO HAS INVALID VALUE:%15.5G%15.5G\n", FLOAT_common.RSO, FLOAT_common.ASO);
    IRET = 0;
//
//     CHECK  SP  FOR EITHER REAL OR IMAG SPIN ORBIT.
//
L34:
    if (JBLOCK.JSP != INTRNL.NOTDEF || !WAVCOM.SOSWS[INTRNL.ICHANW]) goto L40;
    JBLOCK.JSP = 1;
    if (TENSW) JBLOCK.JSP = 2;
    std::printf("0**** WARNING:  SP WAS NOT DEFINED; "
        "IT IS ASSUMED TO BE %1d/2 FOR THE %.8s FORCE.\n",
        (int)JBLOCK.JSP, SOTNAM[(int)JBLOCK.JSP - 1]);
L40:
    if (FLOAT_common.E > 0 && FLOAT_common.AM > 0)
        goto L45;
    std::printf("0**** ENERGY OR M (REDUCED MASS) HAS INVALID"
        " VALUE:%15.5G%15.5G%15.5G\n", FLOAT_common.E, FLOAT_common.AM, 0.0);
    IRET = 0;
L45:
    if (FLOAT_common.AMP != INTRNL.UNDEF && FLOAT_common.AMT != INTRNL.UNDEF) goto L50;
    std::printf("0**** BOTH MP AND MT MUST BE DEFINED.\n");
    IRET = 0;
L50:
    if (FLOAT_common.RC < FLOAT_common.ASYMPT) goto L60;
    std::printf("0**** RC MUST BE LESS THAN ASYMPTOPIA:%15.5G%15.5G\n",
        FLOAT_common.RC, FLOAT_common.ASYMPT);
    IRET = 0;
L60:
    {
        double RCC = FLOAT_common.RC;
        if (FLOAT_common.RC == INTRNL.UNDEF) INTRNL.WASSET[9] = 1;
        if (FLOAT_common.RC == INTRNL.UNDEF) RCC = 1;
        if (INTGER.IZP * INTGER.IZT == 0) goto L70;
        if (FLOAT_common.RC > 0) goto L70;
        std::printf("0**** RC IS INVALID:%15.5G\n", FLOAT_common.RC);
        IRET = 0;
    L70:
        FLOAT_common.RC = RCC;
    }
    if (FLOAT_common.VSI == 0) goto L80;
    if (LNKBLK.LNKADR[2][3] != 0) goto L80;
    if (FLOAT_common.RSI > 0 && FLOAT_common.ASI > 0) goto L80;
    std::printf("0**** RSI OR ASI IS INVALID:%15.5G%15.5G\n", FLOAT_common.RSI, FLOAT_common.ASI);
    IRET = 0;
L80:
    if (!SOISW) goto L90;
    if (LNKBLK.LNKADR[4][3] != 0) goto L90;
    if (FLOAT_common.RSOI > 0 && FLOAT_common.ASOI > 0) goto L90;
    std::printf("0**** RSOI OR ASOI HAS INVALID VALUE:%15.5G%15.5G\n",
        FLOAT_common.RSOI, FLOAT_common.ASOI);
//
L90:
    if (!TENSW) goto L100;
    if ((int)JBLOCK.JSP == 2) goto L93;
    std::printf("0**** SP MUST BE 1 FOR TENSOR POTENTIALS, SP =%8d/2\n"
        " **** TENSOR POTENTIALS ARE SET TO 0\n", (int)JBLOCK.JSP);
    for (int I = 1; I <= 6; I++) {
        VTEN()[I] = 0;
    } // 92
    TENSW = FALSE_F;
    WAVCOM.TCSWS[INTRNL.ICHANW] = FALSE_F;
    WAVCOM.SOSWS[INTRNL.ICHANW] = WAVCOM.SOSWT;
L93:
    for (int I = 1; I <= 6; I++) {
        if (VTEN()[I] == 0) goto L95;
        if (LNKBLK.LNKADR[I + 4][3] != 0) goto L95;
        if (RTEN()[I] > 0 && ATEN()[I] > 0) goto L95;
        std::printf("0**** R%.3s OR A%.3s HAS INVALID VALUE:%15.5G%15.5G\n",
            TENNAM[I - 1], TENNAM[I - 1], RTEN()[I], ATEN()[I]);
        IRET = 0;
    L95:;
    } // 95
L100:
    if (IRET == 0) return;
//
    WAVCOM.NBAKCM = INTGER.NBACK;
//
//
//  SOME DERIVED DATA
//
//  K VALUE
//
    UK = std::sqrt(2 * FLOAT_common.AM * FLOAT_common.E) / CNSTNT.HBARC;
//
    ALAMBD = 2 * CNSTNT.PI / UK;
//
//  COULOMB PARAMETER AND SIMPLE TERMS USING IT
//     INITIALIZE THE COULOMB POTENTIAL
//
    ETA = (INTGER.IZP * INTGER.IZT) * std::sqrt(FLOAT_common.AM / (2.0 * FLOAT_common.E)) / CNSTNT.AFINE;
    // SETVSQ is now called AFTER MAKPOT (below) where RCT/RCP have valid values.
    // Calling it here with potentially UNDEF values (reset by CLRCHN) would corrupt the arrays.
//
//     COUPLED CHANNELS REQUIRE THE SAME STEPSIZE IN ALL CHANNELS
//
    if (SWITCH.PROBLM < 23 || INTRNL.ICHANW == 1) goto L110;
    RMAX = WAVCOM.ASYMPS[1];
    FLOAT_common.STEPSZ = WAVCOM.RSTEPS[1];
    goto L120;
//
//
L110:
    if (FLOAT_common.STPSPR != INTRNL.UNDEF)
        FLOAT_common.STEPSZ = DMIN1(ALAMBD, 1.0) / FLOAT_common.STPSPR;
//
    // Ensure UNDEF has the correct value (may have been corrupted by raw array writes)
    {
        uint64_t undef_pattern = 0xF0F0F0F0F0F0F0F0ull;
        double undef_val;
        std::memcpy(&undef_val, &undef_pattern, sizeof(double));
        if (INTRNL.UNDEF != undef_val) INTRNL.UNDEF = undef_val;
    }
    RMAX = FLOAT_common.ASYMPT;
//
//     ADJUST ASYMPTOPIA IF THERE WILL BE VERY LARGE L'S.
//     MATCHING MUST BE AT OR OUTSIDE THE CLASSICAL TURNING POINT.
//     HOWEVER, FITTER DOES NOT ALLOW VARIABLE ASYMPTOPIA RIGHT NOW
//     ALSO IF SCATASYM WAS READ IN, WE DO NOT OVERRIDE IT
//
    if (WAVCOM.NUMFIT == 0 && INTGER.LMAX != NOTDEF_INT && FLOAT_common.SCTASY < 0)
        RMAX = DMAX1(RMAX, (ETA + std::sqrt(ETA * ETA + (double)INTGER.LMAX * (INTGER.LMAX + 1))) / UK);
//
//     NOTE - NSTEP IS COMPUTED HERE SEPARATLY FOR EACH CHANNEL
//     AND IS STORED IN /WAVCOM/.  IT NEED NOT BE RELATED TO THE
//     "NSTEPS" USED IN THE BOUND STATE CALCULATION.  HOWEVER, IF
//     ASYMPT AND STEPSZ  ARE THE SAME FOR THE ENTIRE JOB, THEN
//        NSTEP = NSTEPS-1.
//
L120:
    NSTEP = (int)(RMAX / FLOAT_common.STEPSZ + 0.5);
    RMAX = NSTEP * FLOAT_common.STEPSZ;
    H = UK * FLOAT_common.STEPSZ;
//
//     TEST FOR IDENTICAL PARTICLES
//
    LSKIP = 1;
    ISTAT = 3;
    if (INTGER.IZP != INTGER.IZT || FLOAT_common.AMP != FLOAT_common.AMT || JBLOCK.JSP != JBLOCK.JST
        || FLOAT_common.EXSPT[1] != FLOAT_common.EXSPT[2]) goto L150;
    if (JBLOCK.JSP != INTRNL.NOTDEF) goto L135;
    std::printf("0*** WARNING:  THE PARTICLES APPEAR TO BE IDENTICAL"
        " BUT THEIR SPINS HAVE NOT BEEN DEFINED;\n"
        "%15s%s\n", "", "NON-IDENTICAL PARTICLE SCATTERING WILL BE ASSUMED.");
    goto L150;
L135:
    if ((int)JBLOCK.JSP == 0) LSKIP = 2;
    ISTAT = ((int)JBLOCK.JSP & 1) + 1;
//
//
//     SAVE THINGS FOR WAVPOT AND WAVELJ
//
L150:
    WAVCOM.HS[INTRNL.ICHANW] = H;
    WAVCOM.V0RS[INTRNL.ICHANW] = TVR_l;
    WAVCOM.V0IS[INTRNL.ICHANW] = TVI_l;
    WAVCOM.V0SORS[INTRNL.ICHANW] = TVSOR_l;
    WAVCOM.V0SOIS[INTRNL.ICHANW] = TVSOI_l;
    WAVCOM.V0SIS[INTRNL.ICHANW] = FLOAT_common.VSI;
    KANDM.RSCTS[INTRNL.ICHANW] = FLOAT_common.R;
    KANDM.ASCTS[INTRNL.ICHANW] = TAR_l;
    KANDM.RCSCTS[INTRNL.ICHANW] = FLOAT_common.RC;
    KANDM.RCSCTP[INTRNL.ICHANW] = FLOAT_common.RCP;
    KANDM.RCSCTT[INTRNL.ICHANW] = FLOAT_common.RCT;
    WAVCOM.RIS[INTRNL.ICHANW] = FLOAT_common.RI;
    WAVCOM.AIS[INTRNL.ICHANW] = TAI_l;
    WAVCOM.RSOS[INTRNL.ICHANW] = FLOAT_common.RSO;
    WAVCOM.ASOS[INTRNL.ICHANW] = FLOAT_common.ASO;
    WAVCOM.RSOIS[INTRNL.ICHANW] = FLOAT_common.RSOI;
    WAVCOM.ASOIS[INTRNL.ICHANW] = FLOAT_common.ASOI;
    WAVCOM.RSIS[INTRNL.ICHANW] = FLOAT_common.RSI;
    WAVCOM.ASIS[INTRNL.ICHANW] = FLOAT_common.ASI;
    for (int I = 1; I <= 6; I++) {
        WAVCOM.VTENS[INTRNL.ICHANW][I] = VTEN()[I-1];  // VTEN() is 0-based; Fortran VTEN(I) = C++ VTEN()[I-1]
        WAVCOM.RTENS[INTRNL.ICHANW][I] = RTEN()[I-1];
        WAVCOM.ATENS[INTRNL.ICHANW][I] = ATEN()[I-1];
    } // 159
    KANDM.ETAS[INTRNL.ICHANW] = ETA;
    WAVCOM.IZPS[INTRNL.ICHANW] = INTGER.IZP;
    WAVCOM.IZTS[INTRNL.ICHANW] = INTGER.IZT;
    WAVCOM.ASYMPS[INTRNL.ICHANW] = RMAX;
    WAVCOM.ES[INTRNL.ICHANW] = FLOAT_common.E;
    WAVCOM.RSTEPS[INTRNL.ICHANW] = FLOAT_common.STEPSZ;
    WAVCOM.NSTPSS[INTRNL.ICHANW] = NSTEP;
    WAVCOM.NSTP2S[INTRNL.ICHANW] = NSTEP;
    AKS()[INTRNL.ICHANW] = UK;
    REDMS()[INTRNL.ICHANW] = FLOAT_common.AM;
    KANDM.RATMSS[INTRNL.ICHANW] = INTRNL.RATMAS;
    WAVCOM.JSPS[INTRNL.ICHANW] = (int)JBLOCK.JSP;
    WAVCOM.STEPR = FLOAT_common.STEP1R;
    WAVCOM.STEPI = FLOAT_common.STEP1I;
    WAVCOM.NUMAIT = INTGER.NAITKN;
    FTIME.NUMLOP = 0;
    FTIME.NUMWAV = 0;
    WAVCOM.LSKIPS[INTRNL.ICHANW] = LSKIP;
    WAVCOM.ISTATS[INTRNL.ICHANW] = ISTAT;
    WAVCOM.POWRLS[INTRNL.ICHANW] = FLOAT_common.POWRL;
    WAVCOM.POWIMS[INTRNL.ICHANW] = FLOAT_common.POWIM;
    WAVCOM.IWAVRS[INTRNL.ICHANW] = 0;
    WAVCOM.IWAVIS[INTRNL.ICHANW] = 0;
//
//     QUANTITIES FOR INTEGRATION PAST ASYMPTOPIA
//
    WAVCOM.XFACS[INTRNL.ICHANW][1] = 1 + H * H / 12;
    WAVCOM.XFACS[INTRNL.ICHANW][2] = -(H * H / 12) * INTGER.IZP * INTGER.IZT * CNSTNT.HBARC / (FLOAT_common.E * CNSTNT.AFINE);
    WAVCOM.XFACS[INTRNL.ICHANW][3] = -(H * H / 12) / (UK * UK);
    FTIME.TIMES[1] = 0;
    FTIME.TIMES[2] = 0;
    FTIME.TIMES[3] = 0;
//
//     INITIALIZE LINKULES IF NEEDED
//
    for (int I = 1; I <= LNKBLK.NUMLNK; I++) {
        if (LNKBLK.LNKADR[I][3] == 0) goto L259;
//
//     GENERATE THE SPECIAL UNIQUE NAME
//
        LNKBLK.IUNIQU = LNKBLK.IUNIQU + 1;
        std::snprintf(linkfmtout, sizeof(linkfmtout), "*%03d", LNKBLK.IUNIQU);
        std::memcpy(LINKID, linkfmtout, 4);
//
//     MAKE INITIALIZING CALL
//
        LINKUL(LNKBLK.LNKADR[I][3], *(char8*)&LNKBLK.LNKADR[I][1], &LNKBLK.LNKADR[I][5], I, 1,
            IRET, INTGER.L, JBLOCK.J, 0.0, FLOAT_common.STEPSZ, NSTEP + 1, &DUMMY, &DUMMY, LINKID);
        if (IRET >= 0) goto L240;
        IRET = 0;
        return;
    L240:
        LNKBLK.LNKADR[I][4] = IRET;
    L259:;
    } // 259
//
    if (IPRNT < 1) goto L600;
//
//     PRINT OUT INPUT NICELY
//
    std::printf("\n0        OPTICAL MODEL SCATTERING FOR THE %.8s CHANNEL\n",
        CHANAM[INTRNL.ICHANW - 1]);
    std::printf("0E LAB =%9.3f MEV,    E CM =%9.3f MEV,     K =", EELAB, FLOAT_common.E);
    print_G(12, 5, UK);
    std::printf("     WAVELENGTH =%8.4f FM\n"
        " PROJECTILE MASS =%7.2f AMU,%5sTARGET MASS =%7.2f AMU,%5sREDUCED MASS =%10.2f MEV/C**2\n",
        ALAMBD, FLOAT_common.AMP, "", FLOAT_common.AMT, "", FLOAT_common.AM);
    if (ISTAT == 3) std::printf(" THIS IS NON-IDENTICAL PARTICLE SCATTERING\n");
    if (ISTAT != 3) std::printf(" THIS IS SCATTERING OF IDENTICAL %.8s\n", STATWD[ISTAT - 1]);
    if (JBLOCK.JSP != INTRNL.NOTDEF && JBLOCK.JST != INTRNL.NOTDEF)
        std::printf(" PROJECTILE SPIN =%3d/2%5sTARGET SPIN =%3d/2\n", (int)JBLOCK.JSP, "", (int)JBLOCK.JST);
    std::printf(" Z PROJECTILE =%4d%5sZ TARGET =%4d\n"
        "0POTENTIAL         COUPLING CONS.    RADIUS    DIFFUSENESS    RADIUS PARAMETER\n",
        INTGER.IZP, "", INTGER.IZT);
//
    IRET = 1;
    {
        double R0M = INTRNL.R0MASS;
        if (R0M == INTRNL.UNDEF || R0M == 0.0) R0M = 1.0;
        PARPRT(1,  "REAL CENTRAL       ", TVR_l,              FLOAT_common.R,    TAR_l,             FLOAT_common.R    / R0M, FLOAT_common.POWRL, "   POWER", IRET);
        PARPRT(2,  "VOLUME ABSORPTION  ", TVI_l,              FLOAT_common.RI,   TAI_l,             FLOAT_common.RI   / R0M, FLOAT_common.POWIM, "   POWER", IRET);
        PARPRT(13, "SURFACE ABSORPTION ", FLOAT_common.VSI,   FLOAT_common.RSI,  FLOAT_common.ASI,  FLOAT_common.RSI  / R0M, 0.0,                "        ", IRET);
        PARPRT(3,  "REAL SPIN-ORBIT    ", FLOAT_common.VSO,   FLOAT_common.RSO,  FLOAT_common.ASO,  FLOAT_common.RSO  / R0M, FLOAT_common.TAU,   "     TAU", IRET);
        PARPRT(4,  "IMAG. SPIN-ORBIT   ", FLOAT_common.VSOI,  FLOAT_common.RSOI, FLOAT_common.ASOI, FLOAT_common.RSOI / R0M, FLOAT_common.TAUI,  "    TAUI", IRET);
        PARPRT(7,  "REAL TR            ", FLOAT_common.VTR,   FLOAT_common.RTR,  FLOAT_common.ATR,  FLOAT_common.RTR  / R0M, 0.0,                "        ", IRET);
        PARPRT(8,  "IMAGINARY TR       ", FLOAT_common.VTRI,  FLOAT_common.RTRI, FLOAT_common.ATRI, FLOAT_common.RTRI / R0M, 0.0,                "        ", IRET);
        PARPRT(9,  "REAL TL            ", FLOAT_common.VTL,   FLOAT_common.RTL,  FLOAT_common.ATL,  FLOAT_common.RTL  / R0M, 0.0,                "        ", IRET);
        PARPRT(10, "IMAG. TL           ", FLOAT_common.VTLI,  FLOAT_common.RTLI, FLOAT_common.ATLI, FLOAT_common.RTLI / R0M, 0.0,                "        ", IRET);
        PARPRT(11, "REAL TPR           ", FLOAT_common.VTP,   FLOAT_common.RTP,  FLOAT_common.ATP,  FLOAT_common.RTP  / R0M, 0.0,                "        ", IRET);
        PARPRT(12, "IMAG. TPR          ", FLOAT_common.VTPI,  FLOAT_common.RTPI, FLOAT_common.ATPI, FLOAT_common.RTPI / R0M, 0.0,                "        ", IRET);
        PARPRT(5,  "PT&SPHERE  COULOMB ", ETA,                FLOAT_common.RC,   0.0,               FLOAT_common.RC   / R0M, 0.0,                "        ", IRET);
    }

    if ((FLOAT_common.RCP == 0) || (FLOAT_common.RCT == 0) ||
        (FLOAT_common.RCP == INTRNL.UNDEF) || (FLOAT_common.RCT == INTRNL.UNDEF)) goto L500;
    TEMP2_l = FLOAT_common.RCP / std::pow(FLOAT_common.AMP, 0.333333333333333333333330);
    TEMP3_l = FLOAT_common.RCT / std::pow(FLOAT_common.AMT, 0.333333333333333333333330);
    std::printf(" FOLDED COULOMB POTENTIALS - RCP = %7.4f  RCT = %7.4f  RC0P = %7.4f  RC0T = %7.4f\n",
        FLOAT_common.RCP, FLOAT_common.RCT, TEMP2_l, TEMP3_l);
//
L500:
    TEMP = DMIN1(ALAMBD, 1.0) / FLOAT_common.STEPSZ;
    std::printf("0ASYMPTOPIA = %8.3f FM\n"
        " STEP SIZE =  %8.3f FM     %8.3f STEPS PER \"WAVELENGTH\"\n\n",
        RMAX, FLOAT_common.STEPSZ, TEMP);
//
//     WILL THE WAVE FUNCTIONS BE FOUND BY A LINKULE
//
    if (LNKBLK.LNKADR[6][3] == 0) goto L600;
    std::printf("0 THE SCATTERING WAVEFUNCTION IS BEING"
        " COMPUTED BY THE %.8s LINKULE:\n", (char*)&LNKBLK.LNKADR[6][1]);
    {
        int I = 6;
//
//     CALL LINKULE TO PRINT STUFF ABOUT ITSELF
//
        LINKUL(LNKBLK.LNKADR[I][3], *(char8*)&LNKBLK.LNKADR[I][1], &LNKBLK.LNKADR[I][5],
            I, 2, IRTN,
            INTGER.L, JBLOCK.J, 0.0, FLOAT_common.STEPSZ, NSTEP + 1, &DUMMY, &DUMMY, (char*)&INTRNL.ICHANW);
        std::printf(" \n");
        IRET = MIN0(IRTN, IRET);
        if (IRET >= 0) goto L600;
    }
    // L570:
    IRET = 0;
    return;
//
//
//
//     SET ASIDE AREAS FOR THE POTENTIAL AND FOR V' USED BY WAVELJ
//     THESE ARE STORED ONLY TO ASYMPTOPIA
//
//     THIS IS DONE AS FOLLOWS:
//
//     1) FOR STAND-ALONE:  ON EVERY CALL TO WAVSET.
//     2) FOR ELASTIC FITS: MAXIMUM NUMBER OF STEPS IS ACCUMULATED IN
//        NSTPSS(2); AS THE SIZE GROWS, THEY ARE ALLOCATED.
//     3) FOR DWBA: POTENTIAL AREAS ARE ALLOCATED AT BOTH CALLS;
//        WAVEFUNCTION AREAS ARE ALLOCATED ONLY ON 2ND CALL.
//
L600:
    NSTPP1 = NSTEP + 1;
//
    if (WAVCOM.NUMFIT == 0) goto L650;
    if (WAVCOM.NSTPSS[2] >= NSTEP) goto L680;
    WAVCOM.NSTPSS[2] = NSTEP;
//
L650:
    WAVCOM.IVWORK = NALLOC("VWORK   ", NSTPP1);
    WAVCOM.IRLVS[INTRNL.ICHANW] = NALLOC(NAMES[INTRNL.ICHANW - 1][0], NSTPP1);
    WAVCOM.IIMVS[INTRNL.ICHANW] = NALLOC(NAMES[INTRNL.ICHANW - 1][1], NSTPP1 + 1);
    WAVCOM.ICENTR[INTRNL.ICHANW] = NALLOC(NAMES[INTRNL.ICHANW - 1][4], NSTPP1);
    if (INTRNL.ICHANW == 1) {
        extern double g_watch_val1, g_watch_val2;
    }
    WAVCOM.ISORS[INTRNL.ICHANW] = 0;
    WAVCOM.ISOIS[INTRNL.ICHANW] = 0;
    if (WAVCOM.V0SORS[INTRNL.ICHANW] != 0)
        WAVCOM.ISORS[INTRNL.ICHANW] = NALLOC(NAMES[INTRNL.ICHANW - 1][2], NSTPP1);
    if (WAVCOM.V0SOIS[INTRNL.ICHANW] != 0)
        WAVCOM.ISOIS[INTRNL.ICHANW] = NALLOC(NAMES[INTRNL.ICHANW - 1][3], NSTPP1);
    for (int II = 1; II <= 6; II++) {
        WAVCOM.ITENS[INTRNL.ICHANW][II] = 0;
        if (WAVCOM.VTENS[INTRNL.ICHANW][II] != 0)
            WAVCOM.ITENS[INTRNL.ICHANW][II] = NALLOC(NAMES[INTRNL.ICHANW - 1][II + 4], NSTPP1);
    } // 659
//
//     AND COMPUTE THE POTENTIAL
//     LX IS SET TO 0 FOR POSSIBLE USE BY DEFORMED POTENTIALS
//
L680:
    LXSAVE = INTGER.LX;
    INTGER.LX = -1;
    // Initialize VCSQ12 Coulomb coefficients for this channel BEFORE MAKPOT.
    // FLOAT_common.RC was computed by SETPOT. RCT/RCP may be UNDEF if SETFIT is not
    // fully translated, so fall back to RCT=RC (full Coulomb radius), RCP=0 (point projectile).
    {
        double RCT_val = FLOAT_common.RCT;
        double RCP_val = FLOAT_common.RCP;
        if (!std::isfinite(RCT_val) || RCT_val <= 0.0 ||
            !std::isfinite(RCP_val) || RCP_val < 0.0) {
            RCT_val = FLOAT_common.RC;
            RCP_val = 0.0;
        }
        SETVSQ(RCT_val, RCP_val, INTGER.IZT, INTGER.IZP, INTRNL.ICHANW);
        // Also fix up RCSCTP/RCSCTT which were saved above with UNDEF values
        if (!std::isfinite(KANDM.RCSCTP[INTRNL.ICHANW]) || KANDM.RCSCTP[INTRNL.ICHANW] < 0.0)
            KANDM.RCSCTP[INTRNL.ICHANW] = RCP_val;
        if (!std::isfinite(KANDM.RCSCTT[INTRNL.ICHANW]) || KANDM.RCSCTT[INTRNL.ICHANW] <= 0.0)
            KANDM.RCSCTT[INTRNL.ICHANW] = RCT_val;
    }
    MAKPOT(INTRNL.ICHANW, IRET);
    INTGER.LX = LXSAVE;
//
//     FOR DWBA MUST SAVE TWO SETS OF LINKULE INDICATORS
//
    if (WAVCOM.NUMFIT != 0) goto L690;
    for (int I = 1; I <= LNKBLK.NUMLNK; I++) {
        for (int II = 1; II <= 6; II++) {
            LOCFIT.linkule.LNKAD2[I][INTRNL.ICHANW][II] = LNKBLK.LNKADR[I][II];
        }
    } // 685
//
//     WORK AREA NOT NEEDED AFTER MAKPOT
//
L690:
    if (WAVCOM.NUMFIT == 0) LOCPTRS.Z[WAVCOM.IVWORK] = -LOCPTRS.Z[WAVCOM.IVWORK];
//
    if (IRET >= 0) goto L700;
    IRET = 0;
    return;
//
//     DETERMINE LCRITICAL.
//     IF BOTH LMIN AND LMAX ARE DEFINED WE BYPASS THIS AND JUST USE
//     THE AVERAGE.
//
L700:
    if (INTGER.LMIN == NOTDEF_INT || INTGER.LMAX == NOTDEF_INT) goto L720;
    KANDM.LCRITS[INTRNL.ICHANW] = (INTGER.LMIN + INTGER.LMAX) / 2;
    goto L750;
L720:
    LCRITL(FLOAT_common.E, UK, ETA, FLOAT_common.RC, NSTEP + 1,
        0.0, FLOAT_common.STEPSZ, WAVCOM.IRLVS[INTRNL.ICHANW], WAVCOM.IIMVS[INTRNL.ICHANW],
        48, (INTGER.IPRINT / 10) % 10,
        LC1, LC2, KANDM.LCRITS[INTRNL.ICHANW]);
    if (IPRNT != 0)
        std::printf(" ESTIMATED CRITICAL L:%4d%5sBY |S|:%4d%5sBY DEFLECTION FUNC.:%4d\n",
            KANDM.LCRITS[INTRNL.ICHANW], "", LC1, "", LC2);
//
//     IF LC1=0  THEN |S|=1 ALWAYS AND WE MUST FORCE LMIN=0.
//
    if (LC1 == 0 && INTGER.LMIN == NOTDEF_INT) INTGER.LMIN = 0;
//
//     ON 2ND CALL ALLOCATE AREAS FOR THE WAVE FUNCTIONS USING LARGEST
//     NSTEP.
//
L750:
    if (!WAVCOM.STANSW && INTRNL.ICHANW == 1) goto L800;
    if (!WAVCOM.STANSW) NSTPP1 = MAX0(WAVCOM.NSTPSS[1] + 1, NSTPP1);
    {
        int I = NSTPP1 + 6;
        WAVCOM.NWAVEF = I;
//
//     WE NEED ROOM FOR 4 WAVE FUNCTIONS IF TENSOR COUPLING IS PRESENT.
        if (WAVCOM.TCSWS[INTRNL.ICHANW] || (WAVCOM.TCSWS[3 - INTRNL.ICHANW] && !WAVCOM.STANSW))
            I = 4 * I;
        WAVCOM.IWAVR = NALLOC("WAVER   ", I);
        WAVCOM.IWAVI = NALLOC("WAVEI   ", I);
    }
//
L800:
    IRET = 1;
    return;
//
//
//
}
