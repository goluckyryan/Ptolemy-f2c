// ingrst_translated.cpp — INGRST subroutine
// Translated from source.f lines 21287-22124
// Initializes the grid for the inelastic scattering radial integral.
//
// This routine is called by MAIN. It sets up the one-dimensional grid and
// potential values for INRDIN. For DWBA it allocates Gauss quadrature points
// and weights and evaluates the nuclear and Coulomb potential form factors on
// the grid. For coupled-channels it builds nuclear and Coulomb form-factor
// radial arrays indexed by radial form number.

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include "ptolemy_alloc.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_forward.h"
#include "ptolemy_io.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <algorithm>

void INGRST(int& IRET) {
    //
    // 1/77 - FIRST VERSION - S.P.
    // 11/22/77 - NEW SCTASYOPIA - S.P.
    // 12/4/77 - LINKULES - S.P.
    // 5/7/78 - STUFF TO ALLOW 2 COUPLED STATES - S.P.
    // 5/23/78 - FIX BUG IN 11/22/77 - ALWAYS DEFINE VWORK - S.P.
    //           BIGGER WAVR, WAVI ARRAYS FOR NEW WAVELJ - S.P.
    // 3/12/79 - CLEAN UP DEF OF SUMMAX, NSTEPS; DEFINE MORE FOR LINKULES
    // 12/20/79 - INELASTIC NAMES INPLACE OF TRANSFER NAMES - S.P.
    //            START OF FULL COUPLED CHANNELS CHANGES - S.P.
    //  1/2/80 PUT IN FULL COUPLED CHANNELS FOR ROTATIONAL MODEL
    //         UP TO 2 ORDER M.R.B.
    // 3/18/80 - FIX BUG IN FREEING ITEMP - S.P.
    // 4/27/80 - GET NPTMIN, FREE UNNEEDED ARRAYS FOR C.C. - S.P.
    // 4/28/80 - REASSIGN LRPTS AFTER GIVALL - S.P.
    // 6/24/80 - SAVE RC IN RCUEFF - S.P.
    // 7/9/80 - FIX BUG IN FREEING COULH, DON'T OVERFLOW AREA
    //          RETURNED BY GIVEAL - S.P.
    // 8/25/80 - FIX TYPING ERROR AT STMNT 820 - S.P.
    // 12/9/80 - MAKPOT NOLONGER COMPUTES COULOMB; DON'T
    //           DEFINE COULOMB COUPLING FOR BETAC=0 - S.P.
    // 12/30/80 - FIX ERRORS IN 2ND ORDER NUCLEAR - S.P.
    // 1/11/81 - ALWAYS 1/R**(LX+1); NO MORE COUPL2 - S.P.
    // 7/13/81 - FIX BUG FROM 1/11/81 (JUMPSZ) - S.P.
    // 7/16/81 - SAVE LNKBLK TO ALLOW FORM-FACTOR CHANGES - S.P.
    // 7/29/81 - FIRST CRUDE STUFF FOR DEFORMED POTENTIALS - S.P.
    // 8/28/81 - MORE FOR DEFORMED POTENTIALS - S.P.
    // 9/21/81 - DON'T LOOK WHEN THERE ARE NO POTENTIALS - S.P.
    // 10/16/81 - DONT' CONDISDER LXTOT FOR MUTUAL - S.P.
    // 12/14/18 - LET GETSCT GET WAVEFUNC ALSO FOR C.C. - S.P.
    // 4/14/82 - FORMASIS ALSO WORKS FOR DWBA - S.P.
    // 6/10/83 - FIX BUGS AT STMNT 235, 408 - S.P.
    // 7/24/84 - HOTFUDGE - S.P.
    // 10/16/84 - USE ABS(SCTASY) - S.P.
    // 12/19/84 - USE PARPRT; DEFORMED SI - S.P.
    // 4/29/86 - FIX BUG FOR IMAG+DEFORMED+HOTFUDGE - S.P.
    // 11/16/01 - initialize IFIRST, ISECND at very start
    //

    // Convenient references to COMMON blocks
    double& UNDEF   = INTRNL.UNDEF;
    double& NOTDEF  = INTRNL.NOTDEF;
    int&    ICHANW  = INTRNL.ICHANW;
    int&    ISTRIP  = INTRNL.ISTRIP;
    int&    IHSAVE  = INTRNL.IHSAVE;

    int&    IPRINT  = INTGER.IPRINT;
    int&    LX      = INTGER.LX;
    int&    NPSUM   = INTGER.NPSUM;
    int&    MAPSUM  = INTGER.MAPSUM;
    int&    IZP     = INTGER.IZP;
    int&    IZT     = INTGER.IZT;

    int&    L_int   = INTGER.L;       // INTGER.L
    double& J_dbl   = JBLOCK.J;       // JBLOCK.J

    int&    PROBLM  = SWITCH.PROBLM;
    int&    KFRMTP  = SWITCH.KFRMTP;

    double& STEPSZ  = FLOAT_common.STEPSZ;
    double& SUMMAX  = FLOAT_common.SUMMAX;
    double& SUMMID  = FLOAT_common.SUMMID;
    double& SUMMIN  = FLOAT_common.SUMMIN;
    double& GAMSUM  = FLOAT_common.GAMSUM;
    double& SUMPTS  = FLOAT_common.SUMPTS;
    double& SCTASY  = FLOAT_common.SCTASY;
    double& R0MASS  = INTRNL.R0MASS;
    double& R       = FLOAT_common.R;
    double& RC      = FLOAT_common.RC;
    double& RI      = FLOAT_common.RI;
    double& RSO     = FLOAT_common.RSO;
    double& RSOI    = FLOAT_common.RSOI;
    double& RSI     = FLOAT_common.RSI;
    double& TVR     = TEMPVS.TVR;
    double& TVI     = TEMPVS.TVI;
    double& TAR     = TEMPVS.TAR;
    double& TAI     = TEMPVS.TAI;
    double& E       = FLOAT_common.E;
    double& VSI     = FLOAT_common.VSI;
    double& VSO     = FLOAT_common.VSO;
    double& VSOI    = FLOAT_common.VSOI;
    double& POWRL   = FLOAT_common.POWRL;
    double& POWIM   = FLOAT_common.POWIM;
    double& TAU     = FLOAT_common.TAU;
    double& TAUI    = FLOAT_common.TAUI;
    double& RCP     = FLOAT_common.RCP;
    double& RCT     = FLOAT_common.RCT;
    double& RC0P    = FLOAT_common.RC0P;
    double& RC0T    = FLOAT_common.RC0T;
    double& AMA     = FLOAT_common.AMA;
    double& AMBIGA  = FLOAT_common.AMBIGA;
    double& AM      = FLOAT_common.AM;
    double& AMP     = FLOAT_common.AMP;
    double& AMT     = FLOAT_common.AMT;

    double& REDMI   = KANDM.REDMI;
    // RCUEFF = GRIDCM.JACOB  (INGRST view of /GRIDCM/ maps RCUEFF to position 1 = JACOB)
    double& RCUEFF  = GRIDCM.JACOB;

    double& HBARC   = CNSTNT.HBARC;
    double& PI      = CNSTNT.PI;

    int&    LXMIN   = INELCM.LXMIN;
    // INGRST view of /GRIDCM/ (RCUEFF,GRDCM1(2),IBINDX,...,NPTMIN,NOFLO,NUMPT,NCRIT,...):
    //   NPTMIN = 17th integer after the 3 doubles = GRIDCM.IWIO
    //   NOFLO  = 18th integer                    = GRIDCM.NOFLO  (same name, same position)
    //   NUMPT  = 19th integer                    = GRIDCM.NRIROI
    //   NCRIT  = 20th integer                    = GRIDCM.NCRIT  (same name, same position)
    //   INUCH  = 11th integer (IRI in GRDSET view) = GRIDCM.IRI
    //   ICOULH = 12th integer (IRO in GRDSET view) = GRIDCM.IRO
    //   IRPTS4 = 26th integer = GRIDCM.ISMHPT
    //   IRPTS  = 27th integer = GRIDCM.ISMIPT
    //   IRWTS  = GRIDCM.ISMIVL  (after GRDCM2/RVRLIM double and two more ints)
    int&    NPTMIN  = GRIDCM.IWIO;
    int&    NUMPT   = GRIDCM.NRIROI;
    int&    NOFLO   = GRIDCM.NOFLO;

    // WAVCOM aliases
    double* HS      = WAVCOM.HS;       // 1-based
    double* RSTEPS  = WAVCOM.RSTEPS;   // 1-based
    int*    NSTPSS  = WAVCOM.NSTPSS;   // 1-based
    int*    LASTNF  = WAVCOM.LASTNF;   // 1-based
    double* ASYMPS  = WAVCOM.ASYMPS;   // 1-based
    double* ES      = WAVCOM.ES;       // 1-based
    int*    IRLVS   = WAVCOM.IRLVS;    // 1-based
    int*    IIMVS   = WAVCOM.IIMVS;    // 1-based
    int*    ICENTR  = WAVCOM.ICENTR;   // 1-based
    int*    NSTP2S  = WAVCOM.NSTP2S;   // 1-based
    int*    IZS     = INTGER.IZS;      // 1-based

    int&    IVWORK  = WAVCOM.IVWORK;
    int&    IWAVR   = WAVCOM.IWAVR;
    int&    IWAVI   = WAVCOM.IWAVI;

    // GRIDCM allocator slot indices — mapped from INGRST view to GRDSET/INELDC struct names:
    //   IRPTS4 (26th int) = GRIDCM.ISMHPT, IRPTS (27th int) = GRIDCM.ISMIPT
    //   INUCH (11th int)  = GRIDCM.IRI,    ICOULH (12th int) = GRIDCM.IRO
    //   IRWTS = GRIDCM.ISMIVL (after GRDCM2 double + 2 ints)
    int&    IRPTS   = GRIDCM.ISMIPT;
    int&    IRWTS   = GRIDCM.ISMIVL;
    int&    IRPTS4  = GRIDCM.ISMHPT;
    int&    INUCH   = GRIDCM.IRI;
    int&    ICOULH  = GRIDCM.IRO;

    // CCBLK
    int&    IBASCP  = CCBLK.IBASCP;

    // LNKBLK
    int&    IUNIQU  = LNKBLK.IUNIQU;

    int&    FACFR4  = ALLOCS.FACFR4;

    // --- Local variables ---
    int    IFIRST, ISECND;      // nuclear radial form indices for 1st/2nd deriv
    int    IPRT;
    int    DBUGSW;              // LOGICAL
    int    CCSW;                // LOGICAL (PROBLM == 24)
    int    DRIVSW;              // LOGICAL (KFRMTP == 1)
    int    SOSW;                // LOGICAL
    int    SOISW;               // LOGICAL

    int    I, II, K, KK, IP, IPP, IB1;
    int    NSTEP, NSTEP1, NSTEPS;
    int    ICUBIC;
    int    NNUCFF, NCOUFF;
    int    NMPTSV;
    int    ICHANW_save;

    // Pointers into ALLOC
    int    LRPTS, LRWTS, LRPTS4;
    int    LREALV, LIMAGV, LVWORK, LNUCH, LCOULH;
    int    LR, LB, LC, LD;
    int    LBASCP, NBASCP, MBASCP;
    int    ISCTMN_loc;

    // Temporaries for potential saving/restoring (HOTFUDGE block)
    double T1, T2, T3, T4, T5, T6;
    int    IT1, IT2;
    int    INDXP;
    int    MODL;

    // Printing temporaries
    double TEMP1_v, TEMP2_v, TEMP3_v, TEMP4_v, TEMP5_v, TEMP6_v;
    double ETA;

    // Loop variables
    int    RVAL_n, N_val;
    double RVAL, WT, X, XX;
    double TERM;
    double VC;

    // Coulomb loop
    int    LABEL, IORDER, IMULT, IPOINT, IMT;
    int    IPOW;
    double TMPFAC;
    int    I1;

    // For ASSIGN-goto emulation (IGOTO)
    // Fortran uses   ASSIGN 652 TO IGOTO / ASSIGN 654 TO IGOTO
    // We use an int flag instead.
    int    IGOTO_real;   // 0 = deriv (652/732), 1 = direct (654/734)
    int    IGOTO_imag;

    // =========================================================================
    // Initialise at start (11/16/01 bug fix)
    // =========================================================================
    IFIRST = 0;
    ISECND = 0;
    IRET   = 0;

    IPRT   = MOD(IPRINT, 10);
    DBUGSW = (IPRT >= 4);
    CCSW   = (PROBLM == 24);

    I = 1;
    if (CCSW) goto L150;

    //
    // FIND WHERE LMIN SCATTERING WAVE EXCEEDS 1E-15
    //
L110:
    NSTEPS   = NSTPSS[1];
    ISCTMN_loc = LOCPTRS.Z[KANDM.ISCTMN];
    for (I = 1; I <= NSTEPS; I++) {
        if (DABS(ALLOC(ISCTMN_loc + I)) > 1.0e-15) goto L140;
    }
    std::printf("0**** ERROR:  U(LMIN) < 1E-15 EVERYWHERE.\n");
    return;

    //
    // USE THIS AS THE LOWER POINT OF THE INTEGRAL, EXCEPT FOR C.C.
    //
L140:
    LOCPTRS.Z[KANDM.ISCTMN] = -LOCPTRS.Z[KANDM.ISCTMN];
    LOCPTRS.Z[KANDM.ISCTCR] = -LOCPTRS.Z[KANDM.ISCTCR];
    goto L160;

    //
    // FOR C.C. WE START AT NFIRST OF THE LOWEST L.
    //
L150:
    I = MAX0(1, LASTNF[1] - 5);

L160:
    if (SUMMIN == UNDEF) SUMMIN = RSTEPS[1] * (I - 1);
    NPTMIN = (int)(SUMMIN / RSTEPS[1] + 0.5);
    if (SUMMAX == UNDEF) SUMMAX = DABS(SCTASY);
    if (CCSW)            SUMMAX = ASYMPS[1];
    if (SUMMID == UNDEF) SUMMID = 0.5 * (SUMMIN + SUMMAX);

    //
    // NUMBER OF GAUSS POINTS IS DETERMINED FROM SUMPTSPER,
    // EXCEPT FOR C.C. IT IS NUMBER OF INCOMING STEPS
    //
    NUMPT  = (int)((SUMMAX - SUMMIN) * (SUMPTS * (KANDM.AKI + KANDM.AKO) / (4.0 * PI)));
    NUMPT  = MAX0(NUMPT, NPSUM);
    NMPTSV = NUMPT;
    if (CCSW) NUMPT  = NSTPSS[1] + 1;
    if (CCSW) NMPTSV = NUMPT - NPTMIN;

    // RMAXS(1) = SUMMAX  (EQUIVALENCE: RMAXS(1)=RIMAX, RMAXS(2)=ROMAX)
    KANDM.RIMAX = SUMMAX;
    KANDM.ROMAX = SUMMAX;

    //
    // MAY HAVE TO REDEFINE WAVELJ WORK AREAS
    //
    I = (int)(SUMMAX / RSTEPS[1] + 0.5);
    NSTP2S[1] = MAX0(I, NSTP2S[1]);
    I = (int)(SUMMAX / RSTEPS[2] + 0.5);
    NSTP2S[2] = MAX0(I, NSTP2S[2]);
    I = MAX0(NSTP2S[1], NSTP2S[2]);
    if (I <= MAX0(NSTPSS[1], NSTPSS[2])) goto L200;
    IWAVR = NALLOC("WAVER   ", I + 6);
    IWAVI = NALLOC("WAVEI   ", I + 6);

    //
    // THE ARRAYS WILL BE:
    //   IRPTS, IRWTS - GAUSS POINTS AND WEIGHTS
    //   IRPTS4 - REAL*4 GAUSS POINTS FOR WAVELJ INTERPOLATOR
    //   INUCH - COMPLEX NUCLEAR POTENTIAL PART OF H WITH R'S AND GAUSS
    //           WEIGHT BUT NO BETA.
    //   ICOULH - COULOMB PART OF H WITH RC AND GAUSS WEIGHT BUT NO BETAS
    //
L200:
    if (CCSW) goto L210;
    IRPTS   = NALLOC("RPTS    ", NUMPT);
    IRWTS   = NALLOC("RWT   S ", NUMPT);
    IRPTS4  = NALLOC("RPTS4 ", (NUMPT + FACFR4 - 1) / FACFR4);
    INUCH   = NALLOC("NUCH    ", 2 * NUMPT);
    ICOULH  = NALLOC("COULH   ", NUMPT);
    IFIRST  = 1;
    goto L300;

    //
    // DEFINE NUCLEAR AND COULOMB H-ARRAYS (RADIAL PART)
    // FOR COUPLED CHANNELS.
    //
    // DETERMIN HOW MANY DIFFERENT RADIAL FORMS WILL BE NEEDED
    // AND SET UP THE POINTERS FOR THE NUCLEAR FORM FACTORS
    //
L210:
    LBASCP = FACFR4 * LOCPTRS.Z[IBASCP] - FACFR4 + 1;
    NBASCP = ILLOC(LBASCP);
    MBASCP = ILLOC(LBASCP + 1);
    LBASCP = ILLOC(LBASCP + 2) + LBASCP;
    NNUCFF = 0;
    IFIRST = 0;
    ISECND = 0;
    for (I = 1; I <= NBASCP; I++) {
        K = LBASCP + (I - 1) * MBASCP;
        if (ALLOC4(K + 4) != 0) goto L215;
        ILLOC(K + 9) = 0;
        goto L229;
L215:
        {
            int IPT = ILLOC(K + 7);
            if (ILLOC(K + 9) != (int)NOTDEF) goto L229;

            //
            // THIS IS A NEW RADIAL FORM
            //
            LX = ILLOC(K + 3);
            NNUCFF = NNUCFF + 1;
            if (IPT == 1) IFIRST = NNUCFF;
            if (IPT == 2) ISECND = NNUCFF;
            for (II = I; II <= NBASCP; II++) {
                KK = LBASCP + (II - 1) * MBASCP;
                if (ILLOC(KK + 7) != IPT) goto L219;

                //
                // FOR DEFORMED MUST ALSO CONSIDER LX EXCEPT FOR MUTUAL
                //
                {
                    int ITEMP = MOD(IPT, 100);
                    if ((ITEMP == 3 || ITEMP == 4) && LX != ILLOC(KK + 3)) goto L219;
                }
                if (ALLOC4(KK + 4) == 0) goto L219;
                ILLOC(KK + 9) = 2 * NNUCFF - 1;
L219:           ;
            }
        }
L229:   ;
    }

    //
    // NOW DETERMIN HOW MANY COULOMB COUPLING RADIAL FORMS
    // ARE NEEDED AND ASSIGN THE COULOMB FORM FACTOR POINTERS
    //
    NCOUFF = 0;
    for (I = 1; I <= NBASCP; I++) {
        K = LBASCP + (I - 1) * MBASCP;
        if (ALLOC4(K + 5) != 0) goto L235;
        ILLOC(K + 10) = 0;
        goto L259;
L235:
        if (ILLOC(K + 10) != (int)NOTDEF) goto L259;
        {
            IORDER = ILLOC(K + 2);
            IMULT  = ILLOC(K + 3);
            NCOUFF = NCOUFF + 1;
            for (II = 1; II <= NBASCP; II++) {
                KK = LBASCP + (II - 1) * MBASCP;
                if (ILLOC(KK + 2) != IORDER ||
                    ILLOC(KK + 3) != IMULT  ||
                    ALLOC4(KK + 5) == 0) goto L249;
                ILLOC(KK + 10) = NCOUFF;
L249:           ;
            }
        }
L259:   ;
    }

    std::printf("0THERE ARE%3d NUCLEAR AND%3d COULOMB FORM-FACTOR RADIAL FUNCTIONS.\n",
                NNUCFF, NCOUFF);

    //
    // NUCLEAR STORAGE FOR DERIVATIVE OF W.S.
    //
    {
        int I_tmp = MAX0((2 * NUMPT * NNUCFF + FACFR4 - 1) / FACFR4, 1);
        INUCH = NALLOC("NUCH    ", I_tmp);
        I_tmp = MAX0((NUMPT * NCOUFF + FACFR4 - 1) / FACFR4, 1);
        ICOULH = NALLOC("COULH   ", I_tmp);
    }

    //
    // NOW SETUP THE POTENTIALS * GAUSS WEIGHTS * R'S
    //
    // FIRST WE GET THE POTENTIALS ON THE ARRAY ( 0, SUMMAX )
    // THEN WE DEFINE CUBICS AND USE THEM TO GET THE DERIVATIVES.
    //
    // TO GET THE POTENTIALS WE SETUP SOME WORK AREAS INPLACE OF
    // THE CHANNEL 1 ARRAYS.
    //

L300:
    ICHANW  = 1;
    NSTEP1  = NSTPSS[1];
    {
        int IREAL1 = IRLVS[1];
        int IIMAG1 = IIMVS[1];
        int ICENT1 = ICENTR[1];

        //
        // WE USE THE STEPSIZE OF THE FIRST CHANNEL
        // WE GO SLIGHTLY FURTHER TO GIVE THE SPLINES TIME TO SETTLE DOWN
        //
        STEPSZ   = RSTEPS[1];
        double H = HS[1];
        E        = ES[1];
        NSTEP    = (int)(SUMMAX / STEPSZ + 20.5);
        NSTEPS   = NSTEP + 1;
        NSTPSS[1] = NSTEP;
        IRLVS[1]  = NALLOC("INELREAL", NSTEPS);
        IIMVS[1]  = NALLOC("INELIMAG", NSTEPS);
        ICENTR[1] = NALLOC("INELCENT", NSTEPS);
        IVWORK    = NALLOC("VWORK   ", NSTEPS);

        //
        // SETUP STUFF AS IF IN FIRST CHANNEL
        //
        AMP = AMA;
        AMT = AMBIGA;
        AM  = REDMI;
        IZP = IZS[1];
        IZT = IZS[3];

        //
        // INITIALIZE LINKULES IF NEEDED
        //
        for (II = 1; II <= 6; II++) {
            I = II;
            if (II == 6) I = 13;
            if (LNKBLK.LNKADR[I][3] == 0) goto L359;

            //
            // GENERATE THE SPECIAL UNIQUE NAME
            //
            IUNIQU = IUNIQU + 1;
            // LNKNAM is a char8* overlay on LNKADR
            // LINKID(4) = '*300' + updated digits
            {
                // Build the 4-character LINKID as in Fortran:
                // character*1 LINKID(4) / '*', '3', '0', '0' /
                // LINKID(3) = char('0' + MOD(IUNIQU/10,10))
                // LINKID(4) = char('0' + MOD(IUNIQU,10))
                char LINKID[5];
                LINKID[0] = '*';
                LINKID[1] = '3';
                LINKID[2] = (char)('0' + MOD(IUNIQU / 10, 10));
                LINKID[3] = (char)('0' + MOD(IUNIQU, 10));
                LINKID[4] = '\0';

                //
                // MAKE INITIALIZING CALL
                //
                int IRETUR_tmp = 0;
                LINKUL(LNKBLK.LNKADR[I][3], LNKNAM()[I - 1],
                       &LNKBLK.LNKADR[I][5], I, 1,
                       IRETUR_tmp, L_int, J_dbl, 0.0, STEPSZ, NSTEP + 1,
                       nullptr, nullptr, LINKID);
                IRET = IRETUR_tmp;
                if (IRET >= 0) goto L340;
                IRET = 0;
                // Restore first-channel pointers before returning
                NSTPSS[1] = NSTEP1;
                IRLVS[1]  = IREAL1;
                IIMVS[1]  = IIMAG1;
                ICENTR[1] = ICENT1;
                return;
L340:
                LNKBLK.LNKADR[I][4] = IRET;
            }
L359:       ;
        }

        //
        //
        // NOW COMPUTE THE POTENTIAL, NWP=3 IS CONVERT TO NWP=1
        // BY MAKPOT EXCEPT THE COULOMB IS NOT COMPUTED.
        //
        // FOR C.C. THE RESULTS OF THIS CALL ARE USED FOR THE
        // FIRST AND SECOND DERIVATIVE FORMS.  LX=-1 INDICATES THIS
        //
        LX   = -1;
        IRET = 1;
        J_dbl = NOTDEF;
        if (IFIRST + ISECND != 0) MAKPOT(3, IRET);
        DRIVSW = (KFRMTP == 1);
        if (IRET >= 0) goto L400;
        IRET = 0;
        // Restore first-channel pointers before returning
        NSTPSS[1] = NSTEP1;
        IRLVS[1]  = IREAL1;
        IIMVS[1]  = IIMAG1;
        ICENTR[1] = ICENT1;
        return;

        //
        // SUMMARIZE THE POTENTIALS TO BE USED
        //
        // FOLLOWING LIFTED FROM WAVSET
        //
L400:
        if (!DRIVSW) std::printf("-POTENTIALS WHOSE DERIVATIVES ARE THE EFFECTIVE INTERACTION OPERATOR:\n");
        if (DRIVSW)  std::printf("-RADIAL FORM FACTORS FOR FIRST-ORDER COUPLINGS.\n");
        if (DRIVSW && ISECND != 0) std::printf(" (DERIVATIVES OF THESE ARE USED FOR SECOND ORDER.)\n");

        SOSW  = false;
        SOISW = false;
        ETA   = KANDM.ETAS[1];

        //
        // IF POSSIBLE, RECOMPUTE R0'S FOR PRINTING
        //
        TEMP1_v = 0.0;
        TEMP2_v = 0.0;
        TEMP3_v = 0.0;
        TEMP4_v = 0.0;
        TEMP5_v = 0.0;
        TEMP6_v = 0.0;
        if (R0MASS == UNDEF) goto L420;
        TEMP1_v = R  / R0MASS;
        if (TVI != 0.0)  TEMP2_v = RI   / R0MASS;
        if (SOSW)        TEMP3_v = RSO  / R0MASS;
        if (SOISW)       TEMP4_v = RSOI / R0MASS;
        if (VSI != 0.0)  TEMP5_v = RSI  / R0MASS;
        TEMP6_v = RC / R0MASS;

L420:
        std::printf("0POTENTIAL           COUPLING CONS.        RADIUS     DIFFUSENESS     RADIUS PARAMETER\n\n");
        PARPRT(1,  "REAL CENTRAL       ", TVR,  R,   TAR,  TEMP1_v, POWRL, "   POWER", IRET);
        PARPRT(2,  "VOLUME ABSORPTION   ", TVI,  RI,  TAI,  TEMP2_v, POWIM, "   POWER   ", IRET);
        PARPRT(13, "SURFACE ABSORPTION   ", VSI,  RSI, 0.0,  TEMP5_v, 0.0,   "        ", IRET);
        PARPRT(3,  "REAL SPIN-ORBIT   ", VSO,  RSO, 0.0,  TEMP3_v, TAU,   "     TAU", IRET);
        PARPRT(4,  "IMAG. SPIN-ORBIT   ", VSOI, RSOI, 0.0, TEMP4_v, TAUI,  "    TAUI", IRET);
        PARPRT(5,  "COULOMB           ", ETA,  RC,  0.0,  TEMP6_v, 0.0,   "        ", IRET);
        if ((RCP == 0.0) || (RCT == 0.0)) goto L500;
        std::printf(" FOLDED COULOMB POTENTIALS - RCP = %7.4f  RCT = %7.4f  RC0P = %7.4f  RC0T = %7.4f\n",
                    RCP, RCT, RC0P, RC0T);

L500:
        std::printf("0\n");
        goto L600;

        //
        // PREPARE TO COMPUTE DERIVATIVE OF POTENTIAL
        // IF THE DERIVATIVE WAS DIRECTLY SUPPLIED, THEN WE
        // DON'T NEED TO COMPUTE IT UNLESS IF THERE IS 2ND ORDER.
        //

L600:
        ICUBIC = 0;
        if (IFIRST + ISECND == 0) goto L605;
        if (ISECND == 0 && DRIVSW && CCSW) goto L605;
        ICUBIC = NALLOC("INELCUBS", 4 * NSTEPS);
        LR = LOCPTRS.Z[ICUBIC];
        LB = LR + NSTEPS;
        LC = LB + NSTEPS;
        LD = LC + NSTEPS;
L605:
        LREALV = LOCPTRS.Z[IRLVS[1]] - 1;
        LIMAGV = LOCPTRS.Z[IIMVS[1]] - 1;
        LVWORK = LOCPTRS.Z[IVWORK];
        LNUCH  = LOCPTRS.Z[INUCH] - 2;    // LNUCH = Z(INUCH) - 2
        LCOULH = LOCPTRS.Z[ICOULH] - 1;   // LCOULH = Z(ICOULH) - 1
        if (IFIRST + ISECND == 0) goto L630;
        if (CCSW) goto L620;

        //
        // GET THE GAUSS POINTS FOR DWBA
        //
        LRPTS  = LOCPTRS.Z[IRPTS] - 1;
        LRWTS  = LOCPTRS.Z[IRWTS] - 1;
        LRPTS4 = FACFR4 * LOCPTRS.Z[IRPTS4] - FACFR4;
        CUBMAP(MAPSUM, SUMMIN, SUMMID, SUMMAX, GAMSUM,
               ALLOC_base(LRPTS + 1), ALLOC_base(LRWTS + 1), NUMPT);

        //
        //
        // FIRST DO THE REAL AND COULOMB PARTS
        //
        // GET THE SPLINES
        // GENERATE INPUT GRID FOR THE SPLINES
        //

L620:
        RVAL = 0.0;
        if (ISECND == 0 && DRIVSW && CCSW) goto L630;
        for (I = 1; I <= NSTEPS; I++) {
            ALLOC(LR - 1 + I) = RVAL;
            RVAL = RVAL + STEPSZ;
        }

        SPLNCB(NSTEPS, &ALLOC(LR), &ALLOC(LREALV + 1),
               &ALLOC(LB), &ALLOC(LC), &ALLOC(LD));

L630:
        VC     = FORMF.RBNDS[1];  // R2S(4) — EQUIVALENCE( R2S(1), RBNDS(1) ), so R2S(4)=RBNDS[4]
        // Note: R2S is a 4-element array equivalenced to RBNDS(1).
        // R2S(4) = RBNDS(4) = FORMF.RBNDS[4]
        VC     = FORMF.RBNDS[4];
        RCUEFF = RC;

        if (CCSW) goto L670;

        //
        // ASSIGN 652 TO IGOTO  (derivative form)
        // ASSIGN 654 TO IGOTO  (direct form, when DRIVSW)
        //
        IGOTO_real = DRIVSW ? 1 : 0;  // 0=>652 (derivative), 1=>654 (direct)
        for (I = 1; I <= NUMPT; I++) {
            RVAL = ALLOC(LRPTS + I);
            WT   = ALLOC(LRWTS + I);
            ALLOC4(LRPTS4 + I) = (float)RVAL;
            N_val = (int)(RVAL / STEPSZ);
            X     = RVAL - N_val * STEPSZ;

            if (IGOTO_real == 0) {
                // 652: derivative of W.S.
                XX = (ALLOC(LB + N_val) + X * (2.0 * ALLOC(LC + N_val) + X * 3.0 * ALLOC(LD + N_val)));
                goto L655;
            } else {
                // 654: direct form (KFRMTP==1)
                XX = ALLOC(LREALV + N_val + 1) + X * (ALLOC(LB + N_val) + X * (ALLOC(LC + N_val) + X * ALLOC(LD + N_val)));
                // fall through to L655
            }
L655:
            // TERM = + (12*E/H**2) * R2S(1) * WT * XX
            // R2S(1) = RBNDS(1) = FORMF.RBNDS[1]
            TERM = (12.0 * E / (HS[1] * HS[1])) * FORMF.RBNDS[1] * WT * XX;
            ALLOC(LNUCH + 2 * I) = TERM;
            if (RVAL < RC) XX = std::pow(RVAL / (RC * RC), (double)LXMIN) / RC;
            if (RVAL >= RC) XX = 1.0 / std::pow(RVAL, (double)(LXMIN + 1));
            ALLOC(LCOULH + I) = -WT * VC * XX;
        }
        goto L700;

        //
        // SET UP RADIAL PARTS OF REAL POT FOR COUPLED CHANNELS
        // WE STORE
        //    - (2*MU/HBARC**2) (STEPSIZE**2/12) DV/DR
        //
        // NOTE THAT MAKPOT HAS GENERATED ESPECIALLY FOR US
        //    - (H**2/12) (V/E)
        //
L670:
        LNUCH = FACFR4 * LOCPTRS.Z[INUCH] - FACFR4 + 1;
        WT    = -(2.0 * REDMI / (HBARC * HBARC)) * (STEPSZ * STEPSZ / 12.0) * (12.0 * E / (HS[1] * HS[1]));
        if (IFIRST == 0) goto L672;
        {
            K = LNUCH + 2 * NUMPT * (IFIRST - 1);
            if (DRIVSW) LB = LREALV + 1;

            for (I = 1; I <= NUMPT; I++) {
                ALLOC4(K + 2 * I - 2) = (float)(-WT * ALLOC(LB + I - 1));
            }
        }

L672:
        if (ISECND == 0) goto L700;

        //
        // NOW GO ON TO SECOND ORDER
        // WE STORE
        //    + (2*MU/HBARC**2) (STEPSIZE**2/12) D**2 V / DR**2
        //
        SPLNCB(NSTEPS, &ALLOC(LR), &ALLOC(LB),
               &ALLOC(LC), &ALLOC(LD), &ALLOC(LVWORK));

        {
            K = LNUCH + 2 * NUMPT * (ISECND - 1);
            for (I = 1; I <= NUMPT; I++) {
                ALLOC4(K + 2 * I - 2) = (float)(+WT * ALLOC(LC + I - 1));
            }
        }

        //
        // NOW GO ON TO THE IMAGINARY PART
        //
L700:
        if (IFIRST + ISECND == 0) goto L770;
        if (ISECND == 0 && DRIVSW && CCSW) goto L750;
        SPLNCB(NSTEPS, &ALLOC(LR),
               &ALLOC(LIMAGV + 1), &ALLOC(LB), &ALLOC(LC), &ALLOC(LD));

        if (CCSW) goto L750;

        //
        // ASSIGN 732 TO IGOTO  (derivative)
        // ASSIGN 734 TO IGOTO  (direct, when DRIVSW)
        //
        IGOTO_imag = DRIVSW ? 1 : 0;  // 0=>732, 1=>734
        for (I = 1; I <= NUMPT; I++) {
            RVAL = ALLOC(LRPTS + I);
            WT   = ALLOC(LRWTS + I);
            N_val = (int)(RVAL / STEPSZ);
            X     = RVAL - N_val * STEPSZ;

            if (IGOTO_imag == 0) {
                // 732: derivative
                XX = (ALLOC(LB + N_val) + X * (2.0 * ALLOC(LC + N_val) + X * 3.0 * ALLOC(LD + N_val)));
                goto L735;
            } else {
                // 734: direct
                XX = ALLOC(LIMAGV + N_val + 1) + X * (ALLOC(LB + N_val) + X * (ALLOC(LC + N_val) + X * ALLOC(LD + N_val)));
            }
L735:
            // R2S(2) = RBNDS(2) = FORMF.RBNDS[2]
            TERM = (12.0 * E / (HS[1] * HS[1])) * FORMF.RBNDS[2] * WT * XX;
            ALLOC(LNUCH + 2 * I + 1) = TERM;
        }
        goto L950;

        //
        // SET UP RADIAL PARTS OF IMAG POT FOR COUPLED CHANNELS
        //
L750:
        if (DRIVSW) LB = LIMAGV + 1;
        if (IFIRST == 0) goto L760;
        {
            K = LNUCH + 2 * NUMPT * (IFIRST - 1);
            for (I = 1; I <= NUMPT; I++) {
                ALLOC4(K + 2 * I - 1) = (float)(-WT * ALLOC(LB + I - 1));
            }
        }

L760:
        if (ISECND == 0) goto L770;

        //
        //  NOW GO ON TO SECOND ORDER
        //
        SPLNCB(NSTEPS, &ALLOC(LR), &ALLOC(LB),
               &ALLOC(LC), &ALLOC(LD), &ALLOC(LVWORK));

        {
            K = LNUCH + 2 * NUMPT * (ISECND - 1);
            for (I = 1; I <= NUMPT; I++) {
                ALLOC4(K + 2 * I - 1) = (float)(+WT * ALLOC(LC + I - 1));
            }
        }

        //
        // NOW PROCESS DEFORMED REQUESTS FOR C.C.
        //
L770:
        // Re-derive LBASCP (pointer may have walked in the nuclear-form loops)
        LBASCP = FACFR4 * LOCPTRS.Z[IBASCP] - FACFR4 + 1;
        LBASCP = ILLOC(LBASCP + 2) + LBASCP;
        IB1 = 1;
        if (NNUCFF == 0) goto L800;
        for (IP = 1; IP <= NNUCFF; IP++) {
            if (IP == IFIRST || IP == ISECND) goto L799;

            //
            // FIND AN EXAMPLE OF THIS RADIAL FORM
            //
L775:
            IPP = (ILLOC(LBASCP + 9) + 1) / 2;
            if (IPP == IP) goto L780;
            LBASCP = LBASCP + MBASCP;
            IB1 = IB1 + 1;
            if (IB1 <= NBASCP) goto L775;
            std::printf("0*** INTERNAL ERROR AT 778 IN INGRST%8d%8d\n", IP, IPP);
            IRET = 0;
            // Restore first-channel pointers before returning
            NSTPSS[1] = NSTEP1;
            IRLVS[1]  = IREAL1;
            IIMVS[1]  = IIMAG1;
            ICENTR[1] = ICENT1;
            return;

L780:
            L_int = ILLOC(LBASCP + 2);
            LX    = ILLOC(LBASCP + 3) / 2;
            J_dbl = ILLOC(LBASCP + 7);

            //
            // FOR VALUES OF MODL SEE IPT1 IN BASCPL JUST AFTER 700
            // HERE FOR MODEL 4 (HOTFUDGE) WE MUST SET UP THE WOODS SAXON
            // PARAMETERS FROM THE PARAMS.  ALSO THERE IS NO IMAG COUPLING
            // IN HOTFUDGE.
            //
            MODL = MOD((int)J_dbl, 100);
            T1   = TVR;
            T2   = R;
            T3   = TAR;
            T4   = POWRL;
            T5   = TVI;
            T6   = VSI;
            IT1  = LNKBLK.LNKADR[1][3];
            IT2  = LNKBLK.LNKADR[2][3];
            if (MODL != 6) goto L782;
            TVR   = 1.0;
            INDXP = (int)J_dbl / 100;
            R     = R0MASS * PARAMS()[INDXP - 1];      // PARAMS(INDXP) — 1-based
            TAR   = PARAMS()[INDXP];                   // PARAMS(INDXP+1) — 1-based
            POWRL = 1.0;
            TVI   = 0.0;
            VSI   = 0.0;
            LNKBLK.LNKADR[1][3] = 0;
            LNKBLK.LNKADR[2][3] = 0;
            std::printf(" HOTFUDGE:%4d     R =%7.3f     A =%7.3f\n", INDXP, R, TAR);

L782:
            if (DBUGSW) std::printf(" GENERATING NUCLEAR FF:%10d%10d%10d%10d%10d\n",
                                    IP, IB1, L_int, LX, (int)J_dbl);
            MAKPOT(3, IRET);
            TVR   = T1;
            R     = T2;
            TAR   = T3;
            POWRL = T4;
            TVI   = T5;
            VSI   = T6;
            LNKBLK.LNKADR[1][3] = IT1;
            LNKBLK.LNKADR[2][3] = IT2;

            //
            // IN THIS CASE WE JUST COPY THE RESULTS OF MAKPOT
            //
            K = LNUCH + 2 * NUMPT * (IP - 1);
            for (I = 1; I <= NUMPT; I++) {
                ALLOC4(K + 2 * I - 2) = (float)(-WT * ALLOC(LREALV + I));
                ALLOC4(K + 2 * I - 1) = (float)(-WT * ALLOC(LIMAGV + I));
            }
L799:       ;
        }

        //
        // SET UP RADIAL PARTS OF COULOMB POT FOR COUPLED CHANNELS
        //
L800:
        if (NCOUFF == 0) goto L950;
        {
            int LCOULH_loc = LOCPTRS.Z[ICOULH];
            LCOULH_loc = FACFR4 * (LCOULH_loc - 1) + 1;
            LBASCP = FACFR4 * LOCPTRS.Z[IBASCP] - FACFR4 + 1;
            LBASCP = ILLOC(LBASCP + 2) + LBASCP;
            IP = 1;

            //
            // LOOP OVER ALL COUPLINGS FINDING 1ST INSTANCE
            // OF EACH RADIAL FORM.
            //
            for (int J_idx = 1; J_idx <= NBASCP; J_idx++) {
                LABEL  = (J_idx - 1) * MBASCP + LBASCP;
                IORDER = ILLOC(LABEL + 2);
                IMULT  = ILLOC(LABEL + 3);
                IPOINT = ILLOC(LABEL + 10);
                if (IPOINT != IP) goto L869;
                IMT = IMULT / 2;

                //
                // THE RADIAL FACTORS IN INGRST ARE COMPUTED SUCH THAT
                // WHEN MULTIPLIED BY THE X FACTORS IN BASCPL THE CORRECT
                // ANSWER RESULTS
                //
                // THEY ARE
                //
                //   - (2*MU/HBARC**2) (STEPSIZE**2/12) (3 ZP ZT ALPHA HBARC) F(R)
                // WHERE
                //   FIRST ORDER:
                //     F(R) = R**LX / RC**(2LX+1)    R < RC
                //          = 1/R**(LX+1)            R > RC
                //   SECOND ORDER:
                //     F(R) = (1-LX)/(2+LX) R**LX / R**(2LX+1)   R < RC
                //          = 1 / R**(LX+1)                      R > RC
                //
                // THE FOLLOWING FORMULAS ARE CORRECT FOR A POINT AND SINGLE
                // DEFORMED SPHERE OF DEFORMATION RADIUS RC.  IN THIS CASE WE
                // WOULD ALSO HAVE RT=RC IN THE COUPLING CONSTANTS.  WE USE
                // THESE FORMULAS EVEN WHEN RP IS NOT ZERO.  THIS GIVES THE
                // CORRECT RESULT FOR R > RC  BUT IS INCORRECT FOR R < RC.
                // HOWEVER, R < RC IS NOT IMPORTANT.
                //
                RVAL   = -STEPSZ;
                WT     = -(2.0 * REDMI / (HBARC * HBARC)) * (STEPSZ * STEPSZ / 12.0) * VC;
                IPOW   = IMT;
                TMPFAC = 1.0 / std::pow(RC, (double)(IMULT + 1));
                if (IORDER == 2) TMPFAC = (1.0 - (double)IMT) / (2.0 + (double)IMT) * TMPFAC;

                // ASSIGN 830 TO IGOTO — start in "inside RC" mode
                // We use a bool to track whether we've crossed RC
                {
                    int crossed_RC = 0;
                    if (DBUGSW) std::printf(" GENERATING COULOMB FF:%10d%10d%10d%10d/2%10d\n",
                                            IP, J_idx, IORDER, IMULT, ILLOC(LABEL + 7));
                    for (I1 = 1; I1 <= NUMPT; I1++) {
                        RVAL = RVAL + STEPSZ;

                        //
                        // CHECK TO SEE IF  R < RC  STILL
                        //
                        if (!crossed_RC) {
                            // L830: test
                            if (RVAL >= RC) {
                                crossed_RC = 1;
                                IPOW   = -(IMT + 1);
                                TMPFAC = 1.0;
                            }
                        }
                        // L835:
                        XX = TMPFAC * std::pow(RVAL, (double)IPOW);
                        ALLOC4(LCOULH_loc + I1 - 1 + (IP - 1) * NUMPT) = (float)(WT * XX);
                    }
                }
                IP = IP + 1;
                if (IP > NCOUFF) goto L950;
L869:           ;
            }
        }

        //
        //
        // FREE WORK AREAS
        //
L950:
        if (ICUBIC != 0) LOCPTRS.Z[ICUBIC] = -LOCPTRS.Z[ICUBIC];
        LOCPTRS.Z[IRLVS[1]]  = -LOCPTRS.Z[IRLVS[1]];
        LOCPTRS.Z[IIMVS[1]]  = -LOCPTRS.Z[IIMVS[1]];
        LOCPTRS.Z[ICENTR[1]] = -LOCPTRS.Z[ICENTR[1]];
        LOCPTRS.Z[IVWORK]    = -LOCPTRS.Z[IVWORK];

        //
        // RESTORE FIRST CHANNEL
        //
        NSTPSS[1] = NSTEP1;
        IRLVS[1]  = IREAL1;
        IIMVS[1]  = IIMAG1;
        ICENTR[1] = ICENT1;

        if (IPRT == 0) goto L990;

        //
        // SUMMARIZE THINGS NICELY
        //
        std::printf("-ONE-DIMENSIONAL INTEGRATION GRID:\n");
        std::printf("0NUM. PTS.   MAP TYPE   GAMMA    MINIMUM   \"MID. PT.\"   MAXIMUM\n");
        std::printf("%6d%11d%11.2f%10.2f%12.2f%14.2f\n",
                    NMPTSV, MAPSUM, GAMSUM, SUMMIN, SUMMID, SUMMAX);

L990:
        IRET = 1;
        return;
    }  // end of local-variable scope block (IREAL1, IIMAG1, ICENT1)
}
