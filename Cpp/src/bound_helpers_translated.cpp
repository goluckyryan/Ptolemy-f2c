// bound_helpers_translated.cpp -- BSPROD+BSSET, CLINTS, JPTOLX
// Translated from source.f lines 5085-5588, 8490-9002, 23065-23147
// Verbatim transliteration from Fortran

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include "ptolemy_alloc.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_forward.h"
#include "ptolemy_io.h"
#include <cstdio>
#include <cmath>

// Forward declarations for routines called here
extern void VCSQ12(double RVAL, double& X, int NWP);
extern void SETVSQ(double RR1, double RR2, int IZ1, int IZ2, int K);
extern int NALLOC(const char* NAME, int IWRDS);
extern double AITLAG(double X, double STPINV, double* TABLE, int NTABLE, int NAIT);
extern void GAUSSL(int N, double* X, double* W);
extern void RCWFN(double RHO, double ETA, int MINL, int MAXL,
                  double* FC, double* FCP, double* GC, double* GCP,
                  double ACCUR, int& IRET);
extern void RCASYM(int L, double ETA, double RHO, int IPRN, double SIGL,
                   double* ZETA, double* PHI, double* DZETA,
                   double* F, double* FP, double* G, double* GP,
                   double* Z, double* DZSQ, double* S, double* ZINV,
                   double EPS, int& NMAX, int& NTZ, int& ISIG);
extern double RACAH(int A_r, int B_r, int Y_r, int X_r, int C_r, int Z_r);

// ============================================================================
// BSPROD + BSSET  (source.f lines 5085-5588)
//
// BSPROD: Computes bound-state product PHI(RT)*VEFF*PHI(RP)
// BSSET:  ENTRY point - initializes form-factor calculation
//
// Uses enum dispatch + _impl pattern for ENTRY BSSET
// Returns: 0 = normal return, 1 = alternate return (RP or RT exceeds max)
// ============================================================================

enum class BsprodEntry { BSPROD, BSSET };

static int BSPROD_impl(BsprodEntry entry,
    double& FPFT, int ITYPE, double RA, double RB, double X,
    int ISCAT, int NAIT, double& RP, double& RT)
{
    // SAVE variables (persist across calls, shared between BSPROD and BSSET)
    static int ISC;
    static int IZS1, IZS2, IZC1, IZC2;
    static int ISEE;
    static int IPRNT;
    static double RNSCAT, RNCORE, VOPT, AOPT;
    static double RCCOP, RCCOT, RCPSCT, RCTSCT;

    // Locals
    double FP, FT, FACTOR, DELVNU, RSCAT, RCORE;
    double POT, DELVCO, DELVSC, DELV, XCORE, FCORE, XSCAT, FSCAT;
    double RBND, VEFF, TEMP, RRP, RRT;
    int I, IVPHI, LVPHI, LPOT, LBD1, LBD2;
    int JBD1, JBD2;
    double AMA3, AMB3, AMBGA3, AMBGB3, POW;
    logical ALLSW;

    if (entry == BsprodEntry::BSSET) goto L_BSSET;

    // ---------------------------------------------------------------
    // BSPROD entry point
    // ---------------------------------------------------------------
    //
    // STEP 1 - COMPUTE  PHI V PHI   OR   PHI' V PHI'
    //
    RP = (BSBLK.S2*RA)*(BSBLK.S2*RA) + (BSBLK.T2*RB)*(BSBLK.T2*RB) +
       2*(BSBLK.S2*RA)*(BSBLK.T2*RB)*X;
    RT = (BSBLK.S1*RA)*(BSBLK.S1*RA) + (BSBLK.T1*RB)*(BSBLK.T1*RB) +
       2*(BSBLK.S1*RA)*(BSBLK.T1*RB)*X;
    if ( DABS(RP) < 1.0e-8 )  RP = 0;
    if ( DABS(RT) < 1.0e-8 )  RT = 0;
    RP = DSQRT(RP);
    RT = DSQRT(RT);

    FPFT = 0;
    if ( RP > BSBLK.BNDMXP  ||  RT > BSBLK.BNDMXT )  return 1;

    FP = BSBLK.VPMAX;
    FT = BSBLK.VTMAX;
    ALLSW = MOD( ITYPE, 2 ) != 0;
    DELVNU = 0;

    if ( ALLSW || RT > BSBLK.RLTMAX )
      FT = AITLAG(RT, BSBLK.BNDSPT, ALLOC_base(LOCPTRS.Z[BSBLK.JBDT]), BSBLK.NSPBDT, NAIT);
    if ( ALLSW || RP > BSBLK.RLPMAX )
      FP = AITLAG(RP, BSBLK.BNDSPP, ALLOC_base(LOCPTRS.Z[BSBLK.JBDP]), BSBLK.NSPBDP, NAIT);
    FACTOR = 1;
    FPFT = FP*FT;
    if ( SWITCH.NUCONL != 2  &&  SWITCH.NUCONL != 3 ) goto L200;
    if ( DABS(FPFT) < CNSTNT.SMLNUM )  goto L200;

    //
    // MULTIPLY V*PHI BY (V+DV)/V ; THIS GIVES (V+DV)*PHI
    // DV IS CORE CORRECTION TO VEFF
    //
    RSCAT = RA;
    if ( BSBLK.IOUTSW )  RSCAT = RB;
    RCORE = DSQRT( ((BSBLK.S1-BSBLK.S2)*RA)*((BSBLK.S1-BSBLK.S2)*RA)
      + ((BSBLK.T1-BSBLK.T2)*RB)*((BSBLK.T1-BSBLK.T2)*RB)
      + 2*((BSBLK.S1-BSBLK.S2)*RA)*((BSBLK.T1-BSBLK.T2)*RB)*X );
    VCSQ12(RCORE, POT, 3);
    DELVCO = POT;
    VCSQ12(RSCAT, POT, ISC);
    DELVSC = POT;
    DELV = DELVCO - DELVSC;
    if ((SWITCH.NUCONL != 3) || (NAIT == 2)) goto L8;
    XCORE = (RCORE - RNCORE) / AOPT;
    if (XCORE < CNSTNT.BIGLOG) goto L400;
    FCORE = 0;
    goto L440;
L400:
    if (XCORE > -CNSTNT.BIGLOG) goto L420;
    FCORE = 1;
    goto L440;
L420:
    FCORE = 1 / (1 + DEXP(XCORE));
L440:
    XSCAT = (RSCAT - RNSCAT) / AOPT;
    if (XSCAT < 174) goto L460;
    FSCAT = 0;
    goto L500;
L460:
    if (XSCAT > -175) goto L480;
    FSCAT = 1;
    goto L500;
L480:
    FSCAT = 1 / (1 + DEXP(XSCAT));
L500:
    DELVNU = VOPT * (FCORE - FSCAT);
    DELV = DELV + DELVNU;
L8:
    RBND = RT * (SWITCH.IVRTEX - 1) + RP * (2 - SWITCH.IVRTEX);
    VEFF = AITLAG(RBND, BSBLK.BNDSPV, ALLOC_base(LOCPTRS.Z[BSBLK.JPOT]), BSBLK.NSPBDV, NAIT);
    FACTOR = 1 + DELV / VEFF;
    FPFT = FACTOR * FPFT;
    if ( IPRNT >= 9 )
      std::printf(" BSP%8.4f%8.4f%8.4f%8.4f%8.4f%12.4g%12.4g%12.4g%12.4g%12.4g%12.4g%12.4g\n",
        RA, RB, RP, RT, RCORE, VEFF, DELVSC, DELVCO, DELVNU, DELV, FACTOR, FPFT);

    //
    // NOW WE MULTIPLY IN THE PSI'S IF DESIRED.
    // IN FACT WE MULTIPLY IN R*PSI
    //
L200:
    if ( ITYPE <= 2 )  return 0;

    if ( RA > BSBLK.SCTMAX )  goto L220;
    FPFT = FPFT *
      AITLAG( RA, BSBLK.SCTSP, ALLOC_base(LOCPTRS.Z[ISCAT]), BSBLK.NSPSCT, NAIT );
    // AVOID UNDERFLOWS
    if ( DABS(FPFT) > 1.0e-34 )  goto L250;
    FPFT = 0;
    return 0;
L220:
    FPFT = FPFT * RA;

L250:
    if ( RB > BSBLK.SCTMAX )  goto L270;
    FPFT = FPFT *
      AITLAG( RB, BSBLK.SCTSP, ALLOC_base(LOCPTRS.Z[ISCAT]), BSBLK.NSPSCT, NAIT );
    return 0;
L270:
    FPFT = FPFT * RB;
    return 0;

    // ---------------------------------------------------------------
    // BSSET entry point
    // INITIALIZE THE FORM-FACTOR CALCULATION
    // ---------------------------------------------------------------
L_BSSET:

    //
    // COMPUTE APPROPRIATE  V*PHI  AND STORE IT IN "VPHI".
    //
    IPRNT = MOD(INTGER.IPRINT / 10000, 10);
    BSBLK.JPOT = FORMF.NPOTS[SWITCH.IVRTEX];
    JBD1 = FORMF.IBDS[SWITCH.IVRTEX];
    JBD2 = FORMF.IBDS[3 - SWITCH.IVRTEX];
    BSBLK.NSPBDV = FORMF.NSTPBD[SWITCH.IVRTEX];
    BSBLK.NSPBDT = FORMF.NSTPBD[3 - SWITCH.IVRTEX];
    IVPHI = NALLOC( "VPHI    ", BSBLK.NSPBDV );
    LVPHI = LOCPTRS.Z[IVPHI] - 1;
    LPOT = LOCPTRS.Z[BSBLK.JPOT] - 1;
    LBD1 = LOCPTRS.Z[JBD1] - 1;
    LBD2 = LOCPTRS.Z[JBD2] - 1;
    VMAXS()[SWITCH.IVRTEX - 1] = 0;
    VMAXS()[3 - SWITCH.IVRTEX - 1] = 0;

    RRP = -FORMF.BNDSTP[SWITCH.IVRTEX];
    for (I = 1; I <= BSBLK.NSPBDV; I++) {
        RRP = RRP + FORMF.BNDSTP[SWITCH.IVRTEX];
        TEMP = ALLOC(LBD1 + I) * ALLOC(LPOT + I);
        if ( DABS(TEMP) <= VMAXS()[SWITCH.IVRTEX - 1] )  goto L215;
        RLMAXS()[SWITCH.IVRTEX - 1] = RRP;
        VMAXS()[SWITCH.IVRTEX - 1] = DABS(TEMP);
L215:
        ALLOC(LVPHI + I) = TEMP;
    }
    RRT = -FORMF.BNDSTP[3 - SWITCH.IVRTEX];
    for (I = 1; I <= BSBLK.NSPBDT; I++) {
        TEMP = ALLOC(LBD2 + I);
        RRT = RRT + FORMF.BNDSTP[3 - SWITCH.IVRTEX];
        if ( DABS(TEMP) <= VMAXS()[3 - SWITCH.IVRTEX - 1] ) goto L229;
        RLMAXS()[3 - SWITCH.IVRTEX - 1] = RRT;
        VMAXS()[3 - SWITCH.IVRTEX - 1] = DABS(TEMP);
L229:   ;
    }

    if (BSBLK.PINFSW)  std::printf(" PROJECTILE AND TARGET FORMFACTOR MAXIMA:%15.5g%15.5g\n"
      " LOCATION OF PROJECTILE AND TARGET FORM-FACTOR MAXIMA:%15.5g%15.5g\n",
      VMAXS()[0], VMAXS()[1], RLMAXS()[0], RLMAXS()[1]);

    //
    // FREE POTENTIAL AT VERTEX WHERE INTERACTION ISNT (IF NOT MULTIPASS)
    //
    if ( SWITCH.IBSPAS > 0 )  goto L1620;
    LOCPTRS.Z[FORMF.NPOTS[3 - SWITCH.IVRTEX]] = -LOCPTRS.Z[FORMF.NPOTS[3 - SWITCH.IVRTEX]];
    if (SWITCH.NUCONL >= 2) goto L1620;
    LOCPTRS.Z[BSBLK.JPOT] = -LOCPTRS.Z[BSBLK.JPOT];

L1620:
    BSBLK.JBDP = FORMF.IBDS[1];
    if ( SWITCH.IVRTEX == 1 )  BSBLK.JBDP = IVPHI;
    BSBLK.JBDT = FORMF.IBDS[2];
    if ( SWITCH.IVRTEX == 2 )  BSBLK.JBDT = IVPHI;
    BSBLK.BNDSPP = 1 / FORMF.BNDSTP[1];
    BSBLK.BNDSPT = 1 / FORMF.BNDSTP[2];
    BSBLK.BNDSPV = 1 / FORMF.BNDSTP[SWITCH.IVRTEX];
    BSBLK.NSPBDP = FORMF.NSTPBD[1];
    BSBLK.NSPBDT = FORMF.NSTPBD[2];
    BSBLK.BNDMXP = FORMF.BNDMAX[1];
    BSBLK.BNDMXT = FORMF.BNDMAX[2];

    //
    // EXPONENTIAL TAIL FACTORS FOR H-INTERPOLATION
    //
    BSBLK.ALPHAP = sqrt( -2*FORMF.BNDMAS[1]*INTRNL.EBNDS[1] ) / CNSTNT.HBARC;
    BSBLK.ALPHAT = sqrt( -2*FORMF.BNDMAS[2]*INTRNL.EBNDS[2] ) / CNSTNT.HBARC;
    if ( IPRNT >= 4 )  std::printf("\n EXPONENTIAL DECAYS OF FORM FACTOR:  RP =%13.5g     RT =%13.5g\n",
      BSBLK.ALPHAP, BSBLK.ALPHAT);

    if ( SWITCH.NUCONL < 2 )  return 0;

    //
    // COMPUTE POTENTIAL PARAMETERS FOR CORE TERMS IN VEFF
    //
    POW = 1.0 / 3.0;
    AMA3 = std::pow(FLOAT_common.AMA, POW);
    AMB3 = std::pow(FLOAT_common.AMB, POW);
    AMBGA3 = std::pow(FLOAT_common.AMBIGA, POW);
    AMBGB3 = std::pow(FLOAT_common.AMBIGB, POW);
    FACTOR = 1;

    //
    // FOLDED COULOMB POTENTIALS IN BOTH CHANNELS
    //
    ISEE = 0;
    for (I = 1; I <= 2; I++) {
        if ((KANDM.RCSCTP[I] == 0) || (KANDM.RCSCTT[I] == 0)) ISEE = 1;
    }

    //
    // PHISGN = +1 ; STRIPPING
    // PHISGN = -1 ; PICK-UP
    // IVRTEX =  1 ; VEFF AT PROJECTILE VERTEX
    // IVRTEX =  2 ; VEFF AT TARGET VERTEX
    //
    if (BSBLK.PHISGN == -1) goto L60;
    if (SWITCH.IVRTEX == 2) goto L40;

    //
    // STRIPPING ; PROJECTILE VERTEX
    //
    ISC = 2;
    IZS1 = INTGER.IZS[2];
    IZS2 = INTGER.IZS[4];
    IZC1 = INTGER.IZS[2];
    IZC2 = INTGER.IZS[3];
    if (ISEE == 0) FACTOR = AMBGA3 / AMBGB3;
    if (ISEE == 1) FACTOR = (AMBGA3 + AMB3) / (AMBGB3 + AMB3);
    RCPSCT = KANDM.RCSCTP[2];
    RCTSCT = KANDM.RCSCTT[2];
    RCCOP = RCPSCT;
    RCCOT = RCTSCT * FACTOR;
    BSBLK.IOUTSW = TRUE_F;
    if (SWITCH.NUCONL != 3) goto L100;
    RNSCAT = KANDM.RSCTS[2];
    RNCORE = RNSCAT * ((AMBGA3 + AMB3) / (AMBGB3 + AMB3));
    VOPT = -WAVCOM.V0RS[2];
    AOPT = KANDM.ASCTS[2];
    goto L100;

    //
    // STRIPPING ; TARGET VERTEX
    //
L40:
    ISC = 1;
    IZS1 = INTGER.IZS[1];
    IZS2 = INTGER.IZS[3];
    IZC1 = INTGER.IZS[2];
    IZC2 = INTGER.IZS[3];
    if (ISEE == 0) FACTOR = AMB3 / AMA3;
    if (ISEE == 1) FACTOR = (AMB3 + AMBGA3) / (AMA3 + AMBGA3);
    RCPSCT = KANDM.RCSCTP[1];
    RCTSCT = KANDM.RCSCTT[1];
    RCCOP = RCPSCT * FACTOR;
    RCCOT = RCTSCT;
    BSBLK.IOUTSW = FALSE_F;
    if (SWITCH.NUCONL != 3) goto L100;
    RNSCAT = KANDM.RSCTS[1];
    RNCORE = RNSCAT * ((AMB3 + AMBGA3) / (AMA3 + AMBGA3));
    VOPT = -WAVCOM.V0RS[1];
    AOPT = KANDM.ASCTS[1];
    goto L100;

    //
    // PICK-UP ; PROJECTILE VERTEX
    //
L60:
    if (SWITCH.IVRTEX == 2) goto L80;
    ISC = 1;
    IZS1 = INTGER.IZS[1];
    IZS2 = INTGER.IZS[3];
    IZC1 = INTGER.IZS[1];
    IZC2 = INTGER.IZS[4];
    if (ISEE == 0) FACTOR = AMBGB3 / AMBGA3;
    if (ISEE == 1) FACTOR = (AMBGB3 + AMA3) / (AMBGA3 + AMA3);
    RCPSCT = KANDM.RCSCTP[1];
    RCTSCT = KANDM.RCSCTT[1];
    RCCOP = RCPSCT;
    RCCOT = RCTSCT * FACTOR;
    BSBLK.IOUTSW = FALSE_F;
    if (SWITCH.NUCONL != 3) goto L100;
    RNSCAT = KANDM.RSCTS[1];
    RNCORE = RNSCAT * ((AMBGB3 + AMA3) / (AMBGA3 + AMA3));
    VOPT = -WAVCOM.V0RS[1];
    AOPT = KANDM.ASCTS[1];
    goto L100;

    //
    // PICK-UP ; TARGET VERTEX
    //
L80:
    ISC = 2;
    IZS1 = INTGER.IZS[2];
    IZS2 = INTGER.IZS[4];
    IZC1 = INTGER.IZS[1];
    IZC2 = INTGER.IZS[4];
    if (ISEE == 0) FACTOR = AMA3 / AMB3;
    if (ISEE == 1) FACTOR = (AMA3 + AMBGB3) / (AMB3 + AMBGB3);
    RCPSCT = KANDM.RCSCTP[2];
    RCTSCT = KANDM.RCSCTT[2];
    RCCOP = RCPSCT * FACTOR;
    RCCOT = RCTSCT;
    BSBLK.IOUTSW = TRUE_F;
    if (SWITCH.NUCONL != 3) goto L100;
    RNSCAT = KANDM.RSCTS[2];
    RNCORE = RNSCAT * ((AMA3 + AMBGB3) / (AMB3 + AMBGB3));
    VOPT = -WAVCOM.V0RS[2];
    AOPT = KANDM.ASCTS[2];

    //
    // INITIALIZE COULOMB-POTENTIAL SUBROUTINE
    //
L100:
    SETVSQ(RCCOP, RCCOT, IZC1, IZC2, 3);

    if (IPRNT <= 4) return 0;
    std::printf("\n  ++++ SETRS ++++\n   CHANNEL %2d(Z1,Z2) = ( %2d %2d )\n"
      "     CORE  (Z1,Z2) = (%2d %2d) NUCONL = %2d\n",
      ISC, IZS1, IZS2, IZC1, IZC2, SWITCH.NUCONL);
    std::printf("   SCATTERING :(RT,RP) = %12.6g %12.6g(COULOMB)\n"
      "   CORE :(RT,RP) = %12.6g %12.6g(COULOMB)\n",
      RCPSCT, RCTSCT, RCCOT, RCCOP);
    if (SWITCH.NUCONL != 3) return 0;
    std::printf("\n   (RNSCAT,RNCORE)=%12.6g %12.6g\n"
      "   (-V,A)=(%12.6g %12.6g %12.6g)\n",
      RNSCAT, RNCORE, VOPT, AOPT, AOPT);
    return 0;
}

// Public entry: BSPROD
// Returns 0 = normal, 1 = alternate return (RP or RT out of range)
int BSPROD(double& FPFT, int ITYPE, double RA, double RB, double X,
           int ISCAT, int NAIT, double& RP, double& RT)
{
    return BSPROD_impl(BsprodEntry::BSPROD, FPFT, ITYPE, RA, RB, X,
                       ISCAT, NAIT, RP, RT);
}

// Public entry: BSSET
void BSSET()
{
    double dummy_FPFT = 0, dummy_RP = 0, dummy_RT = 0;
    BSPROD_impl(BsprodEntry::BSSET, dummy_FPFT, 0, 0.0, 0.0, 0.0,
                0, 0, dummy_RP, dummy_RT);
}


// ============================================================================
// CLINTS  (source.f lines 8490-9002)
//
// Evaluates INTEGRAL(RLOWER, INF) DR FCOUL FCOUL / R**N
// Four integrals: FFINT, FGINT, GFINT, GGINT
// ============================================================================

void CLINTS(double RLOWER, double ETAI, double ETAF, double FKI, double FKF,
            double SIGI, double SIGF, double ACCURA_arg, double ASMULT,
            double& FFINT, double& FGINT, double& GFINT, double& GGINT,
            double* PTS, double* WTS, double* FA, double* FPA,
            double* GA, double* GPA, double* WORK,
            int N, int LI, int LF, int NTERMS, int NPTS,
            int& IRET, int IPRINT)
{
    double PI = 3.14159265358979300;
    double ACCUR = 1.0e-14;
    int ODDSW;
    int ASYMSW;

    double ROI, ROF, ZI, ZF, PHSI, PHSF, DZI, DZF;
    double FI, FPI, GI, GPI, FF, FPF, GF, GPF;
    double RHO0I, RHO0F;
    double TI, TF, CHI, B, F, A0, R, EPS, DB;
    double DELTA, RVAL, RTURN;
    double CI, SI, ZFAC, RI, RF, WI, WF, DI, DJ, TERM, TERM2;
    double CUMC, CUMS, CIP, SIP;
    double DCHI, D1, C1, DF, DTI, DTF, D2ZI, D2ZF, D2CHI, D12, D2, D2F;
    double TIMS, PHANG, SN, CSN;
    double FACTOR_cl, EVEN, ODD, SIGN, AM1, AM2, A, BM1;
    double C2;
    double ARAT, CILAST, SILAST;
    double EPSILO, EPSIL;
    int ISYN, IPIECE, ICHI;
    int NMAX, NASYI, NASYF;
    int I, J, JJ, NIT;

    if ( IPRINT >= 3 )  std::printf("\n CLINTS:%3d%4d%4d%3d%3d%15.6g%15.6g%15.6g%15.6g%15.6g%15.6g%15.6g\n ACCURA=%15.5g\n",
      N, LI, LF, NTERMS, NPTS, RLOWER, ETAI, ETAF, FKI, FKF, SIGI, SIGF, ACCURA_arg);

    //
    // SETUP GAUSS POINTS
    //
    GAUSSL( NPTS, PTS, WTS );

    ICHI = 0;
    EPSILO = ASMULT * ACCURA_arg;
    EPSIL = EPSILO;

    //
    // INITIALIZE BELLING METHOD FOR SUM
    //
    ISYN = 1;
    FFINT = 0;
    FGINT = 0;
    GFINT = 0;
    GGINT = 0;
    DELTA = PI / DMAX1(FKI, FKF);
    IPIECE = 0;
    RVAL = RLOWER;

    //
    // GET TURNING POINT
    //
    RTURN = DMAX1( (ETAI + DSQRT(ETAI*ETAI + LI*(LI+1))) / FKI,
      (ETAF + DSQRT(ETAF*ETAF + LF*(LF+1))) / FKF );
    ASYMSW = FALSE_F;
    if ( IPRINT >= 3 )  std::printf(" RTURN, EPSILO, DELTA =%15.6g%15.6g%15.6g\n",
      RTURN, EPSILO, DELTA);

    //
    // GET COULOMB WAVE FUNCTION AND BELLING PARAMETERS AT R = RVAL
    //
L100:
    ROI = FKI * RVAL;
    ROF = FKF * RVAL;

    if ( ASYMSW )  goto L120;
    if ( RVAL < 1.2*RTURN )  goto L110;

    //
    // GET ASYMPTOTIC EXPANSION PARAMETERS
    //
    NMAX = MAX0( LI, 4*NTERMS );
    RCASYM( LI, ETAI, ROI, IPRINT-4, SIGI, &ZI, &PHSI, &DZI,
      &FI, &FPI, &GI, &GPI, FA, GA, GPA, WORK,
      ACCURA_arg, NMAX, NASYI, IRET );
    if ( IRET != 0 )  goto L110;
    RCASYM( LF, ETAF, ROF, IPRINT-4, SIGF, &ZF, &PHSF, &DZF,
      &FF, &FPF, &GF, &GPF, FPA, GA, GPA, WORK,
      ACCURA_arg, NMAX, NASYF, IRET );
    if ( IRET != 0 )  goto L110;

    //
    // FA NOW HAS EXPANSION PARAMETERS FOR INCOMING
    // FPA HAS PARAMETERS FOR FINAL.
    //
    RHO0I = ROI;
    RHO0F = ROF;
    ASYMSW = TRUE_F;
    goto L130;

L110:
    RCWFN(ROI, ETAI, LI, LI, FA, FPA, GA, GPA, ACCUR, IRET);
    if (IRET != 0)  goto L900;
    FI = FA[LI+1];
    GI = GA[LI+1];
    FPI = FPA[LI+1];
    GPI = GPA[LI+1];
    RCWFN(ROF, ETAF, LF, LF, FA, FPA, GA, GPA, ACCUR, IRET);
    if (IRET != 0) goto L900;
    FF = FA[LF+1];
    GF = GA[LF+1];
    FPF = FPA[LF+1];
    GPF = GPA[LF+1];
    ZI = 1 / (FI*FI + GI*GI);
    ZF = 1 / (FF*FF + GF*GF);
    goto L130;

    //
    // FIND STUFF FROM THE ASYMPTOTIC EXPANSION
    //
L120:
    RI = RHO0I / ROI;
    RF = RHO0F / ROF;
    ZI = 1;
    ZF = 1;
    DZI = 0;
    DZF = 0;
    WI = RI;
    WF = RF;
    PHSI = 0;
    PHSF = 0;

    DI = 1;
    TERM2 = 1;
    for (I = 2; I <= NASYI; I++) {
        TERM = WI * FA[I];
        WI = WI * RI;
        ZI = ZI + TERM;
        DZI = DZI - DI * TERM;
        if ( I > 2 )  PHSI = PHSI + TERM / (DI - 1);
        DI = DI + 1;
        if ( DABS(TERM) + TERM2 < ACCURA_arg * ZI )  goto L125;
        TERM2 = DABS(TERM);
    }
    I = NASYI;
L125:
    NASYI = I;
    PHSI = ROI - ETAI*DLOG(2*ROI) + SIGI - .5*PI*LI - ROI*PHSI;
    DZI = DZI / ROI;

    DI = 1;
    TERM2 = 1;
    for (I = 2; I <= NASYF; I++) {
        TERM = WF * FPA[I];
        WF = WF * RF;
        ZF = ZF + TERM;
        DZF = DZF - DI * TERM;
        if ( I > 2 )  PHSF = PHSF + TERM / (DI - 1);
        DI = DI + 1;
        if ( DABS(TERM) + TERM2 < ACCURA_arg * ZF )  goto L128;
        TERM2 = DABS(TERM);
    }
    I = NASYF;
L128:
    NASYF = I;
    PHSF = ROF - ETAF*DLOG(2*ROF) + SIGF - .5*PI*LF - ROF*PHSF;
    DZF = DZF / ROF;
    FI = DSIN(PHSI) / DSQRT(ZI);
    GI = DCOS(PHSI) / DSQRT(ZI);
    FPI = ZI*GI - (.5*DZI/ZI)*FI;
    GPI = -ZI*FI - (.5*DZI/ZI)*GI;
    FF = DSIN(PHSF) / DSQRT(ZF);
    GF = DCOS(PHSF) / DSQRT(ZF);
    FPF = ZF*GF - (.5*DZF/ZF)*FF;
    GPF = -ZF*FF - (.5*DZF/ZF)*GF;
L130:
    TI = FI*FPI + GI*GPI;
    TF = FF*FPF + GF*GPF;
L140:
    CHI = FKI*ZI + ISYN*FKF*ZF;

    //
    // FOR VERY SMALL CHI (DIFFERENCE PIECE FOR INITIAL = FINAL CASE)
    // WE DO THE WHOLE THING NUMERICALLY
    //
    B = 1.0e-8 * FKI * ZI;
    if ( IPRINT >= 4 )  std::printf(" R, Z's, CHI, TINY, T's:%12.4g%20.10g%20.10g%14.4g%12.4g%14.4g%12.4g\n",
      RVAL, ZI, ZF, CHI, B, TI, TF);

    if ( DABS(CHI) < B )  goto L150;
    if ( DABS(GI) > 1000 )  goto L200;

    F = 1 / ( DSQRT(ZI*ZF) * std::pow(RVAL, N) );
    A0 = F / CHI;

    //
    // Precheck: can we use Belling at this RVAL?
    // Same formula as Fortran (source.f L8855-8860).
    //
    R = ( N/RVAL + 2*(FKI*ZI*DABS(TI) + FKF*ZF*DABS(TF)) ) / DABS(CHI);
    EPS = EPSIL;
    if ( ISYN == -1 )  EPS = EPS * ( 1 + DABS(CUMC/A0) );
    DB = R*R*R;
    if ( DB*DB <= EPS )  goto L300;

    //
    // Not yet — integrate one cycle numerically, then retry at larger RVAL.
    //
L150:
    if ( ISYN == -1 )  DELTA = PI / (DABS(CHI) + N/RVAL);
    goto L200;

L200:
    IPIECE = IPIECE + 1;
    CI = 0;
    SI = 0;

    for (I = 1; I <= NPTS; I++) {
        R = RVAL + .5*DELTA + .5*DELTA*PTS[I];
        ROI = FKI * R;
        ROF = FKF * R;
        ZFAC = .5*DELTA*WTS[I] / std::pow(R, N);
        if ( ASYMSW )  goto L210;
        RCWFN(ROI, ETAI, LI, LI, FA, FPA, GA, GPA, ACCUR, IRET);
        if (IRET != 0) goto L900;
        FI = FA[LI+1];
        GI = GA[LI+1];
        RCWFN(ROF, ETAF, LF, LF, FA, FPA, GA, GPA, ACCUR, IRET);
        if (IRET != 0) goto L900;
        FF = FA[LF+1];
        GF = GA[LF+1];
        if ( ISYN == -1 )  goto L250;
        goto L230;

        //
        // HERE USE ASYMPTOTIC EXPANSION
        //
L210:
        RI = RHO0I / ROI;
        ZI = FA[NASYI] * RI;
        PHSI = 0;
        DJ = NASYI - 2;
        for (JJ = 3; JJ <= NASYI; JJ++) {
            J = NASYI + 2 - JJ;
            ZI = RI * (FA[J] + ZI);
            PHSI = RI * (FA[J+1]/DJ + PHSI);
            DJ = DJ - 1;
        }
        ZI = 1 + ZI;
        PHSI = ROI - ETAI*DLOG(2*ROI) + SIGI - .5*PI*LI - RHO0I*PHSI;

        RF = RHO0F / ROF;
        ZF = FPA[NASYF] * RF;
        PHSF = 0;
        DJ = NASYF - 2;
        for (JJ = 3; JJ <= NASYF; JJ++) {
            J = NASYF + 2 - JJ;
            ZF = RF * (FPA[J] + ZF);
            PHSF = RF * (FPA[J+1]/DJ + PHSF);
            DJ = DJ - 1;
        }
        ZF = 1 + ZF;
        PHSF = ROF - ETAF*DLOG(2*ROF) + SIGF - .5*PI*LF - RHO0F*PHSF;

        if ( ISYN == -1 )  goto L255;

        FI = DSIN(PHSI) / DSQRT(ZI);
        GI = DCOS(PHSI) / DSQRT(ZI);
        FF = DSIN(PHSF) / DSQRT(ZF);
        GF = DCOS(PHSF) / DSQRT(ZF);
L230:
        FFINT = FFINT + ZFAC*FF*FI;
        FGINT = FGINT + ZFAC*FF*GI;
        GFINT = GFINT + ZFAC*GF*FI;
        GGINT = GGINT + ZFAC*GF*GI;
        goto L259;

        //
        // FOLLOWING ARE COS(PHII-PHIF) AND SIN(PHII-PHIF)
        //
L250:
        CI = CI + ZFAC*(GI*GF + FI*FF);
        SI = SI + ZFAC*(FI*GF - GI*FF);
        goto L259;
L255:
        CI = CI + (ZFAC/sqrt(ZI*ZF)) * cos(PHSI - PHSF);
        SI = SI + (ZFAC/sqrt(ZI*ZF)) * sin(PHSI - PHSF);

L259:   ;
    }

    if ( IPRINT < 7 )  goto L280;
    if ( ISYN == +1 )  std::printf(" INTS:%4d%6.1f%14.6g%14.6g%14.6g%14.6g\n",
      IPIECE, RVAL, FFINT, FGINT, GFINT, GGINT);
    if ( ISYN == -1 )  std::printf(" NUMERIC C-, S-:%4d%12.4g%14.6g%14.6g\n",
      IPIECE, RVAL, CI, SI);

L280:
    RVAL = RVAL + DELTA;
    if ( ISYN == +1 )  goto L100;
    CUMC = CUMC + CI;
    CUMS = CUMS + SI;

    //
    // HAVE WE GONE SO FAR IN DIFFERENCE THAT THE BELLING'S
    // PIECE WON'T MATTER.
    //
    if ( DABS(CI) < ACCURA_arg*DABS(FFINT + .5*(CUMC - CIP))  &&
         DABS(SI) < ACCURA_arg*DABS(GFINT + .5*(CUMS + SIP)) )  goto L810;
    goto L100;

    //
    // ALL DONE WITH NUMERIC PART OF TOTALS
    //
L300:
    if ( IPRINT < 3  ||  IPIECE == 0 )  goto L350;
    if ( ISYN == +1 )  std::printf("\n NUMERIC VALUES FOR%4d CYCLES:\n%15.8g%15.8g%15.8g%15.8g\n",
      IPIECE, FFINT, FGINT, GFINT, GGINT);
    if ( ISYN == -1 )  std::printf(" NUMERIC C-, S- FROM%4d CYCLES:%15.8g%15.8g\n",
      IPIECE, CUMC, CUMS);

    //
    // FOR SUM WE USE THE FF INTEGRAL AS A GUIDE TO ERRORS
    //
L350:
    if ( ISYN == -1 )  goto L400;
    CUMC = FFINT;
    CUMS = FFINT;

    //
    // BELLINGS METHOD FOR SUM OR DIF
    //
L400:
    PHSI = DATAN2(FI, GI);
    PHSF = DATAN2(FF, GF);
    DZI = -2*FKI*ZI*ZI*TI;
    DZF = -2*FKF*ZF*ZF*TF;
    DCHI = FKI*DZI + ISYN*FKF*DZF;
    D1 = DZI/ZI + DZF/ZF;
    C1 = DFLOAT(N)/RVAL + .5*D1;
    DF = -F*C1;
    WI = 1 - (2*ETAI)/ROI - (LI*(LI+1))/(ROI*ROI);
    WF = 1 - (2*ETAF)/ROF - (LF*(LF+1))/(ROF*ROF);
    DTI = FKI*(FPI*FPI + GPI*GPI - WI/ZI);
    DTF = FKF*(FPF*FPF + GPF*GPF - WF/ZF);
    D2ZI = -2*FKI*ZI*(2*DZI*TI + ZI*DTI);
    D2ZF = -2*FKF*ZF*(2*DZF*TF + ZF*DTF);
    D2CHI = FKI*D2ZI + ISYN*FKF*D2ZF;
    D12 = (DZI/ZI)*(DZI/ZI) + (DZF/ZF)*(DZF/ZF);
    D2 = D2ZI/ZI + D2ZF/ZF;
    D2F = (F/RVAL - DF)*C1 - (F*D1)/(2*RVAL) - .5*F*(D2 - D12);
    if ( IPRINT >= 4 )  std::printf("\n R =%12.4g   RHO's =%12.4g%12.4g   EPSILON =%10.2g\n",
      RVAL, ROI, ROF, R);
    DB = 0;
    TIMS = 1 / (RVAL * CHI);
    PHANG = PHSI + ISYN*PHSF;
    SN = DSIN(PHANG);
    CSN = DCOS(PHANG);

    //
    // EACH ITERATION ON NITS BRINGS IN A HIGHER DERIVATIVE OF BETA.
    //
    if ( IPRINT >= 6 )  std::printf("\n NIT  K%6s B%9s DB%10s A%13s EVEN%13s ODD\n", "", "", "", "", "");

    for (NIT = 1; NIT <= 2; NIT++) {

    FACTOR_cl = F / CHI;

    //
    // BELLING'S EXPANSION
    //
        EVEN = A0;
        AM1 = A0;
        B = 1;
        ODD = 0;
        SIGN = +1.0;
        ODDSW = TRUE_F;
        for (I = 1; I <= NTERMS; I++) {
            FACTOR_cl = FACTOR_cl * TIMS;
            C2 = DF/F - (I*DCHI)/CHI;
            BM1 = B;
            B = RVAL*DB + (C2*RVAL - I + 1)*B;
            if ( NIT > 1 )
              DB = (C2*RVAL - I + 2)*DB + (C2 - RVAL*((DF/F)*(DF/F)
                - I*(DCHI/CHI)*(DCHI/CHI)) + RVAL*(D2F/F - I*D2CHI/CHI))*BM1;
            A = B * FACTOR_cl;

            if ( ODDSW )  goto L440;
            EVEN = EVEN + SIGN*A;
            goto L445;
L440:
            ODD = ODD + SIGN*A;
            SIGN = -SIGN;
L445:
            ODDSW = !ODDSW;

            if ( IPRINT >= 6 )  std::printf("%3d%4d%11.4g%11.4g%11.4g%17.9g%17.9g\n",
              NIT, I, B, DB, A, EVEN, ODD);
            if ( I == 1 )  goto L465;
            if ( DABS(A) > DABS(AM2) )  goto L700;
            if ( DABS(A) < ACCURA_arg*(DABS(A0) + DMIN1(DABS(CUMC),
                    DABS(CUMS))) )   goto L470;

L465:
            AM2 = AM1;
            AM1 = A;
        }
        goto L700;

L470:
        CI = -EVEN*SN - ODD*CSN;
        SI = -ODD*SN + EVEN*CSN;
        ARAT = A / A0;
        if ( IPRINT >= 4 )  std::printf(" NIT =%2d   A(%2d)/A(0) =%9.2g   COS, SIN =%15.8g%15.8g\n",
          NIT, I, ARAT, CI, SI);
        if ( NIT == 1 )  goto L490;
        EVEN = DABS( (CILAST - CI) / (CUMC + CI) );
        ODD  = DABS( (SILAST - SI) / (CUMS + SI) );
        if ( IPRINT >= 4 )  std::printf("           REL. ERR. =%10.2g%10.2g\n\n", EVEN, ODD);
L490:
        CILAST = CI;
        SILAST = SI;
    }

    //
    // WAS THE DERIVATIVE A BIG CORRECTION?
    //
    if ( DMAX1(EVEN, ODD)*DMAX1(EVEN, ODD) > ACCURA_arg )  goto L700;

    ICHI = 0;
    if (ISYN == -1) goto L800;
    CIP = CI;
    SIP = SI;

    //
    // NOW DO  PH(IN)-PH(OUT)
    //
    IPIECE = 0;
    ISYN = -1;
    CUMS = 0;
    CUMC = 0;
    EPSIL = EPSILO;
    goto L140;

    //
    // WE NEED TO MAKE ASMULT SMALLER
    //
L700:
    R = DMIN1( .5*EPSIL, .75*std::pow(R, 6)*(EPSIL/EPS) );
    ASMULT = ASMULT * R / EPSIL;
    EPSIL = R;
    ICHI = ICHI + 1;
    if ( ICHI > 20 )  goto L720;
    if ( IPRINT >= 4 )  std::printf("\n +++ REDOING CALCULATION WITH NEW ASMULT:%12.4g%12.4g\n",
      ASMULT, EPSIL);
    goto L140;

L720:
    std::printf("\n *** COULD NOT GET CONVERGENCE IN CLINTS *** \n");
    IRET = -1;
    goto L50;

    //
    // ALL DONE, COLLECT PIECES AND LEAVE
    //
L800:
    CUMC = CUMC + CI;
    CUMS = CUMS + SI;
L810:
    FFINT = FFINT + .5*( CUMC - CIP );
    GGINT = GGINT + .5*( CUMC + CIP );
    FGINT = FGINT + .5*( SIP - CUMS );
    GFINT = GFINT + .5*( SIP + CUMS );
    if ( IPRINT >= 3 )  std::printf("\n C+ S+ C- S-:%15.8g%15.8g%15.8g%15.8g\n FF FG GF GG:%15.8g%15.8g%15.8g%15.8g\n",
      CIP, SIP, CUMC, CUMS, FFINT, FGINT, GFINT, GGINT);
    goto L50;

    //
    // ERROR IN RCWFN
    //
L900:
    std::printf("\n *** ERROR IRET = %5d IN RCWFN ***\n"
      " L's, ETA's, RHO's =%5d%5d%20.10g%20.10g%20.10g%20.10g\n",
      IRET, LI, LF, ETAI, ETAF, ROI, ROF);

L50:
    return;
}


// ============================================================================
// JPTOLX  (source.f lines 23065-23147)
//
// Translates elastic S-matrix from (L',L,JP) to (L',L,LX) system.
// ============================================================================

void JPTOLX(int L, int LAS, int JP, int NWP, double* SJR_p, double* SJI_p,
            int* INDX, double* SLX)
{
    int I, K, LX, LX2;
    int LXMN, LXMX;
    double COEF, T;

    //
    // IF NO SPIN-ORBIT, JUST TRANSFER IT.
    //
    double SJR = *SJR_p;
    double SJI = *SJI_p;

    if ( WAVCOM.SOSWS[NWP] )  goto L100;
    I = 2*L + 1;
    SLX[I] = SJR;
    SLX[I+1] = SJI;
    return;

    //
    // SPIN-ORBIT EXISTS.  MULTIPLY BY RACAH'S ETC., AND ADD
    // TO EXISTING VALUES FOR EACH LX.
    //
L100:
    LXMN = MAX0( 0, IABS( L - LAS ) );
    LXMX = MIN0( WAVCOM.JSPS[NWP], L + LAS );
    COEF = (WAVCOM.JSPS[NWP] + 1) * (2*LAS + 1);
    COEF = (JP + 1) / DSQRT(COEF);
    I = ( WAVCOM.JSPS[NWP] - JP + L + LAS ) / 2;
    if ( MOD(I, 2) != 0 )  COEF = -COEF;
    for (LX = LXMN; LX <= LXMX; LX++) {
        LX2 = 2*LX;
        T = LX2 + 1;
        T = COEF * DSQRT(T) * RACAH( 2*L, 2*LAS, WAVCOM.JSPS[NWP], WAVCOM.JSPS[NWP],
              LX2, JP );
        K = 1 + LX*( WAVCOM.JSPS[NWP] + 2 );
        // INDX is 0-based pointer (from &ILLOC), accessed as Fortran INDX(3,*)
        // Fortran INDX(I,K) → C++ INDX[(K-1)*3 + (I-1)]
        if ( INDX[(K-1)*3 + 0] <= 0 )  goto L199;
        I = LAS*WAVCOM.NUMJS[NWP] + INDX[(K-1)*3 + 0] + (L - LAS - INDX[(K-1)*3 + 1]) / 2;
        I = 2*I - 1;
        SLX[I] = SLX[I] + T*SJR;
        SLX[I+1] = SLX[I+1] + T*SJI;
L199:   ;
    }
    return;
}
