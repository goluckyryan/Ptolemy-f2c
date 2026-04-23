// grdset_translated.cpp — GRDSET subroutine
// Translated from source.f lines 18212-19535
// Sets up integration grids for DWBA transfer integrals

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include "ptolemy_alloc.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_forward.h"
#include <cstdio>
#include <cmath>
#include <algorithm>

extern double second();
extern void CUBMAP(int MAPTYPE, double XMIN, double XMID, double XMAX,
                   double GAMMA, double* PTS, double* WTS, int NPTS);
extern void LSQPOL(double* X, double* Y, double* W, double* RESID, int NPTS,
                   double* SMSQ, int NCOL, double* AMAT, double* COEFS,
                   int NTERM, int NROW, int NCOLT);

void GRDSET(int& IRTN) {

    // =========================================================================
    // Local variable declarations
    // =========================================================================

    // LOGICAL variables
    int PAS1SW;   // LOGICAL
    int ZEROSW;   // LOGICAL
    int PRNTSW;   // LOGICAL
    int PGRDSW;   // LOGICAL
    int NPLYSW;   // LOGICAL — .TRUE. means V-polynomials are NOT to be used

    // Local scalars
    double TSTART, TINTRP;
    int    NINTRP;
    int    IPRNT;
    double TEMP, TIME1;
    double FIFO, RP, RT;
    double U, UMAX, VVAL, VMAX, VMIN, VMID;
    double D, SYNE;
    double WVWMAX, RVRLIM;
    double USTEP, ULIM;
    double SUM0, SUM1, SUM2;
    double PVPMAX;
    double VOFMAX;
    double WOW;
    int    INADD, IRECT;
    int    N, I, II, IV, IU;
    int    IPLUNK, IMAX, IEND;
    int    JCNT;
    int    JSIZE;
    int    IBSTYP, NNAIT, IXTOPZ;
    int    NBUMPS;
    int    N4RIO, NRIROH, NRIROI, NPOLSV, NVTERM;
    int    INEED;
    double DV, DXV;
    double X, X0, X0MIN, X0AV, PHI0, PHI, DPHI;
    double RI, RO, AREA;
    double RTMAX, RPMAX;
    double VRNGMN, VRNGMX, VMAXMX;
    double VRANGE[4];   // 1-based, [1..3]

    // Allocator pointer locals (return values from NALLOC/NAMLOC)
    int    ISMHPT, ISMHWK, ISMIPT, ISMIVL;
    int    IWIO_loc, IRI_loc, IRO_loc, IDIFPT, IDIFWT;
    int    IVPOLS, IRIOEX, IRIH, IROH, IVWORK;
    int    IVWORK_loc;
    int    LOGIC_loc, IPHIPT, IPHIWT;
    int    LOCTRP_loc, IPHIT_loc, IPHIP_loc, IPHI_loc;

    // Offset-based pointer aliases into the allocator (computed once)
    int    LVPOLS, LVWTS, LVRESD, LVSMSQ, LVAMAT;
    int    LVMIN, LVMID, LVMAX;
    int    LRIH, LROH, LRIOEX;
    int    LDIFPT, LDIFWT, LSUMPT, LSMIVL;
    int    LWIO, LRI, LRO;
    int    LLOGIC;
    int    LPHIPT, LPHIWT;
    int    LPHIT, LPHIP, LPHI, LTRP;
    int    ITEMP;

    // Polynomial/interpolation counts
    int    NRIROH_h;   // renamed from NRIROH for H-search (= NPDIF*NPSUM)
    double VRNGMN_sv, VRNGMX_sv, X0MIN_sv, X0AV_sv;

    // DATA arrays
    // XS(2) — Fortran DATA XS / 1.D0, 0.D0 /
    // 1-based: XS[1]=1.0, XS[2]=0.0
    double XS[3] = { 0.0, 1.0, 0.0 };   // index 0 unused; XS[1]=1, XS[2]=0

    // VS(5) — work array for 5 test V values
    double VS[6] = {};   // 1-based

    // Convenient COMMON references
    double& UNDEF      = INTRNL.UNDEF;
    int&    ISTRIP     = INTRNL.ISTRIP;
    int&    IHSAVE     = INTRNL.IHSAVE;

    double& SUMMAX     = FLOAT_common.SUMMAX;
    double& SUMMID     = FLOAT_common.SUMMID;
    double& SUMMIN     = FLOAT_common.SUMMIN;
    double& GAMSUM     = FLOAT_common.GAMSUM;
    double& GAMDIF     = FLOAT_common.GAMDIF;
    double& DWCUT      = FLOAT_common.DWCUT;
    double& AMDMLT     = FLOAT_common.AMDMLT;
    double& PHIMID     = FLOAT_common.PHIMID;
    double& GAMPHI     = FLOAT_common.GAMPHI;
    double& ALPHAP     = BSBLK.ALPHAP;    // /BSBLK/ asymptotic decay constant for projectile
    double& ALPHAT     = BSBLK.ALPHAT;    // /BSBLK/ asymptotic decay constant for target

    int&    IPRINT     = INTGER.IPRINT;
    int&    NPSUM      = INTGER.NPSUM;
    int&    NPDIF      = INTGER.NPDIF;
    int&    NPPHI      = INTGER.NPPHI;
    int&    NPHIAD     = INTGER.NPHIAD;
    int&    MAPSUM     = INTGER.MAPSUM;
    int&    MAPDIF     = INTGER.MAPDIF;
    int&    LOOKST     = INTGER.LOOKST;
    int&    NVPOLY     = INTGER.NVPOLY;
    int&    NAITKN     = INTGER.NAITKN;
    int&    MAPPHI     = INTGER.MAPPHI;

    int&    NUCONL     = SWITCH.NUCONL;
    int&    ISCRFL     = SWITCH.ISCRFL;
    int&    IBSPAS     = SWITCH.IBSPAS;
    int&    ISW3       = SWITCH.ISW3;

    int&    ISCTMN     = KANDM.ISCTMN;
    int&    ISCTCR     = KANDM.ISCTCR;
    double& RIMAX_k    = KANDM.RIMAX;
    double& ROMAX_k    = KANDM.ROMAX;

    // RLMAXS(1)=RLPMAX (projectile bound-state max R),
    // RLMAXS(2)=RLTMAX (target bound-state max R)
    // inline double* RLMAXS() { return &BSBLK.RLPMAX; }
    // RLMAXS[1]=BSBLK.RLPMAX, RLMAXS[2]=BSBLK.RLTMAX
    double& RLMAXS1    = BSBLK.RLPMAX;
    double& RLMAXS2    = BSBLK.RLTMAX;
    double& BNDMXP     = BSBLK.BNDMXP;
    double& BNDMXT     = BSBLK.BNDMXT;
    double& ROFMAX     = BSBLK.ROFMAX;
    int&    PINFSW     = BSBLK.PINFSW;
    float&  PHISGN     = BSBLK.PHISGN;
    double& S1         = BSBLK.S1;
    double& S2         = BSBLK.S2;
    double& T1         = BSBLK.T1;
    double& T2         = BSBLK.T2;

    double& JACOB      = GRIDCM.JACOB;
    int&    NPSUMI     = GRIDCM.NPSUMI;
    double& RVRLIM_g   = GRIDCM.RVRLIM;
    int&    ISMHPT_g   = GRIDCM.ISMHPT;
    int&    ISMIPT_g   = GRIDCM.ISMIPT;
    int&    ISMIVL_g   = GRIDCM.ISMIVL;
    int&    ISMHWK_g   = GRIDCM.ISMHWK;
    int&    IRIOEX_g   = GRIDCM.IRIOEX;
    int&    IRI_g      = GRIDCM.IRI;
    int&    IRO_g      = GRIDCM.IRO;
    int&    IWIO_g     = GRIDCM.IWIO;
    int&    NRIROI_g   = GRIDCM.NRIROI;
    int&    MAXCNT     = GRIDCM.MAXCNT;
    int&    NUMHS      = GRIDCM.NUMHS;
    int&    LOGIC_g    = GRIDCM.LOGIC;
    int&    LOCTRP_g   = GRIDCM.LOCTRP;
    int&    IPHIT_g    = GRIDCM.IPHIT;
    int&    IPHIP_g    = GRIDCM.IPHIP;
    int&    IPHI_g     = GRIDCM.IPHI;
    int&    IVWORK_g   = WAVCOM.IVWORK;

    // WAVCOM aliases
    int*    NSTP2S     = WAVCOM.NSTP2S;   // 1-based
    double* RSTEPS     = WAVCOM.RSTEPS;   // 1-based
    int&    NWAVEF     = WAVCOM.NWAVEF;
    int&    IWAVR      = WAVCOM.IWAVR;
    int&    IWAVI      = WAVCOM.IWAVI;
    int*    TCSWS      = WAVCOM.TCSWS;    // 1-based, LOGICAL

    // FORMF
    int*    IBDS       = FORMF.IBDS;      // 1-based
    double* BRATMS     = FORMF.BRATMS;    // 1-based

    // CNSTNT
    double& PI         = CNSTNT.PI;
    double& BIGNUM     = CNSTNT.BIGNUM;

    // KANDM
    double* AKS_ptr    = &KANDM.AKI;     // AKS(1)=AKI, AKS(2)=AKO

    // INTGER for JPOT (from BSBLK)
    int&    JPOT       = BSBLK.JPOT;

    // ALLOCS facfr4
    int&    FACFR4     = ALLOCS.FACFR4;

    // ST1, ST2 — local temporaries for 2*S1*T1, 2*S2*T2
    double ST1, ST2;

    // =========================================================================
    //
    //     TSTART = second()
    //     TINTRP = 0
    //     NINTRP = 0
    //     IRTN = 0
    //
    // =========================================================================
    TSTART  = second();
    TINTRP  = 0.0;
    NINTRP  = 0;
    IRTN    = 0;

    //      IPRNT = MOD( IPRINT, 10 )
    //      PGRDSW = IPRNT .GE. 5
    //      PINFSW = IPRNT .GE. 4
    //      PRNTSW = IPRNT .GE. 1
    IPRNT  = IPRINT % 10;
    PGRDSW = (IPRNT >= 5);
    PINFSW = (IPRNT >= 4);
    PRNTSW = (IPRNT >= 1);

    //      IF ( PINFSW )  write (6, 3)
    //  3   FORMAT ( / '0SETUP OF THE INTEGRATION GRIDS:' )
    // FORMAT '0' = blank line before content
    if (PINFSW) {
        std::printf("\n\n0SETUP OF THE INTEGRATION GRIDS:\n");
    }

    //  MAKE THINGS EASY FOR INELDC WHEN  USEHS  IS IN EFFECT BY
    //  ALLOWING DUMMY  Z(I...)'S.
    //
    //      INTEGER  IZDUMY(45)
    //      EQUIVALENCE ( IZDUMY(1), IMSVAL )
    //
    //      DO 39  I = 1, 45
    // 39   IZDUMY(I) = 1
    //
    // IZDUMY(1) is EQUIVALENCED to IMSVAL which is the first int in GRIDCM.
    // The 45 ints starting at GRIDCM.IMSVAL are set to 1.
    {
        int* IZDUMY = &GRIDCM.IMSVAL;   // EQUIVALENCE ( IZDUMY(1), IMSVAL )
        for (int i = 0; i < 45; i++) {
            IZDUMY[i] = 1;
        }
    }

    // =========================================================================
    //  RESTORE STUFF FROM H'S FILE IF NEEDED
    //
    //      IF ( IHSAVE .NE. 2 )  GO TO 50
    //      READ ( 1, END=9950, ERR=9950 )  NPSUM, NPDIF, NPPHI,
    //   2    MAXCNT, JSIZE, GAMSUM, SUMMIN, SUMMID, SUMMAX,
    //   3    GAMDIF, VRNGMN, VRNGMX, X0MIN, X0AV, RPMAX, RTMAX,
    //   4    NRIROH, ALPHAP, ALPHAT, AREA, WVWMAX, NPLYSW
    // =========================================================================
    if (IHSAVE == 2) {
        // Read from binary H's file (Fortran unit 1 = ichanb or channel 1)
        // In C++ we use fread on the appropriate file handle.
        // For now emit a placeholder call using the same pattern as other translated files.
        // The Fortran reads a single unformatted record with all these items.
        FILE* hfile = nullptr;
        // Attempt to read — if EOF/ERR jump to L9950
        // We replicate the Fortran logic: on error go to 9950.
        // Since actual file I/O infrastructure is handled elsewhere in the C++ port,
        // we call RDHSFL (not yet translated) or inline the read here.
        // Marking as a direct translation of the Fortran binary read:
        //   (items read into local/COMMON variables)
        if (hfile == nullptr) goto L9950;
        // Actual read would be:
        // fread(&NPSUM, sizeof(int), 1, hfile); etc.
        // For now this block is a structural placeholder.
        // The detailed file I/O translation is deferred.
        goto L9950;   // error path until implemented
    }

L50:
    //  50   INADD=2
    //       IF(NUCONL.GE.2) INADD=3
    INADD = 2;
    if (NUCONL >= 2) INADD = 3;

    // ****  TEMPORARY  ****
    //       IRECT = ISW3
    IRECT = ISW3;
    // ****  END OF TEMPORARIES  ****

    // =========================================================================
    //  THE CODE IS BASED ON THE ASSUMPTION THAT THE INCOMING PROJECTILE
    //  IS COMPOSITE ( A = B+X ).  ...
    //
    //      TEMP = 1 / ( BRATMS(1) + BRATMS(2) * (1+BRATMS(1)) )
    //      S1 = ( 1 + BRATMS(1) ) * ( (1 + BRATMS(2) ) * TEMP )
    //      T1 = - ( 1 + BRATMS(2) ) * TEMP
    //      S2 = ( 1 + BRATMS(1) ) * TEMP
    //      T2 = -S1
    //      JACOB = S1**3
    //      PHISGN = 1
    //      IF ( ISTRIP .EQ. 1 )  GO TO 170
    // =========================================================================
    TEMP = 1.0 / (BRATMS[1] + BRATMS[2] * (1.0 + BRATMS[1]));
    S1   = (1.0 + BRATMS[1]) * ((1.0 + BRATMS[2]) * TEMP);
    T1   = -(1.0 + BRATMS[2]) * TEMP;
    S2   = (1.0 + BRATMS[1]) * TEMP;
    T2   = -S1;
    JACOB = S1 * S1 * S1;
    PHISGN = 1.0f;
    if (ISTRIP == 1) goto L170;

    // =========================================================================
    //  B = A+X, FLIP THINGS
    //  IN THE FOLLOWING NOTE THAT BRATMS(1) NOW HAS X/A
    //
    //      T2 = S2
    //      S2 = -S1
    //      S1 = T1
    //      T1 = -S2
    //      JACOB = T1**3
    //      PHISGN = -1
    // =========================================================================
    T2 =  S2;
    S2 = -S1;
    S1 =  T1;
    T1 = -S2;
    JACOB  = T1 * T1 * T1;
    PHISGN = -1.0f;

L170:
    // 170 IF (PINFSW)  write (6, 173) S1, T1, S2, T2, JACOB
    // 173  FORMAT ('0R TO RX MAPPING PARAMETERS (S1, T1, S2, T2) =',
    //     1  4G15.5, 5X, 'JACOB =', G15.5 )
    if (PINFSW) {
        std::printf("\n0R TO RX MAPPING PARAMETERS (S1, T1, S2, T2) =%15.5G%15.5G%15.5G%15.5G     JACOB =%15.5G\n",
                    S1, T1, S2, T2, JACOB);
    }

    //      ST1=2.D0*S1*T1
    //      ST2=2.D0*S2*T2
    ST1 = 2.0 * S1 * T1;
    ST2 = 2.0 * S2 * T2;

    //      IF ( IHSAVE .EQ. 2 )  GO TO 400
    if (IHSAVE == 2) goto L400;

    // =========================================================================
    //  INITIALIZE THE COMPUTATION OF  (PHI V PHI)
    //
    //      CALL BSSET
    // =========================================================================
    BSSET();

    //      IF(LOOKST.GT.254) LOOKST=254
    if (LOOKST > 254) LOOKST = 254;

    //       DXV=2.D0/LOOKST**2
    //      XS(2) = 1 - DXV
    DXV   = 2.0 / ((double)LOOKST * (double)LOOKST);
    XS[2] = 1.0 - DXV;

    // =========================================================================
    //  FIND THE MAXIMUM  ( R PSI PHI V PHI PSI R )   THAT WE WILL
    //  ENCOUNTER AND USE THIS TO SCALE DWCUT
    // =========================================================================

    //      WVWMAX = 0
    //      U = .5*ROFMAX
    //      N = 1.5*ROFMAX/.20D0 + 1
    //      VS(3) = 0
    //      TIME1 = second()
    WVWMAX = 0.0;
    U      = 0.5 * ROFMAX;
    N      = (int)(1.5 * ROFMAX / 0.20 + 1.0);
    VS[3]  = 0.0;
    TIME1  = second();

    //      DO 269  IU = 1, N
    for (IU = 1; IU <= N; IU++) {
        //  WE LOOK AT 5 DIFFERENT V'S FOR EACH U
        //         D = 2*(RLMAXS(2)-(S1+T1)*U)/(S1-T1)
        //         IF ( DABS(D) .GT. 2*U )  D = DSIGN(2*U, D)
        //         VS(1) = D
        //         VS(2) = .5*D
        //         D = 2*(RLMAXS(1)-(S2+T2)*U)/(S2-T2)
        //         IF ( DABS(D) .GT. 2*U )  D = DSIGN(2*U, D)
        //         VS(4) = .5*D
        //         VS(5) = D
        D = 2.0 * (RLMAXS2 - (S1 + T1) * U) / (S1 - T1);
        if (std::fabs(D) > 2.0 * U) D = std::copysign(2.0 * U, D);
        VS[1] = D;
        VS[2] = 0.5 * D;
        D = 2.0 * (RLMAXS1 - (S2 + T2) * U) / (S2 - T2);
        if (std::fabs(D) > 2.0 * U) D = std::copysign(2.0 * U, D);
        VS[4] = 0.5 * D;
        VS[5] = D;

        //         DO 265  IV = 1, 5
        //            VVAL = VS(IV)
        //            DO 265  I = 1, 2
        //               CALL BSPROD ( FIFO, 3, U+.5*VVAL, U-.5*VVAL, XS(I),
        //     1             ISCTMN, 1, RP, RT, *264 )
        //               IF ( DABS(FIFO) .LE. WVWMAX )  GO TO 264
        //               WVWMAX = DABS(FIFO)
        //               UMAX = U
        //               VMAX = VVAL
        //               RPMAX = RP
        //               RTMAX = RT
        //  264        CONTINUE
        //  265     CONTINUE
        for (IV = 1; IV <= 5; IV++) {
            VVAL = VS[IV];
            for (I = 1; I <= 2; I++) {
                if (BSPROD(FIFO, 3, U + 0.5 * VVAL, U - 0.5 * VVAL, XS[I],
                       ISCTMN, 1, RP, RT) != 0) goto L264;
                if (std::fabs(FIFO) <= WVWMAX) goto L264;
                WVWMAX = std::fabs(FIFO);
                UMAX   = U;
                VMAX   = VVAL;
                RPMAX  = RP;
                RTMAX  = RT;
L264:           ;
            }
        }
        //         U = U + .20D0
        U = U + 0.20;
    }   // end DO 269

    TINTRP = TINTRP + second() - TIME1;
    NINTRP = NINTRP + 10 * N * (INADD + 2);

    //      IF ( PINFSW )  write (6, 273) UMAX, VMAX, RPMAX, RTMAX, WVWMAX
    // 273  FORMAT ( '0MAX(RI PSI PHI V PHI PSI RO) OCCURS AT',
    //     1    ' SUM, DIF, RP, RT =', 4F8.2, 5X, 'VALUE =', G13.3 )
    if (PINFSW) {
        std::printf("\n0MAX(RI PSI PHI V PHI PSI RO) OCCURS AT SUM, DIF, RP, RT =%8.2f%8.2f%8.2f%8.2f     VALUE =%13.3G\n",
                    UMAX, VMAX, RPMAX, RTMAX, WVWMAX);
    }

    //      RVRLIM = DWCUT*WVWMAX
    RVRLIM = DWCUT * WVWMAX;

    // =========================================================================
    //  STEP OUT FROM  U = 0  LOOKING FOR SUMMIN.
    //  WE LOOK FOR U SUCH THAT
    //    PSI(LMIN) PHI' V PHI' PSI(LMIN)
    //  EXCEEDS VZERO.
    // =========================================================================

    //      USTEP = .20D0
    //      U = 0
    //      TIME1 = second()
    USTEP  = 0.20;
    U      = 0.0;
    TIME1  = second();

L320:
    //  320     DO 324  I = 1, 2
    //             CALL BSPROD ( FIFO, 4, U, U, XS(I), ISCTMN, 1,
    //     1         RP, RT, *325 )
    //             IF ( DABS(FIFO) .GE. RVRLIM )  GO TO 330
    //  324     CONTINUE
    for (I = 1; I <= 2; I++) {
        if (BSPROD(FIFO, 4, U, U, XS[I], ISCTMN, 1, RP, RT) != 0) goto L325;
        if (std::fabs(FIFO) >= RVRLIM) goto L330;
    }
    //         NINTRP = NINTRP + 2*(INADD+2)
    //         U = U + USTEP
    //         GO TO 320
    NINTRP = NINTRP + 2 * (INADD + 2);
    U      = U + USTEP;
    goto L320;

L325:
    // 325  write (6, 326)
    // 326  FORMAT ('0**** ASYMPTOPIA IS FAR TOO SMALL - CANNOT FIND',
    //     1  ' SUMMIN.' )
    std::printf("\n0**** ASYMPTOPIA IS FAR TOO SMALL - CANNOT FIND SUMMIN.\n");
    U = U - USTEP;
    //      IF ( SUMMIN .NE. UNDEF )  GO TO 330
    if (SUMMIN != UNDEF) goto L330;
    IRTN = 0;
    return;

L330:
    // debug removed
    //  330  TINTRP = TINTRP + second() - TIME1
    //       U = DMAX1( 0.D0, U-USTEP )
    //       IF ( SUMMIN .EQ. UNDEF )  SUMMIN = U
    //       IF ( PINFSW )  write (6, 333) U, SUMMIN
    // 333  FORMAT ( '0SUMMIN:  PROPOSED =', G14.5, 5X, 'USED =', G14.5 )
    TINTRP = TINTRP + second() - TIME1;
    U = std::fmax(0.0, U - USTEP);
    if (SUMMIN == UNDEF) SUMMIN = U;
    if (PINFSW) {
        std::printf("\n0SUMMIN:  PROPOSED =%14.5G     USED =%14.5G\n", U, SUMMIN);
    }

    // =========================================================================
    //  NOW FIND SUMMAX USING
    //    PSI(L CRIT) PHI V PHI PSI(L CRIT)
    //  AND COMPUTE MOMENTS OF THE WEIGHT FUNCTION.
    // =========================================================================

    //      SUM0 = 0
    //      SUM1 = 0
    //      SUM2 = 0
    //      PVPMAX = 0
    SUM0   = 0.0;
    SUM1   = 0.0;
    SUM2   = 0.0;
    PVPMAX = 0.0;

    //      IF ( IPRNT .GE. 7 )  write (6, 347)
    // 347  FORMAT ( '0', 6X, 'U', 9X, 'SUM(F''S)', 6X, 'MAX(F)',
    //     1  5X, 'V OF MAX', 4X, 'RP OF MAX', 4X, 'RT OF MAX' / )
    if (IPRNT >= 7) {
        std::printf("\n0      U         SUM(F'S)      MAX(F)     V OF MAX    RP OF MAX    RT OF MAX\n\n");
    }

    //      U = SUMMIN
    //      ULIM = DMAX1( RSCTS(1)+10*ASCTS(1), ROFMAX+5 )
    U    = SUMMIN;
    ULIM = std::fmax(KANDM.RSCTS[1] + 10.0 * KANDM.ASCTS[1], ROFMAX + 5.0);
    TIME1 = second();

L350:
    //  350     TEMP = 0
    //          WVWMAX = 0
    //          ZEROSW = .TRUE.
    TEMP   = 0.0;
    WVWMAX = 0.0;
    ZEROSW = true;

    //  WE LOOK AT 5 DIFFERENT V'S FOR EACH U
    //  AT EACH U, V PAIR WE LOOK AT 2 X'S
    //  FOR THIS SEARCH, VIOLATIONS OF BOUND-STATE ASYMPTOPIA ARE
    //  IGNORED (TREATED AS ZERO) UNLESS ALL TEST POINTS ARE BAD.
    {
        //         D = 2*(RLMAXS(2)-(S1+T1)*U)/(S1-T1)
        //         IF ( DABS(D) .GT. 2*U )  D = DSIGN(2*U, D)
        //         VS(1) = D
        //         VS(2) = .5*D
        //         D = 2*(RLMAXS(1)-(S2+T2)*U)/(S2-T2)
        //         IF ( DABS(D) .GT. 2*U )  D = DSIGN(2*U, D)
        //         VS(4) = .5*D
        //         VS(5) = D
        D = 2.0 * (RLMAXS2 - (S1 + T1) * U) / (S1 - T1);
        if (std::fabs(D) > 2.0 * U) D = std::copysign(2.0 * U, D);
        VS[1] = D;
        VS[2] = 0.5 * D;
        D = 2.0 * (RLMAXS1 - (S2 + T2) * U) / (S2 - T2);
        if (std::fabs(D) > 2.0 * U) D = std::copysign(2.0 * U, D);
        VS[4] = 0.5 * D;
        VS[5] = D;

        //         DO 359  IV = 1, 5
        //            VVAL = VS(IV)
        //            DO 359  I = 1, 2
        //               CALL BSPROD ( FIFO, 3, U+.5*VVAL, U-.5*VVAL, XS(I),
        //     1             ISCTCR, 1, RP, RT, *359 )
        //               TEMP = TEMP + DABS(FIFO)
        //               IF ( WVWMAX .GT. DABS(FIFO) )  GO TO 359
        //               WVWMAX = DABS(FIFO)
        //               VMAX = VVAL
        //               RPMAX = RP
        //               RTMAX = RT
        //               ZEROSW = .FALSE.
        //  359     CONTINUE
        for (IV = 1; IV <= 5; IV++) {
            VVAL = VS[IV];
            for (I = 1; I <= 2; I++) {
                if (BSPROD(FIFO, 3, U + 0.5 * VVAL, U - 0.5 * VVAL, XS[I],
                       ISCTCR, 1, RP, RT) != 0) goto L359;
                TEMP = TEMP + std::fabs(FIFO);
                if (WVWMAX > std::fabs(FIFO)) goto L359;
                WVWMAX = std::fabs(FIFO);
                VMAX   = VVAL;
                RPMAX  = RP;
                RTMAX  = RT;
                ZEROSW = false;
L359:           ;
            }
        }
    }

    //         IF ( ZEROSW )  GO TO 375
    if (ZEROSW) goto L375;

    //         NINTRP = NINTRP + 10*(INADD+2)
    //         SUM0 = SUM0 + TEMP
    //         SUM1 = SUM1 + TEMP*U
    //         SUM2 = SUM2 + TEMP*U*U
    NINTRP = NINTRP + 10 * (INADD + 2);
    SUM0   = SUM0 + TEMP;
    SUM1   = SUM1 + TEMP * U;
    SUM2   = SUM2 + TEMP * U * U;

    //         IF ( PVPMAX .GT. TEMP )  GO TO 364
    //         PVPMAX = TEMP
    //         UMAX = U
    if (PVPMAX <= TEMP) {
        PVPMAX = TEMP;
        UMAX   = U;
    }

    // 364     IF (IPRNT .GE. 7)  write (6, 367) U, TEMP, WVWMAX, VMAX, RPMAX, RTMAX
    // 367     FORMAT ( 6G13.3 )
    if (IPRNT >= 7) {
        std::printf("%13.3G%13.3G%13.3G%13.3G%13.3G%13.3G\n",
                    U, TEMP, WVWMAX, VMAX, RPMAX, RTMAX);
    }

    //         IF ( WVWMAX .GE. RVRLIM  .OR.  U .LT. ULIM )  GO TO 365
    //         IF ( RP .GT. RLMAXS(1)  .AND.  RT .GT. RLMAXS(2) )  GO TO 380
    // 365     U = U + USTEP
    //         GO TO 350
    if (WVWMAX >= RVRLIM || U < ULIM) goto L365;
    if (RP > RLMAXS1 && RT > RLMAXS2) goto L380;
L365:
    U = U + USTEP;
    goto L350;

L375:
    //  375  write (6, 378)
    // 378  FORMAT ( '0**** WARNING:  SEARCH FOR SUMMAX WAS STOPPED',
    //     1   'BY THE BOUND STATE ASYMPTOPIA.' /
    //     2   '0**** IN FUTURE RUNS INCREASE ONE OR BOTH ASYMPTOPIA.' )
    std::printf("\n0**** WARNING:  SEARCH FOR SUMMAX WAS STOPPED"
                "BY THE BOUND STATE ASYMPTOPIA.\n"
                "0**** IN FUTURE RUNS INCREASE ONE OR BOTH ASYMPTOPIA.\n");
    U = U - USTEP;

L380:
    //  380  TINTRP = TINTRP + second() - TIME1
    //       IF ( SUMMAX .EQ. UNDEF )  SUMMAX = U
    //       IF ( PINFSW )  write (6, 381) U, SUMMAX
    // 381  FORMAT ( ' SUMMAX:  PROPOSED =', G14.5, 5X, 'USED =', G14.5 )
    TINTRP = TINTRP + second() - TIME1;
    if (SUMMAX == UNDEF) SUMMAX = U;
    if (PINFSW) {
        std::printf(" SUMMAX:  PROPOSED =%14.5G     USED =%14.5G\n", U, SUMMAX);
    }

    // =========================================================================
    //  COMPUTE SUMMID AS THE FIRST MOMENT OF U.
    // =========================================================================
    //      U = SUM1/SUM0
    //      SUM2 = SUM2/SUM0
    //      SUM2 = DSQRT(SUM2)
    //      IF ( SUMMID .EQ. UNDEF )  SUMMID = U*AMDMLT
    U    = SUM1 / SUM0;
    SUM2 = SUM2 / SUM0;
    SUM2 = std::sqrt(SUM2);
    if (SUMMID == UNDEF) SUMMID = U * AMDMLT;

    //      IF ( PINFSW )  write (6, 388) U, SUM2, UMAX, U, SUMMID
    // 388  FORMAT ( '0<SUM> =', G14.5, 5X, 'SQRT<SUM**2> =', G14.5,
    //     1  5X, 'SUM(MAX TERM) =', G14.5 /
    //     2  ' SUMMID:  PROPOSED =', G14.5, 5X, 'USED =', G14.5 )
    if (PINFSW) {
        std::printf("\n0<SUM> =%14.5G     SQRT<SUM**2> =%14.5G     SUM(MAX TERM) =%14.5G\n"
                    " SUMMID:  PROPOSED =%14.5G     USED =%14.5G\n",
                    U, SUM2, UMAX, U, SUMMID);
    }

    // =========================================================================
    //  ADJUST SUMMIN, SUMMID TO REQUIREMENTS OF CUBIC MAP.
    //
    //      SUMMIN = DMIN1( SUMMIN, 7.D0*(SUMMID-SUMMAX/7.D0)/6.D0 )
    //      SUMMID = DMIN1( SUMMID, .5*(SUMMIN+SUMMAX) )
    //      IF (PINFSW)  write (6,  393) SUMMIN, SUMMID, SUMMAX
    //  393 FORMAT ( ' FINAL CHOICE OF SUMMIN, SUMMID, SUMMAX =', 5G15.5 )
    // =========================================================================
    SUMMIN = std::fmin(SUMMIN, 7.0 * (SUMMID - SUMMAX / 7.0) / 6.0);
    SUMMID = std::fmin(SUMMID, 0.5 * (SUMMIN + SUMMAX));
    if (PINFSW) {
        std::printf(" FINAL CHOICE OF SUMMIN, SUMMID, SUMMAX =%15.5G%15.5G%15.5G\n",
                    SUMMIN, SUMMID, SUMMAX);
    }

    // =========================================================================
    //  SUMMIN AND SUMMAX NOW KNOWN; PERFORM SINH MAP OF @SUMMIN,SUMMAX@
    //  TO @-1,1@ ; COMPUTE AND STORE GAUSS-LEGENDRE PTS
    //     FOR BOTH THE H'S AND THE INTERPOLATED INTEGRAL GRID
    // =========================================================================

L400:
    //  400  ISMHPT=NALLOC( NPSUM, 'SUMHPTS ' )
    //       ISMHWK=NALLOC( 3*NPSUM, 'SUMHWORK' )
    ISMHPT    = NALLOC("SUMHPTS ", NPSUM);
    ISMHWK    = NALLOC("SUMHWORK", 3 * NPSUM);
    ISMHPT_g  = ISMHPT;
    ISMHWK_g  = ISMHWK;

    //  THE NUMBER OF POINTS FOR THE INTERPOLATED INTEGRAL GRID
    //  IS BASED ON THE WAVE-LENGTH
    //
    //      NPSUMI = (SUMMAX-SUMMIN)*SUMPTS*(AKI+AKO)/(4*PI)
    //      NPSUMI = MAX0( NPSUMI, NPSUM )
    {
        double& SUMPTS = FLOAT_common.SUMPTS;
        NPSUMI = (int)((SUMMAX - SUMMIN) * SUMPTS * (KANDM.AKI + KANDM.AKO) / (4.0 * PI));
        NPSUMI = std::max(NPSUMI, NPSUM);
    }

    //      IF ( IHSAVE .NE. 2 )  NPLYSW = NPSUMI .EQ. NPSUM
    if (IHSAVE != 2) NPLYSW = (NPSUMI == NPSUM);

    //      IF ( .NOT. NPLYSW  .OR.  NPSUM .EQ. NPSUMI )  GO TO 404
    //      write (6, 403) NPSUMI, NPSUM
    // 403  FORMAT ( '0**** WARNING:  NPSUMI WAS COMPUTED TO BE', I5 /
    //     1  6X, 'BUT DUE TO LIMITATIONS OF THE SAVEHS RUN',
    //     2    ' IT WILL BE LIMITED TO', I5 )
    //      NPSUMI = NPSUM
    if (!NPLYSW || NPSUM == NPSUMI) goto L404;
    std::printf("\n0**** WARNING:  NPSUMI WAS COMPUTED TO BE%5d\n"
                "      BUT DUE TO LIMITATIONS OF THE SAVEHS RUN IT WILL BE LIMITED TO%5d\n",
                NPSUMI, NPSUM);
    NPSUMI = NPSUM;

L404:
    //  404  ISMIPT = NALLOC( NPSUMI, 'SUMIPTS ' )
    //       ISMIVL = NALLOC( NPSUMI, 'SUMIVALS' )
    ISMIPT   = NALLOC("SUMIPTS ", NPSUMI);
    ISMIVL   = NALLOC("SUMIVALS", NPSUMI);
    ISMIPT_g = ISMIPT;
    ISMIVL_g = ISMIVL;

    //  GENERATE GAUSS POINTS FOR U
    //      CALL CUBMAP ( MAPSUM,  SUMMIN, SUMMID, SUMMAX, GAMSUM,
    //     1  ALLOC(Z(ISMIPT)), ALLOC(Z(ISMIVL)), NPSUMI )
    CUBMAP(MAPSUM, SUMMIN, SUMMID, SUMMAX, GAMSUM,
           ALLOC_base(LOCPTRS.Z[ISMIPT]), ALLOC_base(LOCPTRS.Z[ISMIVL]), NPSUMI);

    //      NRIROI = NPSUMI*NPDIF
    //      N4RIO = (NRIROI+FACFR4-1)/FACFR4
    //      IWIO = NALLOC( N4RIO, 'RIROWTS ' )
    //      IRI = NALLOC( N4RIO, 'RIPTS   ' )
    //      IRO = NALLOC( N4RIO, 'ROPTS   ' )
    //      IDIFPT = NALLOC( NPDIF, 'DIFPTS  ' )
    //      IDIFWT = NALLOC( NPDIF, 'DIFWTS  ' )
    //      NPOLSV = 3*NVPOLY+3
    //      IVPOLS = NALLOC( NPOLSV, 'VPOLYS  ' )
    NRIROI  = NPSUMI * NPDIF;
    NRIROI_g = NRIROI;
    N4RIO   = (NRIROI + FACFR4 - 1) / FACFR4;
    IWIO_loc = NALLOC("RIROWTS ", N4RIO);
    IRI_loc  = NALLOC("RIPTS   ", N4RIO);
    IRO_loc  = NALLOC("ROPTS   ", N4RIO);
    IDIFPT   = NALLOC("DIFPTS  ", NPDIF);
    IDIFWT   = NALLOC("DIFWTS  ", NPDIF);
    NPOLSV   = 3 * NVPOLY + 3;
    IVPOLS   = NALLOC("VPOLYS  ", NPOLSV);
    // Store into GRIDCM
    IWIO_g  = IWIO_loc;
    IRI_g   = IRI_loc;
    IRO_g   = IRO_loc;

    // =========================================================================
    //  IF USING H'S, RESTORE GENERATED SUM POINTS AND POLYS
    //
    //      IF ( IHSAVE .NE. 2 )  GO TO 405
    //      READ ( 1, END=9950, ERR=9950 )
    //     1  ( ALLOC(Z(IVPOLS)-1+I), I = 1, NPOLSV ),
    //     2  ( ALLOC(Z(ISMHPT)-1+I), I = 1, NPSUM )
    //      II = 3*NPSUM
    //      IF ( NPLYSW )  READ ( 1, END=9950, ERR=9950 )
    //     1  ( ALLOC(Z(ISMHWK)-1+I), I = 1, II )
    //      GO TO 700
    // =========================================================================
    if (IHSAVE != 2) goto L405;
    // Read IVPOLS data from H's file
    {
        // Binary read of NPOLSV doubles starting at ALLOC(Z(IVPOLS)-1+1)
        int base_vpols = LOCPTRS.Z[IVPOLS] - 1;
        for (I = 1; I <= NPOLSV; I++) {
            // placeholder: actual file read would go here
            (void)ALLOC(base_vpols + I);
        }
        int base_smhpt = LOCPTRS.Z[ISMHPT] - 1;
        for (I = 1; I <= NPSUM; I++) {
            (void)ALLOC(base_smhpt + I);
        }
        II = 3 * NPSUM;
        if (NPLYSW) {
            int base_smhwk = LOCPTRS.Z[ISMHWK] - 1;
            for (I = 1; I <= II; I++) {
                (void)ALLOC(base_smhwk + I);
            }
        }
    }
    goto L700;

L405:
    //  405  CALL CUBMAP ( MAPSUM,  SUMMIN, SUMMID, SUMMAX, GAMSUM,
    //     1  ALLOC(Z(ISMHPT)), ALLOC(Z(ISMHWK)), NPSUM )
    CUBMAP(MAPSUM, SUMMIN, SUMMID, SUMMAX, GAMSUM,
           ALLOC_base(LOCPTRS.Z[ISMHPT]), ALLOC_base(LOCPTRS.Z[ISMHWK]), NPSUM);

    // =========================================================================
    //  SETUP THE DIFFERENCE GRID.
    //
    //  THE STAGES ARE
    //  1)  DETERMINE THE RANGE OF V FOR EACH VALUE OF U.
    //  2)  INSIDE EACH V-RANGE, FIND THE MOMENTS OF V AND THUS
    //      THE "MID-POINT".
    //  3)  MAKE A LOW-DEGREE POLYNOMIAL FIT TO THE MID- AND END-
    //      POINTS.
    //  4)  GENERATE ALL THE RI, RO PAIRS FOR FINDING H'S
    //  5)  GENERATE ALL RI, RO PAIRS FOR THE INTERPOLATED INTEGRATION
    //      GRID.
    //
    //  IN ALL THESE OPERATIONS WE DO NOT INCLUDE THE SCATTERING
    //  WAVEFUNCTIONS SINCE  V << U  SO THAT RI, RO WILL VARY ONLY
    //  A LITTLE FOR FIXED U.
    //
    //  IN FINDING THE EXTREMA OF THE V'S, WE USE
    //    PHI' V PHI'
    // =========================================================================

    //      NRIROH = NPDIF*NPSUM
    //      IRIOEX = NALLOC( NRIROH, 'RIROEXPS' )
    //      IRIH = NALLOC( NRIROH, 'RIH     ' )
    //      IROH = NALLOC( NRIROH, 'ROH     ' )
    NRIROH_h = NPDIF * NPSUM;
    IRIOEX   = NALLOC("RIROEXPS", NRIROH_h);
    IRIH     = NALLOC("RIH     ", NRIROH_h);
    IROH     = NALLOC("ROH     ", NRIROH_h);
    IRIOEX_g = IRIOEX;
    // Store NRIROH into GRIDCM.NRIROI (reused for H-grid count separately)
    // Note: Fortran uses NRIROH local for both H-grid and NRIROH_h.
    // We track NRIROH_h as the H-search grid size.
    // The Fortran NRIROI = NRIROH is assigned after the 50x H-grid is found.

    //  WORK AREAS FOR STAGE 3
    //      IVWORK = NALLOC( 4*NPSUM+3+(NVPOLY+1)**2, 'VPOLYWRK' )
    IVWORK_loc = NALLOC("VPOLYWRK", 4 * NPSUM + 3 + (NVPOLY + 1) * (NVPOLY + 1));
    IVWORK_g   = IVWORK_loc;

    //      LVPOLS = Z(IVPOLS)
    //      LVWTS = Z(IVWORK)
    //      LVRESD = LVWTS+NPSUM
    //      LVSMSQ = LVRESD + 3*NPSUM
    //      LVAMAT = LVSMSQ + 3
    //      LVMIN=Z(ISMHWK)-1
    //      LVMID = LVMIN + NPSUM
    //      LVMAX = LVMID + NPSUM
    //      LRIH = Z(IRIH) - 1
    //      LROH = Z(IROH) - 1
    //      LRIOEX = Z(IRIOEX) - 1
    //      LDIFPT=Z(IDIFPT)-1
    //      LDIFWT=Z(IDIFWT)-1
    //      LSUMPT = Z(ISMHPT) - 1
    LVPOLS = LOCPTRS.Z[IVPOLS];
    LVWTS  = LOCPTRS.Z[IVWORK_loc];
    LVRESD = LVWTS  + NPSUM;
    LVSMSQ = LVRESD + 3 * NPSUM;
    LVAMAT = LVSMSQ + 3;
    LVMIN  = LOCPTRS.Z[ISMHWK] - 1;
    LVMID  = LVMIN  + NPSUM;
    LVMAX  = LVMID  + NPSUM;
    LRIH   = LOCPTRS.Z[IRIH]   - 1;
    LROH   = LOCPTRS.Z[IROH]   - 1;
    LRIOEX = LOCPTRS.Z[IRIOEX] - 1;
    LDIFPT = LOCPTRS.Z[IDIFPT] - 1;
    LDIFWT = LOCPTRS.Z[IDIFWT] - 1;
    LSUMPT = LOCPTRS.Z[ISMHPT] - 1;

    // =========================================================================
    //  STAGE 1)  END POINTS OF EACH V-RANGE.
    // =========================================================================
    //      VMAX = 1
    //      VMIN = 1
    //      NBUMPS = 0
    //      RPMAX = 0
    //      RTMAX = 0
    //      DV = 1.D0/LOOKST
    VMAX   = 1.0;
    VMIN   = 1.0;
    NBUMPS = 0;
    RPMAX  = 0.0;
    RTMAX  = 0.0;
    DV     = 1.0 / (double)LOOKST;


    //      DO 489  IU=1,NPSUM
    for (IU = 1; IU <= NPSUM; IU++) {
        //         U=ALLOC(LSUMPT+IU)
        //         VLEN=2.D0*U
        //         IF ( U .LT. 1 )  GO TO 465
        double VLEN;
        U    = ALLOC(LSUMPT + IU);
        VLEN = 2.0 * U;
        if (U < 1.0) goto L465;

        //         VVAL=VMAX
        //         SYNE= +.5*VLEN
        VVAL = VMAX;
        SYNE = +0.5 * VLEN;

        TIME1 = second();

L420:
        //  420     VVAL = DMIN1( 1.D0, VVAL+3*DV )
        VVAL = std::fmin(1.0, VVAL + 3.0 * DV);

L422:
        //  422        IF ( VVAL .LE. .5*DV )  GO TO 450
        if (VVAL <= 0.5 * DV) goto L450;
        //            RI=U+VVAL*SYNE
        //            RO=U-VVAL*SYNE
        //            ULIM = RVRLIM/DMAX1( 1.D-2, RI*RO )
        RI   = U + VVAL * SYNE;
        RO   = U - VVAL * SYNE;
        ULIM = RVRLIM / std::fmax(1.0e-2, RI * RO);

        //            DO 429  I = 1, 2
        //               CALL BSPROD ( FIFO, 2, RI, RO, XS(I), ISCTMN, 1, RP, RT,
        //     1           *435 )
        //               NINTRP=NINTRP+INADD
        //               IF ( DABS(FIFO) .GT. ULIM )  GO TO 450
        //  429        CONTINUE
        for (I = 1; I <= 2; I++) {
            int bsp_ret = BSPROD(FIFO, 2, RI, RO, XS[I], ISCTMN, 1, RP, RT);
            if (bsp_ret != 0) goto L435;
            NINTRP = NINTRP + INADD;
            if (std::fabs(FIFO) > ULIM) goto L450;
        }

L435:
        //  435        VVAL=VVAL-DV
        //             GO TO 422
        VVAL = VVAL - DV;
        goto L422;

L450:
        //  450     VVAL = DMIN1( 1.D0, VVAL+DV )
        VVAL = std::fmin(1.0, VVAL + DV);

        //  MAKE SURE THIS DESIRED VALUE OF V DOES NOT VIOLATE ASYMPTOPIA
        //  BECAUSE  S1*T1 < 0  AND  S2*T2 < 0,  RP2 > RP1 AND
        //  RT2 > RT1 ALWAYS.
        //
        //         RI=U+VVAL*SYNE
        //         RO=U-VVAL*SYNE
        //         RT=DSQRT((S1*RI)**2+(T1*RO)**2+ST1*RI*RO*(1.D0-DXV))
        //         RP=DSQRT((S2*RI)**2+(T2*RO)**2+ST2*RI*RO*(1.D0-DXV))
        //         RTMAX = DMAX1( RTMAX, RT )
        //         RPMAX = DMAX1( RPMAX, RP )
        //         IF ( RT .LE. BNDMXT  .AND. RP .LE. BNDMXP )  GO TO 455
        RI    = U + VVAL * SYNE;
        RO    = U - VVAL * SYNE;
        RT    = std::sqrt((S1 * RI) * (S1 * RI) + (T1 * RO) * (T1 * RO) + ST1 * RI * RO * (1.0 - DXV));
        RP    = std::sqrt((S2 * RI) * (S2 * RI) + (T2 * RO) * (T2 * RO) + ST2 * RI * RO * (1.0 - DXV));
        RTMAX = std::fmax(RTMAX, RT);
        RPMAX = std::fmax(RPMAX, RP);
        if (RT <= BNDMXT && RP <= BNDMXP) goto L455;

        //  CANNOT USE DESIRED VALUE
        //
        //         NBUMPS = NBUMPS + 1
        //         VVAL = VVAL - DV
        NBUMPS = NBUMPS + 1;
        VVAL   = VVAL - DV;

L455:
        //  455     IF(SYNE .GT. 0) VMAX=VVAL
        //          IF(SYNE .LT. 0 ) VMIN=VVAL
        //          IF(SYNE .LT. 0 ) GO TO 460
        //          VVAL=VMIN
        //          SYNE=-SYNE
        //          GO TO 420
        if (SYNE > 0.0) VMAX = VVAL;
        if (SYNE < 0.0) VMIN = VVAL;
        if (SYNE < 0.0) goto L460;
        VVAL = VMIN;
        SYNE = -SYNE;
        goto L420;

L460:
        TINTRP = TINTRP + second() - TIME1;
        // fall through to L465

L465:
        //  465     IF ( PGRDSW )  write (6, 468 ) U, VMIN, VMAX,
        //     1     RP, RT, FIFO
        // 468     FORMAT ( ' U VMIN VMAX =', 8G13.5 )
        if (PGRDSW) {
            std::printf(" U VMIN VMAX =%13.5G%13.5G%13.5G%13.5G%13.5G%13.5G\n",
                        U, VMIN, VMAX, RP, RT, FIFO);
        }

        //  THESE ARE EXTREMA OF V=DIF FOR EACH POINT IN THE U=SUM GRID.
        //         ALLOC(LVMIN+IU) = -VMIN*VLEN
        //         ALLOC(LVMAX+IU) = VMAX*VLEN
        ALLOC(LVMIN + IU) = -VMIN * VLEN;
        ALLOC(LVMAX + IU) =  VMAX * VLEN;
        //  CHECK FOR RANGES CONTRACTED TO ZERO BY ASYMPTOPIA LIMITATIONS
        //
        //         VVAL = 2*U*(VMIN+VMAX)
        //         IF ( VVAL .LT. .1*DV*U )  write (6, 473) U, VVAL
        // 473     FORMAT ( '0**** ERROR:  THE DIF GRID RANGE FOR SUM =', F8.2, ... )
        {
            double VVAL_check = 2.0 * U * (VMIN + VMAX);
            if (VVAL_check < 0.1 * DV * U) {
                std::printf("\n0**** ERROR:  THE DIF GRID RANGE FOR SUM =%8.2f"
                            " IS ESSENTIALLY ZERO DUE TO BOUND STATE ASYMPTOPIA "
                            "LIMITATIONS; RANGE =%13.3G\n"
                            " **** IN FUTURE RUNS ASYMPTOPIA MUST BE INCREASED.\n",
                            U, VVAL_check);
            }
        }

    }   // end DO 489  IU=1,NPSUM

    //  END OF IU (SUM) LOOP
    //
    //      IF ( NBUMPS .NE. 0 )  write (6, 493) NBUMPS
    // 493  FORMAT ('0**** WARNING:  SEARCH FOR THE DIF GRIDS WAS STOPPED',
    //     1  I4, ' TIMES BY BOUND STATE ASYMPTOPIA.' )
    //      IF ( NBUMPS .NE. 0 )  write (6, 897) RPMAX, RTMAX
    if (NBUMPS != 0) {
        std::printf("\n0**** WARNING:  SEARCH FOR THE DIF GRIDS WAS STOPPED%4d"
                    " TIMES BY BOUND STATE ASYMPTOPIA.\n",
                    NBUMPS);
        std::printf(" **** IN FUTURE RUNS INCREASE ASYMPTOPIA TO MORE THAN:%15.2f"
                    " (PROJECTILE)%15.2f (TARGET)\n",
                    RPMAX, RTMAX);
    }

    // =========================================================================
    // Stage 2) Find V-moments of the weight function for each U
    // =========================================================================
    if (PGRDSW) std::printf("\n IU        U           VMIN        VMID        VMAX"
        "      V(MAX)      PVPMAX         <V>      <V**2>     <V VAR>\n");

    for (IU = 1; IU <= NPSUM; IU++) {
        U = ALLOC(LSUMPT + IU);
        VMIN = ALLOC(LVMIN + IU);
        VMAX = ALLOC(LVMAX + IU);
        ALLOC(LVWTS - 1 + IU) = 1.0 / ((VMAX - VMIN) * (VMAX - VMIN));
        IMAX = 2 * NPDIF;
        DV = (VMAX - VMIN) / (IMAX + 1);
        VVAL = VMIN;
        PVPMAX = 0; SUM0 = 0; SUM1 = 0; SUM2 = 0;
        for (IV = 1; IV <= IMAX; IV++) {
            VVAL = VVAL + DV;
            RI = U + 0.5 * VVAL;
            RO = U - 0.5 * VVAL;
            TEMP = 0;
            bool skip539 = false;
            for (I = 1; I <= 2; I++) {
                if (BSPROD(FIFO, 1, RI, RO, XS[I], ISCTCR, 1, RP, RT) != 0)
                    { skip539 = true; break; }
                NINTRP = NINTRP + INADD;
                TEMP = std::fabs(FIFO) + TEMP;
            }
            if (skip539) continue;
            SUM0 = SUM0 + TEMP;
            SUM1 = SUM1 + TEMP * VVAL;
            SUM2 = SUM2 + TEMP * VVAL * VVAL;
            if (TEMP <= PVPMAX) continue;
            PVPMAX = TEMP;
            VOFMAX = VVAL;
        }
        SUM1 = SUM1 / SUM0;
        SUM2 = SUM2 / SUM0;
        SUM0 = std::sqrt(SUM2 - SUM1 * SUM1);
        SUM2 = std::sqrt(SUM2);
        VMID = VOFMAX;
        TEMP = 0.3 * (VMAX - VMIN);
        VMID = std::min(std::max(VMID, VMIN + TEMP), VMAX - TEMP);
        if (PGRDSW)
            std::printf("%4d%14.4G%14.4G%14.4G%14.4G%14.4G%14.4G%14.4G%14.4G%14.4G\n",
                        IU, U, VMIN, VMID, VMAX, VOFMAX, PVPMAX, SUM1, SUM2, SUM0);
        VMID = (VMID - VMIN) / (VMAX - VMIN);
        ALLOC(LVMID + IU) = VMID;
    }

    // =========================================================================
    // Stage 3) Polynomial fits to V-ranges
    // =========================================================================
    NVTERM = NVPOLY + 1;

    LSQPOL(ALLOC_base(LSUMPT + 1), ALLOC_base(LVMIN + 1), ALLOC_base(LVWTS),
           ALLOC_base(LVRESD), NPSUM, ALLOC_base(LVSMSQ), 3, ALLOC_base(LVAMAT),
           ALLOC_base(LVPOLS), NVTERM, NPSUM, NVTERM);

    if (PINFSW) {
        std::printf("\n\n%2d-DEGREE POLYNOMIAL FITS TO DIF RANGES:\n"
                    "\n       V LO           V MID          V HI\n", NVPOLY);
        for (N = 1; N <= NVTERM; N++)
            std::printf("%15.5G%15.5G%15.5G\n",
                        ALLOC(LVPOLS + N - 1), ALLOC(LVPOLS + N - 1 + NVTERM),
                        ALLOC(LVPOLS + N - 1 + 2 * NVTERM));
    }
    if (PGRDSW) std::printf("\n       U              V LO     RESID"
        "               V MID     RESID               V HI     RESID\n\n");

    VRNGMX = 0; VRNGMN = CNSTNT.BIGNUM; VMAXMX = 0;
    for (IU = 1; IU <= NPSUM; IU++) {
        if (!NPLYSW) {
            for (I = 1; I <= 3; I++)
                ALLOC(LVMIN + (I-1)*NPSUM + IU) += ALLOC(LVRESD - 1 + (I-1)*NPSUM + IU);
        }
        VMIN = ALLOC(LVMIN + IU);
        VMAX = ALLOC(LVMAX + IU);
        VMID = (VMAX - VMIN) * ALLOC(LVMID + IU) + VMIN;
        ALLOC(LVMID + IU) = VMID;
        if (PGRDSW) {
            std::printf("%14.4G", ALLOC(LSUMPT + IU));
            for (I = 1; I <= 3; I++)
                std::printf("   %14.4G%12.2G",
                    ALLOC(LVMIN + (I-1)*NPSUM + IU), ALLOC(LVRESD - 1 + (I-1)*NPSUM + IU));
            std::printf("\n");
        }
        VRNGMX = std::max(VRNGMX, VMAX - VMIN);
        VRNGMN = std::min(VRNGMN, VMAX - VMIN);
        VMAXMX = std::max({VMAXMX, VMAX, -VMIN});
        TEMP = 0.2 * (VMAX - VMIN);
        if (VMIN > VMID - TEMP || VMAX < VMID + TEMP)
            std::printf(" **** INVALID VMIN, VMID, VMAX:%15.5G%15.5G%15.5G\n", VMIN, VMID, VMAX);
    }

    // =========================================================================
    // Stage 4) Generate the RI, RO points
    // =========================================================================
    if (IRECT == 1) { VRNGMX = 2.0 * VMAXMX; VRNGMN = VRNGMX; }

    for (IU = 1; IU <= NPSUM; IU++) {
        U = ALLOC(LSUMPT + IU);
        VMIN = ALLOC(LVMIN + IU);
        VMAX = ALLOC(LVMAX + IU);
        VMID = ALLOC(LVMID + IU);
        if (IRECT == 1) { VMAX = 0.5*VRNGMX; VMIN = -VMAX; VMID = 0; }

        if (VMIN < -2*U || VMAX > 2*U) {
            std::printf("0*** ILLEGAL VMIN OR VMAX:%3d", IU);
            print_G(15, 5, U); print_G(15, 5, VMIN); print_G(15, 5, VMAX);
            std::printf("\n *** VALUES SET TO MAXIMUM ALLOWED.\n");
        }
        VMIN = std::max(VMIN, -2.0*U);
        VMAX = std::min(VMAX, 2.0*U);
        TEMP = 0.3*(VMAX - VMIN);
        VMID = std::min(std::max(VMID, VMIN + TEMP), VMAX - TEMP);

        SYNE = 1;
        if (VMID > 0.5*(VMAX + VMIN)) {
            VMID = -VMID; TEMP = VMAX; VMAX = -VMIN; VMIN = -TEMP; SYNE = -1;
        }

        CUBMAP(MAPDIF, VMIN, VMID, VMAX, GAMDIF, ALLOC_base(LDIFPT + 1),
               ALLOC_base(LDIFWT + 1), NPDIF);

        for (IV = 1; IV <= NPDIF; IV++) {
            IPLUNK = IV - 1;
            if (SYNE < 0) IPLUNK = NPDIF - IV;
            IPLUNK = NPSUM * IPLUNK + IU;
            VVAL = ALLOC(LDIFPT + IV);
            RI = U + 0.5*VVAL*SYNE;
            RO = U - 0.5*VVAL*SYNE;
            ALLOC(LRIH + IPLUNK) = RI;
            ALLOC(LROH + IPLUNK) = RO;
            RP = std::sqrt(1 + (S1*RI + T1*RO)*(S1*RI + T1*RO));
            RT = std::sqrt(1 + (S2*RI + T2*RO)*(S2*RI + T2*RO));
            ALLOC(LRIOEX + IPLUNK) = std::exp(ALPHAP*RP + ALPHAT*RT);
        }
    }

    // Check ISMIVL before freeing work arrays
    // Free work arrays (use IVWORK_loc, not uninitialized local IVWORK)
    LOCPTRS.Z[IVWORK_loc] = -LOCPTRS.Z[IVWORK_loc];
    LOCPTRS.Z[ISCTMN] = -LOCPTRS.Z[ISCTMN];
    LOCPTRS.Z[ISCTCR] = -LOCPTRS.Z[ISCTCR];

    // =========================================================================
    // 5) Generate points for the interpolated integral grid (label 700)
    // =========================================================================
L700:
    NRIROH = NRIROH_h;  // H-grid size used for phi-grid passes
    AREA = 0;
    LVPOLS = LOCPTRS.Z[IVPOLS] - 1;
    LSMIVL = LOCPTRS.Z[ISMIVL_g] - 1;
    LSUMPT = LOCPTRS.Z[ISMIPT_g] - 1;
    LWIO = ALLOCS.FACFR4 * LOCPTRS.Z[IWIO_g] - ALLOCS.FACFR4;
    LRI  = ALLOCS.FACFR4 * LOCPTRS.Z[IRI_g]  - ALLOCS.FACFR4;
    LRO  = ALLOCS.FACFR4 * LOCPTRS.Z[IRO_g]  - ALLOCS.FACFR4;
    LVMIN = LOCPTRS.Z[ISMHWK_g] - 1;

    RIMAX_k = 0; ROMAX_k = 0;
    LDIFPT = LOCPTRS.Z[IDIFPT] - 1;
    LDIFWT = LOCPTRS.Z[IDIFWT] - 1;

    for (IU = 1; IU <= NPSUMI; IU++) {
        WOW = ALLOC(LSMIVL + IU);
        U = ALLOC(LSUMPT + IU);

        // Evaluate polynomial-determined ranges
        ITEMP = LVPOLS;
        for (I = 1; I <= 3; I++) {
            VRANGE[I] = ALLOC(LVMIN + (I-1)*NPSUM + IU);
            if (NPLYSW) continue;
            ITEMP = ITEMP + NVPOLY + 1;
            VRANGE[I] = ALLOC(ITEMP);
            for (N = 1; N <= NVPOLY; N++)
                VRANGE[I] = ALLOC(ITEMP - N) + U * VRANGE[I];
        }
        VMIN = VRANGE[1]; VMID = VRANGE[2]; VMAX = VRANGE[3];
        if (!NPLYSW) VMID = (VMAX - VMIN) * VMID + VMIN;
        if (IRECT == 1) { VMAX = 0.5*VRNGMX; VMIN = -VMAX; VMID = 0; }

        if (VMIN < -2*U || VMAX > 2*U) {
            std::printf("0*** ILLEGAL INTERPOLATED VMIN OR VMAX:%4d", IU);
            print_G(15, 5, U); print_G(15, 5, VMIN); print_G(15, 5, VMAX);
            std::printf("\n *** VALUES SET TO MAXIMUM ALLOWED.\n");
        }

        AREA = AREA + WOW * (VMAX - VMIN);
        VMIN = std::max(VMIN, -2.0*U);
        VMAX = std::min(VMAX, 2.0*U);
        TEMP = 0.3*(VMAX - VMIN);
        VMID = std::min(std::max(VMID, VMIN + TEMP), VMAX - TEMP);
        SYNE = 1;
        if (VMID > 0.5*(VMAX + VMIN)) {
            VMID = -VMID; TEMP = VMAX; VMAX = -VMIN; VMIN = -TEMP; SYNE = -1;
        }

        CUBMAP(MAPDIF, VMIN, VMID, VMAX, GAMDIF, ALLOC_base(LDIFPT + 1),
               ALLOC_base(LDIFWT + 1), NPDIF);

        for (IV = 1; IV <= NPDIF; IV++) {
            IPLUNK = IV - 1;
            if (SYNE < 0) IPLUNK = NPDIF - IV;
            IPLUNK = NPSUMI * IPLUNK + IU;
            VVAL = ALLOC(LDIFPT + IV);
            RI = U + 0.5*VVAL*SYNE;
            RO = U - 0.5*VVAL*SYNE;
            RIMAX_k = std::max(RIMAX_k, RI);
            ROMAX_k = std::max(ROMAX_k, RO);
            if (IPLUNK > NRIROI_g) {
                std::printf(" *** IPLUNK OVERFLOW: %d > NRIROI=%d IU=%d IV=%d NPSUMI=%d NPDIF=%d\n",
                    IPLUNK, NRIROI_g, IU, IV, NPSUMI, NPDIF);
                break;
            }
            ALLOC4(LRI + IPLUNK) = (float)RI;
            ALLOC4(LRO + IPLUNK) = (float)RO;

            RP = std::sqrt(1 + (S1*RI + T1*RO)*(S1*RI + T1*RO));
            RT = std::sqrt(1 + (S2*RI + T2*RO)*(S2*RI + T2*RO));
            TEMP = std::exp(-(ALPHAP*RP + ALPHAT*RT));
            ALLOC4(LWIO + IPLUNK) = (float)(JACOB * RI * RO * WOW * ALLOC(LDIFWT + IV) * TEMP);
        }
    }

    LOCPTRS.Z[IDIFPT] = -LOCPTRS.Z[IDIFPT];
    LOCPTRS.Z[IDIFWT] = -LOCPTRS.Z[IDIFWT];
    LOCPTRS.Z[ISMIVL_g] = -LOCPTRS.Z[ISMIVL_g];

    if (PINFSW)
        std::printf(" MAXIMUM RI AND RO:%15.5G%15.5G\n", RIMAX_k, ROMAX_k);

    // Expand wavefunction work area if needed
    NSTP2S[1] = (int)(RIMAX_k / RSTEPS[1] + 3.5);
    NSTP2S[2] = (int)(ROMAX_k / RSTEPS[2] + 3.5);
    I = std::max(NSTP2S[1], NSTP2S[2]) + 6;
    if (I > NWAVEF) {
        NWAVEF = I;
        if (TCSWS[1] || TCSWS[2]) I = 4 * I;
        IWAVR = NALLOC("WAVER   ", I);
        IWAVI = NALLOC("WAVEI   ", I);
    }

    // =========================================================================
    // L800: Phi grid setup (two passes)
    // =========================================================================
L800:
    if (IHSAVE == 2) goto L900;

    INEED = NRIROH / ALLOCS.FACFR4 + 1;
    LOGIC_g = NALLOC("LOGIC   ", INEED);

    IXTOPZ = LOOKST + 1;
    DXV = 2.0 / ((double)LOOKST * (double)LOOKST);
    NNAIT = 2;
    IBSTYP = 2;
    PAS1SW = true;
    X0MIN = 1; X0AV = 0;

    IPHIPT = NALLOC("PHIPTS  ", NPPHI);
    IPHIWT = NALLOC("PHIWTS  ", NPPHI);
    PHIMID = std::max(PHIMID, 0.10);
    PHIMID = std::min(PHIMID, 0.90);
    GAMPHI = std::max(GAMPHI, 1.0e-6);
    CUBMAP(MAPPHI, 0.0, PHIMID, 1.0, GAMPHI,
           ALLOC_base(LOCPTRS.Z[IPHIPT]),
           ALLOC_base(LOCPTRS.Z[IPHIWT]), NPPHI);

    // Come back here to start pass 2
L810:
    JCNT = 0; RTMAX = 0; RPMAX = 0; NBUMPS = 0;
    LLOGIC = ALLOCS.FACFR4 * LOCPTRS.Z[LOGIC_g] - ALLOCS.FACFR4;
    LPHIPT = LOCPTRS.Z[IPHIPT] - 1;
    LPHIWT = LOCPTRS.Z[IPHIWT] - 1;
    LRIH = LOCPTRS.Z[IRIH] - 1;
    LROH = LOCPTRS.Z[IROH] - 1;

    for (IPLUNK = 1; IPLUNK <= NRIROH; IPLUNK++) {
        RI = ALLOC(LRIH + IPLUNK);
        RO = ALLOC(LROH + IPLUNK);
        ULIM = RVRLIM / (std::max(1.0, RI) * std::max(1.0, RO));
        if (PAS1SW) goto L820;

        // Pass 2: set up end point
        IEND = ILLOC(LLOGIC + IPLUNK);
        IEND = std::max(IEND, 2);
        X0 = 1.0 - DXV * (IEND-1) * (IEND-1);
        PHI0 = std::acos(X0);

L820:
        if (PAS1SW) TIME1 = second();

        WOW = 0;
        for (II = 1; II <= IXTOPZ; II++) {
            X = 1.0 - DXV * (II-1) * (II-1);
            if (!PAS1SW) {
                PHI = PHI0 * ALLOC(LPHIPT + II);
                DPHI = PHI0 * ALLOC(LPHIWT + II) * std::sin(PHI);
                X = std::cos(PHI);
            }

            // Calculate bound-state product
            NBUMPS++;
            {
                int bret = BSPROD(FIFO, IBSTYP, RI, RO, X, ISCTMN, NNAIT, RP, RT);
                NBUMPS--;
                if (bret == 0) NINTRP += INADD;
            }
            RTMAX = std::max(RTMAX, RT);
            RPMAX = std::max(RPMAX, RP);
            if (!PAS1SW) {
                // Pass 2: store phi angles and weighted BS product
                JCNT++;
                TEMP = (T1*RO + S1*RI*X) / RT;
                ALLOC4(LPHIT + JCNT) = (float)(PHISGN * std::acos(TEMP));
                TEMP = (T2*RO + S2*RI*X) / RP;
                ALLOC4(LPHIP + JCNT) = (float)(PHISGN * std::acos(TEMP));
                ALLOC4(LPHI + JCNT) = (float)PHI;
                ALLOC4(LTRP + JCNT) = (float)(DPHI * FIFO);
                goto L859;
            }

            // Pass 1: finding X0
            FIFO = std::fabs(FIFO);
            if (II != 1) {
                if (FIFO < ULIM && WOW < ULIM) goto L865;
            }
            WOW = FIFO;
L859:       ;
        } // end II loop

        if (PAS1SW) goto L864;

        // Pass 2: dump if needed
        if (IPRNT >= 9) {
            N = JCNT - NPPHI + 1;
            std::printf(" %7.2f%7.2f%7.4f", RI, RO, PHI0);
            for (I = N; I <= JCNT; I++) std::printf("%9.1E", (double)ALLOC4(LTRP + I));
            std::printf("\n");
        }
        if (ISCRFL != 0) JCNT = 0;  // scratch file write stub
        goto L869;

L864:
        II = LOOKST + 2;
L865:
        TINTRP += second() - TIME1;
        IEND = std::min(II - 1 + NPHIAD, IXTOPZ);
        ILLOC(LLOGIC + IPLUNK) = IEND;
        if (IEND != 1) JCNT += NPPHI;
        X0 = 1.0 - DXV * (IEND-1) * (IEND-1);
        X0MIN = std::min(X0MIN, X0);
        X0AV = X0AV + X0;

L869:   ;
    } // end IPLUNK loop

    if (!PAS1SW) goto L890;

    // End of pass 1 — allocate and setup for pass 2
    if (NBUMPS != 0) {
        std::printf(" **** WARNING:  IN THE PHI-GRID SEARCHES%8d"
                    " BOUND STATE WAVEFUNCTIONS WERE NEEDED BEYOND THE"
                    " BOUND STATE ASYMPTOPIA.\n", NBUMPS);
        std::printf(" **** IN FUTURE RUNS INCREASE ASYMPTOPIA TO MORE THAN:"
                    "%15.2f (PROJECTILE)%15.2f (TARGET)\n", RPMAX, RTMAX);
    }

    goto L876;

L874:
    std::printf("\n**** ERROR IN OPENING OR WRITING SCRATCH FILE.\n"
                " **** ARRAYS WILL BE KEPT IN MEMORY IF POSSIBLE.\n");
    ISCRFL = 0;

L876:
    // debug removed
    JSIZE = (JCNT + ALLOCS.FACFR4 - 1) / ALLOCS.FACFR4;
    if (ISCRFL != 0) JSIZE = (NPPHI + ALLOCS.FACFR4 - 1) / ALLOCS.FACFR4;
    // debug removed
    LOCTRP_g = NALLOC("FIFODX  ", JSIZE);
    IPHIT_g  = NALLOC("PHIANG1 ", JSIZE);
    IPHIP_g  = NALLOC("PHIANG2 ", JSIZE);
    IPHI_g   = NALLOC("PHIANG  ", JSIZE);
    MAXCNT = JCNT;
    X0AV = X0AV / NRIROH;

    if (IPRNT >= 9) std::printf("     RI     RO               PHI V PHI  AT EACH X\n\n");

    LPHIT = ALLOCS.FACFR4 * LOCPTRS.Z[IPHIT_g]  - ALLOCS.FACFR4;
    LPHIP = ALLOCS.FACFR4 * LOCPTRS.Z[IPHIP_g]  - ALLOCS.FACFR4;
    LPHI  = ALLOCS.FACFR4 * LOCPTRS.Z[IPHI_g]   - ALLOCS.FACFR4;
    LTRP  = ALLOCS.FACFR4 * LOCPTRS.Z[LOCTRP_g] - ALLOCS.FACFR4;
    IXTOPZ = NPPHI;
    NNAIT = NAITKN;
    IBSTYP = 1;
    INADD = 0;
    PAS1SW = false;
    goto L810;  // start pass 2

    // =========================================================================
    // End of pass 2 — check asymptopia overruns
    // =========================================================================
L890:
    // debug removed
    if (NBUMPS != 0) {
        std::printf(" **** ERROR:%8d BOUND STATE WAVEFUNCTIONS WERE NEEDED"
                    " BEYOND THE BOUND STATE ASYMPTOPIA.\n", NBUMPS);
        std::printf(" **** IN FUTURE RUNS INCREASE ASYMPTOPIA TO MORE THAN:"
                    "%15.2f (PROJECTILE)%15.2f (TARGET)\n", RPMAX, RTMAX);
    }

    // Free bound state products and phi arrays
    if (IBSPAS == 0) LOCPTRS.Z[FORMF.IBDS[1]] = -LOCPTRS.Z[FORMF.IBDS[1]];
    if (IBSPAS == 0) LOCPTRS.Z[FORMF.IBDS[2]] = -LOCPTRS.Z[FORMF.IBDS[2]];
    LOCPTRS.Z[IPHIPT] = -LOCPTRS.Z[IPHIPT];
    LOCPTRS.Z[IPHIWT] = -LOCPTRS.Z[IPHIWT];
    if (NUCONL >= 2 && IBSPAS == 0) LOCPTRS.Z[BSBLK.JPOT] = -LOCPTRS.Z[BSBLK.JPOT];
    { int Ivphi = NAMLOC("VPHI    "); LOCPTRS.Z[Ivphi] = -LOCPTRS.Z[Ivphi]; }
    LOCPTRS.Z[LOGIC_g] = -LOCPTRS.Z[LOGIC_g];

    // =========================================================================
    // Print summary
    // =========================================================================
L900:
    if (!PRNTSW) goto L950;

    std::printf("\n\n0             SUMMARY OF THE INTEGRATION GRIDS\n"
        "0         NUM. PTS.   MAP TYPE   GAMMA    MINIMUM   \"MID. PT.\"   MAXIMUM\n"
        " (RI+RO)/2:%5d%11d%11.2f%10.2f%11.2f%11.2f"
        "     (H'S FOUND AT%4d (RI+RO)/2 VALUES)\n"
        "0                                        MIN. WIDTH   MAX. WIDTH\n"
        " RI-RO:%9d%11d%11.2f%11.3f%13.3f\n"
        "0                                          NPHIADD    MIN. START   AV. START     PHIMID\n"
        " COS(PHI):%6d%11d%11.2f%9d%16.5f%12.5f%12.4f\n",
        NPSUMI, MAPSUM, GAMSUM, SUMMIN, SUMMID, SUMMAX, NPSUM,
        NPDIF, MAPDIF, GAMDIF, VRNGMN, VRNGMX,
        NPPHI, MAPPHI, GAMPHI, NPHIAD, X0MIN, X0AV, PHIMID);

    if (IRECT == 1) std::printf("\n\"SW3ON\" WAS SPECIFIED:  (RI,RO) GRID IS RECTANGULAR.\n");

    std::printf("0DWCUTOFF =%12.2E\n"
        " MAXIMUM R'S USED FOR BOUND STATES:  PROJECTILE =%5.1f     TARGET =%5.1f\n"
        " MAXIMUM R'S USED FOR SCATTERING STATES:  INCOMING =%5.1f     OUTGOING =%5.1f\n",
        DWCUT, RPMAX, RTMAX, RIMAX_k, ROMAX_k);

    std::printf("0THERE ARE%6d POINTS IN THE (RI, RO) GRID ON WHICH THE H'S ARE FOUND,"
        "  AND%8d POINTS IN THE COMPLETE 3-D GRID.\n"
        " THERE ARE%6d POINTS IN THE (RI, RO) INTERPOLATED INTEGRATION GRID.\n"
        " AREA OF  SUM X DIF  INTEGRATION REGION =%9.3f FM**2.\n\n",
        NRIROH, MAXCNT, NRIROI_g, AREA);

    // =========================================================================
    // Save H's if requested (stub — binary I/O not implemented)
    // =========================================================================
L950:
    if (IHSAVE == 1) {
        // WRITE(1,...) binary H's file — stub for modern Ptolemy
        // goto L9950 on error
    }

    IRTN = 1;
    LOCPTRS.Z[IVPOLS] = -LOCPTRS.Z[IVPOLS];
    if (IHSAVE == 2) goto L1100;
    LOCPTRS.Z[IRIH] = -LOCPTRS.Z[IRIH];
    LOCPTRS.Z[IROH] = -LOCPTRS.Z[IROH];

L1100:
    TEMP = DWCUT;
    if (IHSAVE == 1) TEMP = std::min(1.0e-30, (double)DWCUT);
    RVRLIM_g = WVWMAX * TEMP;

    // Timing
    if (IHSAVE == 2 || !PINFSW) return;
    WOW = second() - TSTART;
    TEMP = WOW - TINTRP;
    U = 1000.0 * TINTRP / NINTRP;
    std::printf("\n        END OF GRID INITIALIZATION\n"
        "\nROUGH INTERPOLATION TIME FOR GRID SEARCHING:%54s%7.3f SECONDS\n"
        " OTHER TIME (INCLUDES FINAL INTERPOLATIONS):%54s%7.3f SECONDS\n"
        " TOTAL TIME:%54s%7.3f SECONDS\n"
        "\n%8d ROUGH INTERPOLATIONS MADE; AVERAGE TIME =%54s%7.3f MILLISECONDS\n\n",
        "", TINTRP, "", TEMP, "", WOW, NINTRP, "", U);
    return;

L9950:
    std::printf("\n**** ERROR IN READING OR WRITING H'S FILE.\n");
    return;

} // end GRDSET
