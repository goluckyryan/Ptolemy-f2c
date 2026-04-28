// inrdin_translated.cpp — INRDIN subroutine
// Translated from source.f lines 22212-22814
// Computes inelastic excitation radial integrals I(LI,LO,LX).
//
// INRDIN integrates U(R)*H(R)*U(R) where H is the full coupling Hamiltonian
// (nuclear + Coulomb), over the distorted scattering wavefunctions.
// Then adds the asymptotic correction from SUMMAX to infinity and
// subtracts the pure Coulomb piece to give the nuclear S-matrix element.
//
// 12/15/76 - BASED ON INELDC - S. P.
// 4/24/78 - INTERCHANGE LI AND LO - S.P.
// ...etc. (see Fortran header)

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include "ptolemy_alloc.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_forward.h"
#include "ptolemy_io.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <complex>
#include <algorithm>

extern double second();

// INDXJS(4,2) <=> JPMIN (INELCM)
// Fortran: INDXJS(I,J) = JPMIN + 4*(J-1) + (I-1)
// C++: INDXJS_f(I,J) = (&INELCM.JPMIN)[ 4*(J-1) + (I-1) ]
inline int INDXJS_f(int I, int J) {
    return (&INELCM.JPMIN)[ 4*(J-1) + (I-1) ];
}

void INRDIN(int& IRTN)
{
    //
    // COMPUTES THE INELASTIC EXCITATION RADIAL INTEGRALLS  I(LI,LO,LX)
    //
    // COMPUTES THE INELASTIC EXCITATION INTEGRALS
    //   INTEGRAL(0 TO SUMMAX) DR U(R) H U(R)
    // WHERE H IS THE FULL H
    // AND THEN COMPUTES
    //   I(NUCLEAR) = I(0 TO SUMMAX) + I(SUMMAX TO INF) - I(COULOMB)
    //

    // WARNING:  VAX VERSION MAY BE INACCURATE (NO COMPLEX*16).
    // IFACTR = (0, 1) = i
    std::complex<double> IFACTR(0.0, 1.0);

    // FOLLOWING ARE  -(I)**LX
    // PHASES(4):  -(i)^1 = -i = (0,-1)
    //             -(i)^2 = -(-1) = 1... wait — Fortran says:
    //   PHASES(4) / (0.D0,-1.D0), (-1.D0,0.D0), (0.D0,1.D0), (+1.D0,0.D0) /
    // These are indexed 1..4 (1-based), so PHASES[0]=(0,-1), [1]=(-1,0), [2]=(0,1), [3]=(1,0)
    static const std::complex<double> PHASES[4] = {
        {0.0, -1.0},   // PHASES(1)
        {-1.0, 0.0},   // PHASES(2)
        {0.0,  1.0},   // PHASES(3)
        {+1.0, 0.0}    // PHASES(4)
    };

    // --- COMMON block references (shorthand aliases) ---
    auto* Z      = LOCPTRS.Z;
    int&  FACFR4 = ALLOCS.FACFR4;

    // /FLOAT/
    double& ELAB   = FLOAT_common.ELAB;

    // /INTGER/
    int& LMIN    = INTGER.LMIN;
    int& LMAX    = INTGER.LMAX;
    int& IPRINT  = INTGER.IPRINT;

    // /INTRNL/
    int& IEXCIT  = INTRNL.IEXCIT;

    // /INELCM/
    int& LXMIN   = INELCM.LXMIN;
    int& LXMAX   = INELCM.LXMAX;
    int& NMLOLX  = INELCM.NMLOLX;
    int& NUMLX   = INELCM.NUMLX;
    int& NUMLIS  = INELCM.NUMLIS;
    int& ILIS    = INELCM.ILIS;
    int& NLX     = INELCM.NLX;
    int& NSPLI   = INELCM.NSPLI;
    int& ISMATR  = INELCM.ISMATR;
    int& ISMATI  = INELCM.ISMATI;
    int& ICL2FF  = INELCM.ICL2FF;
    int* IBETASp = INELCM.IBETAS;   // 1-based: IBETAS[1],IBETAS[2],IBETAS[3]
    int& IINDXS  = INELCM.IINDXS;
    int& IDENSW  = INELCM.IDENSW;

    // /KANDM/
    double& AKI    = KANDM.AKI;
    double& AKO    = KANDM.AKO;
    int*    ISMATS = KANDM.ISMATS;   // 1-based: ISMATS[1], ISMATS[2]
    int&    LOMOST = KANDM.LOMOST;

    // /WAVCOM/
    double* ES     = WAVCOM.ES;      // 1-based: ES[1], ES[2]
    int*    IRLVS  = WAVCOM.IRLVS;  // 1-based
    int*    IIMVS  = WAVCOM.IIMVS;  // 1-based
    int*    ICENTR = WAVCOM.ICENTR; // 1-based
    int*    NFS    = WAVCOM.NFS;    // 1-based
    int*    NGS    = WAVCOM.NGS;    // 1-based
    int*    NF1S   = WAVCOM.NF1S;   // 1-based
    int*    NG1S   = WAVCOM.NG1S;   // 1-based
    double* V0SORS = WAVCOM.V0SORS; // 1-based
    double* V0SOIS = WAVCOM.V0SOIS; // 1-based
    int*    ISORS  = WAVCOM.ISORS;  // 1-based
    int*    ISOIS  = WAVCOM.ISOIS;  // 1-based
    int&    IWAVR  = WAVCOM.IWAVR;
    int&    IWAVI  = WAVCOM.IWAVI;

    // /CNSTNT/
    double& PI     = CNSTNT.PI;

    // /FTIME/
    float*  TIMES  = FTIME.TIMES;   // 1-based: TIMES[1]..TIMES[8]

    // /HEDCOM/
    char*  REACT  = HEDCOM.REACT;   // 1-based: REACT[1]..REACT[45]
    char*  HEADER = HEDCOM.HEADER;  // 1-based: HEADER[1]..HEADER[65]

    // --- GRIDCM inelastic overlay (gridcmf) ---
    // The inelastic GRIDCM starts with: RCUEFF (real*8), GRDCM1(2) (real*8 x2),
    // then a sequence of integers.  In C++ the struct starts with
    // {double JACOB, STPINV, COSSTP; int IMSVAL, ...}
    // So:  GRIDCM.JACOB  = RCUEFF
    //      GRIDCM.STPINV = GRDCM1(1)
    //      GRIDCM.COSSTP = GRDCM1(2)
    // Then GRIDCM.IMSVAL onward = the integer fields below:
    //   IBINDX[0], MBINDX[1], NBINDX[2], NMFFAC[3], MXLXGS[4], NMBFAC[5],
    //   IINHR[6], IINHI[7], IRHSR[8], IRHSI[9],
    //   INUCH[10], ICOULH[11], ILOR[12], ILOI[13], ILIR[14], ILII[15],
    //   NPTMIN[16], NOFLO[17], NUMPT[18], NCRIT[19],
    //   ISMX[20], ISMSQ[21], IPADE[22], IHOMO[23],
    //   ISMTSV[24], IRPTS4[25], IRPTS[26],
    //   IGRDC4(1)[27], IGRDC4(2)[28], IGRDC4(3)[29], IGRDC4(4)[30], IGRDC4(5)[31],
    //   IDW[32], IGRDC5(1)[33], IGRDC5(2)[34],
    //   IABS1[35], IGRDC6(1)[36], IGRDC6(2)[37],
    //   GRDCM2 (real*8 -> occupies int slots [38] and [39]),
    //   IGRDC7(1)[40], IGRDC7(2)[41], IRWTS[42],
    //   IGRDC8(1)[43], IGRDC8(2)[44],
    //   ICL1FF[45], ICL1FG[46], ICL1GF[47], ICL1GG[48]
    double& RCUEFF = GRIDCM.JACOB;         // = RCUEFF in inelastic GRIDCM
    int* GI = reinterpret_cast<int*>(&GRIDCM.IMSVAL);  // base of integer region
    //
    // Named aliases for INRDIN-specific inelastic GRIDCM fields:
    int& INUCH   = GI[10];   // pointer to nuclear H array
    int& ICOULH  = GI[11];   // pointer to Coulomb H array
    // ILOR, ILOI, ILIR, ILII are also in GRIDCM struct directly at same positions
    // GI[12]=ILOR, GI[13]=ILOI, GI[14]=ILIR, GI[15]=ILII
    // But we can also access them via GRIDCM.ILOR etc. since the mapping is same
    int& ILOR_i  = GRIDCM.ILOR;
    int& ILOI_i  = GRIDCM.ILOI;
    int& ILIR_i  = GRIDCM.ILIR;
    int& ILII_i  = GRIDCM.ILII;
    int& NOFLO_i = GRIDCM.NOFLO;   // GI[17]
    int& NUMPT   = GRIDCM.NRIROI;  // GI[18] -- inelastic NUMPT maps here
    int& NCRIT_i = GRIDCM.NCRIT;   // GI[19]
    int& IRPTS4_i= GI[25];         // IRPTS4 (1D grid as REAL*4)
    int& IRPTS_i = GI[26];         // IRPTS  (1D grid as REAL*8)
    int& IDW_i   = GRIDCM.IDW;     // GI[32]
    int& IABS1_i = GRIDCM.IABS1;   // GI[35]
    int& IRWTS_i = GI[42];         // IRWTS
    int& ICL1FF_i= GI[45];         // ICL1FF
    int& ICL1FG_i= GI[46];         // ICL1FG
    int& ICL1GF_i= GI[47];         // ICL1GF
    int& ICL1GG_i= GI[48];         // ICL1GG

    // --- Local variables ---
    int PBUGSW, PINFSW;   // LOGICAL (0/1)
    int NXLISW;           // LOGICAL

    float TSTART;

    int IPRNT;
    int ISPUN;
    int LNUCH, LCOULH, LRPTS, LRPTS4;
    int LLIR, LLII, LLOR, LLOI;
    int LDW, LLIS;
    int LSMATR, LSMATI;
    int LCL1FF, LCL1FG, LCL1GF, LCL1GG;
    int LCL2FF;
    int LSIN, LSOUT;
    int LABS1;
    int LBNRAT, LBRAT;
    int LWAVR, LWAVI;
    int LVRL1, LVRL2, LVIM1, LVIM2;
    int LCENT1, LCENT2;
    int LINDXS;

    double FACTOR;

    int LRANGE;
    int LIMIN, LIL, LIPRTY;
    int ILI, LI;
    int LOMIN, LOMAX;
    int ITRK, LOBASE;
    int LO, LX;
    int IH, IU;
    int NSKIP;
    int NLINE;
    int I, I1, I2;
    int K, KBASE, LDELMN, KOFFS;
    int ITRKLO;
    int LOMNMN, LOMXMX;
    int LOMN, LOMX;
    int LDWBAS;

    double RVAL, RAT;
    double FIR, FII, FOR_, FOI;
    double DWR, DWI, DWLOM;
    double DWRIOS, DWRIOL, DWRIOC, DWCONT;
    double DWRIOS_loc, DWRIOL_loc, DWRIOC_loc, DWCONT_loc;
    double HNUCR, HNUCI, HCOUL;
    double HR, HI;
    double TERMR, TERMI;

    double BETRAT;
    double C, FFI;
    double AMPCOM, AMPNUC, AMPTOT;
    double PHSCOM, PHSNUC, PHSTOT, THISR;

    std::complex<double> INUC, ICOMP, IRTOIN, ITOTAL;
    std::complex<double> SIN_c, SOUT_c;   // SIN/SOUT clash with C math; use _c suffix
    std::complex<double> PHASE;

    double TOTTIM, WAVTIM, OTHTIM, RTIM1, RTIM2;

    // =========================================================================
    // Begin executable code
    // =========================================================================

    TSTART = (float)second();
    IRTN = 0;

    // DO 14 I=1, 3; TIMES(I)=0.0
    for (int iii = 1; iii <= 3; iii++) TIMES[iii] = 0.0f;

    IPRNT = IPRINT % 10;
    // PBUGSW = IPRNT .GE. 4
    PBUGSW = (IPRNT >= 4) ? 1 : 0;
    // PINFSW = IPRNT .GE. 2
    PINFSW = (IPRNT >= 2) ? 1 : 0;

    //
    // write (6, 23) REACT, ELAB, HEADER
    // FORMAT ( '1', T60, 'P T O L E M Y' /
    //   49X, 'COMPUTATION OF INELASTIC S-MATRIX ELEMENTS' /
    //   '0', 45A1, 'ELAB =', F7.2, ' MEV', 5X, 65A1 / )
    //
    std::printf("1%58sP T O L E M Y\n", "");
    std::printf("%49sCOMPUTATION OF INELASTIC S-MATRIX ELEMENTS\n", "");
    {
        char react45[46]; react45[45] = '\0';
        char header65[66]; header65[65] = '\0';
        for (int ii = 0; ii < 45; ii++) react45[ii]  = REACT[ii+1];
        for (int ii = 0; ii < 65; ii++) header65[ii] = HEADER[ii+1];
        std::printf("0%.45sELAB =%7.2f MEV     %.65s\n\n",
                    react45, ELAB, header65);
    }

    NLINE = 1000;

    //
    // NOTE -- NO SPIN ORBIT FORCE
    //
    ISPUN = 1;

    //
    // CALL SETLOG ( 2*(LOMOST+LXMAX) )
    //
    SETLOG(2 * (LOMOST + LXMAX));

    //
    // HERE WE START THE REAL GUTS OF THE CALCULATION.
    // WE WILL HAVE NO MORE IALLOC CALLS SO GET ADDRESSES
    //
    // LNUCH  = Z(INUCH) - 2
    LNUCH  = Z[INUCH]  - 2;
    // LCOULH = Z(ICOULH) - 1
    LCOULH = Z[ICOULH] - 1;
    // LRPTS  = Z(IRPTS) - 1
    LRPTS  = Z[IRPTS_i] - 1;
    // LRPTS4 = FACFR4*(Z(IRPTS4)-1) + 1
    LRPTS4 = FACFR4 * (Z[IRPTS4_i] - 1) + 1;
    // LLIR = FACFR4*Z(ILIR) - FACFR4
    LLIR   = FACFR4 * Z[ILIR_i] - FACFR4;
    // LLII = FACFR4*Z(ILII) - FACFR4
    LLII   = FACFR4 * Z[ILII_i] - FACFR4;
    // LLOR = FACFR4*Z(ILOR) - FACFR4
    LLOR   = FACFR4 * Z[ILOR_i] - FACFR4;
    // LLOI = FACFR4*Z(ILOI) - FACFR4
    LLOI   = FACFR4 * Z[ILOI_i] - FACFR4;
    // LDW = Z(IDW)
    LDW    = Z[IDW_i];
    // LLIS = FACFR4*Z(ILIS) - FACFR4
    LLIS   = FACFR4 * Z[ILIS] - FACFR4;
    // LSMATR = Z(ISMATR) - 1
    LSMATR = Z[ISMATR] - 1;
    // LSMATI = Z(ISMATI) - 1
    LSMATI = Z[ISMATI] - 1;
    // LCL1FF = Z(ICL1FF) - 1
    LCL1FF = Z[ICL1FF_i] - 1;
    // LCL1FG = Z(ICL1FG) - 1
    LCL1FG = Z[ICL1FG_i] - 1;
    // LCL1GF = Z(ICL1GF) - 1
    LCL1GF = Z[ICL1GF_i] - 1;
    // LCL1GG = Z(ICL1GG) - 1
    LCL1GG = Z[ICL1GG_i] - 1;
    // LCL2FF = Z(ICL2FF) - 1
    LCL2FF = Z[ICL2FF] - 1;
    // LSIN  = Z(ISMATS(1))
    LSIN   = Z[ISMATS[1]];
    // LSOUT = Z(ISMATS(2))
    LSOUT  = Z[ISMATS[2]];
    // LABS1 = Z(IABS1) - 1
    LABS1  = Z[IABS1_i] - 1;
    // LBNRAT = Z(IBETAS(2)) - LXMIN/2
    LBNRAT = Z[IBETASp[2]] - LXMIN/2;
    // LBRAT  = Z(IBETAS(3)) - LXMIN/2
    LBRAT  = Z[IBETASp[3]] - LXMIN/2;
    // LWAVR = Z(IWAVR)
    LWAVR  = Z[IWAVR];
    // LWAVI = Z(IWAVI)
    LWAVI  = Z[IWAVI];
    // LVRL1 = Z(IRLVS(1))
    LVRL1  = Z[IRLVS[1]];
    // LVRL2 = Z(IRLVS(2))
    LVRL2  = Z[IRLVS[2]];
    // LVIM1 = Z(IIMVS(1))
    LVIM1  = Z[IIMVS[1]];
    // LVIM2 = Z(IIMVS(2))
    LVIM2  = Z[IIMVS[2]];
    // LCENT1 = Z(ICENTR(1))
    LCENT1 = Z[ICENTR[1]];
    // LCENT2 = Z(ICENTR(2))
    LCENT2 = Z[ICENTR[2]];
    // LINDXS = FACFR4*Z(IINDXS) - FACFR4
    LINDXS = FACFR4 * Z[IINDXS] - FACFR4;

    // FACTOR = DSQRT( AKI*AKO/(ES(1)*ES(2)*PI) )
    FACTOR = DSQRT(AKI * AKO / (ES[1] * ES[2] * PI));
    // IF ( IDENSW )  FACTOR = SQRT(2.D0)*FACTOR
    if (IDENSW) FACTOR = DSQRT(2.0) * FACTOR;

    //
    // FORCE THE OVERLAYS TO GET THE THREEJ COMMON BLOCKS IN.
    // CALLED REPETIVELY BELOW.  HOWEVER, NOTHING WILL HAPPEN
    // FROM NOW ON THAT WOULD CAUSE THE OVERLAYS TO CHANGE.
    //
    DUMMY1();

    // LRANGE = NOFLO - 1
    LRANGE = NOFLO_i - 1;

    //
    // THE LOOPS BEGIN NOW
    //
    //
    // LOOP OVER LI
    //
    // TWO PASSES OVER LI LOOP NOW; EVEN LI FIRST, THEN ODD
    //
    // LIMIN = LMIN + 1
    LIMIN = LMIN + 1;
    // LIL   = LMIN
    LIL   = LMIN;
    // IF(MOD(LMIN,2) .EQ. 1 ) GO TO 150
    if (MOD(LMIN, 2) == 1) goto L150;
    // LIMIN = LMIN
    LIMIN = LMIN;
    // LIL   = LMIN+1
    LIL   = LMIN + 1;

    // LIPRTY = 1 FOR LI EVEN,  = 2 FOR LI ODD

L150:
    // DO 989  LIPRTY = 1, 2
    for (LIPRTY = 1; LIPRTY <= 2; LIPRTY++) {

        // IF(LIMIN.GT.LMAX) GO TO 980
        if (LIMIN > LMAX) goto L980;
        // NXLISW = .FALSE.
        NXLISW = 0;
        // ILI = 1
        ILI = 1;

        // DO 959 LI=LIMIN,LMAX,2
        for (LI = LIMIN; LI <= LMAX; LI += 2) {

            //
            // GET THE REQUIRED SCATTERING WAVEFUNCTIONS.
            // WE ALWAYS DO THIS EVEN IF LI WON'T BE USED SO AS TO GET THE
            // ELASTIC S-MATRIX ELEMENTS.
            //
            LOMIN = LI - LRANGE;
            LOMAX = LI + LRANGE;

            // IF ( NXLISW )  GO TO 250
            if (NXLISW) goto L250;

            //
            // FOR FIRST LI, FILL UP THE OUT-CHANNEL WAVEFUNCTION ARRAY
            //
            ITRK = 0;

            // DO 239  LO=LOMIN,LOMAX,2
            for (LO = LOMIN; LO <= LOMAX; LO += 2) {
                WAVELJ(LO, ISPUN, 2, NUMPT, ALLOC4_base(LRPTS4),
                       ALLOC4_base(LLOR + ITRK + 1), ALLOC4_base(LLOI + ITRK + 1),
                       ALLOC_base(LWAVR), ALLOC_base(LWAVI),
                       ALLOC_base(LVRL2), ALLOC_base(LVIM2), ALLOC_base(LCENT2));
                //  239  ITRK=ITRK+NUMPT
                ITRK = ITRK + NUMPT;
            }

            ITRK   = 0;
            LOBASE = LOMIN;
            // NXLISW = .TRUE.
            NXLISW = 1;
            goto L280;

            //
            // ON SUBSEQUENT LI'S WE NEED ONE NEW OUT-CHANNEL WAVEFUNC.
            //
L250:
            LO = LI + LRANGE;
            WAVELJ(LO, ISPUN, 2, NUMPT, ALLOC4_base(LRPTS4),
                   ALLOC4_base(LLOR + ITRK + 1), ALLOC4_base(LLOI + ITRK + 1),
                   ALLOC_base(LWAVR), ALLOC_base(LWAVI),
                   ALLOC_base(LVRL2), ALLOC_base(LVIM2), ALLOC_base(LCENT2));
            // IF ( ITRK .EQ. 0 )  LOBASE = LO
            if (ITRK == 0) LOBASE = LO;
            ITRK = ITRK + NUMPT;
            // IF(ITRK.GT.NCRIT) ITRK=0
            if (ITRK > NCRIT_i) ITRK = 0;

            //
            // FOR EACH LI, WE MUST GET THE O-CHANNEL WAVE FUNCTION
            //
L280:
            WAVELJ(LI, ISPUN, 1, NUMPT, ALLOC4_base(LRPTS4),
                   ALLOC4_base(LLIR + 1), ALLOC4_base(LLII + 1),
                   ALLOC_base(LWAVR), ALLOC_base(LWAVI),
                   ALLOC_base(LVRL1), ALLOC_base(LVIM1), ALLOC_base(LCENT1));

            //
            // WAVEFUNCTIONS ALL FOUND, WILL WE DO THIS VALUE OF LI
            //
L290:
            {
                int diff = LI - ILLOC(LLIS + ILI);
                if (diff < 0) goto L959;
                if (diff == 0) goto L300;
                // diff > 0
                if (ILI >= NUMLIS) goto L959;
                ILI = ILI + 1;
                goto L290;
            }

            //
            // THIS IS AN LI TO DO
            //
L300:
            // LOMNMN = MAX0( MOD(LOMAX,2) , LOMIN )
            LOMNMN = MAX0(MOD(LOMAX, 2), LOMIN);
            LOMXMX = LOMAX;
            NSKIP  = 0;
            DWRIOS_loc = 10.0;
            DWRIOL_loc = 0.0;
            DWRIOC_loc = 0.0;
            DWCONT_loc = 0.0;

            //
            // IN ABS1 WE ACCUMULATE INTEGRAL |INTEGRAND| FOR THIS LI.
            //
            // DO 329  IH = 1, NMLOLX
            for (IH = 1; IH <= NMLOLX; IH++) {
                ALLOC(LABS1 + IH) = 0.0;
            }
            // 329  CONTINUE

            //
            // LOOP OVER R-VALUES AND COMPUTE THE H'S
            //
            // DO 489  IU = 1, NUMPT
            for (IU = 1; IU <= NUMPT; IU++) {
                // RVAL = ALLOC(LRPTS+IU)
                RVAL = ALLOC(LRPTS + IU);
                // RAT = 1/RVAL
                RAT = 1.0 / RVAL;
                // IF ( RVAL .LT. RCUEFF )  RAT = RVAL/RCUEFF**2
                if (RVAL < RCUEFF) RAT = RVAL / (RCUEFF * RCUEFF);
                // RAT = RAT**2
                RAT = RAT * RAT;

                //
                // COMPUTE IMAG AND REAL PARTS OF PRODUCT OF THE DISTORTED
                // SCATTERING WAVES.
                //
                // FIR=ALLOC4(LLIR+IU)
                FIR = (double)ALLOC4(LLIR + IU);
                // FII=ALLOC4(LLII+IU)
                FII = (double)ALLOC4(LLII + IU);
                // IF ( DABS(FIR)+DABS(FII) .LT. 1.D-20 )  GO TO 480
                if (DABS(FIR) + DABS(FII) < 1.0e-20) goto L480;

                //
                // UNDO THE WRAP AROUND STORAGE AND COMPUTE ALL PSI(A)*PSI(B)'S
                // FOR THIS (RI, RO).
                //
                // DO 439  LO = LOMNMN, LOMXMX, 2
                for (LO = LOMNMN; LO <= LOMXMX; LO += 2) {
                    // ITRKLO = NUMPT*( (LO-LOBASE)/2 )
                    ITRKLO = NUMPT * ((LO - LOBASE) / 2);
                    // IF ( ITRKLO .LT. 0 )  ITRKLO = ITRKLO + NUMPT+NCRIT
                    if (ITRKLO < 0) ITRKLO = ITRKLO + NUMPT + NCRIT_i;
                    // FOR = ALLOC4(LLOR+ITRKLO+IU)
                    FOR_ = (double)ALLOC4(LLOR + ITRKLO + IU);
                    // FOI = ALLOC4(LLOI+ITRKLO+IU)
                    FOI  = (double)ALLOC4(LLOI + ITRKLO + IU);
                    // IF ( DABS(FOR) .LT. 1.D-20 )  FOR = 0
                    if (DABS(FOR_) < 1.0e-20) FOR_ = 0.0;
                    // IF ( DABS(FOI) .LT. 1.D-20 )  FOI = 0
                    if (DABS(FOI)  < 1.0e-20) FOI  = 0.0;
                    // DWR=FOR*FIR-FOI*FII
                    DWR  = FOR_ * FIR - FOI * FII;
                    // DWI=FOR*FII+FOI*FIR
                    DWI  = FOR_ * FII + FOI * FIR;
                    // DWLOM=DMAX1(DABS(DWR),DABS(DWI))
                    DWLOM = DMAX1(DABS(DWR), DABS(DWI));
                    // IF(DWLOM.GT.DWRIOL) DWRIOL=DWLOM
                    if (DWLOM > DWRIOL_loc) DWRIOL_loc = DWLOM;
                    // IF((DWLOM.LT.DWRIOS).AND.(DWLOM.NE.0)) DWRIOS=DWLOM
                    if ((DWLOM < DWRIOS_loc) && (DWLOM != 0.0)) DWRIOS_loc = DWLOM;
                    // DWRIOC=DWRIOC+DWLOM
                    DWRIOC_loc = DWRIOC_loc + DWLOM;
                    // ALLOC(LDW+LO-LOMNMN) = DWR
                    ALLOC(LDW + LO - LOMNMN) = DWR;
                    // ALLOC(LDW+LO-LOMNMN+1) = DWI
                    ALLOC(LDW + LO - LOMNMN + 1) = DWI;
                    // DWCONT = DWCONT + 1
                    DWCONT_loc = DWCONT_loc + 1.0;
                }
                //  439  CONTINUE

                // HNUCR = ALLOC(LNUCH+2*IU)
                HNUCR = ALLOC(LNUCH + 2*IU);
                // HNUCI = ALLOC(LNUCH+2*IU+1)
                HNUCI = ALLOC(LNUCH + 2*IU + 1);
                // HCOUL = ALLOC(LCOULH+IU)
                HCOUL = ALLOC(LCOULH + IU);

                //
                // NOW ADD THESE H'S INTO THE APPROPRIATE I(LI,LO,LX)
                // ALONG WITH THE PRODUCT OF THE DISTORTED WAVES
                // NOTE THAT THE BETA ARRAYS ALSO CONTAIN
                //      CLEBSCH( JIN LX K 0 : JOUT K )
                //
                IH = 0;
                // DO 479  LX = LXMIN, LXMAX, 2
                for (LX = LXMIN; LX <= LXMAX; LX += 2) {
                    // HR = HNUCR*ALLOC(LBNRAT+LX/2) + HCOUL*ALLOC(LBRAT+LX/2)
                    HR = HNUCR * ALLOC(LBNRAT + LX/2) + HCOUL * ALLOC(LBRAT + LX/2);
                    // HI = HNUCI*ALLOC(LBNRAT+LX/2)
                    HI = HNUCI * ALLOC(LBNRAT + LX/2);
                    // HCOUL = RAT * HCOUL
                    HCOUL = RAT * HCOUL;
                    // LOMN = IABS(LI-LX)
                    LOMN = IABS(LI - LX);
                    // LOMX = LI+LX
                    LOMX = LI + LX;
                    // K = LX+1 + NLX*INDXJS(4,3-IEXCIT)*(LX-INDXJS(3,IEXCIT)/2)
                    K = LX + 1 + NLX * INDXJS_f(4, 3-IEXCIT) * (LX - INDXJS_f(3, IEXCIT)/2);
                    // KBASE = ILLOC(LINDXS+3*K-2)
                    KBASE = ILLOC(LINDXS + 3*K - 2);
                    // LDELMN = ILLOC(LINDXS+3*K-1)
                    LDELMN = ILLOC(LINDXS + 3*K - 1);
                    // DO 479  LO = LOMN, LOMX, 2
                    for (LO = LOMN; LO <= LOMX; LO += 2) {
                        // KOFFS = KBASE + (LO-LI-LDELMN)/2
                        KOFFS = KBASE + (LO - LI - LDELMN) / 2;
                        // I = (ILI-1)*NSPLI + KOFFS
                        I = (ILI - 1) * NSPLI + KOFFS;
                        // LDWBAS = LDW-LOMNMN+LO
                        LDWBAS = LDW - LOMNMN + LO;
                        // IH = IH + 1
                        IH = IH + 1;
                        // TERMR = HR*ALLOC(LDWBAS) - HI*ALLOC(LDWBAS+1)
                        TERMR = HR * ALLOC(LDWBAS) - HI * ALLOC(LDWBAS + 1);
                        // ALLOC(LSMATR+I) = ALLOC(LSMATR+I) + TERMR
                        ALLOC(LSMATR + I) = ALLOC(LSMATR + I) + TERMR;
                        // TERMI = HR*ALLOC(LDWBAS+1) + HI*ALLOC(LDWBAS)
                        TERMI = HR * ALLOC(LDWBAS + 1) + HI * ALLOC(LDWBAS);
                        // ALLOC(LSMATI+I) = ALLOC(LSMATI+I) + TERMI
                        ALLOC(LSMATI + I) = ALLOC(LSMATI + I) + TERMI;

                        //
                        // COMPUTE INDICATION OF LOSS OF SIGNIFICANCE:
                        //   ABS1 = INTEGRAL(R) |INTEGRAND|
                        //
                        ALLOC(LABS1 + IH) = ALLOC(LABS1 + IH) + DABS(TERMR) + DABS(TERMI);
                    }
                    //  479  CONTINUE (inner LO loop shared with LX loop below)
                }

                goto L489;

                //
                // WE COME HERE WHEN SKIPPING H CALCULATION DUE TO THE PRODUCT
                // OF THE DISTORTED WAVES BEING SMALL.
                //
L480:
                NSKIP = NSKIP + 1;

L489:           ;
            }
            //  489  CONTINUE

            //
            // END OF INTEGRATION
            //
            // DWRIOC=DWRIOC/DWCONT
            if (DWCONT_loc > 0.0) DWRIOC_loc = DWRIOC_loc / DWCONT_loc;
            // IF (PBUGSW)  WRITE(6, 533)DWRIOS,DWRIOL,DWRIOC
            // 533  FORMAT(' SMALLEST, LARGEST AND AVERAGE DABS(DWRIRO)=', 3G15.5 / )
            if (PBUGSW) {
                std::printf(" SMALLEST, LARGEST AND AVERAGE DABS(DWRIRO)= %15.5G%15.5G%15.5G\n\n",
                            DWRIOS_loc, DWRIOL_loc, DWRIOC_loc);
            }

            // IF ( .NOT. PINFSW )  GO TO 550
            if (!PINFSW) goto L550;
            NLINE = NLINE + NMLOLX + 1;
            // IF ( NLINE .LE. 58 )  GO TO 550
            if (NLINE <= 58) goto L550;
            // IF ( NLINE .LT. 1000 )  write (6, 23) REACT, ELAB, HEADER
            if (NLINE < 1000) {
                char react45[46]; react45[45] = '\0';
                char header65[66]; header65[65] = '\0';
                for (int ii = 0; ii < 45; ii++) react45[ii]  = REACT[ii+1];
                for (int ii = 0; ii < 65; ii++) header65[ii] = HEADER[ii+1];
                std::printf("1%58sP T O L E M Y\n", "");
                std::printf("%49sCOMPUTATION OF INELASTIC S-MATRIX ELEMENTS\n", "");
                std::printf("0%.45sELAB =%7.2f MEV     %.65s\n\n",
                            react45, ELAB, header65);
            }
            // write (6, 538)
            // 538  FORMAT ( '  L   L  LX', T25, 'INTEGRAL(0, SUMMAX)',
            //   T55, 'COULOMB', T70, 'NUCLEAR PART', T94, 'TOTAL AMPLITUDE' /
            //   '  IN OUT', 11X, 'MAG.', 6X, 'PHASE', 4X, 'CANCELLATION',
            //   T54, 'AMPLITUDE',
            //   T68, 2( 'MAG.', 6X, 'PHASE', 9X ) / )
            std::printf("  L   L  LX%14sINTEGRAL(0, SUMMAX)%29sCOULOMB%12sNUCLEAR PART%21sTOTAL AMPLITUDE\n",
                        "", "", "", "");
            std::printf("  IN OUT%11sMAG.      PHASE    CANCELLATION%13sAMPLITUDE%7sMAG.      PHASE         MAG.      PHASE\n",
                        "", "", "");
            NLINE = 7 + NMLOLX;

            //
            // ADD ON PIECE FROM SUMMAX TO INFINITY,
            // MULTIPLY IN CLEBSH GORDEN FACTOR,
            // COMPUTE CANCELLATIONS,
            // AND SEPARATE INTO NUCLEAR AND COULOMB PIECES
            // ALSO PUT IN THE I**(-LX-1)
            //
L550:
            // SIN = DCMPLX( ALLOC(LSIN+2*LI), ALLOC(LSIN+2*LI+1) )
            SIN_c = std::complex<double>(ALLOC(LSIN + 2*LI), ALLOC(LSIN + 2*LI + 1));
            IH = 0;
            // DO 559  LX = LXMIN, LXMAX, 2
            for (LX = LXMIN; LX <= LXMAX; LX += 2) {
                // LOMN = IABS(LI-LX)
                LOMN = IABS(LI - LX);
                // LOMX = LI+LX
                LOMX = LI + LX;
                // BETRAT = ALLOC(LBRAT+LX/2)
                BETRAT = ALLOC(LBRAT + LX/2);
                // PHASE = PHASES(iand(LX,3)+1)
                // Fortran: iand(LX,3) gives LX mod 4 (0-based), +1 gives 1-based index
                PHASE = PHASES[IAND(LX, 3)];  // C++ 0-based: PHASES[0..3]
                // K = LX+1 + NLX*INDXJS(4,3-IEXCIT)*(LX-INDXJS(3,IEXCIT)/2)
                K = LX + 1 + NLX * INDXJS_f(4, 3-IEXCIT) * (LX - INDXJS_f(3, IEXCIT)/2);
                // KBASE = ILLOC(LINDXS+3*K-2)
                KBASE = ILLOC(LINDXS + 3*K - 2);
                // LDELMN = ILLOC(LINDXS+3*K-1)
                LDELMN = ILLOC(LINDXS + 3*K - 1);
                // DO 559  LO = LOMN, LOMX, 2
                for (LO = LOMN; LO <= LOMX; LO += 2) {
                    // KOFFS = KBASE + (LO-LI-LDELMN)/2
                    KOFFS = KBASE + (LO - LI - LDELMN) / 2;
                    // IH = IH + 1
                    IH = IH + 1;
                    // I1 = (ILI-1)*NSPLI + KOFFS
                    I1 = (ILI - 1) * NSPLI + KOFFS;
                    // I2 = (LI-LMIN)*NSPLI + KOFFS
                    I2 = (LI - LMIN) * NSPLI + KOFFS;
                    // C = FACTOR *  DABS( CLEBSH( 2*LI, 2*LX, 0, 0, 2*LO, 0 ) )
                    C = FACTOR * DABS(CLEBSH(2*LI, 2*LX, 0, 0, 2*LO, 0));
                    // ICOMP = C*DCMPLX( ALLOC(LSMATR+I1), ALLOC(LSMATI+I1) )
                    ICOMP = C * std::complex<double>(ALLOC(LSMATR + I1), ALLOC(LSMATI + I1));
                    // SOUT = DCMPLX( ALLOC(LSOUT+2*LO), ALLOC(LSOUT+2*LO+1) )
                    SOUT_c = std::complex<double>(ALLOC(LSOUT + 2*LO), ALLOC(LSOUT + 2*LO + 1));
                    // IRTOIN = .25*C*BETRAT*( (1+SOUT)*(1+SIN)*ALLOC(LCL1FF+I1)
                    //   - (1-SOUT)*(1-SIN)*ALLOC(LCL1GG+I1)
                    //  + IFACTR*( (1+SOUT)*(1-SIN)*ALLOC(LCL1FG+I1)
                    //           + (1-SOUT)*(1+SIN)*ALLOC(LCL1GF+I1) ) )
                    IRTOIN = 0.25 * C * BETRAT * (
                        (1.0 + SOUT_c) * (1.0 + SIN_c) * ALLOC(LCL1FF + I1)
                        - (1.0 - SOUT_c) * (1.0 - SIN_c) * ALLOC(LCL1GG + I1)
                        + IFACTR * ((1.0 + SOUT_c) * (1.0 - SIN_c) * ALLOC(LCL1FG + I1)
                                   + (1.0 - SOUT_c) * (1.0 + SIN_c) * ALLOC(LCL1GF + I1))
                    );
                    // FFI = C*BETRAT*ALLOC(LCL2FF+I2)
                    FFI = C * BETRAT * ALLOC(LCL2FF + I2);
                    // ITOTAL = ICOMP + IRTOIN
                    ITOTAL = ICOMP + IRTOIN;
                    // INUC = ITOTAL - FFI
                    INUC = ITOTAL - FFI;
                    // INUC = PHASE*INUC
                    INUC = PHASE * INUC;
                    // ALLOC(LSMATR+I1) = DREAL(INUC)
                    ALLOC(LSMATR + I1) = DREAL(INUC);
                    // ALLOC(LSMATI+I1) = DIMAG(INUC)
                    ALLOC(LSMATI + I1) = DIMAG(INUC);
                    // IF ( .NOT. PINFSW )  GO TO 559
                    if (!PINFSW) goto L559;
                    // ICOMP = PHASE*ICOMP
                    ICOMP  = PHASE * ICOMP;
                    // ITOTAL = PHASE*ITOTAL
                    ITOTAL = PHASE * ITOTAL;
                    // AMPCOM = CDABS(ICOMP)
                    AMPCOM = std::abs(ICOMP);
                    // THISR = C*ALLOC(LABS1+IH)/AMPCOM
                    THISR  = (AMPCOM > 0.0) ? C * ALLOC(LABS1 + IH) / AMPCOM : 0.0;
                    // AMPNUC = CDABS(INUC)
                    AMPNUC = std::abs(INUC);
                    // AMPTOT = CDABS(ITOTAL)
                    AMPTOT = std::abs(ITOTAL);
                    // PHSCOM = ATAN2( DIMAG(ICOMP), DREAL(ICOMP) )
                    PHSCOM = ATAN2(DIMAG(ICOMP), DREAL(ICOMP));
                    // PHSNUC = ATAN2( DIMAG(INUC), DREAL(INUC) )
                    PHSNUC = ATAN2(DIMAG(INUC), DREAL(INUC));
                    // PHSTOT = ATAN2( DIMAG(ITOTAL), DREAL(ITOTAL) )
                    PHSTOT = ATAN2(DIMAG(ITOTAL), DREAL(ITOTAL));

                    // write (6, 553)  LI, LO, LX, AMPCOM, PHSCOM,
                    //   THISR, FFI, AMPNUC, PHSNUC, AMPTOT, PHSTOT
                    // 553 FORMAT ( 2I4, I3, G16.4, G11.4, F9.2, G16.4,
                    //   2( G13.4, G11.4 ) )
                    std::printf("%4d%4d%3d%16.4G%11.4G%9.2f%16.4G%13.4G%11.4G%13.4G%11.4G\n",
                                LI, LO, LX, AMPCOM, PHSCOM,
                                THISR, FFI, AMPNUC, PHSNUC, AMPTOT, PHSTOT);

L559:               ;
                }
            }
            // 559  CONTINUE

            // IF ( PINFSW )  write (6, 563)
            // 563  FORMAT ( ' ' )
            if (PINFSW) std::printf(" \n");

L959:       ;
        }  // end DO 959 LI=LIMIN,LMAX,2

        //
        // END OF THE LI LOOP FOR A GIVEN PARITY OF LI
        // IF LI IS EVEN HAVE ANOTHER BASH AT LI LOOP WITH ODD LI
        //
L980:
        // LIMIN=LIL
        LIMIN = LIL;

    }  // end DO 989 LIPRTY=1,2

    //
    // END OF THE 2 LI LOOPS (EVEN AND ODD PARITY OF LI)
    //
    // FREE WORK AREAS
    //
    // IF ( MOD(IPRINT/100,10) .GE. 1 )  GO TO 1500
    if (MOD(IPRINT / 100, 10) >= 1) goto L1500;

    Z[INUCH]    = -Z[INUCH];
    Z[ICOULH]   = -Z[ICOULH];
    Z[ILIR_i]   = -Z[ILIR_i];
    Z[ILII_i]   = -Z[ILII_i];
    Z[ILOR_i]   = -Z[ILOR_i];
    Z[ILOI_i]   = -Z[ILOI_i];
    Z[IDW_i]    = -Z[IDW_i];
    Z[IRPTS_i]  = -Z[IRPTS_i];
    Z[IRWTS_i]  = -Z[IRWTS_i];
    Z[IRPTS4_i] = -Z[IRPTS4_i];
    Z[IABS1_i]  = -Z[IABS1_i];
    Z[ICL1FF_i] = -Z[ICL1FF_i];
    Z[ICL1FG_i] = -Z[ICL1FG_i];
    Z[ICL1GF_i] = -Z[ICL1GF_i];
    Z[ICL1GG_i] = -Z[ICL1GG_i];

    //
    // THE AREAS FOR WAVELJ - ONE FOR INCOMING AND ONE FOR OUTGOING
    //
    // DO 1099  I = 1, 2
    for (I = 1; I <= 2; I++) {
        Z[NFS[I]]    = -Z[NFS[I]];
        Z[NGS[I]]    = -Z[NGS[I]];
        Z[NF1S[I]]   = -Z[NF1S[I]];
        Z[NG1S[I]]   = -Z[NG1S[I]];
        Z[IRLVS[I]]  = -Z[IRLVS[I]];
        Z[IIMVS[I]]  = -Z[IIMVS[I]];
        Z[ICENTR[I]] = -Z[ICENTR[I]];
        // IF ( V0SORS(I) .NE. 0 )  Z(ISORS(I)) = -Z(ISORS(I))
        if (V0SORS[I] != 0.0) Z[ISORS[I]] = -Z[ISORS[I]];
        // IF ( V0SOIS(I) .NE. 0 )  Z(ISOIS(I)) = -Z(ISOIS(I))
        if (V0SOIS[I] != 0.0) Z[ISOIS[I]] = -Z[ISOIS[I]];
    }
    // 1099  CONTINUE
    Z[IWAVR] = -Z[IWAVR];
    Z[IWAVI] = -Z[IWAVI];

    //
    // TIMING INFO
    //
L1500:
    // TOTTIM = second() - TSTART
    TOTTIM = second() - (double)TSTART;
    // WAVTIM = TIMES(1) + TIMES(2) + TIMES(3)
    WAVTIM = (double)TIMES[1] + (double)TIMES[2] + (double)TIMES[3];
    // RTIM1 = TIMES(7)-TIMES(5)
    RTIM1 = (double)TIMES[7] - (double)TIMES[5];
    // RTIM2 = TIMES(8)-TIMES(6)
    RTIM2 = (double)TIMES[8] - (double)TIMES[6];
    // OTHTIM = TOTTIM - WAVTIM
    OTHTIM = TOTTIM - WAVTIM;
    // TOTTIM = TOTTIM + TIMES(7) + TIMES(8)
    TOTTIM = TOTTIM + (double)TIMES[7] + (double)TIMES[8];

    // write (6, 1503) ( TIMES(I), I = 1, 3 ), WAVTIM, OTHTIM, TIMES(5),
    //   RTIM1, TIMES(7), TIMES(6), RTIM2,
    //   TIMES(8), TOTTIM
    // 1503 FORMAT ( / '0RADIAL INTEGRAL COMPUTATION TIMES (SECONDS):' /
    //   '0SCATTERING WAVE INITIALIZATIONS:', T35, F9.3 /
    //   ' SCATTERING WAVE LOOP:', T35, F9.3 /
    //   ' SCATTERING WAVE INTERPOLATIONS:', T35, F9.3 /
    //   ' TOTAL SCATTERING WAVE TIME', T45, F9.3 /
    //   ' INTEGRATIONS AND CLEBSCH''S:', T45, F9.3 /
    //   ' BELLINGS (SUMMAX, INFINITY):', T35, F9.3 /
    //   ' RECURSION (SUMMAX, INFINITY):', T35, F9.3 /
    //   ' TOTAL COULOMB (SUMMAX, INFINITY):', T45, F9.3 /
    //   ' BELLINGS (0, INFINITY):', T35, F9.3 /
    //   ' RECURSION (0, INFINITY):', T35, F9.3 /
    //   ' TOTAL PURE COULOMB (0, INFINITY):', T45, F9.3 /
    //   ' TOTAL:', T55, F9.3, ' SEC.' / )
    std::printf("\n0RADIAL INTEGRAL COMPUTATION TIMES (SECONDS):\n");
    std::printf("0SCATTERING WAVE INITIALIZATIONS:%26s%9.3f\n", "", (double)TIMES[1]);
    std::printf(" SCATTERING WAVE LOOP:%33s%9.3f\n",            "", (double)TIMES[2]);
    std::printf(" SCATTERING WAVE INTERPOLATIONS:%28s%9.3f\n",  "", (double)TIMES[3]);
    std::printf(" TOTAL SCATTERING WAVE TIME%38s%9.3f\n",       "", WAVTIM);
    std::printf(" INTEGRATIONS AND CLEBSCH'S:%38s%9.3f\n",      "", OTHTIM);
    std::printf(" BELLINGS (SUMMAX, INFINITY):%26s%9.3f\n",     "", (double)TIMES[5]);
    std::printf(" RECURSION (SUMMAX, INFINITY):%25s%9.3f\n",    "", RTIM1);
    std::printf(" TOTAL COULOMB (SUMMAX, INFINITY):%21s%9.3f\n","", (double)TIMES[7]);
    std::printf(" BELLINGS (0, INFINITY):%31s%9.3f\n",          "", (double)TIMES[6]);
    std::printf(" RECURSION (0, INFINITY):%30s%9.3f\n",         "", RTIM2);
    std::printf(" TOTAL PURE COULOMB (0, INFINITY):%21s%9.3f\n","", (double)TIMES[8]);
    std::printf(" TOTAL:%49s%9.3f SEC.\n\n",                    "", TOTTIM);

    //
    // A L L    D O N E
    //
    IRTN = 1;

    return;
}
