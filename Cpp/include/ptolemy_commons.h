#ifndef PTOLEMY_COMMONS_H
#define PTOLEMY_COMMONS_H

#include "ptolemy_types.h"

//
// All Fortran COMMON blocks translated as C++ global structs.
// Memory layout must match Fortran exactly for interoperability.
//
// Parameters defining dimensions
//
constexpr int NUMSWITCHES  = 34;
constexpr int NUMFLOAT     = 153;
constexpr int NUMINTEGER   = 53;
constexpr int NUMJWORD     = 12;
constexpr int NUMLINKULES  = 13;
constexpr int NAMDIM       = 500;
constexpr int COM_MEM_SIZE = 6250000;

// ============================================================================
// COMMON /SWITCH/  — 34 integer control switches
// ============================================================================
struct SwitchCommon {
    int IASYMP, ICHECK, IECHO, IELAST, IFIT;
    int IVRTEX, IRLWAV, ITSO, LINEAR, NUCONL;
    int NEXT, IFITER, IEPOW, ISAVMX, ISAVB;
    int LABANG, PROBLM, MASTYP, R0TYPE, IBSPAS;
    int ISW1, ISW2, ISW3, ISCRFL, NOBFAC;
    int CONVRG, IEXTYP, IBRNSB, MEORD, KFRMTP;
    int IWRTWV, LORNTZ, IRELPT, ISAVSM;
};
extern SwitchCommon SWITCH;

// Array alias for SWITCH (EQUIVALENCE ( SWITCH(1), IASYMP ))
inline int* SWITCH_arr() { return &SWITCH.IASYMP; }

// ============================================================================
// COMMON /INTRNL/  — internal status
// ============================================================================
struct IntrnlCommon {
    double UNDEF;
    double NOTDEF;
    int    ICHANB, ICHANW, ICNTT;
    double EBNDS[3];     // 1-based: EBNDS[1], EBNDS[2]
    int    IDONE, ISTRIP, IHSAVE, NEXTSV, IPOT, IWFN, LIMOST, IEXCIT;
    double R0MASS;
    logical1 WASSET[41]; // 1-based: WASSET[1]..WASSET[40]
    double RATMAS;
    int    LSPECS[5];    // 1-based
    int    NODESP[5];    // 1-based
    int    LSPCPT[3];    // 1-based
    int    NODEPT[3];    // 1-based
    int    ICD3I;
};
extern IntrnlCommon INTRNL;

// EQUIVALENCE (undef4, notdef) — real*4 overlay on notdef
inline float& undef4() { return *reinterpret_cast<float*>(&INTRNL.NOTDEF); }

// Integer NOTDEF value: the Fortran stores 0xF0F0F0F0 as an int bit pattern in NOTDEF
// (int)INTRNL.NOTDEF gives WRONG result (casts tiny double to 0)
// Use this constant for integer comparisons instead
constexpr int NOTDEF_INT = static_cast<int>(0xF0F0F0F0u); // = -252645136

// ============================================================================
// COMMON /FLOAT/  — 153 floating-point parameters
// ============================================================================
struct FloatCommon {
    double A, ACCURA, AI, ANGSTP, ANGMIN;
    double ANGMAX, ASO, ASOI, ASYMPT, DELTVK;
    double DWCUT, E, ECM, ELAB;
    double EXS[6];   // 1-based: EXS[1]..EXS[5], EXS[0] is padding
    double GAMDIF;
    double GAMSUM, TAU, TAUI, AM, AMA, AMB, AMBIGA, AMBIGB, AMX;
    double AMXCS[6];  // 1-based
    double AMXCGS[6]; // 1-based
    double AMP;
    double AMT, Q, R, R0, RC;
    double RC0, RI, RI0, RSO, RSO0;
    double RSOI, RSOI0, SPAM, SPAMP, SPAMT;
    double STEPSZ, STEP1R, STEP1I, SUMMAX, SUMMID;
    double SUMMIN, V, VI, VSO, VSOI;
    double FITMUL, WRITES, ASI, RSI, RSI0;
    double VSI, AMDMLT, STPSPR, R0E, RI0E;
    double AE, AIE, VE, VIE, FITACC;
    double ALMNMT, SPFACP, SPFACT, R0ESQ, RI0ESQ;
    double AESQ, AIESQ, VESQ, VIESQ, DERIVS;
    double FITRAT, RCP, RCT, RC0P, RC0T;
    double SUMPTS, ACCINE, COULML, ALMXMT;
    double EXSPT[3];  // 1-based
    double POWRL, POWIM, BNDASY, SCTASY;
    double VTR, VTRI, VTL, VTLI, VTP, VTPI;
    double PARAM1, PARAM2, PARAM3, PARAM4, PARAM5;
    double PAR620[16]; // 1-based
    double RTR, RTRI, RTL, RTLI, RTP, RTPI;
    double RTR0, RTRI0, RTL0, RTLI0, RTP0, RTPI0;
    double ATR, ATRI, ATL, ATLI, ATP, ATPI;
    double AMXGPT[3]; // 1-based
    double PHIMID, GAMPHI;
};
extern FloatCommon FLOAT_common;

// Fortran-compatible 1-based FLOAT array access (accounts for struct padding)
// Struct has 1-based padding arrays (EXS[6] etc.) that shift field offsets.
// This accessor maps Fortran position I to the correct C++ struct offset.
inline double& FLOAT_arr(int I) {
    int pad = 0;
    if (I >= 15) pad++;   // EXS[0] padding
    if (I >= 30) pad++;   // AMXCS[0] padding
    if (I >= 35) pad++;   // AMXCGS[0] padding
    if (I >= 100) pad++;  // EXSPT[0] padding
    if (I >= 117) pad++;  // PAR620[0] padding
    if (I >= 150) pad++;  // AMXGPT[0] padding
    return (&FLOAT_common.A)[I - 1 + pad];
}
inline double* FLOAT_arr_raw() { return &FLOAT_common.A; }
inline double* FLOAT_arr() { return &FLOAT_common.A; }  // raw pointer, caller must handle padding

// EQUIVALENCE aliases for FLOAT (use FLOAT_arr for correct padding offset)
inline double& AMS(int I) { return FLOAT_arr(24 + I); }     // AMS(1)=AMA=pos25
inline double* AMS()     { return &FLOAT_arr(25) - 1; }     // AMS()[I] 1-based
inline double* PARAMS()  { return &FLOAT_common.PARAM1; }   // PARAMS(1)..PARAMS(20) — TODO: needs padding fix
inline double* VTEN()    { return &FLOAT_common.VTR; }      // TODO: needs padding fix
inline double* RTEN()    { return &FLOAT_common.RTR; }
inline double* RTEN0()   { return &FLOAT_common.RTR0; }
inline double* ATEN()    { return &FLOAT_common.ATR; }

// ============================================================================
// COMMON /INTGER/  — 53 integer parameters
// ============================================================================
struct IntgerCommon {
    int L, LMAX, LMIN, LOOKST, LX;
    int MAPSUM, MAPDIF, MAXITR, NAITKN, NBACK;
    int NCOSIN, NODES, NPDIF, NPHIAD, NPPHI;
    int NPSUM, IPRINT, ISIZE;
    int IZS[6];   // 1-based
    int IZP, IZT;
    int LSTEP, LBACK, MAXLEX, LMAXAD, MAXFUN;
    int IREINI, MODEFT, NUMRAN, NVPOLY, NPCOUL;
    int MXCOUL, LMINSB;
    int PARITS[6]; // 1-based
    int PARIPT[3]; // 1-based
    int PARITY;
    int IPARM1, IPARM2, IPARM3, IPARM4, IPARM5;
    int NBLKSZ, MAPPHI, NFIROF, LEBACK;
};
extern IntgerCommon INTGER;

// Raw array view (has padding gaps — DO NOT use for Fortran-indexed access)
inline int* INTGER_arr() { return &INTGER.L; }

// Fortran-compatible 1-based INTGER array access (accounts for struct padding)
// Padding: IZS[0] at pos 19, PARITS[0] at pos 38, PARIPT[0] at pos 43
inline int& INTGER_arr_f(int I) {
    int pad = 0;
    if (I >= 19) pad++;   // IZS[0]
    if (I >= 38) pad++;   // PARITS[0]
    if (I >= 43) pad++;   // PARIPT[0]
    return (&INTGER.L)[I - 1 + pad];
}

// ============================================================================
// COMMON /JBLOCK/  — 12 words of angular momentum
// ============================================================================
struct JblockCommon {
    double J;
    double JS[6];  // 1-based
    double JP, JSP, JST;
    double JPARM1, JPARM2, JPARM3;
};
extern JblockCommon JBLOCK;

// Raw array view (has padding — DO NOT use for Fortran-indexed access)
inline double* JBLOCK_arr() { return &JBLOCK.J; }

// Fortran-compatible 1-based JBLOCK array access
// Padding: JS[0] at pos 2
inline double& JBLOCK_arr_f(int I) {
    int pad = 0;
    if (I >= 2) pad++;   // JS[0]
    return (&JBLOCK.J)[I - 1 + pad];
}

// ============================================================================
// COMMON /KANDM/  — kinematic and matching parameters
// ============================================================================
struct KandmCommon {
    double AKI, AKO, REDMI, REDMO, RIMAX, ROMAX;
    int    ISIGS[3];    // 1-based
    double RSCTS[3];    // 1-based
    double ASCTS[3];    // 1-based
    double RCSCTS[3];   // 1-based
    int    ISMATS[3];   // 1-based
    int    LOMOST, LCRIT;
    int    ISCTMN, ISCTCR;
    double ETAS[3];     // 1-based
    double RCSCTP[3];   // 1-based
    double RCSCTT[3];   // 1-based
    int    LCRITS[3];   // 1-based
    double RATMSS[3];   // 1-based
    int    QMIN;
    double ABAR;
};
extern KandmCommon KANDM;

// EQUIVALENCE aliases — 1-based: AKS()[1]=AKI, AKS()[2]=AKO
inline double* AKS()    { return &KANDM.AKI - 1; }    // AKS(1)..AKS(2)
inline double* REDMS()  { return &KANDM.REDMI - 1; }  // REDMS(1)..REDMS(2)
inline double* RMAXS()  { return &KANDM.RIMAX - 1; }  // RMAXS(1)..RMAXS(2)

// ============================================================================
// COMMON /CNSTNT/  — 13 physical constants
// ============================================================================
struct CnstntCommon {
    double PI, RT4PI, PIINV, RADIAN, DEGREE;
    double HBARC, AMUMEV, AFINE;
    double BIGEST, SMLEST, BIGNUM, SMLNUM, BIGLOG;
};
extern CnstntCommon CNSTNT;

// ============================================================================
// COMMON /WAVCOM/  — wave function and scattering parameters
// ============================================================================
struct WavcomCommon {
    double HS[3];       // 1-based
    double RSTEPS[3];   // 1-based
    int    NSTPSS[3], NFS[3], NGS[3], NF1S[3], NG1S[3]; // 1-based
    double STEPR, STEPI;
    int    IWAVR, IWAVI, NUMAIT;
    int    NUMFIT, IOFFIT;
    int    PWAVSW;      // LOGICAL
    int    NBAKCM;
    int    STANSW;      // LOGICAL
    int    SOSWS[3];    // 1-based, LOGICAL
    int    NUMJS[3];    // 1-based
    double RIS[3], AIS[3];   // 1-based
    double RSOS[3], ASOS[3], RSOIS[3], ASOIS[3]; // 1-based
    int    IRLVS[3], IIMVS[3];  // 1-based
    int    IZPS[3], IZTS[3];    // 1-based
    double ASYMPS[3];   // 1-based
    double V0RS[3], V0IS[3], V0SORS[3], V0SOIS[3]; // 1-based
    double ES[3];       // 1-based
    int    JSPS[3];     // 1-based
    double RSIS[3], ASIS[3], V0SIS[3]; // 1-based
    int    ISORS[3], ISOIS[3]; // 1-based
    double ALPHAR, ALPHAI;
    int    NFIRST;
    int    PWBGSW;      // LOGICAL
    int    LASTL[3], LASTNF[3]; // 1-based
    int    LSKIPS[3], ISTATS[3]; // 1-based
    double POWRLS[3], POWIMS[3]; // 1-based
    int    NSTP2S[3];   // 1-based
    double XFACS[4][3]; // 1-based: XFACS(3,2) -> [4][3]
    int    IVWORK;
    int    ICENTR[3];   // 1-based
    int    TCSWS[3];    // 1-based, LOGICAL
    int    NWAVEF;
    double VTENS[3][7]; // 1-based: VTENS(2,6) -> [3][7]
    double RTENS[3][7]; // 1-based
    double ATENS[3][7]; // 1-based
    int    ITENS[3][7]; // 1-based
    int    IWAVRS[3], IWAVIS[3]; // 1-based
    int    LASTZR;
    int    SOSWT;       // LOGICAL
    int    IINDXE[3], ITOCE[3]; // 1-based
};
extern WavcomCommon WAVCOM;

// ============================================================================
// COMMON /FORMF/  — form factor block
// ============================================================================
struct FormfCommon {
    int    IBDS[3];     // 1-based
    int    NPOTS[3];    // 1-based
    int    LBP, LBT, NBP, NBT, JBP, JBT;
    double BNDMAX[3];   // 1-based
    double BNDSTP[3];   // 1-based
    int    NSTPBD[3];   // 1-based
    double RBNDS[3];    // 1-based
    double ABNDS[3];    // 1-based
    double RCBNDS[3];   // 1-based
    double BNDMAS[3];   // 1-based
    double BRATMS[3];   // 1-based
    double RSOBDS[3];   // 1-based
    double ASOBDS[3];   // 1-based
    double VBNDS[3];    // 1-based
    double VSOBDS[3];   // 1-based
    double TAUBDS[3];   // 1-based
    int    LNKADB[13][3][7]; // 1-based: LNKADB(6,2,12) -> [13][3][7]
};
extern FormfCommon FORMF;

// EQUIVALENCE: R2S(4) <=> RBNDS(1)
inline double* R2S() { return &FORMF.RBNDS[1]; }

// ============================================================================
// COMMON /GRIDCM/  — grid and integration control
// ============================================================================
struct GridcmCommon {
    double JACOB, STPINV, COSSTP;
    int    IMSVAL, IINTS, INTOFF, IA12VL, ICOSIN;
    int    IPHIT, IPHIP, IPHI, LOGIC, LOCTRP;
    int    IRI, IRO, ILOR, ILOI, ILIR, ILII, IWIO, NOFLO, NRIROI, NCRIT;
    int    IXLAM, INLAM, IA12N, MSKCOS, MAXCNT, ISMHPT, ISMIPT;
    int    ICOSS4, IJA12S, IHS, IIINDX, IA12TM, IDW, IHINT, IHABS;
    int    IABS1, IGRDC1, NPSUMI;
    double RVRLIM;
    int    ISMHWK, ISMHVL, ISMIVL;
    int    IRIOEX, NUMHS;
    int    IGRDC2[5];   // 1-based
    int    IWFII, IWFIO, NWFI, NWFO, IDWI;
};
extern GridcmCommon GRIDCM;

// ============================================================================
// COMMON /INELCM/  — inelastic scattering control
// ============================================================================
struct InelcmCommon {
    int    JBPF, JBTF, LXMIN, LXMAX, ILILOR, ILILOI, MSTOP;
    int    IA12M, NMLOLX;
    int    ONELSW;      // LOGICAL
    int    NUMLX;
    int    IRDINT, IELCUP, LILOSZ, NUMLIS;
    int    ILIS;
    int    IBETAS[4];   // 1-based
    int    ICL2FF, LXSTEP;
    int    IDENSW;      // LOGICAL
    int    LSKIP;
    int    JPMIN, JPMAX, JPBASE, NJP;
    int    JTMIN, JTMAX, JTBASE, NJT;
    int    IINDXS, ISMATR, ISMATI, NSMAT, NLX, NSPLI, ITOCS;
    int    NASPLI, IUNITR, LIFIT, ILIFIT, NOLFIT;
    int    IDELSR, IDELSI;
};
extern InelcmCommon INELCM;

// EQUIVALENCE: INDXJS(4,2) <=> JPMIN
inline int* INDXJS() { return &INELCM.JPMIN; }

// ============================================================================
// COMMON /CCBLK/  — coupled channels block
// ============================================================================
struct CcblkCommon {
    int ICHNDF, ICHNCP, ITRANS, IBASDF, IBASCP;
    int ICCBLK[3];      // 1-based
    int IORDX, MCHNVL, ICHNVL;
    int IHOMOA, IHOMOB;
    int IINHR8, IINHI8;
};
extern CcblkCommon CCBLK;

// ============================================================================
// COMMON /HEDCOM/  — header/reaction information
// ============================================================================
struct HedcomCommon {
    char REACT[46];     // 1-based: REACT(45) -> [46]
    char HEADER[66];    // 1-based: HEADER(65) -> [66]
};
extern HedcomCommon HEDCOM;

// ============================================================================
// COMMON /FTIME/  — timing
// REAL*4 TBEGIN, TIMES
// ============================================================================
struct FtimeCommon {
    float  TBEGIN;
    int    NUMLOP, NUMWAV;
    float  TIMES[9];    // 1-based: TIMES(6) or TIMES(8)
    int    NUMCAL, NUMGRA;
};
extern FtimeCommon FTIME;

// ============================================================================
// COMMON /FITCOM/  — fit parameters
// ============================================================================
struct FitcomCommon {
    float  T0;
    int    IP, IPL, IEPARM, IFS, MAXCAL;
    int    NROC, ISCSIG;
};
extern FitcomCommon FITCOM;

// ============================================================================
// COMMON /LOCFIT/  — fit parameter locations (union of different views)
// ============================================================================
struct LocfitCommon {
    union {
        struct {
            int IDATA, LPL, LSIGMA, LSMAT, IESETS, IMASK;
            int IEQUAL, NPARA, IPARA, ITEMP, LSCSIG, LEPARM, LFS, NPTS, LEQNAM;
            int MAXEQL, MAXVS, LVDERV, LSDELS, LDVTYP, LDSDP, LDFDX, LFGRAD;
            int LTGRAD, LMASK, LEQUAL, MNDRIV, MXDRIV, IESIZE, LENPRE;
            int LENKEY, LENHED, ILNKBL;
        } full;
        struct {
            int LFITSP[27]; // 1-based: LFITSP(26) -> [27]
            int MNDRIV, MXDRIV;
        } compact;
        struct {
            int LNKAD2[NUMLINKULES + 1][3][7]; // 1-based: LNKAD2(6,2,NUMLINKULES)
        } linkule;
    };
};
extern LocfitCommon LOCFIT;

// ============================================================================
// COMMON /INPBUF/  — input buffer
// ============================================================================
constexpr int IBUFSIZE = 200;

struct InpbufCommon {
    int  NOCH, INCH;
    char IBUF[IBUFSIZE + 1]; // 1-based: IBUF[1]..IBUF[200]
};
extern InpbufCommon INPBUF;

// ============================================================================
// COMMON /TEMPVS/  — temporary potential values
// ============================================================================
struct TempvsCommon {
    double TVR, TVI, TVSOR, TVSOI, TAR, TAI;
};
extern TempvsCommon TEMPVS;

// Array view of TEMPVS
inline double* TEMPVS_arr() { return &TEMPVS.TVR; }

// ============================================================================
// COMMON /LOCptrs/  — memory location pointers
// ============================================================================
struct LocptrsCommon {
    integer4 Z[NAMDIM + 1]; // 1-based: Z(500) -> [501]
};
extern LocptrsCommon LOCPTRS;

// ============================================================================
// COMMON /LENGTH/  — object lengths
// ============================================================================
struct LengthCommon {
    int LENG[NAMDIM + 1]; // 1-based
};
extern LengthCommon LENGTH;

// ============================================================================
// COMMON /NAMCOM/  — named object registry
// ============================================================================
struct NamcomCommon {
    char8 NAMES[NAMDIM + 1]; // 1-based
};
extern NamcomCommon NAMCOM;

// ============================================================================
// COMMON /USAGE/  — allocator statistics and control
// ============================================================================
struct UsageCommon {
    int    IBASE, NWORDS, IMAX, INAME, IGOT, NLOC;
    int    IPEAK, NAMPEK, NUMCOM, ICHEC;
    char8  THENAM;
    int    IWASAT, NAMLC2;
    int    NSPTSW;      // LOGICAL
};
extern UsageCommon USAGE;

// ============================================================================
// COMMON /ALLOCS/  — the main memory pool
// ============================================================================
struct AllocsCommon {
    int    FACFR4, FACFR2, FACFR1;
    int    intjunk;     // padding to align ALLOC on 8-byte boundary
    double ALLOC_pool[COM_MEM_SIZE + 1]; // 1-based
};
extern AllocsCommon ALLOCS;

// ============================================================================
// COMMON /LNKBLK/  — linkule address block
// ============================================================================
struct LnkblkCommon {
    int NUMLNK, IUNIQU;
    int LNKADR[NUMLINKULES + 1][7]; // 1-based: LNKADR(6,NUMLINKULES)
};
extern LnkblkCommon LNKBLK;

// EQUIVALENCE: LNKNAM(3,NUMLINKULES) <=> LNKADR(1,1) — character*8 name overlay
inline char8* LNKNAM() { return reinterpret_cast<char8*>(&LNKBLK.LNKADR[1][1]); }

// ============================================================================
// COMMON /SCRTCH/  — scratch work area
// ============================================================================
struct ScrtchCommon {
    int LNKAD3[NUMLINKULES + 1][7]; // 1-based: LNKAD3(6,NUMLINKULES)
};
extern ScrtchCommon SCRTCH;

// ============================================================================
// COMMON /FACTRL/  — factorial table
// ============================================================================
struct FactrlCommon {
    int    MAXFAC, ISPACE;
    double FACTBL[98];  // 1-based: FACTBL[1]..FACTBL[97] = sqrt(0!)..sqrt(96!)
};
extern FactrlCommon FACTRL;

// ============================================================================
// COMMON /LOGFAC/  — log factorial table
// ============================================================================
struct LogfacCommon {
    int    MAXLF, LFBIAS;
    double LF[1002];    // 1-based: LF[1]..LF[1001]
};
extern LogfacCommon LOGFAC;

// ============================================================================
// COMMON /CKBLOK/  — Coulomb parameter block
// ============================================================================
struct CkblokCommon {
    double CPARAM[5];   // 1-based
    double R1, FACTOR;
};
extern CkblokCommon CKBLOK;

// ============================================================================
// COMMON /VCSBLK/  — Coulomb correction factors
// ============================================================================
struct VcsblkCommon {
    double VCSQ12[12][4]; // 1-based: VCSQ12(3,11)
};
extern VcsblkCommon VCSBLK;

// ============================================================================
// COMMON /BSBLK/  — bound state parameters
// ============================================================================
struct BsblkCommon {
    double S1, S2, T1, T2, BNDMXP, BNDMXT, RLPMAX, RLTMAX;
    double VPMAX, VTMAX, SCTMAX, SCTSP;
    double BNDSPT, BNDSPP, BNDSPV;
    int    NSPBDT, NSPBDP, NSPBDV;
    int    JBDT, JBDP, JPOT;
    int    IOUTSW;       // LOGICAL
    int    NSPSCT;
    double ALPHAP, ALPHAT;
    int    PINFSW;       // LOGICAL
    float  PHISGN;
    double ROFMAX;
};
extern BsblkCommon BSBLK;

// EQUIVALENCE: RLMAXS(2) <=> RLPMAX, VMAXS(2) <=> VPMAX
inline double* RLMAXS() { return &BSBLK.RLPMAX; }  // RLMAXS(1)=RLPMAX, RLMAXS(2)=RLTMAX
inline double* VMAXS()  { return &BSBLK.VPMAX; }    // VMAXS(1)=VPMAX, VMAXS(2)=VTMAX

// ============================================================================
// COMMON /GRABCOM/  — GRAB state
// ============================================================================
struct GrabcomCommon {
    int LASTHLM;
};
extern GrabcomCommon GRABCOM;

// ============================================================================
// COMMON /ROOTIS/  — roots of orthogonal polynomials
// ============================================================================
struct RootisCommon {
    int    IMXDIM, ISPACE;
    double ROOTI[603];  // 1-based: ROOTI(602)
};
extern RootisCommon ROOTIS;

// ============================================================================
// COMMON /ROOTIX/  — roots index/control
// ============================================================================
struct RootixCommon {
    double X;
};
extern RootixCommon ROOTIX;

// ============================================================================
// COMMON /F152/  — numerical workspace (in linkulesfitters)
// ============================================================================
struct F152Common {
    double P[26];       // 1-based: P(25)
    double Q_arr[26];   // 1-based: Q(25)
    int    PIVOT[26];   // 1-based: PIVOT(25)
};
extern F152Common F152;

// ============================================================================
// COMMON /F402/  — numerical workspace (in fortlib)
// ============================================================================
struct F402Common {
    double PIVOT[51];   // 1-based: PIVOT(50)
    int    INDEX[51];   // 1-based: INDEX(50)
};
extern F402Common F402;

// ============================================================================
// COMMON /INTWRK/  — integer workspace
// ============================================================================
struct IntwrkCommon {
    int MAXINT, INBIAS, WKBIAS;
    int INTTBL[1];      // dynamic size
};
extern IntwrkCommon INTWRK;

// ============================================================================
// COMMON /VMM/, /VMM2/  — variable metric minimization (in linkulesfitters)
// ============================================================================
struct VmmCommon {
    double F, GS, EL, SL, FP, GSP, T0, Z, Q_val, A_val;
    double GSS, F0, GTP, FB, GTT, GSB, F1;
    int    IY;
};
extern VmmCommon VMM;

struct Vmm2Common {
    double dummy[1]; // placeholder - populated during linkulesfitters translation
};
extern Vmm2Common VMM2;

#endif // PTOLEMY_COMMONS_H
