#ifndef PTOLEMY_FORWARD_H
#define PTOLEMY_FORWARD_H

//
// Forward declarations of all Ptolemy subroutines/functions.
// This file is populated incrementally as files are translated.
//

// --- gfortran_stuff.f ---
void CMPUTR(char* mach, char* hostname);
double DGAMMA(double X);
double DGAMAP(double X);
double DLGAMA(double X);
double DLNGAM(double X);
double DIGAMA(double X);
void GRAB(double* S, int I4OR8, int& IBASE, int& NWORDS);

// --- dtime.cpp ---
double dtime_();

// --- srread.cpp ---
void srread(int* idcb, int* icode, char* line, int* numchr, int* maxchr, int* rtype);
void swrite(int* idcb, int* icode, char* line, int* numchr);
void snwrte(int* idcb, int* icode, char* line, int* numchr);
void sropen(int* idcb, int* icode, char* fname);
void srclse(int* idcb);
void swopen(int* idcb, int* icode, char* fname, int* oversw, int* append);
void swclse(int* idcb);
void swdlte(int* idcb, int* icode, char* fname);

// --- rcwfn.f ---
void RCWFN(double RHO, double ETA, int MINL, int MAXL,
           double* FC, double* FCP, double* GC, double* GCP,
           double ACCUR, int& IRET);

// --- fortlib.f ---
// Angular momentum coefficients
double CLEBSH(int A, int B, int X, int Y, int CIN, int ZIN);
double THRJ(int A, int B, int CIN, int X, int Y, int ZIN);
double THREEJ(int A, int B, int CIN, int X, int Y, int ZIN);
double SIXJ(int A, int B, int C, int X, int Y, int Z);
double RACAH(int A_r, int B_r, int Y_r, int X_r, int C_r, int Z_r);
double WIG9J(int J1, int J2, int J3, int J4, int J5, int J6,
             int J7, int J8, int J9);

// Complex continued fraction
void CCNFRC(int NUMPTS, complex16* XS, complex16* YS);
void CCONTF(int NUMPTS, complex16* XS, complex16* YS, complex16 X, complex16& Y);

// Coulomb phases
void DSGMAL(double ETA, int NO, double* DSG);

// Legendre and spherical harmonics
void PLMSUB(int LMAX, int MMAXin, double X, double* PRAY);
void PLSUB(int LMAX, double X, double* PRAY);
void YLMSUB(int LMAX, int MMAX, double Z, double* YLMRAY);

// Interpolation and splines
void INTRPC(int NCUBIC, double* XCUBES, double* AS, double* BS, double* CS,
            double* DS, int NPTS, double* XS, double* YS);
void SPLNCB(int N, double* X, double* Y, double* B, double* C, double* D);

// Linear algebra
void LIN(double* A, int NMAX, int NDIM, double* B, int MMAX,
         double POWER, double& DET);
void LINBSV(double* A, int NMAX, int NDIM, double* BSV, int MMAX,
            double* B, double POWER, double& DET);
void LINSAV(double* ASV, int NMAX, int NDIM, double* BSV, int MMAX,
            double* A, double* B, double POWER, double& DET);
void LINTRI(double* A, int NMAX, int NDIM, double* BSV, int MMAX,
            double* B, double POWER, double& DET);
void MATINV(double* A, int N, double* B, int M, double& DETERM, int NMAX);
void LSQPOL(double* X, double* Y, double* W, double* RESID, int NSUB,
            double* SUM, int LSUB, double* A, double* B, int MSUB,
            int NMAX, int MMAX);

// Quadrature
void LAGBC(int NN, double ALF, double* A, double* B, double* C);
void LGRECR(double& PN, double& DPN, double& PN1, double X,
            int NN, double ALF, double* A, double* B, double* C);
void LGROOT(double& X, int NN, double ALF, double& DPN, double& PN1,
            double* A, double* B, double* C, double EPS);
void LAGUER(int NN, double* X, double* W, double ALF, double EPS,
            double& CSX, double& CSW, double& TSX, double& TSW,
            double* A, double* B, double* C);
void GAUSSL(int N, double* X, double* W);

// Date and error
void get_date(char* date);
void SYSERR();

// --- masstable.f ---
double EXCESS(int KZ, int KA, int& NOTTAB);
void MASEXX(int IZ, int IA, double& EX, double& ER, int& IERR);
void GETAS(int IZ, int& NUMAS, int* IAS);
void AZCODE(char* SIN, int& IZ, int& IA, char* SYMBOL, int& IRTN);
void CODEZ(int IZ, char* ZSYM);

// --- numbered_store.f ---
int IALLOC(int IWRDS);
int NALLOC(const char* NAME, int IWRDS);
int NAMLOC(const char* NAME);
void IREDEF(int N1, int N2);
int IREDEF_F(int IWRDS, int NAMLCA);
void NSDUMP(int IDUMP, int NSKEY);
void ISTART(int N1);
void ICLEAR(int N1);
void NSSTAT(int N1);
void NSCAT(int N1);

// --- keep.f ---
void KEEPIT();

// --- keepsub.f ---
void KEEP(char* OBJ, int KIND, int IFAM, int NDIM, int* IDIM,
          char* LIB, char* NAME, int& IER);
void KEEPCK(char* OBJ, int KIND, int IFAM, int NDIM, int* IDIM,
            char* LIB, char* NAME, int& IER);

// --- keptsub.f ---
void KEPT(char* OBJ, int& KIND, int& IFAM, int& NDIMS, int* IDIM,
          char* LIB, char* NAME, int& IER1, int MBYTE1);
void KEPTCN(char* OBJ, int& IER2, int MBYTE2);
void KEPTCL();

// --- linkule.f ---
void GETLNK(char8 AKEY, int& IRET);
void LINKUL(int LOC, char8 ALIAS, int* MYINTS, int IPOTYP, int IREQUE,
            int& IRETUR, int L, double J, double RSTART, double STEPSZ,
            int NUMPTS, double* ARRAY1, double* ARRAY2, char* ID);

// --- av18.f ---
void av18(char8 alias, real4* savint, int ipotyp, int ireque, int& iretur,
          int iunit, int& numout, int l, double& jp, double rstart,
          double stepsz, int numpts, double* array1, double* array2,
          double* flt, double* tempvs, double* param, int* intger,
          double* jblock, int* iswtch, double* intrnl, double* fintrn,
          double* cnstnt, double* wavbl, double* fkandm,
          char* id, double* alloc, int* illoc, int facfr4,
          integer4* loc, int* length,
          int (*nalloc)(const char*, int), int (*namloc)(const char*),
          int* iparam);

// --- phiffer.f ---
void phiffer(char8 alias, real4* savint, int ipotyp, int ireque, int& iretur,
             int iunit, int& numout, int l, double& jp, double rstart,
             double stepsz, int numpts, double* array1, double* array2,
             double* flt, double* tempvs, double* param, int* intger,
             double* jblock, int* iswtch, double* intrnl, double* fintrn,
             double* cnstnt, double* wavbl, double* fkandm,
             char* id, double* alloc, int* illoc, int facfr4,
             integer4* loc, int* length,
             int (*nalloc)(const char*, int), int (*namloc)(const char*),
             int* iparam);
void PHIFER(int itype, int ievopt, int nodes, int l, int jtwo,
            double h2o2mu, double* params, bool prntsw, int npts,
            double rmax, double* rgrid, double* phi_by_rl,
            int nvgrid, double* vongrid, double& vrgrid);
void STOP_F(const char* msg, int icode);

// --- linkulesfitters.f — potential linkules ---
void AITKEN(int NDEGRE, double DELTA, double EPSLON, int NIN,
            double* XIN, double* FIN,
            int NOUT, double* XOUT, double* FOUT,
            int& NFAIL, double& WORST, double* WORK);
void BKGPTElp(char8 ALIAS, int* MYINTS, int IPOTYP, int IREQUE, int& IRETUR,
              int IUNIT, int& NUMOUT, int L, double J, double RSTART,
              double STEPSZ, int NUMPTS, double* ARRAY1, double* ARRAY2,
              char* ID, int* IPARAM);
void DEFORMEd(char8 ALIAS, int* MYINTS, int IPOTYP, int IREQUE, int& IRETUR,
              int IUNIT, int& NUMOUT, int L, double J, double RSTART,
              double STEPSZ, int NUMPTS, double* ARRAY1, double* ARRAY2,
              char* ID, int* IPARAM);
void FIXEDWOo(char8 ALIAS, int* MYINTS, int IPOTYP, int IREQUE, int& IRETUR,
              int IUNIT, int& NUMOUT, int L, double J, double RSTART,
              double STEPSZ, int NUMPTS, double* ARRAY1, double* ARRAY2,
              char* ID);
void GAUSSIAn(char8 ALIAS, int* MYINTS, int IPOTYP, int IREQUE, int& IRETUR,
              int IUNIT, int& NUMOUT, int L, double J, double RSTART,
              double STEPSZ, int NUMPTS, double* ARRAY1, double* ARRAY2,
              char* ID);
void JDEPEN(char8 ALIAS, int* MYINTS, int IPOTYP, int IREQUE, int& IRETUR,
            int IUNIT, int& NUMOUT, int L, double J, double RSTART,
            double STEPSZ, int NUMPTS, double* ARRAY1, double* ARRAY2,
            char* ID);
void JDEPENWS(char8 ALIAS, int* MYINTS, int IPOTYP, int IREQUE, int& IRETUR,
              int IUNIT, int& NUMOUT, int L, double J, double RSTART,
              double STEPSZ, int NUMPTS, double* ARRAY1, double* ARRAY2,
              char* ID);
void LAGRANGE(char8 ALIAS, int* MYINTS, int IPOTYP, int IREQUE, int& IRETUR,
              int IUNIT, int& NUMOUT, int L, double J, double RSTART,
              double STEPSZ, int NUMPTS, double* ARRAY1, double* ARRAY2,
              char* ID);
void LTSTELp(char8 ALIAS, int* MYINTS, int IPOTYP, int IREQUE, int& IRETUR,
             int IUNIT, int& NUMOUT, int L, double J, double RSTART,
             double STEPSZ, int NUMPTS, double* ARRAY1, double* ARRAY2,
             char* ID, int* IPARAM);
void PARITWOO(char8 ALIAS, int* MYINTS, int IPOTYP, int IREQUE, int& IRETUR,
              int IUNIT, int& NUMOUT, int L, double J, double RSTART,
              double STEPSZ, int NUMPTS, double* ARRAY1, double* ARRAY2,
              char* ID);
void SHAPE(char8 ALIAS, int* MYINTS, int IPOTYP, int IREQUE, int& IRETUR,
           int IUNIT, int& NUMOUT, int L, double J, double RSTART,
           double STEPSZ, int NUMPTS, double* ARRAY1, double* ARRAY2,
           char* ID);
void SPLINE_linkule(char8 ALIAS, int* MYINTS, int IPOTYP, int IREQUE, int& IRETUR,
                    int IUNIT, int& NUMOUT, int L, double J, double RSTART,
                    double STEPSZ, int NUMPTS, double* ARRAY1, double* ARRAY2,
                    char* ID);
void TWOSHApe(char8 ALIAS, int* MYINTS, int IPOTYP, int IREQUE, int& IRETUR,
              int IUNIT, int& NUMOUT, int L, double J, double RSTART,
              double STEPSZ, int NUMPTS, double* ARRAY1, double* ARRAY2,
              char* ID);

// --- linkulesfitters.f — wavefunction linkules ---
complex16 CK_func(double R);
void CKSET(double RR1, double STEPSZ, double A1, double B1, double DC1, double D1,
           double A2, double B2, double DC2, double D2);
void OHTA(char8 ALIAS, int* MYINTS, int IPOTYP, int IREQUE, int& IRETUR,
          int IUNIT, int& NUMOUT, int L, double J, double RSTART,
          double STEPSZ, int NUMPTS, double* WAVR, double* WAVI,
          char* ID, int* IPARAM);
void RAWITSch(char8 ALIAS, int* MYINTS, int IPOTYP, int IREQUE, int& IRETUR,
              int IUNIT, int& NUMOUT, int L, double J, double RSTART,
              double STEPSZ, int NUMPTS, double* WAVR, double* WAVI,
              char* ID, int* IPARAM);
void REID(char8 ALIAS, real4* SAVINT, int IPOTYP, int IREQUE, int& IRETUR,
          int IUNIT, int& NUMOUT, int L, double& JP, double RSTART,
          double STEPSZ, int NUMPTS, double* ARRAY1, double* ARRAY2,
          double* FLT, double* TEMPVS, double* PARAM, int* INTGER_arr,
          double* JBLOCK_arr, int* ISWTCH, double* INTRNL_arr, double* FINTRN,
          double* CNSTNT_arr, double* WAVBL, double* FKANDM,
          char* ID, double* ALLOC_arr, int* ILLOC_arr, int FACFR4,
          integer4* LOC_arr, int* LENGTH_arr,
          int (*NALLOC_fn)(const char*, int), int (*NAMLOC_fn)(const char*),
          int* IPARAM);

// --- linkulesfitters.f — fitter routines ---
void CMTRIX(int N, int M, double* RJ, double* F, double* RJTJ, double* RJTF);
double FCUBIC(double F0, double G0, double F1, double G1, double XLAM);
void GENINV(int M, int N, double* A, int IA, double* W);
void MINMON(int ITCNT, int NFCALL, int NJCALL, int NPARAM, int NPTS,
            double* PARAM, double* F, double& FS);
void MOVE(double* X0, double* X, int N);
void MTMM(double* A1, double* A2, double* RESULT, int N1, int N2, int N3);
void MATMPY(int M, int N, double* H, double* G, double* S);
void RANDU(int& IX, double& YFL);

// --- source.f subroutines (Phase 8) ---
// I/O and parsing
void NXINT(int& INTARG, int& iret);
int NXWORD(char* CVARG);
int NXVAL(double& VALARG);
int NXVALF(double& VALARG);
int NXHINT(int& INTARG);
int FITKEY(char* CVARG);
void NEWCD();
void MSCAN(int KEYEX, int& ITOKEN, char8& CVALU, double& VALU, int& IVALU,
           char* MESS, int& IMESS, char& STOP);
void DATAIN(int& IRET);
void GETIT(int& IRET);
int GETNUM(int& NUM);
void DEFINE(char8 ANAME, int& ILEN, int& ISIZE);
void QVSCAN(char8* GUY, double* ESTARS, int* NODVLS, int* LVALS,
            int* JVALS, int* IPARYS, int& IRTN);
void LSTKEY(int VALSW, int NUMFLT, char8* FLTWRD, int NUMINT, char8* INTWRD,
            int NUMJ, char8* JWRDS, int NUMSWW, char8* SWTWRD, int* SWTNUM,
            int* SWTSET, int NUMKEY_arg, char8* KEYWRD, int NUMALI, char8* ALIAI);
void ALSORT(char8* A, int* IX, int N);

// Control
void CONTRL(int& IGOTO, int& JGOTO, int& IPRM1, int& IPRM2, int& IPRM3, double& PRM1);
void DEFALT(int NUMFLT, int NUMINT, int NUMJ, int NUMSW, int RESTSW);
void PARAM(int& IRET);
void REACTN(int& IRET);
void GSINFO(int KZ, int KA, int& JVAL, int& PARITY, char* JPI, int& NOTTAB);

// Channels
void BASCPL(int& IRET);
void BASLBL();
void CCHAN(int& IRET);
void CCOUP(int& IRET);
void CHANEL(int& IRET);
void CHNSCN(char8* GUY, double* ESTARS, int* NODVLS, int* LVALS, int* JVALS,
            int* IPARYS, int& IRTN);
void CLRCHN(int ICHAN);
void CTRAN(int& IRET);
void SETCHN(int ICHAN, int& IRET);
void SETFG(int L, int NWP, double* FG, int NSTEP, int NBAK);

// Potentials
void WOODSX(int NUMPTS, double RSTRT, double STEPSZ, double* VRAY,
            int& N1, int& N2, int ITYPE, double V, double R, double A, double POWER);
void DEFWOO(int NUMPTS, double RSTRT, double STEPSZ, double* VRAY,
            int& N1, int& N2, int ITYPE, double V, double R, double A, double POWER,
            int IPORT, int LXMIN, int* LDEFMX, double* RBETAS,
            int* LXS, int NUMPT);
void MAKPOT(int NWPIN, int& IRET);
void SETPOT(int& IRET);
void WAVPOT(int& IRET);
void FFACST(int& IRTN);
void SETSPT(int& ISIZE, int LOC, int IPASS, int JBGBIN, int JBIN,
            int PBUGSW, int* BASDEF, int MBASDF, int NBASDF);
void SFROMI(int LI, int ILI, double* SMATR, double* SMATI, int* INDXS,
            double* XIREAL, double* XIIMAG, int* IINDEX, int NI,
            int* INDXDW, int* IDWFI, int* IDWFO, double* RIROAB,
            double* ATERM, double FACTOR, int PINFSW, int& NLINE);
void BEBETA(int MODEL, int IORD, int I12, int& IBETAC, int MEORD,
            double* XCS, int LXX2, int* LXX2S, int IZ, double AMVAL,
            double RC, double DELTAE, int J2, int J1, int K2, int K1,
            int IBETAN, double* XNS, double RNUC, double SCALE, int IPT,
            double& RMENUC, double& RMECOU, int& IRTN, double* VALS);

// Bound state
void BOUND(int& IRET, int IVOUT, int IWF2);
int BSPROD(double& FPFT, int ITYPE, double RA, double RB, double X,
           int ISCAT, int NAIT, double& RP, double& RT);
void BSSET();
void CLINTS(double RLOWER, double ETAI, double ETAF, double FKI, double FKF,
            double SIGI, double SIGF, double ACCURA_arg, double ASMULT,
            double& FFINT, double& FGINT, double& GFINT, double& GGINT,
            double* PTS, double* WTS, double* FA, double* FPA,
            double* GA, double* GPA, double* WORK,
            int N, int LI, int LF, int NTERMS, int NPTS,
            int& IRET, int IPRINT);
void JPTOLX(int L, int LAS, int JP, int NWP, double* SJR, double* SJI,
            int* INDX, double* SLX);

// Scattering
// COULIN — Coulomb integrals by recursion (source.f L10647)
// SUBROUTINE COULIN(N, MAXDEL, LMIN, LMAX, ETAOUT, AKOUT, SIGOUT,
//   ETAIN, AKIN, SIGIN, R, ALLSW, FF, FG, GF, GG, LDLDIM,
//   ACCURA, ASMULT, NTERMS, NPTS, WORK, FI, FO, GI, GO, STARTS,
//   IPRINT, IRET, CLTIME)
void COULIN(int N, int MAXDEL, int LMIN_arg, int LMAX_arg,
            double ETAOUT, double AKOUT, double* SIGOUT,
            double ETAIN,  double AKIN,  double* SIGIN,
            double R, int ALLSW,
            double* FF, double* FG, double* GF, double* GG,
            int LDLDIM,
            double ACCURA, double ASMULT, int NTERMS, int NPTS,
            double* WORK, double* FI, double* FO, double* GI, double* GO,
            double* STARTS,
            int IPRINT, int& IRET, double& CLTIME);
double COULNG(int L, double ETA, double RHO, double ANORM);
void COULST(int& IRTN);
void GETSCT(int& IRET);
void WAVEF(int& IRET);
void WAVELJ(int L, int JP, int NWP, int NUMPTS, float* RGRID, float* WAVER,
            float* WAVEI, double* WAVR, double* WAVI, double* VREAL,
            double* VIMAG, double* VCENT);
// NOTE: TCOEF with 5 args (the real signature from source.f L34242):
double TCOEF_5(int ITYP, int LAS, int L, int JP, int JSP);
void WAVETC(int L, int LAS, int JP, int NWP, int NWAVF, int NUMPTS,
            float* RGRID, float* WAVER, float* WAVEI,
            double* WAVR_flat, double* WAVI_flat,
            double* VREAL, double* VIMAG, double* VCENT);
void WAVINH(int NFIRST, int NSTEP, double* WAVR, double* WAVI,
            double* RHSR, double* RHSI, double& FLP, int NSTPDM);
void WFGET(int L, int LAS, int JP, int NWP, int NUMPTS, float* RGRID,
           float* WAVER, float* WAVEI,
           double* WAVR, double* WAVI, double* VREAL, double* VIMAG,
           double* VCENT);
void WAVSET(int& IRET, int ALONSW, int FITSW);
void LCRITL(double E, double FK, double ETA, double RC, int NSTEPS,
            double RSTART, double STEPSZ, int IV, int IW, int NGAUSS,
            int IPRINT, int& LC1, int& LC2, int& LC);
void SWKB(double E, double FK, int L, double EPS, int IPRINT,
          double ETA, double RC,
          int NSTEPS, double ROSTRT, double H, double VCONS,
          double* VARRAY, double* VCOEFS_flat, double* WARRAY, double* WCOEFS_flat,
          int NGAUSS, double* PTS, double* WTS,
          double& DELTR, double& DELTIM, double& THETA, double& THETC,
          double& RHOT, double& ROTC, double& DUEFF);
void TMATCH(double* WAVR, double* WAVI, double* CR, double* CI,
            double* SMATS, int JP, int NWP, int NWAVF);

// Angular/cross sections
void AMPCAL(double ANGLE, int ELSW, int NSPL, int LMN, int LMX, int LSKP,
            int LXMX, int IDNPAR, double ETA, double UK, double SIGZRO,
            int LOHISW, int LHI,
            int* JTOCS, float* BETAS, double* ALOWFC,
            double* F, double* FLO, double* FHI, double* FERROR,
            double* FCOUL, double* PLM, double* CONTR, double* CONTI,
            int LEBACK, double* FEPSLO, int& NFLOP);
void ANAPOW(double ANGMIN, double ANGMAX, double ANGSTP, int JA, int JB,
            int JBIGA, int JBIGB, int NSPL, int IPARIT, int PRNTSW,
            int PBUGSW, const char* APNAME, double ELAB, const char* CHNM,
            int IF_arg, int ITOCS, int& IAPOW);
void ANGSET(int& IRTN);
void BETCAL(int ELSW, double UK, int JSP, int NSPL, int LMN, int LMX, int LSKP,
            int LXMN, int LXMX, int LDELMX, int NTEMPS, int ISTAT, int PRNTSW,
            double& SIGTOT, double& SIGREA,
            int* JTOCS, double* S, float* SMAG, float* SPHASE, double* SIGIN, double* SIGOT,
            float* BETAS, double* TOTLX, double* TOTMX, double* ALOWFC, double* TEMPS);
void CALANG(int GRADSW, int LLMIN, int LLMAX, int LSKIP, int MPARA,
            double* FS, double* XS);
double CALFUN(int NPTST, int MPARA, double* FS, double* XS);
void CALGRA(int NPTST, int MPARA, double* FS, double* XS, double* FGRAD, double* TGRAD);
void ELDCS(double UK, double ETA, double ANGMIN_arg, double ANGMAX_arg, double ANGSTP_arg,
           int LMN, int LMX, int LSKP, int ISTAT, int JSP,
           int ISMAT, int ITOC, int NSPL, int ISIG, double ELAB, double TAU,
           int PRNTSW, int& ITORUT, const char* RUTNAM,
           int FKEPSW, int& IFEL, const char* FNAME, double& SIGREA,
           int& IRUTH, const char* RUTHNM);
void XSECTN(int& IRTN);

// Coupled channels
void A12(int LI, int LXMIN_a, int LXMAX_a, int LMIN_a, int LMAX_a, double* XLAM,
         int* NLAM, double* A12VL, double* DMSVAL,
         int* JA12S, int& JA12M, int& JA12N, int& JA12AN, int* IINDEX,
         int& IHMAX, int& LOMNMN, int& LOMXMX, double* DINTS, double* OUTTMP,
         double* XLOTMP, int* LXTMP, int LHSM1, int LA12VL, int IPRINT, int INIT,
         int* INDXDW, int& NDW, int& NI, int* IDWFI, int* IDWFO);
void CCDUMP();
void CCMTCH(int ILI, int J, int IB, double* SMATS, double* CFAC_out,
            int NBASDF_arg, double* WAVR, double* WAVI, int NPTSAV,
            double* HOMO, int MBASDF_arg, int* BASDEF,
            int MBASCP_arg, int* BASCUP, float* RBASCP,
            int MCHNVL, double* CHNVAL,
            float* BFACS, int MBNDX, int* BINDEX,
            float* FG, double VSCALE);
void CCONV(double* CSMATS, int ITR, int& ITEST, int NBASDF_arg,
           double* PADE, int MAXIT1, double* ELAST, double* EXR,
           double* EXI, int* INDXS, int* CHNDEF, int MCHNDF_arg,
           int* BASDEF, int MBASDF_arg, int IJ, int NCHNDF_arg,
           int JT, int& NLINE, int& IHEAD, int& NFLOP);
void COUPLN(int& IRET);
void DRIVE(double* RHSR, double* RHSI, int I5, int I6, int NFIRST,
           int ITR, int MBASCP_arg, int* BASCUP, float* RBASCP,
           float* NUCH, float* COULH, int NPTIN, int NPTSAV,
           double* INHR, double* INHI, double& FLP, double VMULT);
void INELDC(int& IRTN);
void INELD2(float TSTART, int NUMHIN, int NUMTIM, int NUMTRM, int NUMANG,
            int NUM3, int NUM4, int NUM5, int NUM6);
void INGRST(int& IRET);
void INRDIN(int& IRTN);

// Fitting
void FITEL(int& IRET);
void FITEL2(int& IRET);
void FITINP(int& IRET, int NUMFLT, char8* FLTWRD, int NUMALI, char8* ALIAI);
void LINTRP();
void LMCFUN(int ITYPE, int MPARA, int MPTS, double* X, double* FS,
            double* FJACOB, int& IERR);
double LXTRPM(int IEXTYP, double BARA, double B, double BARC, double BARL,
              double DLMAX, int LX, int LDEL, double* ETAS, double WEEBOY);
void LXTRP1(int IEXTYP, int N, int& ISIG, int IPRINT, double* XVALS, double* SVALS,
            double& FLCRIT, double& AVAL, double& WID, double STEPIN, double DELIN,
            double& BARL, double& BARA, double& B, double& BARC, double DLMAX,
            double& CHISQ, int LX, int LDEL, double* ETAS, int MCHN, int JP, int JT);
void LXTRP2(int IEXTYP, double BARA, double B, double BARC, double BARL,
            double DLMAX, int LI, double& SIZE, int LX, int LDEL, double* ETAS);
void MAKDER(int NSTEPS, int NUMVS, int MXVSIN, int MPARA,
            double* DERIV, double* POTRL, double* POTIM, double STEPSZ);
void SDERIV(int& IRTN);
void SECDER(int& IRTN);
void DERCHK(int NPARA, int NPTS, double STEP, double* PARAM, double* FS,
            double* FGRAD, double* TGRAD, double* FSAVE);
void GENBNX(int IPASS, int* NBINDX, int* MBINDX, int* BINDEX, int NBASCP,
            int* BASCP, int MBASCP, int NBASDF, int* BASDF, int MBASDF);
void GETBFC(int ILI, int J, double DELTA, double RASYMP, int NMBFC,
            float* BFACIN, float* BFACOT, int MBNDX, int* BINDEX,
            int MCHNDF, int* CHNDEF, int MCHNVL, double* CHNVAL,
            int MBASDF, int* BASDEF, int NBASDF, float* FG);
void SETFIT(int& IRET);
void SETBFC(int NBINDX, int MBINDX, int* BINDEX, int MCHNDF, int* CHNDF_i,
            float* CHNDF_f, int* NUMJS, int LMIN2, double* R2S, int LMXMX,
            double* SIG1, double* SIG2, double DUM, int FLGSW, int IPASS,
            int NMFFAC, float* CL2FF, double* WRKFF, int IDIM1, int LDLDIM,
            double* WRKWK, double* WRKFI, int IDIM2, double* WRKST,
            int& IRET, double& CLTIME);
void SETBRN(int NBINDX, int MBINDX, int* BINDEX, int MCHNDF, int* CHNDF_i,
            int MBASDF, int* BASDF_i, int MBASCP, int* BASCP_i, float* BASCP_f,
            int NUMLIS, double* LIS, int LMIN2, float* CL2FF,
            double* INDXS, double* DELSR, double* DELSI);

// Misc utilities
double AITLAG(double X, double STPINV, double* TABLE, int NTABLE, int NAIT);
void AKFCR(int L1, int J1, int L2, int J2, int LXX2, int J, double& AKREST);
void CUBMAP(int MAPTYP, double XLO, double XMID, double XHI, double GAMMA,
            double* ARGS, double* WTS, int NPTS);
void CUPAB(int IORDER, int IPT, int* LXX2S, int LXX2,
           double* RMENUC, double* RMECOU, double* RNUC, int* IZPT,
           double& RMEGN, double& RMEGC, int& IPT1);
void CUPSPN(int& IRTN, int LXX2, int IPT, int* SCHN,
            int* LXX2S, int* JNUC, double RMEGN, double RMEGC,
            double& XN, double& XC);
void DUMMY1();
void DUMMY2();
void EPSLON(double* XIN, int N, double* FRET, double& DEPS, double& DERR, int& IER);
void GRDSET(int& IRTN);
void INIT8(int I, double X8);
void INIT4(int I, float x4);
void INIT4i(int I, int i4);
void INIT2(int I, int16_t x2);
void INIT1(int I, int8_t x1);
void INTRCF(int NCHN, int JP, int JT, int LX, int LDEL,
            int NUMLIS, int* LIS, int LSTEP, int LXMAX,
            double* SR, double* SI, int LOFF1, int ISDIM1, int LISW,
            int LSTRT, int LMAX_arg, double* SOUT, float* SOUT4,
            int LOFF2, int ISDIM2, int LBASE2, int FOURSW,
            double& BIGBOY, int& LPEAK, double* CMPLXL, double* CMPLXS);
void IXSORT(int* IX, double* X, int N);
void LINLSQ(int IFUNC, int N, double* XVALS, double* SVALS, double& CVAL,
             double& AVAL, double& B, double& CHI, int PBUGSW);
void MEBDEF(int IORD, int I12, int IBN, double* XNS, int IBC, double* XCS,
            int MEORD, int LXX2, int* LXX2S, int IZ,
            double RNUC, double RC, double DELTAE, int J2, int J1, int K2, int K1,
            double& RMENUC, double& RMECOU, int& IRTN, int IPT);
void MEBROT(int IORD, int I12, int IBN, double* XNS, int IBC, double* XCS,
            int MEORD, int LXX2, int* LXX2S, int IZ,
            double RNUC, double RC, double DELTAE, int J2, int J1, int K2, int K1,
            double& RMENUC, double& RMECOU, int& IRTN);
void MEBVIB(int IORD, int I12, int IBN, double* XNS, int IBC, double* XCS,
            int MEORD, int LXX2, int* LXX2S, int IZ,
            double RNUC, double RC, double DELTAE, int J2, int J1, int K2, int K1,
            double& RMENUC, double& RMECOU, int& IRTN);
void MUELCO(int KA, int KB, int KBIGA, int KBIGB, int JA, int JB,
            int JBIGA, int JBIGB, int NQA, int NQB, int NQBIGA, int NQBIGB,
            int NSPL, int IPARIT, int* ITOC, double* COEF, double* COEP,
            double* COET, double* COEX, double* TERMK, int PBUGSW);
void print_G(int w, int d, double val);
void PARPRT(int I, const char* NAME, double VPAR, double RPAR, double APAR,
            double RPAR0, double P5, const char* P5NAME, int& IRET);
void PHSPRT(bool PRNTSW);
void PRBPRT(int& IRET);
void PRTCHI(int& IRET);
void PRTDER(int& IRET);
void RCASYM(int L, double ETA, double RHO, int IPRN, double SIGL,
            double* ZETA, double* PHI, double* DZETA,
            double* F, double* FP, double* G, double* GP,
            double* Z, double* DZSQ, double* S, double* ZINV,
            double EPS, int& NMAX, int& NTZ, int& ISIG);
void REQUIV(double& RC, double& RC0, int IZ, int IA);
double RTXLNX(double X);
void SAVEHS();
void SAVECM();
void USEHS();
void SETINT(int N);
void SETLOG(int N);
double TCOEF(int ITYP, int LAS, int L, int JP, int JSP);
void TEMP1();
void TEMP2();
void TEMP3();
void VCSQ12(double RVAL, double& X, int NWP);
void SETVSQ(double RR1, double RR2, int IZ1, int IZ2, int K);

#endif // PTOLEMY_FORWARD_H
