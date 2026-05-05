// source_fitting.cpp — Fitting, chi-square, L-extrapolation, derivatives
// Translated from source.f: CALFUN+CALGRA, FITEL, FITEL2, FITINP, LINTRP, LMCFUN,
//   LXTRPM, LXTRP1, LXTRP2, MAKDER, SDERIV, SECDER, DERCHK, GENBNX, GETBFC,
//   SETFIT, SETBFC, SETBRN
// Phase 8 — stub bodies to be filled in from Fortran source

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include "ptolemy_alloc.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_forward.h"
#include "ptolemy_io.h"
#include <cstdio>
#include <cstring>
#include <cmath>

void FITEL(int& IRET)
{ IRET = 0; std::printf(" FITEL: stub — translate from source.f L15647-16315\n"); }

void FITEL2(int& IRET)
{ IRET = 0; std::printf(" FITEL2: stub — translate from source.f L16317-16915\n"); }

void FITINP(int& IRET, int NUMFLT, char8* FLTWRD, int NUMALI, char8* ALIAI)
{ IRET = 0; std::printf(" FITINP: stub — translate from source.f L16917-17283\n"); }

// LINTRP → lintrp_translated.cpp

void LMCFUN(int ITYPE, int MPARA, int MPTS, double* X, double* FS,
            double* FJACOB, int& IERR)
{ IERR = 0; std::printf(" LMCFUN: stub — translate from source.f L24798-24850\n"); }

// 4-argument RTXLNX: finds real solution of  A*X + B*LN(X) + C = 0
extern double RTXLNX(double A, double B, double C, double ACC);

double LXTRPM(int IEXTYP, double BARA, double B, double BARC, double BARL,
              double DLMAX, int LX, int LDEL, double* ETAS, double WEEBOY)
{ return 0.0; }

void LXTRP1(int IEXTYP, int N, int& ISIG, int IPRINT, double* XVALS, double* SVALS,
            double& FLCRIT, double& AVAL, double& WID, double STEPIN, double DELIN,
            double& BARL, double& BARA, double& B, double& BARC, double DLMAX,
            double& CHISQ, int LX, int LDEL, double* ETAS, int MCHN, int JP, int JT)
{ ISIG = 0; BARA = 0; B = 0; BARL = 0; BARC = 1; FLCRIT = 0; AVAL = 0; WID = 1; CHISQ = 0; }

void LXTRP2(int IEXTYP, double BARA, double B, double BARC, double BARL,
            double DLMAX, int LI, double& SIZE, int LX, int LDEL, double* ETAS)
{ SIZE = 0; }

void MAKDER(int NSTEPS, int NUMVS, int MXVSIN, int MPARA,
            double* DERIV, double* POTRL, double* POTIM, double STEPSZ)
{ std::printf(" MAKDER: stub — translate from source.f L25596-25748\n"); }

void SDERIV(int& IRTN)
{ IRTN = 0; std::printf(" SDERIV: stub — translate from source.f L31178-31291\n"); }

void SECDER(int& IRTN)
{ IRTN = 0; std::printf(" SECDER: stub — translate from source.f L31293-31510\n"); }

void DERCHK(int NPARA, int NPTS, double STEP, double* PARAM_arr, double* FS,
            double* FGRAD, double* TGRAD, double* FSAVE)
{ std::printf(" DERCHK: stub — translate from source.f L14564-14633\n"); }

void GENBNX(int IPASS, int* NBINDX, int* MBINDX, int* BINDEX, int NBASCP,
            int* BASCP, int MBASCP, int NBASDF, int* BASDF, int MBASDF)
{ std::printf(" GENBNX: stub — translate from source.f L17285-17556\n"); }

void GETBFC(int ILI, int J, double* DELTA, double RASYMP, int NMBFC,
            double* BFACIN, double* BFACOT, int NBASDF, int* BASDF,
            int MBASDF)
{ std::printf(" GETBFC: stub — translate from source.f L17558-17753\n"); }

void SETFIT(int& IRET)
{ IRET = 0; std::printf(" SETFIT: stub — translate from source.f L32576-32878\n"); }

void SETBFC(int NBINDX, int MBINDX, int* BINDEX, int MCHNDF, int* CHNDF_i,
            float* CHNDF_f, int* NUMJS_arr, int LMIN2, double* R2S, int LMXMX,
            double* SIG1, double* SIG2, double DUM, int FLGSW, int IPASS,
            int NMFFAC, float* CL2FF, double* WRKFF, int IDIM1, int LDLDIM,
            double* WRKWK, double* WRKFI, int IDIM2, double* WRKST,
            int& IRET, double& CLTIME)
{ IRET = 0; CLTIME = 0.0; std::printf(" SETBFC: stub — translate from source.f L31512-31798\n"); }

void SETBRN(int NBINDX, int MBINDX, int* BINDEX, int MCHNDF, int* CHNDF_i,
            int MBASDF, int* BASDF_i, int MBASCP, int* BASCP_i, float* BASCP_f,
            int NUMLIS, double* LIS, int LMIN2, float* CL2FF,
            double* INDXS_arr, double* DELSR, double* DELSI)
{ std::printf(" SETBRN: stub — translate from source.f L31800-32000\n"); }
