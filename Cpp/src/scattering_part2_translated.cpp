// scattering_part2_translated.cpp — Verbatim C++ translation of LCRITL, SWKB, TMATCH
// Translated from source.f lines 23149-23405, 33975-34240, 34575-34778

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include "ptolemy_alloc.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_forward.h"
#include "ptolemy_io.h"
#include <cstdio>
#include <cstring>
#include <cmath>

// ============================================================================
// SUBROUTINE LCRITL — lines 23149-23405
//
// USE CLASSICAL DEFLECTION AND WKB PHASES TO ESTIMATE CRITICAL L
// LC IS TAKEN TO BE GREATER OF TWO ESTIMATES LC1 AND LC2
// LC1 IS THE FIRST L FOR WHICH |S(L)| > (|S|MIN + 1)/2.
// LC2 IS L-VALUE OF ALGEBRAICALLY MINIMUM THETA(L)
// WHERE THETA(L) IS THE CLASSICAL DEFLECTION FUNCTION.
// ============================================================================
void LCRITL(double E, double FK, double ETA, double RC, int NSTEPS,
            double RSTART, double STEPSZ, int IV, int IW, int NGAUSS,
            int IPRINT, int& LC1, int& LC2, int& LC)
{
    // COMMON block aliases
    auto& PROBLM = SWITCH.PROBLM;
    auto* Z = LOCPTRS.Z;  // 1-based

    // Local variables
    int MINSW;  // LOGICAL
    double EPS = 1.0e-12;

    // SET MINIMUM AND MAXIMUM L FOR SEARCH
    // ROW IS THE L WHERE COULOMB CLASSICAL TURNING POINT IS MAX(RC,1.2)

    MINSW = FALSE_F;
    double RCL = DMAX1(RC, 1.20e0);
    double ROW = FK * RCL;
    ROW = ROW * ROW - 2 * ETA * ROW + .25;
    if (ROW < 0) ROW = 0;
    if (ROW > 0) ROW = 0.5 + DSQRT(ROW);
    int LMIN = (int)(.2 * ROW - 4.);
    LMIN = MAX0(LMIN, 0);
    int LMAX = (int)(1.5 * ROW);
    LMAX = MAX0(LMAX, LMIN + 50);
    LC1 = 0;
    LC2 = 0;
    LC = 0;
    double DLSMAT = 0;
    double DLDFLC = 0;
    int NOLS = LMAX - LMIN + 1;
    double AMPMIN = 1.0e20;
    double VCONS = 1;
    if (PROBLM == 24) VCONS = 0;

    // ALLOCATE SPACE FOR GAUSS POINTS AND WKB STUFF
    // SET UP STARTING ADDRESSES

    int IWKBPT = NALLOC("WKBPTS", NGAUSS);
    int IWKBWT = NALLOC("WKBWTS", NGAUSS);
    int NSEMIC = 5 * NOLS;
    int ISEMIC = NALLOC("SCSTUF", NSEMIC);
    int ICUBIC = NALLOC("WKBCUBIC", 7 * NSTEPS);
    int LRVAL = Z[ICUBIC];
    int LVCUB = LRVAL + NSTEPS;
    int LWCUB = LVCUB + 3 * NSTEPS;
    int LV = Z[IV];
    int LW = Z[IW];
    int IAMP = Z[ISEMIC];
    int LAMP = IAMP - 1;
    int IDNUC = IAMP + NOLS;
    int LDNUC = IDNUC - 1;
    int IDEFL = IDNUC + NOLS;
    int LDEFL = IDEFL - 1;
    int ICTP = IDEFL + NOLS;
    int LCTP = ICTP - 1;
    int ICTPC = ICTP + NOLS;
    int LCTPC = ICTPC - 1;

    // DETERMINE UNMAPPED GAUSS POINTS FOR (-1,1)

    GAUSSL(NGAUSS, ALLOC_base(Z[IWKBPT]), ALLOC_base(Z[IWKBWT]));

    // GET CUBIC SPLINES FOR THE POTENTIALS

    double ROSTRT = FK * RSTART;
    double RHO = ROSTRT;
    double H = FK * STEPSZ;
    for (int I = 1; I <= NSTEPS; I++) {
        ALLOC(LRVAL - 1 + I) = RHO;
        RHO = RHO + H;
    }
    double ROEND = RHO - H - EPS;

    SPLNCB(NSTEPS, ALLOC_base(LRVAL), ALLOC_base(LV), ALLOC_base(LVCUB),
           ALLOC_base(LVCUB + NSTEPS), ALLOC_base(LVCUB + 2 * NSTEPS));

    SPLNCB(NSTEPS, ALLOC_base(LRVAL), ALLOC_base(LW), ALLOC_base(LWCUB),
           ALLOC_base(LWCUB + NSTEPS), ALLOC_base(LWCUB + 2 * NSTEPS));

    // CALCULATE AND STORE WKB S-MATRIX ELEMENTS (AMP AND ARG)
    // CLASSICAL DEFLECTION FUNCTIONS(DEFL) AND CLASSICAL TPS(CTP CTPC)

    double DELTR, DELTIM, THETA, THETC, RHOT, ROTC, DUEFF;
    double AMPL;
    int L;

    for (int I = 1; I <= NOLS; I++) {
        L = LMIN + I - 1;
        // SWKB uses 0-based array macros (VARR[N-1], PTS[I-1]), so pass &ALLOC(x)
        // not ALLOC_base(x) (which is 1-based). Off-by-one was causing SWKB to fail.
        SWKB(E, FK, L, EPS, IPRINT, ETA, RC,
             NSTEPS, ROSTRT, H, VCONS, &ALLOC(LV), &ALLOC(LVCUB), &ALLOC(LW),
             &ALLOC(LWCUB), NGAUSS, &ALLOC(Z[IWKBPT]), &ALLOC(Z[IWKBWT]),
             DELTR, DELTIM, THETA, THETC, RHOT, ROTC, DUEFF);
        if ((DELTR != 5.0e0) || (DELTIM != 5.0e0)) goto L250;
        std::printf("0**** ERROR RETURN IN LCRITL FROM SWKB \n");
        goto L640;
    L250:
        ALLOC(LCTP + I) = RHOT;
        ALLOC(LCTPC + I) = ROTC;
        AMPL = 0;
        if (DELTIM < 20) AMPL = DEXP(-2.0e0 * DELTIM);
        ALLOC(LAMP + I) = AMPL;
        AMPMIN = DMIN1(AMPMIN, AMPL);
        if (THETA > THETC) THETA = THETC;
        ALLOC(LDEFL + I) = THETA;
        ALLOC(LDNUC + I) = 2 * DELTR;
        if (RHOT >= ROEND) goto L280;

        // CHECK ON THE PROGRESS OF THE DEFLECTION FUNCTION

        if (AMPL < .05 || I == 1) goto L299;
        if (MINSW) goto L260;
        MINSW = (THETA > ALLOC(LDEFL + I - 1) + EPS) ? TRUE_F : FALSE_F;
        goto L299;
    L260:
        if (THETA > ALLOC(LDEFL + I - 1) - EPS) goto L299;

        // DEFLECTION FUNCTION IS FALLING FROM A MAX AFTER A MIN.
        // IF |S| IS CLOSE TO 1 WE STOP

        if (AMPL < .95 + .05 * AMPMIN) goto L299;
    L280:
        LMAX = L;
        goto L300;
    L299: ;
    }

L300:
    NOLS = LMAX - LMIN + 1;

    // PRINT WKB PHASES AND CLASSICAL DEFLECTION FUNCTION

    if (IPRINT == 0) goto L500;
    std::printf("\n     ---  WKB S-MATRIX ELEMENTS AND CLASSICAL TPS ---\n"
                "  L      AMP(S(L))        ARG(S(L))        DEFL.ANGLE "
                "       CLASSICAL TP      COULOMB TP\n");
    for (int I = 1; I <= NOLS; I++) {
        L = LMIN + I - 1;
        double CTR = ALLOC(LCTP + I) / FK;
        double CTRC = ALLOC(LCTPC + I) / FK;
        double AMP = ALLOC(LAMP + I);
        double DNUC = ALLOC(LDNUC + I);
        double DEFL = ALLOC(LDEFL + I);
        std::printf("  %3d%16.5G%16.5G%16.5G%16.5G%16.5G\n", L, AMP, DNUC, DEFL, CTR, CTRC);
    }

    // ESTIMATE CRITICAL L FROM PHASE-SHIFTS AND DEFLECTION FN

L500:
    MINSW = FALSE_F;
    AMPMIN = 0.5 * (AMPMIN + 1.);

    for (int I = 1; I <= NOLS; I++) {
        int I1 = NOLS - I + 1;
        L = LMAX + 1 - I;
        if (DLSMAT != 0) goto L550;
        if (ALLOC(LAMP + I1) > AMPMIN) goto L550;

        // |S| < .5*(MIN(|S|)+1).  INTERPOLATE UNLESS LAST POINT

        if (AMPMIN > .99) goto L550;
        if (I > 1) goto L520;
        std::printf("\n**** LCRITL PICKED LMAX TOO SMALL TO FIND"
                    " |S| > .5*(MIN(|S|)+1):%5d%15.5G"
                    "  LCRIT SET TO LMAX\n", LMAX, ALLOC(LAMP + I1));
        DLSMAT = LMAX;
        goto L550;

    L520:
        DLSMAT = L + (AMPMIN - ALLOC(LAMP + I1)) / (ALLOC(LAMP + I1 + 1)
                 - ALLOC(LAMP + I1));

        // LOOK FOR OUTERMOST MINIMUM OF DEFLECT FUNCTION

    L550:
        if (DLDFLC != 0 || L == LMAX) goto L599;
        if (MINSW) goto L560;
        MINSW = (ALLOC(LDEFL + I1) < ALLOC(LDEFL + I1 + 1) - EPS) ? TRUE_F : FALSE_F;
        goto L599;
        // HAVE WE PASSED A MINIMUM
    L560:
        if (ALLOC(LDEFL + I1) < ALLOC(LDEFL + I1 + 1) + EPS) goto L599;
        { double C = ALLOC(LDEFL + I1) - 2 * ALLOC(LDEFL + I1 + 1)
                     + ALLOC(LDEFL + I1 + 2);
          DLDFLC = L + 1 + (ALLOC(LDEFL + I1) - ALLOC(LDEFL + I1 + 2)) / (2 * C);
        }

    L599: ;
    }

    LC1 = (int)(DLSMAT + .5);
    LC2 = (int)(DLDFLC + .5);
    LC = MAX0(LC1, LC2);

    if (IPRINT != 0)
        std::printf("\n   --- 1/2 POINT OF |S(L)|  =  %15.5G"
                    "     MINIMUM OF DEFLECTION FUNC =%15.5G"
                    "     L CRIT =%5d\n", DLSMAT, DLDFLC, LC);

L640:
    if (DLSMAT != 0 || DLDFLC != 0) goto L700;
    LC = (int)(ROW + .5);
    std::printf("\n**** NO CRITICAL L FOUND; SET LC =%5d"
                " = L FOR WHICH CLASSICAL COULOMB TURNING POINT =%10.3f\n", LC, RCL);

    // FREE ALLOCATOR AREAS

L700:
    Z[IWKBPT] = -Z[IWKBPT];
    Z[IWKBWT] = -Z[IWKBWT];
    Z[ISEMIC] = -Z[ISEMIC];
    Z[ICUBIC] = -Z[ICUBIC];

    return;
}


// ============================================================================
// SUBROUTINE SWKB — lines 33975-34240
//
// EVALUATES WKB PHASE-SHIFT FOR PARTIAL WAVE L
// AND THE CORRESPONDING CLASSICAL DEFLECTION FUNCTION
// ============================================================================
void SWKB(double E, double FK, int L, double EPS, int IPRINT,
          double ETA, double RC,
          int NSTEPS, double ROSTRT, double H, double VCONS,
          double* VARRAY, double* VCOEFS_flat, double* WARRAY, double* WCOEFS_flat,
          int NGAUSS, double* PTS, double* WTS,
          double& DELTR, double& DELTIM, double& THETA, double& THETC,
          double& RHOT, double& ROTC, double& DUEFF)
{
    // VCOEFS and WCOEFS are (NSTEPS,3) arrays — access as VCOEFS[col*NSTEPS + row]
    // But in Fortran column-major: VCOEFS(N,J) -> VCOEFS_flat[(J-1)*NSTEPS + (N-1)]
    // Using 1-based pointer arithmetic: element (N,J) is at offset (J-1)*NSTEPS + (N-1)
    // We'll define a macro-like lambda for access
    // VCOEFS(N,J) -> VCOEFS_flat[ (J-1)*NSTEPS + N - 1 ]  (1-based N, 1-based J)
    // But VARRAY, WARRAY, VCOEFS_flat, WCOEFS_flat are passed as pointers to ALLOC(LV) etc.
    // which are 1-based. So the pointer already points to element [1].
    // Actually — they are passed as &ALLOC(LV) which is ALLOC(LV), i.e., the address of element LV.
    // In the Fortran, VARRAY(N) means the N-th element starting from 1.
    // The C++ pointer starts at element [0] of the passed array.
    // So VARRAY[N-1] = Fortran VARRAY(N).
    // VCOEFS(N,1) = VCOEFS_flat[N-1], VCOEFS(N,2) = VCOEFS_flat[NSTEPS+N-1],
    // VCOEFS(N,3) = VCOEFS_flat[2*NSTEPS+N-1]

    // Helper macros (0-based pointer, 1-based indices)
    #define VARR(n) VARRAY[(n)-1]
    #define WARR(n) WARRAY[(n)-1]
    #define VCOF(n,j) VCOEFS_flat[((j)-1)*NSTEPS + (n)-1]
    #define WCOF(n,j) WCOEFS_flat[((j)-1)*NSTEPS + (n)-1]

    double PI = 3.14159265358979300e0;
    double RT2 = 1.41421356237309500e0;
    double TP = 1.0 / RT2;
    int NMAX = 50;
    double RHOC = FK * RC;
    double UC0 = (2 * ETA) / RHOC;

    // NOTE QUANTUM-MECHANICAL VALUE.  THIS AVOIDS SOME TROUBLES

    double FLL = L * (L + 1.0e0);
    double ECC = DSQRT(ETA * ETA + FLL);
    DELTR = 5.0e0;
    DELTIM = 5.0e0;

    // OUTER LIMIT OF THE WKB INTEGRALS WILL BE ASYMPTOPIA

    double ROMAX = ROSTRT + H * (NSTEPS - 1);

    // CALCULATE CTP IN PURE COULOMB POTENTIAL

    ROTC = ETA + ECC;
    if (ROTC >= RHOC || ETA == 0) goto L24;
    { double F1 = (3 * ETA - RHOC) / (2 * ETA);
      ROTC = DSQRT(F1 * F1 + FLL / (ETA * RHOC)) + F1;
      if (ROTC < 0) ROTC = 0;
      ROTC = DSQRT(ROTC) * RHOC;
    }

    // CALCULATE CTP IN FULL(NUCLEAR+COULOMB) POTENTIAL

    // FIND INTERVAL (RHO,RHO-H) CONTAINING CTP

L24:
    { double ULAS = 0;
      double RHO = ROMAX;
      int N = NSTEPS;
      int NITS = 0;
      double UEFF;
      double X;

    L2:
      UEFF = ((12.0 / (H * H)) + VCONS) - (12.0 / (H * H)) * VARR(N) + FLL / (RHO * RHO);
      if (UEFF > 1) goto L4;
      ULAS = UEFF;
      RHO = RHO - H;
      N = N - 1;
      if (RHO > 0) goto L2;
      RHOT = 0;
      RHO = 0;
      if (L == 0) goto L11;
      std::printf("\n   --- RHOT < H FOR L = %3d ---\n", L);
      RHOT = H;
      goto L11;
    L4:
      RHOT = RHO;
      if (RHO >= ROMAX) goto L160;
      if (DABS(UEFF - 1.0e0) > EPS) goto L5;
      goto L11;

      // USE REGULA FALSI TO LOCATE CTP WITHIN (RHO,RHO-H)

    L5:
      { double RO1 = RHO;
        double RO2 = RHO + H;
        double U1 = UEFF;
        double U2 = ULAS;
        double RO3, RO, ACC;
    L6:
        RO3 = (RO2 * (U1 - 1) + RO1 * (1 - U2)) / (U1 - U2);
        RO = RO3;
        // ASSIGN 7 TO IK / GO TO 710
        // --- inline label 710: compute UEFF from RO ---
        X = (RO - ROSTRT - (N - 1) * H);
        UEFF = ((12.0 / (H * H)) + VCONS) - (12.0 / (H * H)) * (VARR(N)
            + X * (VCOF(N, 1) + X * (VCOF(N, 2) + X * VCOF(N, 3))))
            + FLL / (RO * RO);
        // --- return to label 7 ---
        if (DABS(UEFF - 1.0e0) > EPS) goto L8;
        RHOT = RO3;
        goto L11;
    L8:
        if (UEFF < 1) goto L9;
        RO1 = RO3;
        U1 = UEFF;
        goto L10;
    L9:
        RO2 = RO3;
        U2 = UEFF;
    L10:
        NITS = NITS + 1;
        if (NITS <= NMAX) goto L6;
        ACC = DABS(UEFF - 1.0e0);
        std::printf("\n**** SWKB: FOR L = %3d AFTER %4d ITERATIONS\n"
                    "   REQUESTED ACCURACY OF %12.4E CANNOT BE ACHIEVED "
                    "   ACCURACY ACHIEVED = %12.4E\n", L, NMAX, EPS, ACC);
        RHOT = RO3;
      }

      // EVALUATE D(UEFF)/D(RHO) AT CTP

    L11:
      X = (RHOT - RHO);
      { double DUN = -(12.0 / (H * H)) * (VCOF(N, 1) + X * (2 * VCOF(N, 2) +
                3 * X * VCOF(N, 3)));
        double DUL = 0;
        if (L != 0) DUL = -(2 * FLL) / (RHOT * RHOT * RHOT);
        DUEFF = DUN + DUL;
        double RT = RHOT / FK;
        if (IPRINT > 1)
            std::printf("\n   --- CTP FOR L = %3d FOUND IN %4d"
                        " ITERATIONS   RHOT = %13.5G RT = %13.5G\n"
                        " DUEFF =%14.5G    DUN+DUC = %13.5G DUL = %13.5G\n",
                        L, NITS, RHOT, RT, DUEFF, DUN, DUL);
        if (DUEFF < 0 || RHOT == 0) goto L160;
        std::printf("\n****** DUEFF=%15.5G ERROR IN SWKB ***\n", DUEFF);
      }

      // INTEGRATE FROM RHOT TO ROMAX TO EVALUATE WKB PHASE-SHIFT
      // AND CLASSICAL DEFLECTION FUNCTION

    L160:
      { double PHC, TEMP, PHI;
        if (ETA != 0) goto L165;
        PHC = 0;
        TEMP = DSQRT(FLL) / ROMAX;
        goto L170;
    L165:
        TEMP = ETA / ECC;
        PHC = std::acos(TEMP);
        TEMP = (ETA + FLL / ROMAX) / ECC;
    L170:
        if (TEMP > 1) TEMP = 1;
        PHI = std::acos(TEMP);

        double TURN;
        if (ROMAX > RHOT) goto L220;
        DELTR = 0;
        DELTIM = 0;
        TURN = PHC;
        goto L230_inline;

    L220:
        { double AP = (ROMAX - RHOT) / 2.0;
          double BP = (ROMAX + RHOT) / 2.0;
          double AD = ROMAX / 2.0;
          double BD = 0;
          if (RHOT == 0) goto L225;
          AD = 1.0 / RHOT - 1.0 / ROMAX;
          AD = .5 * DSQRT(AD);
          BD = 1.0 / RHOT;
    L225:
          double TEMR = 0;
          double TEMI = 0;
          double TEMP2 = 0;

          for (int I = 1; I <= NGAUSS; I++) {

              double WTP = AP * WTS[I - 1];
              double WTD = AD * WTS[I - 1];
              double RO = AP * PTS[I - 1] + BP;
              double XVAL = AD * (1 + PTS[I - 1]);
              // ASSIGN 185 TO IK / GO TO 700
              // --- inline label 700 + 710: compute N, X, UEFF from RO ---
              N = (int)((RO - ROSTRT) / H) + 1;
              X = (RO - ROSTRT - (N - 1) * H);
              UEFF = ((12.0 / (H * H)) + VCONS) - (12.0 / (H * H)) * (VARR(N)
                  + X * (VCOF(N, 1) + X * (VCOF(N, 2) + X * VCOF(N, 3))))
                  + FLL / (RO * RO);
              // --- return to label 185 ---
              double UEFFC = FLL / (RO * RO);
              double Y = RO / RHOC;
              if (RO < RHOC) UEFFC = UEFFC + UC0 * (1.5 - .5 * (Y * Y));
              if (RO >= RHOC) UEFFC = UEFFC + UC0 / Y;
              double W = -(12.0 / (H * H)) * (WARR(N) +
                  X * (WCOF(N, 1) + X * (WCOF(N, 2) + X * WCOF(N, 3))));
              double UE1 = 1 - UEFF;
              if (UE1 < 0) UE1 = 0;
              double UE2 = DSQRT(UE1 * UE1 + W * W);
              if ((UE2 - UE1) < 0) UE2 = UE1;
              double FINTI = TP * DSQRT(UE2 - UE1);
              double FINTR = TP * DSQRT(UE2 + UE1);
              if (RO <= ROTC) goto L195;
              { double UEC = 1 - UEFFC;
                if (UEC < 0) UEC = 0;
                FINTR = FINTR - DSQRT(UEC);
              }
          L195:
              if (RHOT == 0) goto L190;
              RO = BD - XVAL * XVAL;
              RO = 1.0 / RO;
              goto L197;
          L190:
              RO = .5 * AD * (1 + PTS[I - 1]) * (1 + PTS[I - 1]);
          L197:
              // ASSIGN 200 TO IK / GO TO 700
              // --- inline label 700 + 710: compute N, X, UEFF from RO ---
              N = (int)((RO - ROSTRT) / H) + 1;
              X = (RO - ROSTRT - (N - 1) * H);
              UEFF = ((12.0 / (H * H)) + VCONS) - (12.0 / (H * H)) * (VARR(N)
                  + X * (VCOF(N, 1) + X * (VCOF(N, 2) + X * VCOF(N, 3))))
                  + FLL / (RO * RO);
              // --- return to label 200 ---
              { double UER = 1 - UEFF;
                if (UER <= 0) goto L900;
                double FINTD = XVAL / DSQRT(UER);
                TEMI = TEMI + WTP * FINTI;
                TEMR = TEMR + WTP * FINTR;
                TEMP2 = TEMP2 + WTD * FINTD;
              }
              continue;

          L900:
              { double UER = 1 - UEFF;
                double RT_local = RHOT / FK;
                double XVAL_local = XVAL;
                std::printf("\n**** SWKB ERROR AT STMNT 205:%5d%5d\n"
                            "      %15.5G%15.5G%15.5G%15.5G%15.5G%15.5G\n",
                            L, I, RHOT, RT_local, XVAL_local, RO, UEFF, UER);
              }
              return;
          }

          DELTR = TEMR;
          DELTIM = TEMI;
          TURN = 2 * (L + .5) * TEMP2 + PHC - PHI;
        }

    L230_inline:
        THETA = PI - 2 * TURN;
        THETC = PI - 2 * PHC;
      }
    }

    #undef VARR
    #undef WARR
    #undef VCOF
    #undef WCOF

    return;
}


// ============================================================================
// SUBROUTINE TMATCH — lines 34575-34778
//
// CALCULATES S-MATRIX AND RENORMALIZATION FOR TENSOR-COUPLED WF.
// ============================================================================
void TMATCH(double* WAVR, double* WAVI, double* CR, double* CI,
            double* SMATS, int JP, int NWP, int NWAVF)
{
    // WAVR, WAVI: dimension (NWAVF,2,2)  — 1-based
    // CR, CI: dimension (2,2) — 1-based
    // SMATS: dimension (*) — 1-based

    // Helper macros for multi-dimensional 1-based access (Fortran column-major)
    // WAVR(i,j,k) -> WAVR[ (k-1)*NWAVF*2 + (j-1)*NWAVF + (i-1) ]
    #define WAVRm(i,j,k) WAVR[((k)-1)*NWAVF*2 + ((j)-1)*NWAVF + ((i)-1)]
    #define WAVIm(i,j,k) WAVI[((k)-1)*NWAVF*2 + ((j)-1)*NWAVF + ((i)-1)]
    // CR(i,j) -> CR[ (j-1)*2 + (i-1) ]
    #define CRm(i,j) CR[((j)-1)*2 + ((i)-1)]
    #define CIm(i,j) CI[((j)-1)*2 + ((i)-1)]

    // COMMON block aliases
    auto& PWAVSW = WAVCOM.PWAVSW;
    auto& NBAKCM = WAVCOM.NBAKCM;
    auto& STANSW = WAVCOM.STANSW;
    auto& IOFFIT = WAVCOM.IOFFIT;
    auto& NWAVEF = WAVCOM.NWAVEF;
    auto* NFS = WAVCOM.NFS;       // 1-based
    auto* NGS = WAVCOM.NGS;       // 1-based
    auto* NF1S = WAVCOM.NF1S;     // 1-based
    auto* NG1S = WAVCOM.NG1S;     // 1-based
    auto* NSTPSS = WAVCOM.NSTPSS; // 1-based
    auto* IINDXE = WAVCOM.IINDXE; // 1-based
    auto& LOMOST = KANDM.LOMOST;
    auto* Z = LOCPTRS.Z;          // 1-based
    auto& FACFR4 = ALLOCS.FACFR4;

    // Local arrays
    double F1[3], G1[3], F2[3], G2[3], ABI[3];   // 1-based
    double U1R[3][3], U1I[3][3], U2R[3][3], U2I[3][3]; // 1-based: (IL,IS)
    double BYR[3][3], BYI[3][3];  // 1-based: (IL,IS)
    const char* CHANAM[4] = { "", "INCOMING", "OUTGOING", " " }; // 1-based

    // F IS THE REGULAR COULOMB WAVE
    // G IS THE IRREGULAR COULOMB WAVE
    // NOTE THAT F1,G1 ARE REGULAR AND IRREGULAR COULOMB FUNCTIONS AT NBACK

    int I_idx = JP / 2 - 1 + (IOFFIT - 1) * (LOMOST + 1);
    int LF = (Z[NFS[NWP]] + I_idx);
    int LG = (Z[NGS[NWP]] + I_idx);
    int LF1 = (Z[NF1S[NWP]] + I_idx);
    int LG1 = (Z[NG1S[NWP]] + I_idx);
    int NSTEP = NSTPSS[NWP];

    // REDUCE THE 4-BY-4 SYSTEM TO A 2-BY-2 SYSTEM FOR THE NORMALIZATIONS
    // ONLY.  THE RIGHT-HAND SIDE (ABI) IS A DIAGONAL, PURELY IMAGINARY
    // MATRIX.

    for (int IL = 1; IL <= 2; IL++) {
        int I_off = 2 * IL - 2;
        F1[IL] = 0.5 * ALLOC(LF + I_off);
        G1[IL] = 0.5 * ALLOC(LG + I_off);
        F2[IL] = 0.5 * ALLOC(LF1 + I_off);
        G2[IL] = 0.5 * ALLOC(LG1 + I_off);
        ABI[IL] = F2[IL] * G1[IL] - F1[IL] * G2[IL];
        ABI[IL] = ABI[IL] + ABI[IL];
        for (int IS = 1; IS <= 2; IS++) {
            U1R[IL][IS] = WAVRm(NSTEP + 1, IL, IS);
            U1I[IL][IS] = WAVIm(NSTEP + 1, IL, IS);
            U2R[IL][IS] = WAVRm(NSTEP + 1 - NBAKCM, IL, IS);
            U2I[IL][IS] = WAVIm(NSTEP + 1 - NBAKCM, IL, IS);
            BYR[IL][IS] = F2[IL] * U1R[IL][IS] + G2[IL] * U1I[IL][IS]
                        - F1[IL] * U2R[IL][IS] - G1[IL] * U2I[IL][IS];
            BYI[IL][IS] = F2[IL] * U1I[IL][IS] - G2[IL] * U1R[IL][IS]
                        - F1[IL] * U2I[IL][IS] + G1[IL] * U2R[IL][IS];
        }
    }

    // SOLVE THE 2-BY-2 SYSTEM FOR THE NORMALIZATIONS (CR,CI).
    //
    // THE SYSTEM IS SOLVED BY GAUSSIAN ELIMINATION.

    double TEMP, PIVR, PIVI;

    TEMP = BYR[1][1] * BYR[1][1] + BYI[1][1] * BYI[1][1];
    PIVR = BYR[1][1] / TEMP;
    PIVI = -BYI[1][1] / TEMP;
    TEMP = PIVR * BYR[1][2] - PIVI * BYI[1][2];
    BYI[1][2] = PIVR * BYI[1][2] + PIVI * BYR[1][2];
    BYR[1][2] = TEMP;
    CRm(1, 1) = -PIVI * ABI[1];
    CIm(1, 1) = PIVR * ABI[1];
    BYR[2][2] = BYR[2][2] - BYR[2][1] * BYR[1][2] + BYI[2][1] * BYI[1][2];
    BYI[2][2] = BYI[2][2] - BYR[2][1] * BYI[1][2] - BYI[2][1] * BYR[1][2];
    CRm(2, 1) = -BYR[2][1] * CRm(1, 1) + BYI[2][1] * CIm(1, 1);
    CIm(2, 1) = -BYR[2][1] * CIm(1, 1) - BYI[2][1] * CRm(1, 1);
    TEMP = BYR[2][2] * BYR[2][2] + BYI[2][2] * BYI[2][2];
    PIVR = BYR[2][2] / TEMP;
    PIVI = -BYI[2][2] / TEMP;
    TEMP = PIVR * CRm(2, 1) - PIVI * CIm(2, 1);
    CIm(2, 1) = PIVR * CIm(2, 1) + PIVI * CRm(2, 1);
    CRm(2, 1) = TEMP;
    CRm(2, 2) = -PIVI * ABI[2];
    CIm(2, 2) = PIVR * ABI[2];
    CRm(1, 1) = CRm(1, 1) - BYR[1][2] * CRm(2, 1) + BYI[1][2] * CIm(2, 1);
    CIm(1, 1) = CIm(1, 1) - BYR[1][2] * CIm(2, 1) - BYI[1][2] * CRm(2, 1);
    CRm(1, 2) = -BYR[1][2] * CRm(2, 2) + BYI[1][2] * CIm(2, 2);
    CIm(1, 2) = -BYR[1][2] * CIm(2, 2) - BYI[1][2] * CRm(2, 2);

    // NOW SUBSTITUTE BACK TO FIND THE S-MATRIX ELEMENTS.

    int NNWP = NWP;
    if (STANSW) NNWP = 3;
    for (int IL = 1; IL <= 2; IL++) {
        int LL = (JP / 2 - 3) + 2 * IL;
        double DENA = F1[IL] * F1[IL] + G1[IL] * G1[IL];
        double DENR = F1[IL] / DENA;
        double DENI = G1[IL] / DENA;
        for (int ILAS = 1; ILAS <= 2; ILAS++) {
            int LLAS = (JP / 2 - 3) + 2 * ILAS;
            double TR = U1R[IL][1] * CRm(1, ILAS) - U1I[IL][1] * CIm(1, ILAS)
                      + U1R[IL][2] * CRm(2, ILAS) - U1I[IL][2] * CIm(2, ILAS);
            double TI = U1R[IL][1] * CIm(1, ILAS) + U1I[IL][1] * CRm(1, ILAS)
                      + U1R[IL][2] * CIm(2, ILAS) + U1I[IL][2] * CRm(2, ILAS);
            if (IL != ILAS) goto L530;
            TR = TR - F1[IL];
            TI = TI - G1[IL];
        L530:
            double SJR = TR * DENR - TI * DENI;
            double SJI = TR * DENI + TI * DENR;

            // PRINT ELASTIC S-MATRIX IF REQUESTED (1000 DIGIT OF IPRINT)

            if (!PWAVSW) goto L540;
            std::printf(" %s ELASTIC S-MATRIX FOR L =%4d,  LAS =%4d,   JP =%4d/2:%15.5G +%13.5G I\n",
                        CHANAM[NNWP], LL, LLAS, JP, SJR, SJI);

            // TRANSLATE S-MATRICES TO LX SYSTEM AND STORE THEM.

        L540:
            { int LINDXE = FACFR4 * Z[IINDXE[NWP]] - FACFR4 + 1;
              JPTOLX(LL, LLAS, JP, NWP, &SJR, &SJI,
                     &ILLOC(LINDXE), SMATS);
            }
        }
    }

    #undef WAVRm
    #undef WAVIm
    #undef CRm
    #undef CIm

    return;
}
