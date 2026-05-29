// fortlib_part1.cpp — translated from fortlib.f
// Math library: angular momentum coefficients, special functions,
// interpolation, linear algebra, Legendre polynomials, etc.
//
// Part 1: BLOCK DATA initializers, CCNFRC/CCONTF, CLEBSH/THRJ/THREEJ,
//         SIXJ/RACAH, WIG9J, DSGMAL, get_date, GAUSSL

#include "ptolemy_commons.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_forward.h"
#include "ptolemy_io.h"
#include <cstdio>
#include <cmath>
#include <cstring>
#include <complex>
#include <ctime>
#include <algorithm>

// ============================================================================
// BLOCK DATA FACDUM — sqrt(factorial) table
// COMMON /FACTRL/ MAXFAC, ISPACE, FACTBL(97)
// ============================================================================
static struct FacdumInit {
    FacdumInit() {
        FACTRL.MAXFAC = 96;
        FACTRL.ISPACE = 0;
        FACTRL.FACTBL[ 1] =  0.1000000000000000e+01;
        FACTRL.FACTBL[ 2] =  0.1000000000000000e+01;
        FACTRL.FACTBL[ 3] =  0.1414213562373095e+01;
        FACTRL.FACTBL[ 4] =  0.2449489742783178e+01;
        FACTRL.FACTBL[ 5] =  0.4898979485566356e+01;
        FACTRL.FACTBL[ 6] =  0.1095445115010332e+02;
        FACTRL.FACTBL[ 7] =  0.2683281572999748e+02;
        FACTRL.FACTBL[ 8] =  0.7099295739719539e+02;
        FACTRL.FACTBL[ 9] =  0.2007984063681781e+03;
        FACTRL.FACTBL[10] =  0.6023952191045344e+03;
        FACTRL.FACTBL[11] =  0.1904940943966505e+04;
        FACTRL.FACTBL[12] =  0.6317974358922328e+04;
        FACTRL.FACTBL[13] =  0.2188610518114176e+05;
        FACTRL.FACTBL[14] =  0.7891147445080469e+05;
        FACTRL.FACTBL[15] =  0.2952597012800765e+06;
        FACTRL.FACTBL[16] =  0.1143535905863913e+07;
        FACTRL.FACTBL[17] =  0.4574143623455652e+07;
        FACTRL.FACTBL[18] =  0.1885967730625315e+08;
        FACTRL.FACTBL[19] =  0.8001483428544984e+08;
        FACTRL.FACTBL[20] =  0.3487765766344294e+09;
        FACTRL.FACTBL[21] =  0.1559776268628498e+10;
        FACTRL.FACTBL[22] =  0.7147792818185865e+10;
        FACTRL.FACTBL[23] =  0.3352612008237171e+11;
        FACTRL.FACTBL[24] =  0.1607856235454059e+12;
        FACTRL.FACTBL[25] =  0.7876854713229383e+12;
        FACTRL.FACTBL[26] =  0.3938427356614691e+13;
        FACTRL.FACTBL[27] =  0.2008211794424596e+14;
        FACTRL.FACTBL[28] =  0.1043497458090740e+15;
        FACTRL.FACTBL[29] =  0.5521669535672285e+15;
        FACTRL.FACTBL[30] =  0.2973510046012911e+16;
        FACTRL.FACTBL[31] =  0.1628658527169496e+17;
        FACTRL.FACTBL[32] =  0.9067986906793549e+17;
        FACTRL.FACTBL[33] =  0.5129628026803635e+18;
        FACTRL.FACTBL[34] =  0.2946746955341073e+19;
        FACTRL.FACTBL[35] =  0.1718233974287565e+20;
        FACTRL.FACTBL[36] =  0.1016520927791757e+21;
        FACTRL.FACTBL[37] =  0.6099125566750542e+21;
        FACTRL.FACTBL[38] =  0.3709953246501409e+22;
        FACTRL.FACTBL[39] =  0.2286968774309350e+23;
        FACTRL.FACTBL[40] =  0.1428211541796153e+24;
        FACTRL.FACTBL[41] =  0.9032802905233224e+24;
        FACTRL.FACTBL[42] =  0.5783815921445271e+25;
        FACTRL.FACTBL[43] =  0.3748341123420972e+26;
        FACTRL.FACTBL[44] =  0.2457951648494613e+27;
        FACTRL.FACTBL[45] =  0.1630420674178431e+28;
        FACTRL.FACTBL[46] =  0.1093719437815202e+29;
        FACTRL.FACTBL[47] =  0.7417966136220958e+29;
        FACTRL.FACTBL[48] =  0.5085501366740237e+30;
        FACTRL.FACTBL[49] =  0.3523338699662023e+31;
        FACTRL.FACTBL[50] =  0.2466337089763416e+32;
        FACTRL.FACTBL[51] =  0.1743963680863606e+33;
        FACTRL.FACTBL[52] =  0.1245439180886559e+34;
        FACTRL.FACTBL[53] =  0.8980989654316716e+34;
        FACTRL.FACTBL[54] =  0.6538259159791714e+35;
        FACTRL.FACTBL[55] =  0.4804619624270389e+36;
        FACTRL.FACTBL[56] =  0.3563201278858420e+37;
        FACTRL.FACTBL[57] =  0.2666455677120592e+38;
        FACTRL.FACTBL[58] =  0.2013129889124823e+39;
        FACTRL.FACTBL[59] =  0.1533154046820762e+40;
        FACTRL.FACTBL[60] =  0.1177637968756484e+41;
        FACTRL.FACTBL[61] =  0.9121944481710788e+41;
        FACTRL.FACTBL[62] =  0.7124466393192018e+42;
        FACTRL.FACTBL[63] =  0.5609810447812647e+43;
        FACTRL.FACTBL[64] =  0.4452649004137245e+44;
        FACTRL.FACTBL[65] =  0.3562119203309796e+45;
        FACTRL.FACTBL[66] =  0.2871872314724746e+46;
        FACTRL.FACTBL[67] =  0.2333120097803461e+47;
        FACTRL.FACTBL[68] =  0.1909741105966688e+48;
        FACTRL.FACTBL[69] =  0.1574812859496909e+49;
        FACTRL.FACTBL[70] =  0.1308137807832727e+50;
        FACTRL.FACTBL[71] =  0.1094466613011557e+51;
        FACTRL.FACTBL[72] =  0.9222139602976428e+51;
        FACTRL.FACTBL[73] =  0.7825244940376377e+52;
        FACTRL.FACTBL[74] =  0.6685892207860282e+53;
        FACTRL.FACTBL[75] =  0.5751421947239992e+54;
        FACTRL.FACTBL[76] =  0.4980877514193197e+55;
        FACTRL.FACTBL[77] =  0.4342228346904444e+56;
        FACTRL.FACTBL[78] =  0.3810289910601106e+57;
        FACTRL.FACTBL[79] =  0.3365156932181068e+58;
        FACTRL.FACTBL[80] =  0.2991016905800262e+59;
        FACTRL.FACTBL[81] =  0.2675246849288189e+60;
        FACTRL.FACTBL[82] =  0.2407722164359370e+61;
        FACTRL.FACTBL[83] =  0.2180285150390389e+62;
        FACTRL.FACTBL[84] =  0.1986334304622628e+63;
        FACTRL.FACTBL[85] =  0.1820505461284133e+64;
        FACTRL.FACTBL[86] =  0.1678423103505356e+65;
        FACTRL.FACTBL[87] =  0.1556505553593457e+66;
        FACTRL.FACTBL[88] =  0.1451811729660402e+67;
        FACTRL.FACTBL[89] =  0.1361920123419132e+68;
        FACTRL.FACTBL[90] =  0.1284832874770429e+69;
        FACTRL.FACTBL[91] =  0.1218899489080934e+70;
        FACTRL.FACTBL[92] =  0.1162756005221389e+71;
        FACTRL.FACTBL[93] =  0.1115276380752381e+72;
        FACTRL.FACTBL[94] =  0.1075533591796017e+73;
        FACTRL.FACTBL[95] =  0.1042768505784838e+74;
        FACTRL.FACTBL[96] =  0.1016365017512855e+75;
        FACTRL.FACTBL[97] =  0.9958302741285533e+75;
    }
} facdum_init_;

// ============================================================================
// BLOCK DATA LOGDUM — log(factorial) table
// COMMON /LOGFAC/ MAXLF, lfbias, LF(1001)
// ============================================================================
static struct LogdumInit {
    LogdumInit() {
        LOGFAC.MAXLF = 0;
        LOGFAC.LFBIAS = 0;
        LOGFAC.LF[1] = 0.0;
    }
} logdum_init_;

// ============================================================================
// BLOCK DATA ROOTIS — sqrt(integer) table
// COMMON /ROOTIS/ IMXDIM, ISPACE, ROOTI(602)
// Compute at static init: ROOTI[k] = sqrt(k-1) for k=1..602
// ============================================================================
static struct RootisInit {
    RootisInit() {
        ROOTIS.IMXDIM = 602;
        ROOTIS.ISPACE = 0;
        for (int k = 1; k <= 602; k++) {
            ROOTIS.ROOTI[k] = std::sqrt(static_cast<double>(k - 1));
        }
    }
} rootis_init_;

// ============================================================================
// CCNFRC / CCONTF — complex continued fraction interpolation
// Fortran lines 1-174
// ============================================================================

// File-static saved variable (Fortran SAVE NMAX)
static int ccnfrc_nmax = 0;

void CCNFRC(int NUMPTS, complex16* XS, complex16* YS) {
    // Fortran uses 1-based arrays; C++ uses 0-based.
    // All indices shifted by -1 from the original Fortran.

    complex16 D, XJ, YJ;
    double TINY = 1.0e-14;
    double COMP, TEMP;
    int K, I, J;

    ccnfrc_nmax = NUMPTS - 1;

    COMP = -1.0;
    for (I = 0; I < NUMPTS; I++) {
        D = YS[I];
        TEMP = std::norm(D);
        if (TEMP <= COMP) continue;
        COMP = TEMP;
        K = I;
    }
    if (ccnfrc_nmax <= 0) return;

    YJ = YS[K];
    if (K != 0) {
        YS[K] = YS[0];
        YS[0] = YJ;
        XJ = XS[K];
        XS[K] = XS[0];
        XS[0] = XJ;
    }

    COMP = -1.0;
    for (I = 1; I < NUMPTS; I++) {
        D = 1.0 - YJ / YS[I];
        YS[I] = D;
        TEMP = std::norm(D);
        if (TEMP <= COMP) continue;
        COMP = TEMP;
        K = I;
    }

    for (J = 1; J < NUMPTS; J++) {
        XJ = XS[K];
        YJ = YS[K];
        if (K != J) {
            YS[K] = YS[J];
            XS[K] = XS[J];
            XS[J] = XJ;
        }
        YJ = YJ / (XS[J-1] - XJ);
        YS[J] = YJ;
        if (J == NUMPTS - 1) goto label_450;

        if (COMP < TINY) goto label_430;

        {
            int JP1 = J + 1;
            COMP = -1.0;
            XJ = XS[J-1];
            for (I = JP1; I < NUMPTS; I++) {
                D = 1.0 + YJ * (XS[I] - XJ) / YS[I];
                YS[I] = D;
                TEMP = std::norm(D);
                if (TEMP <= COMP) continue;
                COMP = TEMP;
                K = I;
            }
        }
    }

    return;

label_430:
    ccnfrc_nmax = J;
    std::printf("\n**** WARNING:  CONTINUED FRACTION USED ONLY%3d OUT OF%3d POINTS.\n",
               ccnfrc_nmax, NUMPTS);

label_450:
    return;
}

void CCONTF(int NUMPTS, complex16* XS, complex16* YS, complex16 X, complex16& Y) {
    // 0-based indexing (shifted from Fortran 1-based)
    Y = complex16(0.0, 0.0);
    if (ccnfrc_nmax < 1) goto label_520;
    for (int J = 1; J <= ccnfrc_nmax; J++) {
        int K = ccnfrc_nmax - J;
        Y = YS[K+1] * (X - XS[K]) / (1.0 + Y);
    }
label_520:
    Y = YS[0] / (1.0 + Y);
    return;
}


// ============================================================================
// Begin CLEBSH/THRJ/THREEJ and SIXJ/RACAH translations
// ============================================================================


// ============================================================================
// CLEBSH / THRJ / THREEJ
// ============================================================================
//
// IMPLICIT INTEGER*4 (A-C, I-Z), REAL*8 (D-H)
//   Variables A-C, I-Z are int; D-H are double.
//   SUM, ANSWER are explicitly REAL*8.
//
// ENTRY points:
//   CLEBSH(A,B,X,Y,CIN,ZIN)  — Clebsch-Gordan coefficient
//   THRJ(A,B,CIN,X,Y,ZIN)    — 3-J symbol (Z = -ZIN)
//   THREEJ same as THRJ
//
// Handled via static helper with THRESW flag.
// DBLSW is always true (double precision only).

static double clebsh_impl(int A, int B, int X, int Y, int CIN, int ZIN,
                           bool THRESW)
{
    // ---- All local variable declarations before any gotos ----

    // IMPLICIT INTEGER*4 (A-C, I-Z)
    int C, Z;
    int X1, Y1, Z1, T30;
    int X2, Y2, Z2;
    int R;
    int A1, A2, A3;
    int B1, B2;
    int NMIN, NMAX, NMAXM2;
    int N;
    int JSUM;
    int TEMP;
    int I;

    // IMPLICIT REAL*8 (D-H)
    double DB1, DB2;
    double DX, DN, DC;
    double E, E1, E2, E3, E4, E5, E6;
    double F1, F2, F3, F4;
    double G1, G4;
    double H1, H4;
    double DJSUM;

    // REAL*8 (explicit)
    double SUM, ANSWER;

    // DIMENSION DRAY(9) with EQUIVALENCE:
    //   DRAY(1)=DA1, DRAY(2)=DA2, DRAY(3)=DA3,
    //   DRAY(4)=DY1, DRAY(5)=DX2, DRAY(6)=DY2,
    //   DRAY(7)=DX1, DRAY(8)=DZ1, DRAY(9)=DZ2
    double DRAY[10]; // 1-based: DRAY[1]..DRAY[9]
    double& DA1 = DRAY[1];
    double& DA2 = DRAY[2];
    double& DA3 = DRAY[3];
    double& DY1 = DRAY[4];
    double& DX2 = DRAY[5];
    double& DY2 = DRAY[6];
    double& DX1 = DRAY[7];
    double& DZ1 = DRAY[8];
    double& DZ2 = DRAY[9];

    // ---- Begin translated code ----

    // CLEBSH entry: THRESW = .FALSE., DBLSW = .TRUE., GO TO 10
    // THRJ/THREEJ entry: Z = -ZIN, THRESW = .TRUE., DBLSW = .TRUE., GO TO 20

    if (THRESW) {
        // THRJ / THREEJ entry
        Z = -ZIN;
        goto L20;
    }

    // CLEBSH entry
    // GO TO 10

// 10   Z = ZIN
L10:
    Z = ZIN;

// 20   C = CIN
L20:
    C = CIN;

//      IF (Z .EQ. X+Y)  GO TO 120
    if (Z == X+Y) goto L120;

// 110  CLEBSH=0.0
//      RETURN
L110:
    return 0.0;

// 120  X1 = ISHFT(A+X, -1)
L120:
    X1 = ISHFT(A+X, -1);
    Y1 = ISHFT(B+Y, -1);
    Z1 = ISHFT(C+Z, -1);

//      T30 = 30
    T30 = 30;

//      IF (BTEST(X1,T30).OR.BTEST(Y1,T30).OR.BTEST(Z1,T30)) GO TO 110
    if (BTEST(X1,T30) || BTEST(Y1,T30) || BTEST(Z1,T30)) goto L110;

//      X2 = X1 - X
    X2 = X1 - X;
    Y2 = Y1 - Y;
    Z2 = Z1 - Z;

//      R = X1 + Y1 + Z2 + 1
    R = X1 + Y1 + Z2 + 1;

//      A1 = R-1 - A
    A1 = R-1 - A;
    A2 = R-1 - B;
    A3 = R-1 - C;

//      B1 = A1 - Y1
    B1 = A1 - Y1;
    B2 = A2 - X2;

//      DA1 = A1
    DA1 = A1;
    DA2 = A2;
    DA3 = A3;
    DY1 = Y1;
    DX2 = X2;

//      DB1 = DA1 - DY1
    DB1 = DA1 - DY1;
    DB2 = DA2 - DX2;

//      SUM = 1
    SUM = 1;

//      NMIN = MAX0 (0, -B1, -B2)
    NMIN = MAX0(0, -B1, -B2);
//      NMAX = MIN0 (A3, X2, Y1)
    NMAX = MIN0(A3, X2, Y1);

//      IF (NMAX - NMIN)  110, 300, 200
    if (NMAX - NMIN < 0) goto L110;
    if (NMAX - NMIN == 0) goto L300;
    // fall through to 200

// 200  IF (IABS(X)+IABS(Y) .EQ. 0)  GO TO 500
L200:
    if (IABS(X)+IABS(Y) == 0) goto L500;

//      DX = NMAX
    DX = NMAX;
    E1 = -DX;
    E2 = -DB1-DX;
    E3 = -DB2 - DX;
    E4 = DA3 - (DX-1);
    E5 = DX2 - (DX-1);
    E6 = DY1 - (DX-1);

//      F1 = E1*E2*E3
    F1 = E1*E2*E3;
    F4 = E4*E5*E6;

//      NMAXM2 = NMAX - 2
    NMAXM2 = NMAX - 2;
    SUM = 1 + F4/F1;
    if (NMAXM2 < NMIN) goto L300;

//      G1 = (E1+1)*(E2+1)*(E3+1) - F1
    G1 = (E1+1)*(E2+1)*(E3+1) - F1;
    G4 = (E4+1)*(E5+1)*(E6+1) - F4;
    H1 = (E1+2)*(E2+2)*(E3+2) - G1 - G1 - F1;
    H4 = (E4+2)*(E5+2)*(E6+2) - G4 - G4 - F4;

//      DO  259  N = NMIN, NMAXM2
    for (N = NMIN; N <= NMAXM2; N++) {
        F1 = F1 + G1;
        F4 = F4 + G4;
        G1 = G1 + H1;
        G4 = G4 + H4;
        H1 = H1 + 6;
        H4 = H4 + 6;
        SUM = 1 + SUM*(F4/F1);
    } // 259

// 300  IF (R .LE. MAXFAC)  GO TO 400
L300:
    if (R <= FACTRL.MAXFAC) goto L400;

//      IF (R .GT. MAXLF)  GO TO 350
    if (R > LOGFAC.MAXLF) goto L350;

//      E = lf(lfbias+1+A3) - lf(lfbias+1+NMIN) - lf(lfbias+1+A3-NMIN)
//        + lf(lfbias+1+A1) - lf(lfbias+1+B1+NMIN) - lf(lfbias+1+Y1-NMIN)
//        + lf(lfbias+1+A2) - lf(lfbias+1+B2+NMIN) - lf(lfbias+1+X2-NMIN)
    E = LOGFAC.LF[LOGFAC.LFBIAS+1+A3] - LOGFAC.LF[LOGFAC.LFBIAS+1+NMIN]
          - LOGFAC.LF[LOGFAC.LFBIAS+1+A3-NMIN]
      + LOGFAC.LF[LOGFAC.LFBIAS+1+A1] - LOGFAC.LF[LOGFAC.LFBIAS+1+B1+NMIN]
          - LOGFAC.LF[LOGFAC.LFBIAS+1+Y1-NMIN]
      + LOGFAC.LF[LOGFAC.LFBIAS+1+A2] - LOGFAC.LF[LOGFAC.LFBIAS+1+B2+NMIN]
          - LOGFAC.LF[LOGFAC.LFBIAS+1+X2-NMIN];

//      E2 = .5*(-lf(lfbias+1+A1) - lf(lfbias+1+A2) - lf(lfbias+1+A3)
//         + lf(lfbias+1+X1) + lf(lfbias+1+X2)
//         + lf(lfbias+1+Y1) + lf(lfbias+1+Y2) + lf(lfbias+1+Z1)
//         + lf(lfbias+1+Z2) - lf(lfbias+1+R)
//         + lf(lfbias+2+C) - lf(lfbias+1+C) )
    E2 = .5*(-LOGFAC.LF[LOGFAC.LFBIAS+1+A1] - LOGFAC.LF[LOGFAC.LFBIAS+1+A2]
             - LOGFAC.LF[LOGFAC.LFBIAS+1+A3]
       + LOGFAC.LF[LOGFAC.LFBIAS+1+X1] + LOGFAC.LF[LOGFAC.LFBIAS+1+X2]
       + LOGFAC.LF[LOGFAC.LFBIAS+1+Y1] + LOGFAC.LF[LOGFAC.LFBIAS+1+Y2]
             + LOGFAC.LF[LOGFAC.LFBIAS+1+Z1]
       + LOGFAC.LF[LOGFAC.LFBIAS+1+Z2] - LOGFAC.LF[LOGFAC.LFBIAS+1+R]
       + LOGFAC.LF[LOGFAC.LFBIAS+2+C] - LOGFAC.LF[LOGFAC.LFBIAS+1+C] );

//      IF (DABS(E) .GT. 40)  GO TO 330
    if (DABS(E) > 40) goto L330;

//      SUM = ANINT(SUM*DEXP(E))
    SUM = ANINT(SUM*DEXP(E));

//      ANSWER = DEXP(E2)*SUM
    ANSWER = DEXP(E2)*SUM;
    goto L450;

// 330  E = E + E2
L330:
    E = E + E2;
    ANSWER = DEXP(E)*SUM;
    goto L450;

// 350  (large J — use DLGAMA)
L350:
    DX1 = DA2 + DA3 - DX2;
    DC = C;
    DN = NMIN;
    DZ1 = DX1 + DY1 - DA3;
    DY2 = DA1 + DA3 - DY1;
    DZ2 = DC - DZ1;

//      E = DLOG(DC+1) - DLGAMA(2+DX1+DY1+DZ2)
    E = DLOG(DC+1) - DLGAMA(2+DX1+DY1+DZ2);
//      DO 369  I = 1, 9
    for (I = 1; I <= 9; I++) {
        E = E + DLGAMA(1+DRAY[I]);
    } // 369
//      E = .5*E
    E = .5*E;

//      DA1 = DB1 + (DN+1)
    DA1 = DB1 + (DN+1);
//      DA2 = DB2 + (DN+1)
    DA2 = DB2 + (DN+1);
//      DA3 = DA3 - (DN-1)
    DA3 = DA3 - (DN-1);
//      DY1 = DY1 - (DN-1)
    DY1 = DY1 - (DN-1);
//      DX2 = DX2 - (DN-1)
    DX2 = DX2 - (DN-1);
//      DY2 = (DN+1)
    DY2 = (DN+1);

//      DO  379  I = 1, 6
    for (I = 1; I <= 6; I++) {
        E = E - DLGAMA(DRAY[I]);
    } // 379

//      ANSWER = DEXP(E) * SUM
    ANSWER = DEXP(E) * SUM;

    goto L450;

// 400  (factorial tables — small J case)
L400:
    F1 = FACTRL.FACTBL[1+A1];
    F2 = FACTRL.FACTBL[1+A2];
    F3 = FACTRL.FACTBL[1+A3];

//      E = (F3/(FACTBL(1+NMIN)*FACTBL(1+A3-NMIN)))
//        * (F2/(FACTBL(1+X2-NMIN)*FACTBL(1+B2+NMIN)))
//        * (F1/(FACTBL(1+Y1-NMIN)*FACTBL(1+B1+NMIN)))
    E = (F3/(FACTRL.FACTBL[1+NMIN]*FACTRL.FACTBL[1+A3-NMIN]))
      * (F2/(FACTRL.FACTBL[1+X2-NMIN]*FACTRL.FACTBL[1+B2+NMIN]))
      * (F1/(FACTRL.FACTBL[1+Y1-NMIN]*FACTRL.FACTBL[1+B1+NMIN]));
//      SUM = ANINT(E*SUM*E)
    SUM = ANINT(E*SUM*E);

//      ANSWER = SUM * (FACTBL(1+X1)*FACTBL(1+X2)/(F3*F2))
//        * (FACTBL(1+Y1)*FACTBL(1+Y2)/F1)
//        * (FACTBL(2+C)/FACTBL(1+C))
//        * (FACTBL(1+Z1)*FACTBL(1+Z2)/FACTBL(1+R))
    ANSWER = SUM * (FACTRL.FACTBL[1+X1]*FACTRL.FACTBL[1+X2]/(F3*F2))
      * (FACTRL.FACTBL[1+Y1]*FACTRL.FACTBL[1+Y2]/F1)
      * (FACTRL.FACTBL[2+C]/FACTRL.FACTBL[1+C])
      * (FACTRL.FACTBL[1+Z1]*FACTRL.FACTBL[1+Z2]/FACTRL.FACTBL[1+R]);

// 450  IF (BTEST(NMIN,0))  ANSWER = -ANSWER
L450:
    if (BTEST(NMIN,0)) ANSWER = -ANSWER;
    goto L800;

// 500  (special case: M1 = M2 = M3 = 0)
//      RESULT IS 0 IF J1+J2+J3 IS ODD
L500:
    if (!BTEST(R,0)) goto L110;
    JSUM = R-1;
    if (R > FACTRL.MAXFAC) goto L550;

//      ANSWER = ((FACTBL(1+A1)*FACTBL(1+A2))/FACTBL(1+R))
//        * ((FACTBL(2+C)*FACTBL(1+A3))/FACTBL(1+C))
//        * (FACTBL(1+JSUM/2)/
//           (FACTBL(1+A1/2)*FACTBL(1+A2/2)*FACTBL(1+A3/2)) )**2
    ANSWER = ((FACTRL.FACTBL[1+A1]*FACTRL.FACTBL[1+A2])/FACTRL.FACTBL[1+R])
      * ((FACTRL.FACTBL[2+C]*FACTRL.FACTBL[1+A3])/FACTRL.FACTBL[1+C])
      * (FACTRL.FACTBL[1+JSUM/2]
         / (FACTRL.FACTBL[1+A1/2]*FACTRL.FACTBL[1+A2/2]*FACTRL.FACTBL[1+A3/2]))
      * (FACTRL.FACTBL[1+JSUM/2]
         / (FACTRL.FACTBL[1+A1/2]*FACTRL.FACTBL[1+A2/2]*FACTRL.FACTBL[1+A3/2]));
    goto L570;

// 550
L550:
    if (R > LOGFAC.MAXLF) goto L560;

//      E = .5 * (lf(lfbias+1+A1) + lf(lfbias+1+A2) + lf(lfbias+1+A3)
//          + lf(lfbias+2+C) - lf(lfbias+1+C) - lf(lfbias+1+R) )
//        - lf(lfbias+1+ISHFT(A1,-1)) - lf(lfbias+1+ISHFT(A2,-1))
//        - lf(lfbias+1+ISHFT(A3,-1)) + lf(lfbias+1+ISHFT(JSUM,-1))
    E = .5 * (LOGFAC.LF[LOGFAC.LFBIAS+1+A1] + LOGFAC.LF[LOGFAC.LFBIAS+1+A2]
              + LOGFAC.LF[LOGFAC.LFBIAS+1+A3]
        + LOGFAC.LF[LOGFAC.LFBIAS+2+C] - LOGFAC.LF[LOGFAC.LFBIAS+1+C]
              - LOGFAC.LF[LOGFAC.LFBIAS+1+R] )
      - LOGFAC.LF[LOGFAC.LFBIAS+1+ISHFT(A1,-1)]
              - LOGFAC.LF[LOGFAC.LFBIAS+1+ISHFT(A2,-1)]
      - LOGFAC.LF[LOGFAC.LFBIAS+1+ISHFT(A3,-1)]
              + LOGFAC.LF[LOGFAC.LFBIAS+1+ISHFT(JSUM,-1)];

    goto L565;

// 560
L560:
    DJSUM = DA1 + DA2 + DA3;
    E = DLOG(DFLOAT(C+1)) - DLGAMA(DJSUM+2);
    E2 = -DLGAMA(.5*DJSUM+1);
//      DO 564  I = 1, 3
    for (I = 1; I <= 3; I++) {
        E = E + DLGAMA(1+DRAY[I]);
        E2 = E2 + DLGAMA(1 + .5*DRAY[I]);
    } // 564
//      E = .5*E - E2
    E = .5*E - E2;

// 565  ANSWER = DEXP(E)
L565:
    ANSWER = DEXP(E);

// 570  IF (BTEST(A3,1))  ANSWER = -ANSWER
L570:
    if (BTEST(A3,1)) ANSWER = -ANSWER;

// 800  IF (THRESW) GO TO 950
L800:
    if (THRESW) goto L950;

// 810  (DBLSW always true => label 850)
L810:
// 850  CLEBSH = ANSWER
//      RETURN
    return ANSWER;

// 950  (3-J factor)
L950:
    TEMP = A2-Z2;
    if (BTEST(TEMP, 0)) ANSWER = -ANSWER;
    if (C > FACTRL.MAXFAC-1) goto L970;
    ANSWER = ANSWER * (FACTRL.FACTBL[1+C]/FACTRL.FACTBL[2+C]);
    goto L810;
L970:
    ANSWER = ANSWER / DSQRT(C+1.0);
    goto L810;
}

// ---- Public entry points ----

// CLEBSH(A,B,X,Y,CIN,ZIN)
double CLEBSH(int A, int B, int X, int Y, int CIN, int ZIN)
{
    return clebsh_impl(A, B, X, Y, CIN, ZIN, false);
}

// THRJ(A,B,CIN,X,Y,ZIN)
// Fortran ENTRY THRJ(A,B,CIN,X,Y,ZIN) reorders the args vs CLEBSH:
//   THRJ's 3rd arg is CIN (= J3), 4th is X (= M1), 5th is Y (= M2), 6th is ZIN (= M3)
// The impl expects (A, B, X, Y, CIN, ZIN, thresw) matching CLEBSH's arg order.
// The impl's THRESW path does Z = -ZIN internally, so we pass ZIN as-is.
double THRJ(int A, int B, int CIN, int X, int Y, int ZIN)
{
    return clebsh_impl(A, B, X, Y, CIN, ZIN, true);
}

// THREEJ — same signature and behavior as THRJ
double THREEJ(int A, int B, int CIN, int X, int Y, int ZIN)
{
    return clebsh_impl(A, B, X, Y, CIN, ZIN, true);
}


// ============================================================================
// SIXJ / RACAH
// ============================================================================
//
// IMPLICIT INTEGER*4 (A-C, I-Z), REAL*8 (D-H)
//   Variables A-C, I-Z are int; D-H are double.
//   SUM, ANSWER, SIXJ, RACAH are explicitly REAL*8.
//
// ENTRY points:
//   SIXJ(A,B,C,X,Y,Z)    — 6-J coefficient
//   RACAH(A,B,Y,X,C,Z)   — Racah W coefficient
//     Fortran: ENTRY RACAH(A,B,Y,X,C,Z)
//     The parameter names in the ENTRY match the FUNCTION's names,
//     so RACAH's 1st arg -> A, 2nd -> B, 3rd -> Y, 4th -> X, 5th -> C, 6th -> Z
//
// DBLSW is always true (double precision only).

static double sixj_impl(int A, int B, int C, int X, int Y, int Z,
                          bool RACHSW)
{
    // ---- All local variable declarations before any gotos ----

    // IMPLICIT INTEGER*4 (A-C, I-Z)
    int A1, A2, A3, A4;
    int B1, B2;
    int R;
    int NMIN, NMAX, NMAXM2;
    int N;
    int TEMP;
    int I;

    // IMPLICIT REAL*8 (D-H)
    double DB1, DB2;
    double DR, DN, DX, DY;
    double D1, D2, D3, D4, D5, D6, D7, D8;
    double E, E1, E2, E5;
    double F1, F5;
    double G1, G5;
    double H1, H5;
    double DENOM;

    // REAL*8 (explicit)
    double SUM, ANSWER;

    // DIMENSION DRAY(8) with EQUIVALENCE:
    //   DRAY(1)=DA1, DRAY(2)=DA2, DRAY(3)=DA3, DRAY(4)=DA4,
    //   DRAY(5)=DC,  DRAY(7)=DZ
    //   DRAY(6) and DRAY(8) are not named but used directly
    double DRAY[9]; // 1-based: DRAY[1]..DRAY[8]
    double& DA1 = DRAY[1];
    double& DA2 = DRAY[2];
    double& DA3 = DRAY[3];
    double& DA4 = DRAY[4];
    double& DC  = DRAY[5];
    double& DZ  = DRAY[7];
    // DRAY[6] and DRAY[8] are used directly (not named)

    // ---- Begin translated code ----

    // SIXJ entry: DBLSW = .TRUE., RACHSW = .FALSE., GO TO 100
    // RACAH entry: RACHSW = .TRUE., DBLSW = .TRUE., fall through to 100

// 100  SUM=1
    SUM = 1;

//      A2 = ISHFT(X+Y-C, -1)
    A2 = ISHFT(X+Y-C, -1);
    A3 = ISHFT(A+Y-Z, -1);
    A4 = ISHFT(B+X-Z, -1);
    A1 = ISHFT(A+B-C, -1);

//      IF (BTEST(A1,30) .OR. BTEST(A2,30) .OR.
//          BTEST(A3,30) .OR. BTEST(A4,30)) GO TO 900
    if (BTEST(A1,30) || BTEST(A2,30) ||
        BTEST(A3,30) || BTEST(A4,30)) goto L900;

//      B1 = Y - A2 - A3
    B1 = Y - A2 - A3;
//      B2 = A3 - A1 + Z - Y
    B2 = A3 - A1 + Z - Y;

//      R = A1 + A2 + C + 1
    R = A1 + A2 + C + 1;

//      NMAX = MIN0 (A1, A2, A3, A4)
    NMAX = std::min({A1, A2, A3, A4});
//      NMIN = MAX0 (0, -B1, -B2)
    NMIN = MAX0(0, -B1, -B2);

//      DA1 = A1
    DA1 = A1;
    DA2 = A2;
    DA3 = A3;
    DA4 = A4;
    DY = Y;
    DZ = Z;
    DC = C;
    DB1 = DY - DA2 - DA3;
    DB2 = DA3 - DA1 + DZ - DY;
    DR = DA1 + DA2 + DC + 1;
    DN = NMIN;
    DX = NMAX;

//      IF (NMAX - NMIN)  900, 300, 200
    if (NMAX - NMIN < 0) goto L900;
    if (NMAX - NMIN == 0) goto L300;
    // fall through to 200

// 200
    D1 = DA1 - (DX-1);
    D2 = DA2 - (DX-1);
    D3 = DA3 - (DX-1);
    D4 = DA4 - (DX-1);
    D5 = -DX;
    D6 = -DB1 - DX;
    D7 = -DB2 - DX;
    D8 = DR - (DX-1);

//      E1 = D1*D2*D3*D4
    E1 = D1*D2*D3*D4;
    E5 = D5*D6*D7*D8;

//      NMAXM2 = NMAX-2
    NMAXM2 = NMAX-2;
    SUM = 1 + E1/E5;
    if (NMAXM2 < NMIN) goto L300;

//      F1 = (D1+1)*(D2+1)*(D3+1)*(D4+1) - E1
    F1 = (D1+1)*(D2+1)*(D3+1)*(D4+1) - E1;
    F5 = (D5+1)*(D6+1)*(D7+1)*(D8+1) - E5;

//      G1 = (D1+2)*(D2+2)*(D3+2)*(D4+2) - F1 - F1 - E1
    G1 = (D1+2)*(D2+2)*(D3+2)*(D4+2) - F1 - F1 - E1;
    G5 = (D5+2)*(D6+2)*(D7+2)*(D8+2) - F5 - F5 - E5;

//      H1 = (D1+3)*(D2+3)*(D3+3)*(D4+3) - 3*(G1+F1) - E1
    H1 = (D1+3)*(D2+3)*(D3+3)*(D4+3) - 3*(G1+F1) - E1;
    H5 = (D5+3)*(D6+3)*(D7+3)*(D8+3) - 3*(G5+F5) - E5;

//      DO 259  N = NMIN, NMAXM2
    for (N = NMIN; N <= NMAXM2; N++) {
        E1 = E1 + F1;
        E5 = E5 + F5;
        F1 = F1 + G1;
        F5 = F5 + G5;
        G1 = G1 + H1;
        G5 = G5 + H5;
        H1 = H1 + 24;
        H5 = H5 + 24;
        SUM = 1 + SUM*(E1/E5);
    } // 259

// 300  IF (R .GT. MAXFAC)  GO TO 400
L300:
    if (R > FACTRL.MAXFAC) goto L400;

//      E = FACTBL(1+R-NMIN) / (FACTBL(1+NMIN) * ...)
    E = FACTRL.FACTBL[1+R-NMIN] / (FACTRL.FACTBL[1+NMIN]
      * FACTRL.FACTBL[1+B1+NMIN] * FACTRL.FACTBL[1+B2+NMIN]
      * FACTRL.FACTBL[1+A1-NMIN]
      * FACTRL.FACTBL[1+A2-NMIN] * FACTRL.FACTBL[1+A3-NMIN]
      * FACTRL.FACTBL[1+A4-NMIN]);

//      SUM = ANINT(E*SUM*E)
    SUM = ANINT(E*SUM*E);

//      DENOM = (FACTBL(2+A1+C)/(FACTBL(1+A1)*FACTBL(1+A4+B1)
//                              *FACTBL(1+A3+B2)))
//        * (FACTBL(2+A2+C)/(FACTBL(1+A2)*FACTBL(1+A3+B1)
//                              *FACTBL(1+A4+B2)))
//        * (FACTBL(2+A3+Z)/(FACTBL(1+A3)*FACTBL(1+A2+B1)
//                              *FACTBL(1+A1+B2)))
//        * (FACTBL(2+A4+Z)/(FACTBL(1+A4)*FACTBL(1+A1+B1)
//                              *FACTBL(1+A2+B2)))
    DENOM = (FACTRL.FACTBL[2+A1+C]/(FACTRL.FACTBL[1+A1]*FACTRL.FACTBL[1+A4+B1]
                                    *FACTRL.FACTBL[1+A3+B2]))
      * (FACTRL.FACTBL[2+A2+C]/(FACTRL.FACTBL[1+A2]*FACTRL.FACTBL[1+A3+B1]
                                    *FACTRL.FACTBL[1+A4+B2]))
      * (FACTRL.FACTBL[2+A3+Z]/(FACTRL.FACTBL[1+A3]*FACTRL.FACTBL[1+A2+B1]
                                    *FACTRL.FACTBL[1+A1+B2]))
      * (FACTRL.FACTBL[2+A4+Z]/(FACTRL.FACTBL[1+A4]*FACTRL.FACTBL[1+A1+B1]
                                    *FACTRL.FACTBL[1+A2+B2]));
//      ANSWER = SUM/(DENOM)
    ANSWER = SUM/(DENOM);
    goto L600;

// 400  IF (R .GT. MAXLF)  GO TO 500
L400:
    if (R > LOGFAC.MAXLF) goto L500;

//      E = lf(lfbias+1+R-NMIN) - lf(lfbias+1+NMIN)
//        - lf(lfbias+1+B1+NMIN) - lf(lfbias+1+B2+NMIN) - lf(lfbias+1+A1-NMIN)
//        - lf(lfbias+1+A2-NMIN) - lf(lfbias+1+A3-NMIN) - lf(lfbias+1+A4-NMIN)
    E = LOGFAC.LF[LOGFAC.LFBIAS+1+R-NMIN] - LOGFAC.LF[LOGFAC.LFBIAS+1+NMIN]
      - LOGFAC.LF[LOGFAC.LFBIAS+1+B1+NMIN] - LOGFAC.LF[LOGFAC.LFBIAS+1+B2+NMIN]
      - LOGFAC.LF[LOGFAC.LFBIAS+1+A1-NMIN]
      - LOGFAC.LF[LOGFAC.LFBIAS+1+A2-NMIN] - LOGFAC.LF[LOGFAC.LFBIAS+1+A3-NMIN]
      - LOGFAC.LF[LOGFAC.LFBIAS+1+A4-NMIN];

//      E2 = .5 * ( lf(lfbias+1+A1) + lf(lfbias+1+A4+B1) + lf(lfbias+1+A3+B2)
//          - lf(lfbias+2+A1+C)
//        + lf(lfbias+1+A2) + lf(lfbias+1+A3+B1) + lf(lfbias+1+A4+B2)
//          - lf(lfbias+2+A2+C)
//        + lf(lfbias+1+A3) + lf(lfbias+1+A2+B1) + lf(lfbias+1+A1+B2)
//          - lf(lfbias+2+A3+Z)
//        + lf(lfbias+1+A4) + lf(lfbias+1+A1+B1) + lf(lfbias+1+A2+B2)
//          - lf(lfbias+2+A4+Z) )
    E2 = .5 * ( LOGFAC.LF[LOGFAC.LFBIAS+1+A1] + LOGFAC.LF[LOGFAC.LFBIAS+1+A4+B1]
      + LOGFAC.LF[LOGFAC.LFBIAS+1+A3+B2]
        - LOGFAC.LF[LOGFAC.LFBIAS+2+A1+C]
      + LOGFAC.LF[LOGFAC.LFBIAS+1+A2] + LOGFAC.LF[LOGFAC.LFBIAS+1+A3+B1]
      + LOGFAC.LF[LOGFAC.LFBIAS+1+A4+B2]
        - LOGFAC.LF[LOGFAC.LFBIAS+2+A2+C]
      + LOGFAC.LF[LOGFAC.LFBIAS+1+A3] + LOGFAC.LF[LOGFAC.LFBIAS+1+A2+B1]
      + LOGFAC.LF[LOGFAC.LFBIAS+1+A1+B2]
        - LOGFAC.LF[LOGFAC.LFBIAS+2+A3+Z]
      + LOGFAC.LF[LOGFAC.LFBIAS+1+A4] + LOGFAC.LF[LOGFAC.LFBIAS+1+A1+B1]
      + LOGFAC.LF[LOGFAC.LFBIAS+1+A2+B2]
        - LOGFAC.LF[LOGFAC.LFBIAS+2+A4+Z] );

// 450  ANSWER = DEXP(E+E2) * SUM
    ANSWER = DEXP(E+E2) * SUM;
    goto L600;

// 500  (large J — use DLGAMA)
L500:
//      DRAY(6) = DC
    DRAY[6] = DC;
//      DRAY(8) = DZ
    DRAY[8] = DZ;
    E2 = 0;
//      E = DLGAMA(1-DN+DR) - DLGAMA(1+DN)
//        - DLGAMA(1+DN+DB1) - DLGAMA(1+DN+DB2)
    E = DLGAMA(1-DN+DR) - DLGAMA(1+DN)
      - DLGAMA(1+DN+DB1) - DLGAMA(1+DN+DB2);
//      DO 549  I = 1, 4
    for (I = 1; I <= 4; I++) {
        E2 = E2 + DLGAMA(1+DRAY[I]) + DLGAMA(1+DB1+DRAY[I])
          + DLGAMA(1+DB2+DRAY[I]) - DLGAMA(2+DRAY[I]+DRAY[4+I]);
        E = E - DLGAMA(1-DN+DRAY[I]);
    } // 549
//      E = E + .5*E2
    E = E + .5*E2;

//      ANSWER = DEXP(E) * SUM
    ANSWER = DEXP(E) * SUM;

// 600  TEMP = NMIN + R
L600:
    TEMP = NMIN + R;
//      IF (.NOT. BTEST(TEMP, 0))  ANSWER = -ANSWER
    if (!BTEST(TEMP, 0)) ANSWER = -ANSWER;

// 700  IF (.NOT. RACHSW)  GO TO 800
    if (!RACHSW) goto L800;

//      IF (.NOT. BTEST(R, 0))  ANSWER = -ANSWER
    if (!BTEST(R, 0)) ANSWER = -ANSWER;

// 800  (DBLSW always true => label 850)
L800:
// 850  SIXJ = ANSWER
//      RETURN
    return ANSWER;

// 900  (triangle rules not satisfied)
L900:
    return 0.0;
}

// ---- Public entry points ----

// SIXJ(A,B,C,X,Y,Z)
double SIXJ(int A, int B, int C, int X, int Y, int Z)
{
    return sixj_impl(A, B, C, X, Y, Z, false);
}

// RACAH(A,B,Y,X,C,Z)
// Fortran: ENTRY RACAH(A,B,Y,X,C,Z)
// The ENTRY statement reuses the FUNCTION's parameter names.
// So RACAH's 1st arg maps to A, 2nd to B, 3rd to Y, 4th to X, 5th to C, 6th to Z.
// The caller calls RACAH(j1, j2, j5, j4, j3, j6) to compute W(j1,j2,j5,j4;j3,j6).
double RACAH(int A_r, int B_r, int Y_r, int X_r, int C_r, int Z_r)
{
    return sixj_impl(A_r, B_r, C_r, X_r, Y_r, Z_r, true);
}

// ============================================================================
// Begin WIG9J, DSGMAL, get_date, GAUSSL translations
// ============================================================================


// ============================================================================
// FUNCTION WIG9J — Wigner 9-J coefficient (lines 1070-1518)
//
// COMPUTES 9-J COEFFICIENTS.
// Input: 9 integers that are TWICE the corresponding J's.
// Returns double precision 9-J coefficient.
//
// S. PIEPER
// 9/6/74 - FIRST VERSION
// 4/11/88 - VAX/VMS VERSION - R.OSBORN AND G.L.GOODMAN
// 7/3/91 - RS6000 - s.p.
// ============================================================================
double WIG9J(int J1, int J2, int J3, int J4, int J5, int J6,
             int J7, int J8, int J9)
{
    // IMPLICIT INTEGER*4 (A-C, I-R, T-Z), REAL*8 (D-H, S)
    // Vars A-C, I-R, T-Z are int; D-H, S are double.

    // Statement function: DLF(N) = DLGAMA(DFLOAT(N))
    auto DLF = [](int N) -> double { return DLGAMA(DFLOAT(N)); };

    // DATA ZERONE(2) /0.D0, 1.D0/  (1-based)
    double ZERONE[3] = {0, 0.0, 1.0};

    // Local integer variables
    int J1S[6];   // DIMENSION J1S(5), 1-based
    int J2S[6];   // DIMENSION J2S(5)
    int J3S[6];   // DIMENSION J3S(5)
    int AS[5][4]; // DIMENSION AS(4,3) -> AS[5][4], 1-based
    int BS[3][4]; // DIMENSION BS(2,3) -> BS[3][4], 1-based
    int MX2S[4];  // DIMENSION MX2S(3), 1-based
    int MN2S[4];  // DIMENSION MN2S(3), 1-based
    int MN1S[4];  // DIMENSION MN1S(3), 1-based
    int RS[4];    // DIMENSION RS(3), 1-based
    int IDS[3][4]; // DIMENSION IDS(2,3) -> IDS[3][4], 1-based

    // Local double variables
    double DAS[5][4]; // DIMENSION DAS(4,3) -> DAS[5][4], 1-based
    double DBS[3][4]; // DIMENSION DBS(2,3) -> DBS[3][4], 1-based
    double DRS[4];    // DIMENSION DRS(3), 1-based
    double DMN1S[4];  // DIMENSION DMN1S(3), 1-based
    double DMN2S[4];  // DIMENSION DMN2S(3), 1-based
    double DMX2S[4];  // DIMENSION DMX2S(3), 1-based
    double SUMS[4];   // DIMENSION SUMS(3), 1-based

    // Scalars — integer (A-C, I-R, T-Z)
    int I, J, N, X;
    int XSTART, XEND, XSUM, XHALF, XPIECE;
    int MX2, MN2, NMAX, NMIN, NMAXM2;
    int ISIGN;

    // Scalars — double (D-H, S)
    double DX, DXHALF;
    double SUMOUT, SUMLOG;
    double SUM, SUMP;
    double D1, D2, D3, D4, D5, D6, D7, D8;
    double E1, E5;
    double F1, F5;
    double G1, G5;
    double H1, H5;
    double DNMAX, DNMIN;

    double wig9j_result;

    //
    // PICK UP THE J'S IN THE APPROPRIATE ORDER
    // NOTE THAT THE JIS ARRAYS REALLY CONTAIN 2*J
    //
    J1S[1] = J1;
    J1S[2] = J6;
    J1S[3] = J8;
    J2S[1] = J2;
    J2S[2] = J4;
    J2S[3] = J9;
    J3S[1] = J3;
    J3S[2] = J5;
    J3S[3] = J7;
    //
    // EXPAND THE ARRAYS TO ALLOW CYCLIC INDICES
    //
    for (I = 1; I <= 2; I++) {
        J1S[I+3] = J1S[I];
        J2S[I+3] = J2S[I];
        J3S[I+3] = J3S[I];
    } // 59

    //
    // DETERMINE THE RANGE OF 2*X
    //
    XSTART = IABS(J2S[1] - J1S[2]);
    XEND = J2S[1] + J1S[2];
    XSUM = XEND;
    for (I = 2; I <= 3; I++) {
        XSTART = MAX0(XSTART, IABS(J2S[I]-J1S[I+1]));
        XPIECE = J2S[I] + J1S[I+1];
        XEND = MIN0(XEND, XPIECE);
        XSUM = XSUM + XPIECE;
    } // 129

    if (XEND < XSTART) goto label_900;

    XSUM = XSUM - XEND;

    //
    // EXTRACT 1/2 INTEGER PART OF X AND CONVERT 2*X TO INTPART(X)
    //
    XHALF = 0;
    if (BTEST(XSTART, 0)) XHALF = 1;
    XSTART = ISHFT(XSTART, -1);
    XEND = ISHFT(XEND, -1);
    XSUM = ISHFT(XSUM, -1);
    DXHALF = ZERONE[XHALF+1];
    DX = XSTART - 1;
    XSUM = XSUM + 1;

    //
    // SETUP THE J1+J2-J3, ETC
    //
    for (I = 1; I <= 3; I++) {
        //
        // THE AS, FOR EACH 6-J THEY ARE
        //   A1 = J1+J2-J3
        //   A2 = J4+J5-J3
        //   A3 = J1+J5-J6
        //   A4 = J4+J2-J6
        // IN A3 AND A4 THE J6 = X IS LEFT OUT.
        //
        AS[1][I] = ISHFT(J1S[I]+J2S[I]-J3S[I], -1);
        AS[2][I] = ISHFT(J1S[I+1]+J2S[I+2]-J3S[I], -1);
        AS[3][I] = ISHFT(J1S[I]+J2S[I+2], -1);
        AS[4][I] = ISHFT(J2S[I]+J1S[I+1], -1);
        //
        // NOW MAKE SURE NONE OF THEM ARE NEGATIVE WHICH WOULD INDICATE
        // A VIOLATION OF THE TRIANGLE RULES.
        //
        if (BTEST(AS[1][I],30) || BTEST(AS[2][I],30)) goto label_900;
        //
        // THE BS, FOR EACH 6-J THEY ARE
        //    B1 = J3 + J6 - J1 - J4
        //    B2 = J3 + J6 - J2 - J5
        // THE J6 = X IS LEFT OUT
        //
        BS[1][I] = J2S[I+2] - AS[2][I] - AS[3][I];
        BS[2][I] = J1S[I] - AS[1][I] - AS[3][I];
    } // 159

    //
    // SETUP THE QUANTITIES FOR EACH OF THE 3 6-J SUMS.
    //
    for (I = 1; I <= 3; I++) {
        MN1S[I] = MIN0(AS[1][I], AS[2][I]);
        MN2S[I] = MIN0(AS[3][I], AS[4][I]);
        MX2S[I] = -MIN0(BS[1][I], BS[2][I]);
        RS[I] = 1 + XHALF + AS[3][I] + AS[4][I];
        for (J = 1; J <= 4; J++) {
            DAS[J][I] = AS[J][I];
        } // 239
        DBS[1][I] = BS[1][I];
        DBS[2][I] = BS[2][I];
        DRS[I] = 1 + DXHALF + DAS[3][I] + DAS[4][I];
        DMN1S[I] = MN1S[I];
        DMN2S[I] = MN2S[I];
        DMX2S[I] = MX2S[I];
        //
        // THESE ARE INDICES FOR THE 4 FACTORIALS THAT DEPEND ON X
        //
        IDS[1][I] = AS[2][I] + BS[2][I];
        IDS[2][I] = AS[1][I] + BS[1][I];
    } // 299

    //
    // FIND THE LOG OF THE OUTSIDE FACTORS NOW
    //
    SUMOUT = 0;

    //
    // XSUM IS >= ANY OF THE REQUIRED FACTORIAL ARGUMENTS
    //
    if (XSUM > LOGFAC.MAXLF) goto label_350;

    //
    // WE CAN USE THE LOG(FACTORIAL) TABLE
    //
    for (I = 1; I <= 3; I++) {
        for (J = 1; J <= 2; J++) {
            SUMOUT = SUMOUT +
                LOGFAC.LF[LOGFAC.LFBIAS+1+AS[J][I]] - LOGFAC.LF[LOGFAC.LFBIAS+2+J3S[I]+AS[J][I]]
                + LOGFAC.LF[LOGFAC.LFBIAS+1-AS[J][I]+J1S[I-1+J]]
                + LOGFAC.LF[LOGFAC.LFBIAS+1-AS[J][I]+J2S[I-2+2*J]];
        }
    } // 339
    goto label_400;

    //
    // WE MUST USE LOG GAMMA
    //
label_350:
    for (I = 1; I <= 3; I++) {
        for (J = 1; J <= 2; J++) {
            SUMOUT = SUMOUT +
                DLGAMA(1+DAS[J][I]) - DLF(2+J3S[I]+AS[J][I])
                + DLF(1-AS[J][I]+J1S[I-1+J])
                + DLF(1-AS[J][I]+J2S[I-2+2*J]);
        }
    } // 359

label_400:
    SUMOUT = .5*SUMOUT;

    //
    // WE ARE NOW READY TO DO THE 9-J LOOP
    //
    wig9j_result = 0;

    for (X = XSTART; X <= XEND; X++) {
        DX = DX + 1;
        //
        // DO THE PARTS OF EACH 6-J THAT DEPEND ON X
        //
        SUMLOG = SUMOUT;
        ISIGN = 0;

        for (I = 1; I <= 3; I++) {
            //
            // GET THE COEFFICIENTS FOR THIS 6-J
            //
            MX2 = MX2S[I] - X;
            MN2 = MN2S[I] - X;
            NMAX = MN1S[I];
            DNMAX = DMN1S[I];
            DNMIN = 0;
            NMIN = 0;
            //
            // GET THE RANGE OF THE 6-J SUM
            //
            if (MN2 >= NMAX) goto label_520;
            NMAX = MN2;
            DNMAX = DMN2S[I] - DX;

        label_520:
            if (MX2 <= 0) goto label_540;
            NMIN = MX2;
            DNMIN = DMX2S[I] - DX;

        label_540:
            SUM = 1;
            ISIGN = ISIGN + NMIN;

            //
            // THIS ONE TEST CHECKS ALL 12 INEQUALITIES IMPLIED BY THE
            // 4 TRIANGLE CONDITIONS FOR THE GIVEN 6-J.
            //
            if (NMAX - NMIN < 0) goto label_900;
            if (NMAX - NMIN == 0) goto label_700;
            // else (NMAX - NMIN > 0) fall through to label_600

            // label_600:
            D1 = DAS[1][I] - (DNMAX-1);
            D2 = DAS[2][I] - (DNMAX-1);
            D8 = DRS[I] - (DNMAX-1);
            D3 = DAS[3][I] - (DNMAX-1+DX);
            D4 = DAS[4][I] - (DNMAX-1+DX);
            D5 = -DNMAX;
            D6 = -DBS[1][I] - (DX + DNMAX);
            D7 = -DBS[2][I] - (DX + DNMAX);

            E1 = D1*D2*D3*D4;
            E5 = D5*D6*D7*D8;

            NMAXM2 = NMAX-2;
            SUM = 1 + E1/E5;
            if (NMAXM2 < NMIN) goto label_700;

            F1 = (D1+1)*(D2+1)*(D3+1)*(D4+1) - E1;
            F5 = (D5+1)*(D6+1)*(D7+1)*(D8+1) - E5;

            G1 = (D1+2)*(D2+2)*(D3+2)*(D4+2) - F1 - F1 - E1;
            G5 = (D5+2)*(D6+2)*(D7+2)*(D8+2) - F5 - F5 - E5;

            H1 = (D1+3)*(D2+3)*(D3+3)*(D4+3) - 3*(G1+F1) - E1;
            H5 = (D5+3)*(D6+3)*(D7+3)*(D8+3) - 3*(G5+F5) - E5;

            //
            // COMPUTE THE SUM IN THE RACAH FORMULA
            //
            for (N = NMIN; N <= NMAXM2; N++) {
                E1 = E1 + F1;
                E5 = E5 + F5;
                F1 = F1 + G1;
                F5 = F5 + G5;
                G1 = G1 + H1;
                G5 = G5 + H5;
                H1 = H1 + 24;
                H5 = H5 + 24;
                SUM = 1 + SUM*(E1/E5);
            } // 659

        label_700:
            SUMS[I] = SUM;

            //
            // HERE ARE THE PARTS OF THE FACTORIALS FOR NMIN AND ALSO THOSE
            // FACTORIALS THAT DEPEND ON X
            //
            if (RS[I] > LOGFAC.MAXLF) goto label_730;

            //
            // WE CAN USE THE LOG FACTORIAL TABLE
            //
            SUMLOG = SUMLOG +
                LOGFAC.LF[LOGFAC.LFBIAS+1-X+AS[4][I]] + LOGFAC.LF[LOGFAC.LFBIAS+1+X+IDS[1][I]]
                + LOGFAC.LF[LOGFAC.LFBIAS+1+X+IDS[2][I]] - LOGFAC.LF[LOGFAC.LFBIAS+2+XHALF+X+AS[4][I]]
                + LOGFAC.LF[LOGFAC.LFBIAS+1-NMIN+RS[I]] - LOGFAC.LF[LOGFAC.LFBIAS+1+NMIN]
                - LOGFAC.LF[LOGFAC.LFBIAS+1+NMIN+X+BS[1][I]]
                - LOGFAC.LF[LOGFAC.LFBIAS+1+NMIN+X+BS[2][I]]
                - LOGFAC.LF[LOGFAC.LFBIAS+1-NMIN+AS[1][I]] - LOGFAC.LF[LOGFAC.LFBIAS+1-NMIN+AS[2][I]]
                - LOGFAC.LF[LOGFAC.LFBIAS+1-NMIN-X+AS[3][I]]
                - LOGFAC.LF[LOGFAC.LFBIAS+1-NMIN-X+AS[4][I]];

            goto label_759;

            //
            // WE MUST USE LOG GAMMA
            //
        label_730:
            SUMLOG = SUMLOG +
                DLGAMA(1-DX+DAS[4][I]) + DLGAMA(1+DX+DAS[2][I]+DBS[2][I])
                + DLGAMA(1+DX+DAS[1][I]+DBS[1][I])
                - DLGAMA(2+DXHALF+DX+DAS[4][I])
                + DLGAMA(1-DNMIN+DRS[I]) - DLGAMA(1+DNMIN)
                - DLGAMA(1+DNMIN+DX+DBS[1][I])
                - DLGAMA(1+DNMIN+DX+DBS[2][I])
                - DLGAMA(1-DNMIN+DAS[1][I]) - DLGAMA(1-DNMIN+DAS[2][I])
                - DLGAMA(1-DNMIN-DX+DAS[3][I])
                - DLGAMA(1-DNMIN-DX+DAS[4][I]);

        label_759:;
        } // DO 759 I = 1, 3

        //
        // NOW MULTIPLY THE 3 6-J PARTS TOGETHER AND SUM
        //
        SUMP = SUMS[1] * SUMS[2] * DEXP(SUMLOG) * SUMS[3]
            * (1+DXHALF+(DX+DX));

        if (BTEST(ISIGN, 0)) SUMP = -SUMP;
        wig9j_result = wig9j_result + SUMP;
    } // 799 DO X = XSTART, XEND

    //
    return wig9j_result;

    //
    // THE TRIANGLE RULES ARE NOT SATISFIED
    //
label_900:
    wig9j_result = 0;
    return wig9j_result;
}

// ============================================================================
// SUBROUTINE DSGMAL — Coulomb phases (lines 1945-2124)
//
// COULOMB PHASES
// ?/?/? - ORIGINAL VERSION
// 12/28/79 - CDC TO CNI, TITLE COMMENT - RPG
// 7/9/80 - USE CDM PREFIX - S.P.
// 7/3/91 - RS6000 version - s.p.
// ============================================================================
void DSGMAL(double ETA, int NO, double* DSG)
{
    // IMPLICIT REAL*8 ( A-H, O-Z )
    // I-N are int, rest are double

    // Local integer variables
    int I, J, K, L, M, N;

    // Local double variables
    double X, XSQ, RSQ, QSUM, PSUM, R;

    // DIMENSION P(26), Q(23) — 1-based
    static const double P[27] = {0,
        0.3101657810129948870e-08,
        0.3015947478999109100e-05,
        0.2056442871538789580e-03,
        0.4445588614720562960e-02,
        0.4241372519181223210e-01,
        0.2063788681972964780e+00,
        0.5468922695201207380e+00,
        0.7949480657792201980e+00,
        0.5931315608708119800e+00,
        0.1770600075369600650e+00,
        0.2410878569735943570e-08,
        0.9326569949959554900e-06,
        0.7585813793142706060e-04,
        0.2007067188606569890e-02,
        0.1889244677027970660e-01,
        0.4443159218008918380e-01,
       -0.1596016183850622740e+00,
       -0.6851417732759696180e+00,
       -0.5772019207036128280e+00,
        0.1395410517881128990e+03,
        0.1090522693580003650e+04,
        0.4853585241219512340e+03,
       -0.6587580539038850700e+03,
       -0.4130231383850326670e+02,
        0.9229335782342382280e+01,
       -0.9999999999999999440e+00
    };

    static const double Q[24] = {0,
        0.1626320933946765800e-05,
        0.1808585994795038710e-03,
        0.5771439575139208560e-02,
        0.7876590720775496210e-01,
        0.5451579810646677990e+00,
        0.2085536101102752230e+01,
        0.4573132992231460970e+01,
        0.5697178591946184810e+01,
        0.3737311340573166380e+01,
        0.8418712796550826150e-09,
        0.5251990683341128460e-06,
        0.6387938953567890980e-04,
        0.2628228899692266520e-02,
        0.4482759373310784680e-01,
        0.3453892205334134070e+00,
        0.1223260295405861340e+01,
        0.1881032630842315980e+01,
       -0.1411156525262207530e+03,
       -0.9965362519625687700e+03,
       -0.5358167734466346130e+03,
        0.6544271625884707650e+03,
        0.4215892515370163450e+02,
       -0.9312669115675720870e+01
    };

    static const double EZ  = 0.1805547071605106970e+01;
    static const double EZ1 = 0.1805547714233398440e+01;
    static const double EZ2 = 0.6426282915176236330e-06;

    X = DABS(ETA);
    XSQ = X*X;
    if (X > 4.0) goto label_200;
    if (X > 2.0) goto label_100;

    //
    // ABS(ETA) IN (0,2) RANGE SO=ETA(ETA**2-EZSQ)R(ETA**2)
    //
    I = 7;
    J = 1;
    K = 1;
    M = 1;
    RSQ = XSQ;
    goto label_300;

    //
    // ABS(ETA) IN (2,4) RANGE SO=ETA*R(ETA**2)
    //
label_100:
    I = 6;
    J = 11;
    K = 10;
    M = 2;
    RSQ = XSQ;
    goto label_300;

    //
    // ABS(ETA) IN (4,-) RANGE SO=ATAN(ETA)/2+ETA(LN(1+ETA**2)/2
    //                            +R(1/ETA**2))
    //
label_200:
    I = 4;
    J = 20;
    K = 18;
    M = 3;
    RSQ = 1.0e0/XSQ;

label_300:
    QSUM = Q[K];
    N = K+I;
    for (L = K; L <= N; L++) {
        QSUM = QSUM*RSQ + Q[L+1];
    } // 305
    PSUM = P[J];
    N = J+I+1;
    for (L = J; L <= N; L++) {
        PSUM = PSUM*RSQ + P[L+1];
    } // 310
    R = PSUM/(QSUM*RSQ + 1.0e0);
    if (M == 1) goto label_400;
    if (M == 2) goto label_500;
    if (M == 3) goto label_600;

label_400:
    DSG[0] = X*R*(X+EZ)*((X-EZ1)+EZ2);
    goto label_610;

label_500:
    DSG[0] = X*R;
    goto label_610;

label_600:
    DSG[0] = DATAN(X)/2.0e0 + X*(DLOG(1.0e0+XSQ)/2.0e0 + R);

    //
    // SL(ETA)=SO(ETA)+SUM(I=1TOL)DATAN(ETA/I)
    //
label_610:
    if (NO == 0) goto label_650;
    R = 1.0e0;
    for (I = 1; I <= NO; I++) {
        RSQ = X/R;
        DSG[I] = DSG[I-1] + DATAN(RSQ);
        R = R + 1.0e0;
    } // 620
    if (ETA >= 0.0e0) goto label_660;
    for (I = 1; I <= NO; I++) {
        DSG[I] = -DSG[I];
    } // 640

label_650:
    if (ETA < 0.0e0) DSG[0] = -DSG[0];

label_660:
    return;
}

// ============================================================================
// SUBROUTINE get_date — return date as "DD Mon YY" (lines 2125-2152)
//
// RETURNS THE DATE IN THE FORM  19 Jan 01
// date - a character*9
//
// 11/16/01 - all new
// ============================================================================
void get_date(char* date)
{
    // character*3 months(12)
    static const char* months[13] = {"",
        "Jan", "Feb", "Mar", "Apr", "May", "Jun",
        "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"
    };

    // Use C++ time functions in place of Fortran date_and_time
    std::time_t now = std::time(nullptr);
    struct std::tm* tm_info = std::localtime(&now);

    int day   = tm_info->tm_mday;        // 1-31
    int month = tm_info->tm_mon + 1;     // 1-12
    int year  = tm_info->tm_year % 100;  // 2-digit year

    // Format: "DD Mon YY" (9 characters)
    // date(1:2) = day, date(4:6) = month name, date(8:9) = year
    std::memset(date, ' ', 9);
    date[9] = '\0';

    // Write day (positions 0-1, i.e. Fortran 1:2)
    date[0] = (day / 10 == 0) ? ' ' : ('0' + day / 10);
    date[1] = '0' + day % 10;

    // date(3) = ' '  (space separator, already set)

    // Write month name (positions 3-5, i.e. Fortran 4:6)
    date[3] = months[month][0];
    date[4] = months[month][1];
    date[5] = months[month][2];

    // date(7) = ' '  (space separator, already set)

    // Write year (positions 7-8, i.e. Fortran 8:9)
    date[7] = '0' + year / 10;
    date[8] = '0' + year % 10;

    return;
}

// ============================================================================
// SUBROUTINE GAUSSL — Gauss-Legendre integration points and weights
//                     (lines 2153-2311)
//
// POINTS AND WEIGHTS FOR GAUSS-LEGENDRE INTEGRATION OVER (-1 1)
// COMPUTED BY METHOD GIVEN ON PP.88 AND 89 OF
// 'METHODS OF NUMERICAL INTEGRATION' DAVIS AND RABINOWITZ
// ACADEMIC PRESS (1975)
//
// THIS SUBROUTINE IS AN ADAPTATION OF SUBROUTINE GRULE LISTED ON
// P.369 (APP 2) OF DAVIS AND RABINOWITZ
// MHM  JUNE 29, 1975
// 7/5/91 - RS6000 - s.p.
// ============================================================================
void GAUSSL(int N, double* X, double* W)
{
    // IMPLICIT REAL*8 ( A-H, O-Z )
    // N, I, K, M are int. All others double.

    // Local integer variables
    int I, K, M;

    // Local double variables
    double X0, E1, E2, E3, T, T1;
    double PK, PKM1, PKP1, FK, DEN;
    double D1, DPN, D2PN, D3PN, D4PN;
    double U, V, H, H1, H2;
    double P, DP, FX;

    //
    // TAKE CARE OF SPECIAL CASES (N.LE.1)
    //
    if (N > 0) goto label_200;
    std::printf("\n  0000 NUMBER OF GAUSS POINTS = %5d 0000\n"
                "   0000             ABSURD             0000\n", N);
    goto label_600;

label_200:
    if (N != 1) goto label_300;
    X[1] = 0;
    W[1] = 2;
    goto label_600;

    //
    // BEGIN COMPUTATION
    //
label_300:
    M = (N+1)/2;
    E1 = N*(N+1);
    E2 = 3.141592653589793200e0 / (4*N+2);
    E3 = 1.-(1.-1.0e0/N)/(8.0e0*N*N);
    for (I = 1; I <= M; I++) {
        T = (4*I-1)*E2;
        X0 = E3*DCOS(T);
        //
        // USE RECURSION FORMULAS AND THE 2ND ORDER DIFFERENTIAL EQUATION
        // SATISFIED BY PN TO COMPUTE PK AND ITS FIRST FOUR DERIVATIVES
        // AT X0
        //
        PKM1 = 1;
        PK = X0;
        FK = 1;
        for (K = 2; K <= N; K++) {
            FK = FK + 1;
            T1 = X0*PK;
            PKP1 = T1 - PKM1 - (T1-PKM1)/FK + T1;
            PKM1 = PK;
            PK = PKP1;
        } // 400
        DEN = 1. - X0*X0;
        D1 = N*(PKM1-X0*PK);
        DPN = D1/DEN;
        D2PN = (2.*X0*DPN - E1*PK)/DEN;
        D3PN = (4.*X0*D2PN + (2.-E1)*DPN)/DEN;
        D4PN = (6.*X0*D3PN + (6.-E1)*D2PN)/DEN;
        //
        // USE THIRD-ORDER NEWTON-RAPHSON METHOD
        //
        U = PK/DPN;
        V = D2PN/DPN;
        //
        // Intel compiler workaround: compute H in two steps
        //
        H1 = (1.+.5*U*(V+U*(V*V-U*D3PN/(3.*DPN))));
        H2 = -U*H1;
        H = H2;
        //
        // EVALUATE PN AND D(PN)/DX AT X0+H USING TAYLOR EXPANSION
        //
        P = PK + H*(DPN + .5*H*(D2PN + H/3.*(D3PN + .25*H*D4PN)));
        DP = DPN + H*(D2PN + .5*H*(D3PN + H*D4PN/3.));
        H = H - P/DP;
        X[N+1-I] = X0 + H;
        X[I] = -X[N+1-I];
        //
        // CALCULATE WEIGHT W(I) CORRESPONDING TO I'TH GAUSS POINT X0I
        // W(I)=2*(1-X0I**2)/(N*PN-1(X0I))**2
        //
        FX = D1 - H*E1*(PK + .5*H*(DPN + H/3.*(D2PN + .25*H*(D3PN + .2*H*D4PN))));
        W[N+1-I] = 2.*(1-X[I]*X[I])/(FX*FX);
        W[I] = W[N+1-I];
    } // 500
    if (M+M > N) X[M] = 0;

label_600:
    return;
}
