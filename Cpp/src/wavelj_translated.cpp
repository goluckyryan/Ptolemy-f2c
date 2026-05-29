// wavelj_translated.cpp — WAVELJ subroutine
// Translated from source.f lines 35432-36167
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
#include <cstdlib>

// External timing function
extern double second();

// Correct TCOEF signature (5 integer args) — overrides the 2-arg stub in ptolemy_forward.h
// The real Fortran FUNCTION TCOEF(ITYP, LAS, L, JP, JSP) is at source.f L34242
double TCOEF_5(int ITYP, int LAS, int L, int JP, int JSP);

void WAVELJ(int L, int JP, int NWP, int NUMPTS, float* RGRID, float* WAVER,
            float* WAVEI, double* WAVR, double* WAVI, double* VREAL,
            double* VIMAG, double* VCENT)
{
//
//     COMPUTES SCATTERING WAVE FOR A GIVEN CHANNEL, L, JP.
//
//     THE INPUT IS
//
//     L - L VALUE.
//     JP - J OF PROJECTILE - USED ONLY FOR SPIN ORBIT FORCE.
//     NWP - CHANNEL NUMBER (1 OR 2).
//     NUMPTS - NUMBER OF DESIRED OUTPUT (BY INTERPOLATION) POINTS.
//     RGRID - VALUES OF R AT WHICH OUTPUT IS TO BE FOUND.
//     WAVER, WAVEI - THE COMPUTED WAVE FUNCTION AT RGRID IS STORED HERE.
//     WAVR, WAVI - TWO WORK ARRAYS OF DIMENSION NSTEP+4.
//     VREAL, VIMAG - THE REAL AND IMAGINARY PARTS OF THE POTENTIAL IN NUMEROV FORM.
//     VCENT - 1/RHO**2 IN THE NUMEROV FORM.
//

    // Local variables
    double DL, DLSQ, H, H2;
    double THISR, THISI, THSRSV, AA1, BB1;
    double D, A1, A2, CR, CI, SJR, SJI;
    double F, G, F1, G1;
    double SDOTL, TCO;
    double RT, RINV, P, PS, X1, X2, X3, X4, X5;
    double C1, C2, C3, C4, C5;
    double RBYH, STPINV;
    double CCH, RSTART;
    double TT1, T1, T2;
    int I, II, IR, N, N1, N2, NSTEP, NSTEP2, ISTRT, IRET;
    int LSOR, LSOI, LIMV, LTEN, LINDXE, NNWP, IMX;
    int N1ADD;
    logical BADSW;

    // CHARACTER*8 CHANAM(3) / 'INCOMING', 'OUTGOING', ' ' /
    const char* CHANAM[4]; // 1-based
    CHANAM[1] = "INCOMING";
    CHANAM[2] = "OUTGOING";
    CHANAM[3] = "        ";

    // DATA TWELTH / .08333333333333333333330D0 /
    // DATA OVER24 / .0416666666666666666670D0 /
    double TWELTH = .0833333333333333333333e0;
    double OVER24 = .041666666666666666667e0;

//
//     SETUP THE PARAMETERS FOR THIS CHANNEL
//
    if (NUMPTS < 0) goto L5;
    if (L < 0) return;
    I = IABS(2 * L - WAVCOM.JSPS[NWP]);
    if (WAVCOM.SOSWS[NWP] && JP < I) return;
//
    TT1 = second();
L5:
    DL = L;
    DLSQ = DL * (DL + 1);
    H = WAVCOM.HS[NWP];
    H2 = H * H;
    FLOAT_common.STEPSZ = WAVCOM.RSTEPS[NWP];
    NSTEP = WAVCOM.NSTPSS[NWP];
//     DEFINITION OF COULOMB FUNCTIONS FOR MULTI-CHANNEL
    if (SWITCH.PROBLM != 24) goto L40;
//
//     WE ARE PASSING THE F & G IN IN A REAL*4 ARRAY.
//
    F1 = WAVER[1];
    G1 = WAVER[2];
    F = WAVER[3];
    G = WAVER[4];
    goto L45;
//
//     FOR NON-C.C. CALCULATIONS, STORE J-TOTAL FOR L-DEPENDANT LINKULES
//
L40:
    JBLOCK.J = 2 * L;
//
//     IF LINKULE POTENTIALS OR WAVEFUNCTIONS ARE IN USE,
//     AND IF THIS IS A DWBA, THEN WE MUST SETUP THE PROPER /LNKBLK/
//
L45:
    if (SWITCH.PROBLM < 20 || SWITCH.PROBLM == 21) goto L50;
    for (I = 1; I <= LNKBLK.NUMLNK; I++) {
        for (II = 1; II <= 6; II++) {
            LNKBLK.LNKADR[I][II] = LOCFIT.linkule.LNKAD2[I][NWP][II];
        }
    }
//
//     START THE INTEGRATION AT THE POINT DETERMINED FOR THE PREVIOUS
//     SOLUTION IF WE HAVE NOT RESET TO A NEW, LOWER, L VALUE.
//
L50:
    WAVCOM.NFIRST = -NUMPTS;
    ISTRT = WAVCOM.NFIRST - 2;
    if (NUMPTS < 0) goto L220;
    if (SWITCH.PROBLM == 24)
        std::fprintf(stderr, "WAVELJ entry: L=%d NWP=%d NUMPTS=%d LASTL[%d]=%d LASTNF[%d]=%d\n",
                     L, NWP, NUMPTS, NWP, WAVCOM.LASTL[NWP], NWP, WAVCOM.LASTNF[NWP]);
    WAVCOM.NFIRST = WAVCOM.LASTNF[NWP];
    ISTRT = WAVCOM.NFIRST - 2;
    if (L >= WAVCOM.LASTL[NWP]) goto L160;
    WAVCOM.NFIRST = 2;
    ISTRT = 0;
    WAVCOM.LASTZR = 1;
//
//     STEP1I (=STEPI) IS THE SMALLEST VALUE OF THE WAVEFUNCTION TO
//     ALLOW AND THUS CONTROLS THE CHOICE OF NFIRST.
//
L160:
    if (WAVCOM.STEPI == 1) WAVCOM.STEPI = 1.e-6;
//
//     OUR INITIAL GUESS IS ALWAYS STEPR*STEPI.
//
    THISR = WAVCOM.STEPI * WAVCOM.STEPR;
    THISI = THISR;
//
//     ASSUME   U(R) = R**(L+1)  TO GET RATIO OF FIRST TWO POINTS
//
    AA1 = ISTRT;
//     ZERO THE FIRST FEW POINTS IF NECESSARY.
    if (WAVCOM.NUMFIT != 0 || ISTRT < WAVCOM.LASTZR) goto L180;
    for (I = WAVCOM.LASTZR; I <= ISTRT; I++) {
        WAVR[I] = 0.;
        WAVI[I] = 0.;
    }
L180:
    D = AA1 / (AA1 + 1);
    WAVR[ISTRT + 1] = THISR * std::pow(D, (L + 1));
    WAVI[ISTRT + 1] = THISI * std::pow(D, (L + 1));
    WAVR[ISTRT + 2] = THISR;
    WAVI[ISTRT + 2] = THISI;
    THISR = THISR + THISI;
    THSRSV = THISR;
//
//     COMPUTE THE TOTAL EFFECTIVE POTENTIAL FOR THIS L AND TEMPORARILY
//     STORE IN THE WAVEFUNCTION ARRAYS.
//     MODIFY POTENTIAL FOR MULTI-CHANNELS
    if (SWITCH.PROBLM == 24) goto L220;
//
    for (I = ISTRT; I <= NSTEP; I++) {
        WAVR[I + 4] = VREAL[I + 1] + DLSQ * VCENT[I + 1];
        WAVI[I + 4] = VIMAG[I + 1];
    }
    goto L230;
L220:
    CCH = H2 * TWELTH;
    for (I = ISTRT; I <= NSTEP; I++) {
        WAVR[I + 4] = VREAL[I + 1] + CCH + DLSQ * VCENT[I + 1];
        WAVI[I + 4] = VIMAG[I + 1];
    }
//
//        MAKE ONE EXTRA SAFE
//
L230:
    WAVR[NSTEP + 5] = 1;
    WAVI[NSTEP + 5] = 0;
//
//     ADD IN THE SPIN ORBIT FORCES IF THEY EXIST
//
    if (WAVCOM.ISORS[NWP] + WAVCOM.ISOIS[NWP] == 0) goto L250;
//
//     SDOTL IS REALLY  SIGMA.L
//
    SDOTL = .25 * (JP * (JP + 2) - WAVCOM.JSPS[NWP] * (WAVCOM.JSPS[NWP] + 2)) - DLSQ;
    SDOTL = SDOTL / WAVCOM.JSPS[NWP];
    if (WAVCOM.ISORS[NWP] == 0) goto L240;
    LSOR = LOCPTRS.Z[WAVCOM.ISORS[NWP]];
    for (I = ISTRT; I <= NSTEP; I++) {
        WAVR[I + 4] = WAVR[I + 4] + SDOTL * ALLOC(LSOR + I);
    }
L240:
    if (WAVCOM.ISOIS[NWP] == 0) goto L250;
    LSOI = LOCPTRS.Z[WAVCOM.ISOIS[NWP]];
    LIMV = LOCPTRS.Z[WAVCOM.IIMVS[NWP]];
    for (I = ISTRT; I <= NSTEP; I++) {
        WAVI[I + 4] = WAVI[I + 4] + SDOTL * ALLOC(LSOI + I);
    }
//
//     ADD IN THE TENSOR POTENTIALS IF THEY EXIST.
//
L250:
    for (II = 1; II <= 3; II++) {
        if (WAVCOM.ITENS[NWP][2 * II - 1] + WAVCOM.ITENS[NWP][2 * II] == 0) goto L259;
        TCO = TCOEF_5(II, L, L, JP, WAVCOM.JSPS[NWP]);
        if (TCO == 0) goto L259;
        if (WAVCOM.ITENS[NWP][2 * II - 1] == 0) goto L255;
        LTEN = LOCPTRS.Z[WAVCOM.ITENS[NWP][2 * II - 1]];
        for (I = ISTRT; I <= NSTEP; I++) {
            WAVR[I + 4] = WAVR[I + 4] + TCO * ALLOC(LTEN + I);
        }
L255:
        if (WAVCOM.ITENS[NWP][2 * II] == 0) goto L259;
        LTEN = LOCPTRS.Z[WAVCOM.ITENS[NWP][2 * II]];
        for (I = ISTRT; I <= NSTEP; I++) {
            WAVI[I + 4] = WAVI[I + 4] + TCO * ALLOC(LTEN + I);
        }
L259: ;
    }
//
//     NOW PROCESS L-DEPENDENT POTENTIALS WHICH COME ONLY FROM LINKULES.
//
    A1 = -(FLOAT_common.STEPSZ / CNSTNT.HBARC) * (FLOAT_common.STEPSZ / CNSTNT.HBARC) * REDMS()[NWP] / 6.0;
    N = NSTEP - ISTRT + 1;
    RSTART = ISTRT * FLOAT_common.STEPSZ;
//
    if (LNKBLK.LNKADR[1][3] == 0 || LNKBLK.LNKADR[1][4] != 1) goto L270;
    LINKUL(LNKBLK.LNKADR[1][3], *reinterpret_cast<char8*>(&LNKBLK.LNKADR[1][1]),
           &LNKBLK.LNKADR[1][5],
           1, 4, IRET,
           L, (double)JP, RSTART, FLOAT_common.STEPSZ, N, &WAVR[ISTRT + 4], &A1, (char*)&NWP);
    if (IRET < 0) std::exit(7777);
//
L270:
    if (LNKBLK.LNKADR[2][3] == 0 ||
        LNKBLK.LNKADR[2][4] != 1) goto L300;
    LINKUL(LNKBLK.LNKADR[2][3], *reinterpret_cast<char8*>(&LNKBLK.LNKADR[2][1]),
           &LNKBLK.LNKADR[2][5],
           2, 4, IRET,
           L, (double)JP, RSTART, FLOAT_common.STEPSZ, N, &WAVI[ISTRT + 4], &A1, (char*)&NWP);
    if (IRET < 0) std::exit(7777);
//
//     COMPUTE FIRST TWO XSI'S
//
L300:
    II = ISTRT + 1;
    for (I = ISTRT; I <= II; I++) {
        WAVR[I + 3] =
            TWELTH * (WAVR[I + 1] * WAVR[I + 4] - WAVI[I + 1] * WAVI[I + 4]);
        WAVI[I + 3] =
            TWELTH * (WAVR[I + 1] * WAVI[I + 4] + WAVI[I + 1] * WAVR[I + 4]);
    }
//
    if (NUMPTS < 0) return;
//
//     READY TO DO THE INTEGRATION.
//     IS IT DONE HERE OR IN A LINKULE?
//
    T1 = second();
    FTIME.TIMES[1] = FTIME.TIMES[1] + (float)(T1 - TT1);
//
    if (LNKBLK.LNKADR[6][3] == 0) goto L350;
//
    LINKUL(LNKBLK.LNKADR[6][3], *reinterpret_cast<char8*>(&LNKBLK.LNKADR[6][1]),
           &LNKBLK.LNKADR[6][5],
           WAVCOM.NFIRST, 3, IRET, L, (double)JP, H * ISTRT, H, NSTEP + 1,
           WAVR, WAVI, (char*)&NWP);
    if (IRET >= 0) goto L500;
    std::exit(7777);
//
//     WILL DO INTEGRATION HERE
//
L350: ;
//
//      NUMEROV METHOD (A LA RAYNAL) IS USED IN THE FOLLOWING LOOP
//
    N1 = WAVCOM.NFIRST;
    N2 = NSTEP;
    if (L >= WAVCOM.LASTL[NWP]) goto L420;
    N1ADD = 20;
    if (L > 16) N1ADD = std::exp(50.0 / DL);
//
L410:
    N2 = N1 + N1ADD;
    N2 = MIN0(N2, NSTEP);
//
L420:
    D = (12. / (WAVR[N1 + 4] * WAVR[N1 + 4] + WAVI[N1 + 4] * WAVI[N1 + 4]));
//
    for (I = N1; I <= N2; I++) {
//
        WAVR[I + 3] = ((-10.) * WAVR[I + 2]) + (WAVR[I] - WAVR[I + 1]);
        WAVI[I + 3] = ((-10.) * WAVI[I + 2]) + (WAVI[I] - WAVI[I + 1]);
//
        WAVR[I + 1] = (WAVR[I + 4] * WAVR[I + 3] + WAVI[I + 4] * WAVI[I + 3]) * D;
        WAVI[I + 1] = (WAVR[I + 4] * WAVI[I + 3] - WAVI[I + 4] * WAVR[I + 3]) * D;
//
        D = (12. / (WAVR[I + 5] * WAVR[I + 5] + WAVI[I + 5] * WAVI[I + 5]));
    }
    if (N2 == NSTEP) goto L500;
//
//     CHECK FOR NEARNESS OF OVERFLOW
//
    THISR = DABS(WAVR[N2 + 1]) + DABS(WAVI[N2 + 1]);
    if (WAVCOM.PWBGSW) std::printf(" CHECKING SIZE: NWP, L, N1, N2, MAG =%5d%5d%5d%5d%13.3G\n",
        NWP, L, N1, N2, THISR);
//
    if (THISR > CNSTNT.BIGNUM) goto L470;
    if (THISR < 10000 * THSRSV) N1ADD = N1ADD + 5;
    goto L490;
//
//     MUST SCALE DOWN.  WE WILL MULTIPLY EVERYTHING BY 1/THISR.
//     FIRST WE FIND WHERE THIS WILL RESULT IN |U| = STEPI.
//
L470:
    THISI = WAVCOM.STEPI * THISR;
    N1 = ISTRT;
    for (ISTRT = N1; ISTRT <= N2; ISTRT++) {
        if (DABS(WAVR[ISTRT + 1]) >= THISI) goto L480;
        WAVR[ISTRT + 1] = 0;
        WAVI[ISTRT + 1] = 0.;
    }
//
//
//     SHOULD WE REDUCE N1ADD
//
L480:
    if (THISR * 1.e-10 > CNSTNT.BIGNUM) N1ADD = MAX0(N1ADD / 2, 1);
    if (WAVCOM.PWBGSW) std::printf("     RESCALING: NEW ISTRT =%6d%6d\n", ISTRT, N1ADD);
//
    THISI = 1 / THISR;
    N1 = N2 + 2;
    for (I = ISTRT; I <= N1; I++) {
        WAVR[I + 1] = THISI * WAVR[I + 1];
        WAVI[I + 1] = THISI * WAVI[I + 1];
    }
    THISR = 0;
//
L490:
    N1 = N2 + 1;
    THSRSV = THISR;
    goto L410;
//
//     END OF THE INTEGRATION
//
L500:
    T2 = second();
    FTIME.TIMES[2] = FTIME.TIMES[2] + (float)(T2 - T1);
    FTIME.NUMWAV = FTIME.NUMWAV + 1;
    FTIME.NUMLOP = FTIME.NUMLOP + (NSTEP - WAVCOM.NFIRST + 1);
//
//     EXTRACT S-MATRIX AND NORMALIZATION FROM U AND U AT NBAKCM
//     STEPS BACK FROM END.
//
// F IS THE REGULAR COULOMB WAVE
//  G IS THE IRREGULAR COULOMB WAVE
//  NOTE THAT F1,G1 ARE REGULAR AND IRREGULAR COULOMB FUNTIONS AT NBACK
//
    I = L + (WAVCOM.IOFFIT - 1) * (KANDM.LOMOST + 1);
    if (SWITCH.PROBLM == 24) goto L520;
    F = ALLOC(LOCPTRS.Z[WAVCOM.NFS[NWP]] + I);
    G = ALLOC(LOCPTRS.Z[WAVCOM.NGS[NWP]] + I);
    F1 = ALLOC(LOCPTRS.Z[WAVCOM.NF1S[NWP]] + I);
    G1 = ALLOC(LOCPTRS.Z[WAVCOM.NG1S[NWP]] + I);
//
L520: ;
//  NOW FINDING REAL S AND IMAG S
//
    A1 = WAVR[NSTEP + 1] * F1 + WAVI[NSTEP + 1] * G1
         - WAVR[NSTEP + 1 - WAVCOM.NBAKCM] * F
         - WAVI[NSTEP + 1 - WAVCOM.NBAKCM] * G;
    A2 = -WAVR[NSTEP + 1] * G1 + WAVI[NSTEP + 1] * F1
         + WAVR[NSTEP + 1 - WAVCOM.NBAKCM] * G
         - WAVI[NSTEP + 1 - WAVCOM.NBAKCM] * F;
    CR = -WAVR[NSTEP + 1] * F1 + WAVI[NSTEP + 1] * G1
         + WAVR[NSTEP + 1 - WAVCOM.NBAKCM] * F
         - WAVI[NSTEP + 1 - WAVCOM.NBAKCM] * G;
    CI = -WAVR[NSTEP + 1] * G1 - WAVI[NSTEP + 1] * F1
         + WAVR[NSTEP + 1 - WAVCOM.NBAKCM] * G
         + WAVI[NSTEP + 1 - WAVCOM.NBAKCM] * F;
//
//     TRANSLATE THE ELASTIC S-MATRICES TO THE LX SYSTEM AND SAVE THEM.
//
    SJR = (CR * A1 + CI * A2) / (A1 * A1 + A2 * A2);
    SJI = (CI * A1 - CR * A2) / (A1 * A1 + A2 * A2);
    LINDXE = ALLOCS.FACFR4 * LOCPTRS.Z[WAVCOM.IINDXE[NWP]] - ALLOCS.FACFR4 + 1;
    if ((NWP != 2) || (SWITCH.PROBLM != 24)) {
        JPTOLX(L, L, JP, NWP, &SJR, &SJI,
               &ILLOC(LINDXE), ALLOC_base(LOCPTRS.Z[KANDM.ISMATS[NWP]]));
    }
//
//  NORMALIZATION HERE
//
    A1 = .5 * (F * (1. + SJR) + SJI * G);
    A2 = .5 * (G * (1. - SJR) + SJI * F);
    CR = (WAVR[NSTEP + 1] * A1 + WAVI[NSTEP + 1] * A2) /
         (WAVR[NSTEP + 1] * WAVR[NSTEP + 1] + WAVI[NSTEP + 1] * WAVI[NSTEP + 1]);
    WAVCOM.ALPHAR = CR;
    CI = (WAVR[NSTEP + 1] * A2 - WAVI[NSTEP + 1] * A1) /
         (WAVR[NSTEP + 1] * WAVR[NSTEP + 1] + WAVI[NSTEP + 1] * WAVI[NSTEP + 1]);
    WAVCOM.ALPHAI = CI;
//
//     FIND NEXT STARTING POINT
//     FIND POINT WHERE ALL RENORMALIZED WAVEFUNCTION WILL BE < STEPI
//
    AA1 = DABS(CR) + DABS(CI);
    BB1 = WAVCOM.STEPI / AA1;
    if (SWITCH.PROBLM == 24)
        std::fprintf(stderr, "WAVELJ norm: L=%d NWP=%d CR=%.6E CI=%.6E AA1=%.6E BB1=%.6E STEPI=%.6E ISTRT=%d NSTEP=%d\n",
                     L, NWP, CR, CI, AA1, BB1, WAVCOM.STEPI, ISTRT, NSTEP);
    if (SWITCH.PROBLM == 24)
        std::fprintf(stderr, "  F=%.6E G=%.6E F1=%.6E G1=%.6E SJR=%.6E SJI=%.6E\n",
                     F, G, F1, G1, SJR, SJI);
    if (SWITCH.PROBLM == 24)
        std::fprintf(stderr, "  WAVR[NSTEP+1]=%.6E WAVI[NSTEP+1]=%.6E WAVR[N-NBACK+1]=%.6E WAVI[N-NBACK+1]=%.6E\n",
                     WAVR[NSTEP+1], WAVI[NSTEP+1], WAVR[NSTEP+1-WAVCOM.NBAKCM], WAVI[NSTEP+1-WAVCOM.NBAKCM]);
    for (I = ISTRT; I <= NSTEP; I++) {
        if (DABS(WAVR[I + 1]) + DABS(WAVI[I + 1]) >= BB1) goto L630;
        WAVR[I + 1] = 0.;
        WAVI[I + 1] = 0.;
    }
//
L630:
    ISTRT = I;
    WAVCOM.LASTZR = I + 1;
//
//     SET IT BACK SOME FOR LINKULES
//
    I = MAX0(I + 1, 2);
    if (SWITCH.PROBLM == 24)
        std::fprintf(stderr, "WAVELJ exit: L=%d NWP=%d LASTNF=%d (loop exited at I=%d)\n",
                     L, NWP, I, ISTRT);
    WAVCOM.LASTNF[NWP] = I;
    WAVCOM.LASTL[NWP] = L;
    BADSW = ((AA1 > CNSTNT.BIGNUM * 1.e-10 && WAVCOM.NFIRST > 2) ||
             AA1 < CNSTNT.SMLNUM * 1.e-5) ? TRUE_F : FALSE_F;
    if (BADSW) std::printf("0**** WARNING, LARGE WAVEFUNCTION RENORMALIZATION:\n");
    if (WAVCOM.PWBGSW || BADSW)
        std::printf(" NWP, L =%2d%4d  NFIRST =%4d  NEW NFIRST =%4d  RENORM =%11.3E   NUM PTS:%6d%6d%6d\n",
                    NWP, L, WAVCOM.NFIRST, I, AA1, NSTEP, WAVCOM.NSTP2S[NWP], NUMPTS);
    WAVCOM.NFIRST = ISTRT;
//
    if (WAVCOM.NUMFIT != 0) goto L880;
//
    if (!WAVCOM.PWAVSW) goto L750;
    NNWP = NWP;
    if (WAVCOM.STANSW) NNWP = 3;
    if (WAVCOM.SOSWS[NWP]) goto L730;
    std::printf(" %8s ELASTIC S-MATRIX FOR L =%4d:%15.5G +%13.5G I\n",
                CHANAM[NNWP], L, SJR, SJI);
    goto L750;
L730:
    if (WAVCOM.TCSWS[NWP]) goto L740;
    std::printf(" %8s ELASTIC S-MATRIX FOR L =%4d,   JP =%4d/2:%15.5G +%13.5G I\n",
                CHANAM[NNWP], L, JP, SJR, SJI);
    goto L750;
L740:
    std::printf(" %8s ELASTIC S-MATRIX FOR L =%4d,  LAS =%4d,   JP =%4d/2:%15.5G +%13.5G I\n",
                CHANAM[NNWP], L, L, JP, SJR, SJI);
//
//  STORE NORMALIZED WAVE FUNCTIONS
//
L750:
    N2 = NSTEP + 2;
    for (I = ISTRT; I <= N2; I++) {
        A1 = CR * WAVR[I + 1] - CI * WAVI[I + 1];
        WAVI[I + 1] = CR * WAVI[I + 1] + CI * WAVR[I + 1];
        WAVR[I + 1] = A1;
    }
    II = ISTRT + 2;
    if (WAVCOM.PWBGSW)
        std::printf(" STARTS:%5d%14.5G%14.5G%5d%14.5G%14.5G%5d%14.5G%14.5G\n",
                    ISTRT, WAVR[ISTRT], WAVI[ISTRT],
                    ISTRT + 1, WAVR[ISTRT + 1], WAVI[ISTRT + 1],
                    ISTRT + 2, WAVR[ISTRT + 2], WAVI[ISTRT + 2]);
//
//     IF NECESSARY, CONTINUE THE INTEGRATION USING JUST COULOMB
//
    NSTEP2 = WAVCOM.NSTP2S[NWP];
    if (NSTEP2 <= NSTEP) goto L800;
//
    NSTEP = NSTEP + 1;
    RT = NSTEP * FLOAT_common.STEPSZ;
//
    for (I = NSTEP; I <= NSTEP2; I++) {
        RINV = 1. / RT;
        WAVR[I + 3] = (-10.) * WAVR[I + 2] + (WAVR[I] - WAVR[I + 1]);
        WAVI[I + 3] = (-10.) * WAVI[I + 2] + (WAVI[I] - WAVI[I + 1]);
        D = 12. / (WAVCOM.XFACS[NWP][1] + RINV * (WAVCOM.XFACS[NWP][2] +
            (DLSQ * WAVCOM.XFACS[NWP][3]) * RINV));
        WAVR[I + 1] = D * WAVR[I + 3];
        WAVI[I + 1] = D * WAVI[I + 3];
        RT = RT + FLOAT_common.STEPSZ;
    }
//
L800:
    if (WAVCOM.STANSW || NUMPTS == 0) goto L880;
    // Debug: check max wave function and potential before interpolation
    if (NWP == 2 && L <= 1) {
        // Also check potential at a few points
        static int potdbg = 0;
        if (potdbg++ < 2) {
        }
        double wmax = 0;
        int wmax_idx = 0;
        for (I = 1; I <= NSTEP2; I++) {
            double wm = std::fabs(WAVR[I]) + std::fabs(WAVI[I]);
            if (wm > wmax) { wmax = wm; wmax_idx = I; }
        }
    }
//
//     NOW INTERPOLATE TO THE DESIRED GRID
//
    STPINV = 1 / WAVCOM.RSTEPS[NWP];
    IMX = NSTEP2 - 2;
    if (WAVCOM.PWBGSW) std::printf(" INTERPOLATING:%6d%6d%15.5G\n", NUMPTS, IMX, STPINV);
    for (IR = 1; IR <= NUMPTS; IR++) {
        RT = RGRID[IR];
//
//     FIND THE LOCATION OF R IN THE TABLE
//
        RBYH = RT * STPINV;
        I = (int)(RBYH + 0.5);
        I = MAX0(2, MIN0(I, IMX));
        P = RBYH - I;
//     USE 5-POINT LAGRANGIAN INTERPOLATION (ABRAMOWITZ & STEGUN 25.2.15)
        PS = P * P;
        X1 = P * (PS - 1.) * OVER24;
        X2 = X1 + X1;
        X3 = X1 * P;
        X4 = X2 + X2 - 0.5 * P;
        X5 = X4 * P;
        C1 = X3 - X2;
        C5 = X3 + X2;
        C3 = X5 - X3;
        C2 = X5 - X4;
        C4 = X5 + X4;
        C3 = C3 + C3 + 1.;
        WAVER[IR] = (float)(C1 * WAVR[I - 1] - C2 * WAVR[I] + C3 * WAVR[I + 1]
                            - C4 * WAVR[I + 2] + C5 * WAVR[I + 3]);
        WAVEI[IR] = (float)(C1 * WAVI[I - 1] - C2 * WAVI[I] + C3 * WAVI[I + 1]
                            - C4 * WAVI[I + 2] + C5 * WAVI[I + 3]);
//
    }
//
L880:
    FTIME.TIMES[3] = FTIME.TIMES[3] + (float)(second() - T2);
    return;
}
