// wavetc_translated.cpp — Verbatim C++ translation of SUBROUTINE WAVETC
// Translated from source.f lines 36169-36900

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include "ptolemy_alloc.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_forward.h"
#include "ptolemy_io.h"
#include <cstdio>
#include <cstring>
#include <cmath>

extern double second();

// TCOEF has a different signature than what's in ptolemy_forward.h;
// the tensor coefficient function takes 5 integer arguments.
extern double TCOEF(int ITYP, int LAS, int L, int JP, int JSP);

// ============================================================================
// SUBROUTINE WAVETC — lines 36169-36900
//
// COMPUTES TENSOR-COUPLED SCATTERING WAVE FOR A GIVEN L, LAS, JP.
//
// THE INPUT IS
//   L - L VALUE.
//   LAS - ASYMPTOTIC L VALUE.
//   JP - 2*J OF PROJECTILE.
//   NWP - CHANNEL NUMBER (1 OR 2).
//   NWAVF - DIMENSION OF WAVR, WAVI FIRST INDEX.
//   NUMPTS - NUMBER OF DESIRED OUTPUT POINTS.
//   RGRID - VALUES OF R AT WHICH OUTPUT IS TO BE FOUND.
//   WAVER, WAVEI - THE COMPUTED WAVE FUNCTION AT RGRID IS STORED HERE.
//   WAVR, WAVI - TWO WORK ARRAYS OF DIMENSION 4*NWAVEF.
//   VREAL, VIMAG - THE REAL AND IMAGINARY PARTS OF THE POTENTIAL.
//   VCENT - 1/RHO**2 IN THE NUMEROV FORM.
// ============================================================================
void WAVETC(int L, int LAS, int JP, int NWP, int NWAVF, int NUMPTS,
            float* RGRID, float* WAVER, float* WAVEI,
            double* WAVR_flat, double* WAVI_flat,
            double* VREAL, double* VIMAG, double* VCENT)
{
    // WAVR, WAVI: dimension (NWAVF,2,2) — 1-based Fortran column-major indexing
    // Fortran: WAVR(i,j,k) at offset (i-1) + (j-1)*NWAVF + (k-1)*NWAVF*2
    #define WAVR(i,j,k) WAVR_flat[((i)-1) + ((j)-1)*NWAVF + ((k)-1)*NWAVF*2]
    #define WAVI(i,j,k) WAVI_flat[((i)-1) + ((j)-1)*NWAVF + ((k)-1)*NWAVF*2]

    // CR, CI: dimension (2,2) — Fortran code walks IL=1..4 through CR(IL,1),
    // treating the 2x2 column-major array as flat[4].
    // Use flat 1-based arrays of 5 elements (index 1..4).
    double CR_flat[5]; // 1-based: CR_flat[1]..CR_flat[4]
    double CI_flat[5]; // 1-based
    // Macro for 2D (2,2) column-major 1-based access:
    #define CR(i,j) CR_flat[((i)-1) + ((j)-1)*2 + 1]
    #define CI(i,j) CI_flat[((i)-1) + ((j)-1)*2 + 1]

    // U1R, U1I: dimension (2,2) — same flat treatment, IL=1..4
    double U1R_flat[5]; // 1-based: U1R_flat[1]..U1R_flat[4]
    double U1I_flat[5]; // 1-based
    #define U1R(i,j) U1R_flat[((i)-1) + ((j)-1)*2 + 1]
    #define U1I(i,j) U1I_flat[((i)-1) + ((j)-1)*2 + 1]

    // COMMON block aliases
    auto& FACFR4     = ALLOCS.FACFR4;
    auto* Z          = LOCPTRS.Z;       // 1-based
    auto& STEPSZ     = FLOAT_common.STEPSZ;
    auto& STEPR      = WAVCOM.STEPR;
    auto& STEPI      = WAVCOM.STEPI;
    auto& HBARC      = CNSTNT.HBARC;
    auto& BIGNUM     = CNSTNT.BIGNUM;
    auto& SMLNUM     = CNSTNT.SMLNUM;
    auto& PROBLM     = SWITCH.PROBLM;
    auto& NUMLNK     = LNKBLK.NUMLNK;
    auto& NFIRST     = WAVCOM.NFIRST;
    auto& LASTZR     = WAVCOM.LASTZR;
    auto& ALPHAR     = WAVCOM.ALPHAR;
    auto& ALPHAI     = WAVCOM.ALPHAI;
    auto& NUMFIT     = WAVCOM.NUMFIT;
    auto& NWAVEF     = WAVCOM.NWAVEF;
    auto& NBAKCM     = WAVCOM.NBAKCM;
    auto& PWBGSW     = WAVCOM.PWBGSW;
    auto& STANSW     = WAVCOM.STANSW;
    auto& NUMWAV     = FTIME.NUMWAV;
    auto& NUMLOP     = FTIME.NUMLOP;

    // Local variables — all at function scope for goto compatibility
    double DLSQ[3];   // 1-based: DLSQ[1], DLSQ[2]
    double SDOTL[3];   // 1-based
    double TCO[4];     // 1-based: TCO[1], TCO[2], TCO[3]
    double WIR[3];     // 1-based: WIR[1], WIR[2]
    double WII[3];     // 1-based
    int BADSW; // LOGICAL
    double THISR, THISI, THSRSV;
    double AA1, BB1;
    double DL, H, D, A1, RT;
    double T1, T2;
    double STPINV, RBYH, P, PS;
    double X1, X2, X3, X4, X5;
    double C1, C2, C3, C4, C5;
    double TR, TI, TD, WTR, WTI, RINV;
    int LL, NSTEP, NSTEP2, ISTRT;
    int LWAVRS, LWAVIS;
    int N1, N2, N1ADD;
    int LFORS, LSOR, LSOI, LTEN;
    int I, II, IL, ILP, IS, ILAS, IR, IMX;
    int IRET;

    // DATA
    double TWELTH = .08333333333333333333330e0;
    double OVER24 = .0416666666666666666670e0;

    //
    // SETUP THE PARAMETERS FOR THIS CHANNEL
    //
    LL = JP/2 - 1;
    if ( L < 0  ||  LL < 0 )  return;
    double TT1 = second();
    T2 = TT1;
    LWAVRS = Z[WAVCOM.IWAVRS[NWP]]*FACFR4 - FACFR4;
    LWAVIS = Z[WAVCOM.IWAVIS[NWP]]*FACFR4 - FACFR4;
    //
    // IF THIS WAVE HAS ALREADY BEEN CALCULATED, GO TRANSFER IT.
    //
    if ( ( LL != L  ||  LAS != L )
      &&  L >= WAVCOM.LASTL[NWP] )  goto L900;
    //
    DL = LL;
    DLSQ[1] = DL*(DL+1);
    DLSQ[2] = (DL+2)*(DL+3);
    H = WAVCOM.HS[NWP];
    STEPSZ = WAVCOM.RSTEPS[NWP];
    NSTEP = WAVCOM.NSTPSS[NWP];
    //
    // IF LINKULE POTENTIALS OR WAVEFUNCTIONS ARE IN USE,
    // AND IF THIS IS A DWBA, THEN WE MUST SETUP THE PROPER /LNKBLK/
    //
    if ( PROBLM < 20  ||  PROBLM == 21 )  goto L50;
    for (I = 1; I <= NUMLNK; I++) {
        for (II = 1; II <= 6; II++) {
            LNKBLK.LNKADR[I][II] = LOCFIT.linkule.LNKAD2[I][NWP][II];
        }
    }
    //
    // START THE INTEGRATION AT THE POINT DETERMINED FOR THE PREVIOUS
    // SOLUTION IF WE HAVE NOT RESET TO A NEW, LOWER, L VALUE.
    //
L50:
    NFIRST = WAVCOM.LASTNF[NWP];
    if ( LL >= WAVCOM.LASTL[NWP] )  goto L160;
    NFIRST = 2;
    LASTZR = 1;
    //
    // STEP1I (=STEPI) IS THE SMALLEST VALUE OF THE WAVEFUNCTION TO
    // ALLOW AND THUS CONTROLS THE CHOICE OF NFIRST.
    //
L160:
    if ( STEPI == 1 )  STEPI = 1.0e-6;
    //
    // OUR INITIAL GUESS IS ALWAYS STEPR*STEPI.
    //
    THISR = STEPI*STEPR;
    THISI = THISR;
    //
    // ASSUME   U(R) = R**(L+1)  TO GET RATIO OF FIRST TWO POINTS
    //
    AA1 = NFIRST - 2;
    ISTRT = NFIRST - 2;
    // ZERO THE FIRST FEW POINTS IF NECESSARY.
    if ( NFIRST < LASTZR )  goto L180;
    for (I = LASTZR; I <= NFIRST; I++) {
        for (IL = 1; IL <= 4; IL++) {
            WAVR(I,IL,1) = 0.;
            WAVI(I,IL,1) = 0.;
        }
    }
L180:
    D = AA1/(AA1+1);
    WAVR(ISTRT+1,1,1) = THISR*std::pow(D, LL+1);
    WAVI(ISTRT+1,1,1) = THISI*std::pow(D, LL+1);
    WAVR(ISTRT+2,1,1) = THISR;
    WAVI(ISTRT+2,1,1) = THISI;
    // THE RATIO OF THE L+2 WAVE FUNCTION TO THE L WAVE FUNCTION
    // DEPENDS ON THE LOCAL WAVE NUMBER.
    AA1 = (2*LL+3)*(2*LL+5);
    AA1 = (NFIRST-1)*(NFIRST-1) * 12.*(VREAL[NFIRST]-1.) / AA1;
    WAVR(ISTRT+1,2,2) = THISR*AA1*std::pow(D, LL+3);
    WAVI(ISTRT+1,2,2) = THISI*AA1*std::pow(D, LL+3);
    WAVR(ISTRT+2,2,2) = THISR*AA1;
    WAVI(ISTRT+2,2,2) = THISI*AA1;
    //
    // COMPUTE THE TOTAL EFFECTIVE POTENTIAL FOR THIS L AND TEMPORARILY
    // STORE IN THE WAVEFUNCTION ARRAYS.
    //
    for (I = ISTRT; I <= NSTEP; I++) {
        WAVR(I+4,1,1) = VREAL[I+1] + DLSQ[1]*VCENT[I+1];
        WAVR(I+4,2,1) = VREAL[I+1] + DLSQ[2]*VCENT[I+1];
        WAVI(I+4,1,1) = VIMAG[I+1];
        WAVI(I+4,2,1) = WAVI(I+4,1,1);
        WAVR(I+4,1,2) = 0.;
        WAVI(I+4,1,2) = 0.;
    }
    THSRSV = THISR + THISI;
    //
    // ADD IN THE SPIN ORBIT FORCES IF THEY EXIST
    //
    if ( WAVCOM.ISORS[NWP]+WAVCOM.ISOIS[NWP] == 0 )  goto L250;
    //
    // SDOTL IS REALLY  SIGMA.L.
    //
    AA1 = .25*(JP*(JP+2) - WAVCOM.JSPS[NWP]*(WAVCOM.JSPS[NWP]+2));
    SDOTL[1] = ( AA1 - DLSQ[1] ) / WAVCOM.JSPS[NWP];
    SDOTL[2] = ( AA1 - DLSQ[2] ) / WAVCOM.JSPS[NWP];
    if ( WAVCOM.ISORS[NWP] == 0 )  goto L240;
    LSOR = Z[WAVCOM.ISORS[NWP]];
    for (I = ISTRT; I <= NSTEP; I++) {
        WAVR(I+4,1,1) = WAVR(I+4,1,1) + SDOTL[1]*ALLOC(LSOR+I);
        WAVR(I+4,2,1) = WAVR(I+4,2,1) + SDOTL[2]*ALLOC(LSOR+I);
    }
L240:
    if ( WAVCOM.ISOIS[NWP] == 0 )  goto L250;
    LSOI = Z[WAVCOM.ISOIS[NWP]];
    for (I = ISTRT; I <= NSTEP; I++) {
        WAVI(I+4,1,1) = WAVI(I+4,1,1) + SDOTL[1]*ALLOC(LSOI+I);
        WAVI(I+4,2,1) = WAVI(I+4,2,1) + SDOTL[2]*ALLOC(LSOI+I);
    }
    //
    // ADD IN THE TENSOR POTENTIALS IF THEY EXIST.
    //
L250:
    for (II = 1; II <= 3; II++) {
        if ( WAVCOM.ITENS[NWP][2*II-1]+WAVCOM.ITENS[NWP][2*II] == 0 ) goto L259;
        TCO[1] = TCOEF( II, LL, LL, JP, WAVCOM.JSPS[NWP] );
        TCO[2] = TCOEF( II, LL+2, LL+2, JP, WAVCOM.JSPS[NWP] );
        TCO[3] = TCOEF( II, LL, LL+2, JP, WAVCOM.JSPS[NWP] );
        if ( WAVCOM.ITENS[NWP][2*II-1] == 0 ) goto L255;
        LTEN = Z[WAVCOM.ITENS[NWP][2*II-1]];
        for (IL = 1; IL <= 3; IL++) {
            if ( TCO[IL] == 0 ) goto L254;
            for (I = ISTRT; I <= NSTEP; I++) {
                WAVR(I+4,IL,1) = WAVR(I+4,IL,1) + TCO[IL]*ALLOC(LTEN+I);
            }
L254:       ;
        }
L255:
        if ( WAVCOM.ITENS[NWP][2*II] == 0 ) goto L259;
        LTEN = Z[WAVCOM.ITENS[NWP][2*II]];
        for (IL = 1; IL <= 3; IL++) {
            if ( TCO[IL] == 0 ) goto L258;
            for (I = ISTRT; I <= NSTEP; I++) {
                WAVI(I+4,IL,1) = WAVI(I+4,IL,1) + TCO[IL]*ALLOC(LTEN+I);
            }
L258:       ;
        }
L259:   ;
    }
    //
    // NOW PROCESS L-DEPENDENT POTENTIALS WHICH COME ONLY FROM LINKULES.
    //
    // 260
    A1 = -(STEPSZ/HBARC)*(STEPSZ/HBARC) * REDMS()[NWP] / 6.0;
    //
    if ( LNKBLK.LNKADR[1][3] == 0  ||  LNKBLK.LNKADR[1][4] != 1 )  goto L270;
    {
    double J_dummy = 0.0;
    LINKUL( LNKBLK.LNKADR[1][3], *reinterpret_cast<char8*>(&LNKBLK.LNKADR[1][1]),
      &LNKBLK.LNKADR[1][5],
      1, 4, IRET,
      LL, J_dummy, 0.0, STEPSZ, NSTEP+1, &WAVR(4,1,1), &A1, (char*)&NWP );
    if ( IRET < 0 )  std::exit(7777);
    LINKUL( LNKBLK.LNKADR[1][3], *reinterpret_cast<char8*>(&LNKBLK.LNKADR[1][1]),
      &LNKBLK.LNKADR[1][5],
      1, 4, IRET,
      LL+2, J_dummy, 0.0, STEPSZ, NSTEP+1, &WAVR(4,2,1), &A1, (char*)&NWP );
    if ( IRET < 0 )  std::exit(7777);
    }
    //
L270:
    if ( LNKBLK.LNKADR[2][3] == 0  ||
      LNKBLK.LNKADR[2][4] != 1 )  goto L300;
    {
    double J_dummy = 0.0;
    LINKUL( LNKBLK.LNKADR[2][3], *reinterpret_cast<char8*>(&LNKBLK.LNKADR[2][1]),
      &LNKBLK.LNKADR[2][5],
      2, 4, IRET,
      LL, J_dummy, 0.0, STEPSZ, NSTEP+1, &WAVI(4,1,1), &A1, (char*)&NWP );
    if ( IRET < 0 )  std::exit(7777);
    LINKUL( LNKBLK.LNKADR[2][3], *reinterpret_cast<char8*>(&LNKBLK.LNKADR[2][1]),
      &LNKBLK.LNKADR[2][5],
      2, 4, IRET,
      LL+2, J_dummy, 0.0, STEPSZ, NSTEP+1, &WAVI(4,2,1), &A1, (char*)&NWP );
    if ( IRET < 0 )  std::exit(7777);
    }
    //
    // COMPUTE FIRST TWO XSI'S
    //
L300:
    II = ISTRT+1;
    for (I = ISTRT; I <= II; I++) {
        for (IL = 1; IL <= 2; IL++) {
            ILP = 3 - IL;
            for (IS = 1; IS <= 2; IS++) {
                WAVR(I+3,IL,IS) =
                    TWELTH*( WAVR(I+1,IL,IS)*WAVR(I+4,IL,1)
                    - WAVI(I+1,IL,IS)*WAVI(I+4,IL,1)
                    + WAVR(I+1,ILP,IS)*WAVR(I+4,1,2)
                    - WAVI(I+1,ILP,IS)*WAVI(I+4,1,2) );
                WAVI(I+3,IL,IS) =
                    TWELTH*( WAVR(I+1,IL,IS)*WAVI(I+4,IL,1)
                    + WAVI(I+1,IL,IS)*WAVR(I+4,IL,1)
                    + WAVR(I+1,ILP,IS)*WAVI(I+4,1,2)
                    + WAVI(I+1,ILP,IS)*WAVR(I+4,1,2) );
            }
        }
    }
    //
    // READY TO DO THE INTEGRATION.
    // IS IT DONE HERE OR IN A LINKULE?
    //
    T1 = second();
    FTIME.TIMES[1] = FTIME.TIMES[1] + (float)(T1 - TT1);
    //
    if ( LNKBLK.LNKADR[6][3] == 0 )  goto L350;
    //
    LINKUL( LNKBLK.LNKADR[6][3], *reinterpret_cast<char8*>(&LNKBLK.LNKADR[6][1]),
      &LNKBLK.LNKADR[6][5],
      NFIRST, 3, IRET, L, (double)JP, H*ISTRT, H, NSTEP+1,
      WAVR_flat, WAVI_flat, (char*)&NWP );
    if ( IRET >= 0 )  goto L500;
    std::exit(7777);
    //
    // WILL DO INTEGRATION HERE
    //
L350:
    //
    // NUMEROV METHOD (A LA RAYNAL)
    //
    N1 = NFIRST;
    N2 = NSTEP;
    N1ADD = 20;
    if ( LL >= WAVCOM.LASTL[NWP] )  goto L420;
    N1ADD = 20;
    if ( L > 16 ) N1ADD = (int)std::exp( 50.0 / L );
    //
L410:
    N2 = N1 + N1ADD;
    N2 = MIN0(N2, NSTEP);
    //
    // LOOP THROUGH RADII.
    //
L420:
    for (I = N1; I <= N2; I++) {
        //
        // WI(I) = 12./W(I)  (MATRIX INVERSION)
        //
        TR = WAVR(I+4,1,1)*WAVR(I+4,2,1) - WAVI(I+4,1,1)*WAVI(I+4,2,1)
              - WAVR(I+4,1,2)*WAVR(I+4,1,2) + WAVI(I+4,1,2)*WAVI(I+4,1,2);
        TI = WAVR(I+4,1,2)*WAVI(I+4,1,2);
        TI = TI + TI
             - WAVR(I+4,1,1)*WAVI(I+4,2,1) - WAVI(I+4,1,1)*WAVR(I+4,2,1);
        TD = 12. / ( TR*TR + TI*TI );
        TR = TR*TD;
        TI = TI*TD;
        WIR[1] = TR*WAVR(I+4,2,1) - TI*WAVI(I+4,2,1);
        WII[1] = TR*WAVI(I+4,2,1) + TI*WAVR(I+4,2,1);
        WIR[2] = TR*WAVR(I+4,1,1) - TI*WAVI(I+4,1,1);
        WII[2] = TR*WAVI(I+4,1,1) + TI*WAVR(I+4,1,1);
        WTR = (-TR)*WAVR(I+4,1,2) + TI*WAVI(I+4,1,2);
        WTI = (-TR)*WAVI(I+4,1,2) - TI*WAVR(I+4,1,2);
        //
        // LOOP THROUGH THE 2 INDEPENDENT SOLUTIONS.
        //
        for (IS = 1; IS <= 2; IS++) {
            //
            // XSI(I) = (-10.)*XSI(I-1) + U(I-1) - XSI(I-2)  (FOR EACH L)
            //
            for (IL = 1; IL <= 2; IL++) {
                WAVR(I+3,IL,IS) = (-10.)*WAVR(I+2,IL,IS)
                    + WAVR(I,IL,IS) - WAVR(I+1,IL,IS);
                WAVI(I+3,IL,IS) = (-10.)*WAVI(I+2,IL,IS)
                    + WAVI(I,IL,IS) - WAVI(I+1,IL,IS);
            }
            //
            // U(I) = (12/W(I))*XSI(I)  (FOR EACH L)
            //
            for (IL = 1; IL <= 2; IL++) {
                ILP = 3 - IL;
                WAVR(I+1,IL,IS)
                    = WIR[IL]*WAVR(I+3,IL,IS) - WII[IL]*WAVI(I+3,IL,IS)
                    + WTR*WAVR(I+3,ILP,IS) - WTI*WAVI(I+3,ILP,IS);
                WAVI(I+1,IL,IS)
                    = WIR[IL]*WAVI(I+3,IL,IS) + WII[IL]*WAVR(I+3,IL,IS)
                    + WTR*WAVI(I+3,ILP,IS) + WTI*WAVR(I+3,ILP,IS);
            }
            //
        }
        //
    }
    if ( N2 == NSTEP )  goto L500;
    //
    // CHECK FOR NEARNESS OF OVERFLOW
    //
    THISR = 0.;
    for (IL = 1; IL <= 4; IL++) {
        THISR = THISR + DABS(WAVR(N2+1,IL,1)) + DABS(WAVI(N2+1,IL,1));
    }
    if ( PWBGSW )  std::printf(" CHECKING SIZE: NWP, LL, N1, N2, MAG =%5d%5d%5d%5d%13.3E\n",
        NWP, LL, N1, N2, THISR);
    //
    if ( THISR > BIGNUM )  goto L470;
    if ( THISR < 100*THSRSV ) N1ADD = N1ADD + 5;
    goto L490;
    //
    // MUST SCALE DOWN.
    // FIRST WE FIND WHERE THIS WILL RESULT IN |U| = STEPI.
    //
L470:
    THISI = STEPI*THISR;
    N1 = ISTRT;
    for (ISTRT = N1; ISTRT <= N2; ISTRT++) {
        for (IL = 1; IL <= 4; IL++) {
            if ( DABS(WAVR(ISTRT+1,IL,1)) + DABS(WAVI(ISTRT+1,IL,1))
                >= THISI )  goto L480;
            WAVR(ISTRT+1,IL,1) = 0;
            WAVI(ISTRT+1,IL,1) = 0.;
        }
    }
    //
    // SHOULD WE REDUCE N1ADD?
    //
L480:
    if ( THISR*1.0e-10 > BIGNUM )  N1ADD = MAX0(N1ADD/2, 1 );
    if ( PWBGSW )  std::printf("     RESCALING: NEW ISTRT =%6d\n", ISTRT);
    //
    THISI = 1/THISR;
    N1 = N2+2;
    for (I = ISTRT; I <= N1; I++) {
        for (IL = 1; IL <= 4; IL++) {
            WAVR(I+1,IL,1) = THISI*WAVR(I+1,IL,1);
            WAVI(I+1,IL,1) = THISI*WAVI(I+1,IL,1);
        }
    }
    //
L490:
    N1 = N2+1;
    THSRSV = THISR;
    goto L410;
    //
    // END OF THE INTEGRATION
    //
L500:
    T2 = second();
    FTIME.TIMES[2] = FTIME.TIMES[2] + (float)(T2 - T1);
    NUMWAV = NUMWAV + 1;
    NUMLOP = NUMLOP + (NSTEP-NFIRST+1);
    //
    // EXTRACT S-MATRIX AND NORMALIZATION FROM U AND U AT NBAKCM
    // STEPS BACK FROM END.
    //
    LFORS = Z[KANDM.ISMATS[NWP]];
    TMATCH( WAVR_flat, WAVI_flat, CR_flat+1, CI_flat+1, ALLOC_base(LFORS), JP, NWP, NWAVF );
    // Note: TMATCH expects pointers to the start of the (2,2) arrays.
    // CR_flat+1 points to element [1] which is the first valid 1-based element.
    //
    // FIND NEXT STARTING POINT
    // FIND POINT WHERE ALL RENORMALIZED WAVEFUNCTION WILL BE < STEPI
    //
    AA1 = 0.;
    for (IL = 1; IL <= 4; IL++) {
        // Fortran: CR(IL,1) for IL=1..4 walks flat through the (2,2) array
        AA1 = AA1 + DABS(CR_flat[IL]) + DABS(CI_flat[IL]);
    }
    BB1 = STEPI/AA1;
    for (I = ISTRT; I <= NSTEP; I++) {
        for (IL = 1; IL <= 4; IL++) {
            if ( DABS(WAVR(I+1,IL,1)) + DABS(WAVI(I+1,IL,1)) >= BB1 )
                goto L630;
            WAVR(I+1,IL,1) = 0.;
            WAVI(I+1,IL,1) = 0.;
        }
    }
    //
L630:
    ISTRT = I;
    LASTZR = I+1;
    //
    // SET IT BACK SOME FOR LINKULES
    //
    if ( LNKBLK.LNKADR[6][3] != 0 )  I = I - 10;
    I = MAX0( I+1, 2 );
    WAVCOM.LASTNF[NWP] = I;
    WAVCOM.LASTL[NWP] = LL;
    BADSW = ( ( AA1 > BIGNUM*1.0e-10  &&  NFIRST > 2 )  ||
      AA1 < SMLNUM*1.0e-5 ) ? TRUE_F : FALSE_F;
    if ( BADSW )  std::printf("0**** WARNING, LARGE WAVEFUNCTION RENORMALIZATION:\n");
    if ( PWBGSW  ||  BADSW )  std::printf(" NWP, LL =%2d%4d  NFIRST =%4d  NEW NFIRST =%4d  RENORM =%11.3E\n",
        NWP, LL, NFIRST, I, AA1);
    NFIRST = ISTRT;
    //
    //********* REVIEW THESE 3 LINES IF/WHEN TENSOR FITS ARE ADDED.
    ALPHAR = 1.;
    ALPHAI = 0.;
    //CCC      IF ( NUMFIT .NE. 0 )  GO TO 880
    //
    // STORE NORMALIZED WAVE FUNCTIONS
    //
    // CR AND CI WILL NORMALIZE THE WAVEFUNCTION TO
    //   U(L) ---> ( (1+S)*F(L) + I (1-S)*G(L) ) / 2
    //
    // 750
    N2 = NSTEP+2;
    for (I = ISTRT; I <= N2; I++) {
        for (IL = 1; IL <= 4; IL++) {
            U1R_flat[IL] = WAVR(I+1,IL,1);
            U1I_flat[IL] = WAVI(I+1,IL,1);
        }
        for (IL = 1; IL <= 2; IL++) {
            for (ILAS = 1; ILAS <= 2; ILAS++) {
                WAVR(I+1,IL,ILAS)
                    = U1R(IL,1)*CR(1,ILAS) - U1I(IL,1)*CI(1,ILAS)
                    + U1R(IL,2)*CR(2,ILAS) - U1I(IL,2)*CI(2,ILAS);
                WAVI(I+1,IL,ILAS)
                    = U1R(IL,1)*CI(1,ILAS) + U1I(IL,1)*CR(1,ILAS)
                    + U1R(IL,2)*CI(2,ILAS) + U1I(IL,2)*CR(2,ILAS);
            }
        }
    }
    //
    // IF NECESSARY, CONTINUE THE INTEGRATION USING JUST COULOMB
    //
    NSTEP2 = WAVCOM.NSTP2S[NWP];
    if ( NSTEP2 <= NSTEP )  goto L800;
    //
    NSTEP = NSTEP+1;
    RT = NSTEP*STEPSZ;
    DLSQ[1] = DLSQ[1] * WAVCOM.XFACS[NWP][3];
    DLSQ[2] = DLSQ[2] * WAVCOM.XFACS[NWP][3];
    //
    for (I = NSTEP; I <= NSTEP2; I++) {
        RINV = 1./RT;
        for (IL = 1; IL <= 2; IL++) {
            D = 12. / ( WAVCOM.XFACS[NWP][1] + RINV*( WAVCOM.XFACS[NWP][2]
                + DLSQ[IL]*RINV ) );
            for (ILAS = 1; ILAS <= 2; ILAS++) {
                WAVR(I+3,IL,ILAS) = (-10.)*WAVR(I+2,IL,ILAS)
                    + WAVR(I,IL,ILAS) - WAVR(I+1,IL,ILAS);
                WAVI(I+3,IL,ILAS) = (-10.)*WAVI(I+2,IL,ILAS)
                    + WAVI(I,IL,ILAS) - WAVI(I+1,IL,ILAS);
                WAVR(I+1,IL,ILAS) = D*WAVR(I+3,IL,ILAS);
                WAVI(I+1,IL,ILAS) = D*WAVI(I+3,IL,ILAS);
            }
        }
        RT = RT + STEPSZ;
    }
    //
L800:
    if ( STANSW  ||  NUMPTS == 0 )  goto L950;
    //
    // NOW INTERPOLATE TO THE DESIRED GRID
    //
    STPINV = 1/WAVCOM.RSTEPS[NWP];
    IMX = NSTEP2 - 2;
    for (IR = 1; IR <= NUMPTS; IR++) {
        RT = RGRID[IR];
        //
        // FIND THE LOCATION OF R IN THE TABLE
        //
        RBYH = RT*STPINV;
        I = (int)(RBYH + 0.5);
        I = MAX0( 2, MIN0( I, IMX ) );
        P = RBYH - I;
        // USE 5-POINT LAGRANGIAN INTERPOLATION (ABRAMOWITZ & STEGUN 25.2.15)
        PS = P*P;
        X1 = P*(PS-1.)*OVER24;
        X2 = X1 + X1;
        X3 = X1*P;
        X4 = X2 + X2 - 0.5*P;
        X5 = X4*P;
        C1 = X3 - X2;
        C5 = X3 + X2;
        C3 = X5 - X3;
        C2 = X5 - X4;
        C4 = X5 + X4;
        C3 = C3 + C3 + 1.;
        II = IR;
        for (IL = 1; IL <= 4; IL++) {
            ALLOC4(LWAVRS+II) = (float)(C1*WAVR(I-1,IL,1) - C2*WAVR(I,IL,1)
                + C3*WAVR(I+1,IL,1) - C4*WAVR(I+2,IL,1) + C5*WAVR(I+3,IL,1));
            ALLOC4(LWAVIS+II) = (float)(C1*WAVI(I-1,IL,1) - C2*WAVI(I,IL,1)
                + C3*WAVI(I+1,IL,1) - C4*WAVI(I+2,IL,1) + C5*WAVI(I+3,IL,1));
            II = II + NUMPTS;
        }
        //
    }
    //
    // MOVE THE DESIRED WAVE FUNCTION FROM THE STORAGE ARRAY
    // TO THE OUTPUT ARRAY.
    //
L900:
    if ( STANSW  ||  NUMPTS == 0 )  goto L950;
    I = ( (L-LL)/2 + LAS-LL ) * NUMPTS;
    LWAVRS = LWAVRS + I;
    LWAVIS = LWAVIS + I;
    for (IR = 1; IR <= NUMPTS; IR++) {
        WAVER[IR] = ALLOC4(LWAVRS+IR);
        WAVEI[IR] = ALLOC4(LWAVIS+IR);
    }
    //
L950:
    FTIME.TIMES[3] = FTIME.TIMES[3] + (float)(second() - T2);
    return;

    #undef WAVR
    #undef WAVI
    #undef CR
    #undef CI
    #undef U1R
    #undef U1I
}
