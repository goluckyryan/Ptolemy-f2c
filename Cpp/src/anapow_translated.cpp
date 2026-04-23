// anapow_translated.cpp — Translated from Fortran source.f L582-914 (ANAPOW)
//                          and L27278-27637 (MUELCO)
// ANAPOW: compute and print analyzing powers (reaction or scattering)
// MUELCO: compute angle-independent Mueller matrix coefficients

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
// MUELCO (source.f L27278-27637)
// Computes angle-independent parts of the Mueller matrix for fixed K's.
//
// COEF(NQA,NQB,NQBIGA,NQBIGB,NSPL,NSPL) — 6D Fortran column-major array,
// passed as double* (0-based flattened). Use macro below.
// COEP(NQA,NQB), COET(NQBIGA,NQBIGB) — 2D scratch arrays, column-major.
// COEX(KMX+1), TERMK(KMX+1) — 1D scratch, 0-based index K (Fortran K+1).
// ITOC(4,NSPL) — column-major integer array, passed as int* 0-based.
// ============================================================================

// 6-dimensional COEF access macro: Fortran COEF(IQA,IQB,IQBIGA,IQBIGB,KOFF,KOFFP)
// dims = (NQA,NQB,NQBIGA,NQBIGB,NSPL,NSPL), 1-based all indices
#define COEF6(iq1,iq2,iqbiga,iqbigb,koff,koffp) \
    COEF[ (int)(koffp-1)*(int)NSPL*(int)NQBIGB*(int)NQBIGA*(int)NQB*(int)NQA \
        + (int)(koff-1)*(int)NQBIGB*(int)NQBIGA*(int)NQB*(int)NQA \
        + (int)(iqbigb-1)*(int)NQBIGA*(int)NQB*(int)NQA \
        + (int)(iqbiga-1)*(int)NQB*(int)NQA \
        + (int)(iq2-1)*(int)NQA \
        + (int)(iq1-1) ]

// ITOC(I,J) with Fortran dims (4,NSPL), passed as 0-based int*
#define ITOCF(i,j) ITOC[(int)(j-1)*4 + (int)(i-1)]

void MUELCO(int KA, int KB, int KBIGA, int KBIGB, int JA, int JB,
            int JBIGA, int JBIGB, int NQA, int NQB, int NQBIGA, int NQBIGB,
            int NSPL, int IPARIT, int* ITOC, double* COEF, double* COEP,
            double* COET, double* COEX, double* TERMK, int PBUGSW)
{
    // Fortran: COEP(NQA,NQB), COET(NQBIGA,NQBIGB) column-major
    // COEP[col*NQA + row], COET[col*NQBIGA + row] (0-based)
    #define COEPF(iq1,iq2) COEP[(int)(iq2-1)*(int)NQA + (int)(iq1-1)]
    #define COETF(iqbiga,iqbigb) COET[(int)(iqbigb-1)*(int)NQBIGA + (int)(iqbiga-1)]

    int QAS = KA - NQA + 1;
    int QBS = KB - NQB + 1;
    int QBIGAS = KBIGA - NQBIGA + 1;
    int QBIGBS = KBIGB - NQBIGB + 1;
    double PARITY_d = (double)IPARIT;

    // Clear COEF (N elements)
    int N = NSPL * NSPL * NQA * NQB * NQBIGA * NQBIGB;
    if (PBUGSW)
        std::printf("\nMUELCO... %12d%12d%12d%12d%12d%12d%12d%12d%12d%12d%12d%12d%12d%12d%12d\n",
                    KA, KB, KBIGA, KBIGB, JA, JB, JBIGA, JBIGB,
                    NQA, NQB, NQBIGA, NQBIGB, NSPL, IPARIT, N);
    for (int I = 0; I < N; I++) COEF[I] = 0.0;

    DUMMY2();  // gets factorial tables into core

    // Outer loops over JP and JT groups (grouped by JP,JT in ITOC)
    int JTPRV = -1, JTPPRV = -1, KTPRV = -1, JPPRV = -1, JPPPRV = -1, KPPRV = -1;
    double C9JT = 0.0, C9JP = 0.0;

    int KOFE = 0;
    while (true) {
        // Find next valid KOFS
        int KOFS = KOFE + 1;
        while (true) {
            KOFE = KOFS;
            if (KOFS > NSPL) return;
            if (ITOCF(4,KOFS) >= 0) break;
            KOFS++;
        }
        int JT = ITOCF(4,KOFS);
        int JP = ITOCF(3,KOFS);
        int LXMN = ITOCF(2,KOFS);

        // Find last element with this JP,JT
        KOFE = KOFS;
        while (true) {
            KOFE++;
            if (KOFE > NSPL) break;
            if (ITOCF(4,KOFE) != JT || ITOCF(3,KOFE) != JP) break;
        }
        KOFE--;
        int LXMX = ITOCF(2,KOFE);

        // Inner: loop over JP', JT' groups
        int KOFPE = 0;
        while (true) {
            // Find next valid KOFPS
            int KOFPS = KOFPE + 1;
            while (true) {
                KOFPE = KOFPS;
                if (KOFPS > NSPL) goto next_outer;
                if (ITOCF(4,KOFPS) >= 0) break;
                KOFPS++;
            }
            int JTP = ITOCF(4,KOFPS);
            int JPP = ITOCF(3,KOFPS);
            int LXPMN = ITOCF(2,KOFPS);
            KOFPE = KOFPS;
            while (true) {
                KOFPE++;
                if (KOFPE > NSPL) break;
                if (ITOCF(4,KOFPE) != JTP || ITOCF(3,KOFPE) != JPP) break;
            }
            KOFPE--;
            int LXPMX = ITOCF(2,KOFPE);

            // KT loop
            int KTMN = DMAX1((double)IABS(KBIGA-KBIGB), (double)IABS(JT-JTP)/2);
            int KTMX = DMIN1((double)(KBIGA+KBIGB), (double)(JT+JTP)/2);
            if (KTMN > KTMX) continue;

            for (int KT = KTMN; KT <= KTMX; KT++) {
                int QTMIN = DMAX1(-(double)KT, (double)(QBIGBS-KBIGA));
                int QTMAX = DMIN1((double)KT, (double)(KBIGB-QBIGAS));
                if (QTMIN > QTMAX) continue;

                // Compute COET if KT/JT/JTP changed
                if (KT != KTPRV || JTP != JTPPRV || JT != JTPRV) {
                    double fac = (double)(JBIGA+1)*(JBIGB+1)*(2*KBIGA+1)*(JT+1)*(JTP+1)*(2*KT+1);
                    C9JT = DSQRT(fac) * WIG9J(JBIGA, JT, JBIGB,
                                               2*KBIGA, 2*KT, 2*KBIGB, JBIGA, JTP, JBIGB);
                    if (PBUGSW)
                        std::printf(" 9J AT 203:%8d%8d%8d%16.6g\n", JT, KT, JTP, C9JT);
                    KTPRV = KT; JTPPRV = JTP; JTPRV = JT;
                }
                if (C9JT == 0.0) continue;

                // Fill COET
                for (int IQBIGA = 1; IQBIGA <= NQBIGA; IQBIGA++) {
                    int QBIGA = IQBIGA + QBIGAS - 1;
                    for (int IQBIGB = 1; IQBIGB <= NQBIGB; IQBIGB++) {
                        int QBIGB_v = IQBIGB + QBIGBS - 1;
                        COETF(IQBIGA,IQBIGB) = 0.0;
                        if (IABS(QBIGB_v - QBIGA) <= KT)
                            COETF(IQBIGA,IQBIGB) = C9JT * CLEBSH(2*KBIGA, 2*KT,
                                2*QBIGA, 2*(QBIGB_v-QBIGA), 2*KBIGB, 2*QBIGB_v);
                    }
                }

                // KP loop
                int KPMN = DMAX1((double)IABS(KB-KA), (double)IABS(JP-JPP)/2);
                int KPMX = DMIN1((double)(KB+KA), (double)(JP+JPP)/2);
                if (KPMN > KPMX) continue;

                for (int KP = KPMN; KP <= KPMX; KP++) {
                    int QPMIN = DMAX1(-(double)KP, (double)(QAS-KB));
                    int QPMAX = DMIN1((double)KP, (double)(KA-QBS));
                    if (QPMIN > QPMAX) continue;

                    // Compute COEP if KP/JP/JPP changed
                    if (KP != KPPRV || JPP != JPPPRV || JP != JPPRV) {
                        double fac2 = (double)(JA+1)*(JB+1)*(2*KB+1)*(JP+1)*(JPP+1);
                        C9JP = (double)(2*KP+1) * DSQRT(fac2)
                             * WIG9J(JB, JP, JA, 2*KB, 2*KP, 2*KA, JB, JPP, JA);
                        if (PBUGSW)
                            std::printf(" 9J AT 253:%8d%8d%8d%16.6g\n", JP, KP, JPP, C9JP);
                        KPPRV = KP; JPPPRV = JPP; JPPRV = JP;
                    }
                    if (C9JP == 0.0) continue;

                    // Fill COEP
                    for (int IQA = 1; IQA <= NQA; IQA++) {
                        int QA_v = IQA + QAS - 1;
                        for (int IQB = 1; IQB <= NQB; IQB++) {
                            int QB_v = IQB + QBS - 1;
                            COEPF(IQA,IQB) = 0.0;
                            if (IABS(QA_v - QB_v) <= KP)
                                COEPF(IQA,IQB) = C9JP * CLEBSH(2*KB, 2*KP,
                                    2*QB_v, 2*(QA_v-QB_v), 2*KA, 2*QA_v);
                        }
                    }

                    // LX, LXP loops with KOFFS/KOFFPS tracking
                    int KOFFS = KOFS;
                    for (int LX = LXMN; LX <= LXMX; LX++) {
                        int MXMN = -LX;
                        int MXZ = (ITOCF(1,KOFFS) + LX + 1) / 2;

                        int KOFFPS = KOFPS;
                        for (int LXP = LXPMN; LXP <= LXPMX; LXP++) {
                            int MXPZ = (ITOCF(1,KOFFPS) + LXP + 1) / 2;

                            // Compute 3rd 9-J for all K
                            int KMN = DMAX1(DMAX1((double)IABS(LX-LXP), (double)IABS(KP-KT)),
                                            (double)(QTMIN-QPMAX));
                            int KMX = DMIN1((double)(LX+LXP), (double)(KP+KT));
                            if (KMN > KMX) goto next_lxp;

                            {
                                int QMIN_c = DMAX1(-(double)KMX, (double)(QTMIN-QPMAX));
                                int QMAX_c = DMIN1((double)KMX, (double)(QTMAX-QPMIN));
                                if (QMIN_c > QMAX_c) goto next_lxp;

                                double T_lx = (double)(2*LX + 1);
                                for (int K = KMN; K <= KMX; K++) {
                                    COEX[K] = (double)(2*K+1) * DSQRT(T_lx)
                                        * WIG9J(2*LX, JP, JT, 2*K, 2*KP, 2*KT, 2*LXP, JPP, JTP);
                                    if (PBUGSW)
                                        std::printf(" 9J:%5d%5d%5d%5d%5d%5d%5d%5d%5d%16.6g\n",
                                                    KT, KP, LX, LXP, K, JP, JT, JPP, JTP, COEX[K]);
                                }

                                // MX, MXP loops
                                for (int MX = MXMN; MX <= LX; MX++) {
                                    int KOFF = KOFFS + IABS(MX) - MXZ;
                                    if (KOFF < KOFFS) continue;
                                    double FACT;
                                    if (MX >= 0) {
                                        FACT = 1.0;
                                    } else {
                                        FACT = PARITY_d;
                                        if ((LX + MX) % 2 != 0) FACT = -FACT;
                                    }

                                    int MXPMN = DMAX1(-(double)LXP, (double)(MX+QMIN_c));
                                    int MXPMX = DMIN1((double)LXP, (double)(MX+QMAX_c));
                                    if (MXPMN > MXPMX) continue;

                                    for (int MXP = MXPMN; MXP <= MXPMX; MXP++) {
                                        int KOFFP = KOFFPS + IABS(MXP) - MXPZ;
                                        if (KOFFP < KOFFPS) continue;
                                        double FAC = FACT;
                                        if (MXP < 0) {
                                            FAC *= PARITY_d;
                                            if ((LXP + MXP) % 2 != 0) FAC = -FAC;
                                        }

                                        int Q = MXP - MX;

                                        // Compute TERMK for all K
                                        int KMIN = DMAX1((double)KMN, (double)IABS(Q));
                                        for (int K = KMIN; K <= KMX; K++) {
                                            TERMK[K] = FAC * COEX[K]
                                                * CLEBSH(2*LX, 2*K, 2*MX, 2*Q, 2*LXP, 2*MXP);
                                        }

                                        // QP loop
                                        int QPMN = DMAX1((double)QPMIN, (double)(QTMIN-Q));
                                        int QPMX = DMIN1((double)QPMAX, (double)(QTMAX-Q));
                                        for (int QP = QPMN; QP <= QPMX; QP++) {
                                            int QT = Q + QP;

                                            // Sum over K
                                            double SUMK = 0.0;
                                            for (int K = KMIN; K <= KMX; K++)
                                                SUMK += TERMK[K] * CLEBSH(2*K, 2*KP, 2*Q, 2*QP, 2*KT, 2*QT);

                                            // QBIGA, QB loops
                                            int QAMN = DMAX1((double)QBIGAS, (double)(QBIGBS-QT));
                                            int QAMX = DMIN1((double)KBIGA, (double)(KBIGB-QT));
                                            int QBMN = DMAX1((double)QBS, (double)(QAS-QP));
                                            int QBMX = DMIN1((double)KB, (double)(KA-QP));

                                            for (int QBIGA_v = QAMN; QBIGA_v <= QAMX; QBIGA_v++) {
                                                int IQBIGA = QBIGA_v - QBIGAS + 1;
                                                int QBIGB_v = QBIGA_v + QT;
                                                int IQBIGB = QBIGB_v - QBIGBS + 1;
                                                double TERMT = SUMK * COETF(IQBIGA,IQBIGB);
                                                if (TERMT == 0.0) continue;

                                                for (int QB_v = QBMN; QB_v <= QBMX; QB_v++) {
                                                    int IQB = QB_v - QBS + 1;
                                                    int QA_v = QB_v + QP;
                                                    int IQA = QA_v - QAS + 1;
                                                    double TERM = TERMT * COEPF(IQA,IQB);
                                                    COEF6(IQA,IQB,IQBIGA,IQBIGB,KOFF,KOFFP)
                                                        += TERM;
                                                    if (PBUGSW)
                                                        std::printf(
                                                            " COEF:%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%5d%15.5g%15.5g%15.5g\n",
                                                            QB_v, QBIGA_v, QP, MXP, MX,
                                                            IQA, IQB, KOFFP, KOFF, IQBIGA, IQBIGB,
                                                            SUMK, TERMT, TERM);
                                                }
                                            }
                                        } // QP loop
                                    } // MXP loop
                                } // MX loop
                            } // block for QMIN_c/QMAX_c

                        next_lxp:
                            KOFFPS = KOFFPS + LXP - MXPZ + 1;
                        } // LXP loop

                        KOFFS = KOFFS + LX - MXZ + 1;
                    } // LX loop
                } // KP loop
            } // KT loop
        } // inner JP'/JT' while
    next_outer:;
    } // outer JP/JT while

    #undef COEPF
    #undef COETF
}

#undef ITOCF
#undef COEF6


// ============================================================================
// ANAPOW (source.f L582-914)
// Computes and prints analyzing powers.
// ============================================================================

void ANAPOW(double ANGMIN, double ANGMAX, double ANGSTP, int JA, int JB,
            int JBIGA, int JBIGB, int NSPL, int IPARIT, int PRNTSW,
            int PBUGSW, const char* APNAME, double ELAB, const char* CHNM,
            int IF_arg, int ITOCS, int& IAPOW)
{
    // Common block accesses
    auto& LABANG    = SWITCH.LABANG;
    auto& FACFR4    = ALLOCS.FACFR4;
    double SMLNUM   = CNSTNT.SMLNUM;
    auto* Z_ptr     = LOCPTRS.Z;

    // LABCM: 4-char label
    const char* LABCM = (LABANG != 0) ? "LAB." : "C.M.";

    int NUMANG = (int)((ANGMAX - ANGMIN) / ANGSTP + 1.5);
    int NUMAP  = ((JA + 2) * JA) / 2;

    IAPOW = NALLOC(APNAME, NUMANG * NUMAP);
    INIT8(IAPOW, 0.0);

    // NQ for K=1 (first iteration of outer K loop uses NQ = JA+1-QMIN)
    // Allocate scratch arrays
    int NQ = JA + (JA + 1) % 2;   // = JA + MOD(JA+1,2)
    int N  = NSPL * NSPL * NQ;
    int IAPCO = NALLOC("APCO    ", N);
    N = NQ + 2 * JA + 3;
    int IAPSCR = NALLOC("APSCR   ", N);
    int IFTF  = NALLOC("FTIMESF ", NUMANG);
    int ICROS = NALLOC("TEMPCROS", NUMANG);
    INIT8(ICROS, SMLNUM);
    int KCROS = MIN0(JA, 2);

    // Get allocator pointers (SETLOG equivalent)
    SETLOG(1);

    // Fortran: LCOEP = Z(IAPSCR); LCOET = LCOEP+NQ; etc.
    // In C++: 1-based offsets into ALLOC array
    int LCOEP  = Z_ptr[IAPSCR];
    int LCOET  = LCOEP + NQ;
    int LCOEX  = LCOET + 1;
    int LTERMK = LCOEX + JA + 1;
    int LTOCS  = Z_ptr[ITOCS] * FACFR4 - FACFR4;  // ILLOC offset (1-based)
    int LF     = Z_ptr[IF_arg] - 2 * NUMANG - 2;   // Fortran: LF = Z(IF) - 2*NUMANG - 2
    int LFTF   = Z_ptr[IFTF] - 1;
    int LCROS  = Z_ptr[ICROS] - 1;

    // Outer K loop
    for (int K = 1; K <= JA; K++) {
        int QMIN_k = K % 2;
        NQ = K + 1 - QMIN_k;

        // ASSIGN 180/170 equivalent: use REAL part for even K (QMIN=0), IMAG for odd K (QMIN=1)
        bool USEREAL = (QMIN_k == 0);

        // Fortran: LAPOW = Z(IAPOW) - 1 + NUMANG*((K+1)*(K-1)+QMIN)/2 - QMIN
        //   Actually: LAPOW = Z(IAPOW) - 1 + NUMANG * LAPOW_offset
        // where LAPOW_offset = ((K+1)*(K-1) + QMIN_k)/2 - QMIN_k
        int LAPOW_base = Z_ptr[IAPOW] - 1 + NUMANG * (((K+1)*(K-1) + QMIN_k)/2 - QMIN_k);

        // Call MUELCO to get angle-independent coefficients
        // ITOC → ILLOC starting at LTOCS+1 (1-based), passed as int* 0-based
        // COEF → ALLOC starting at Z(IAPCO) (1-based)
        MUELCO(K, 0, 0, 0, JA, JB, JBIGA, JBIGB,
               NQ, 1, 1, 1, NSPL, IPARIT,
               &ILLOC(LTOCS + 1),      // ITOC: 0-based from &ILLOC(LTOCS+1)
               &ALLOC(Z_ptr[IAPCO]),   // COEF
               &ALLOC(LCOEP),          // COEP
               &ALLOC(LCOET),          // COET
               &ALLOC(LCOEX),          // COEX
               &ALLOC(LTERMK),         // TERMK
               PBUGSW);

        // Loops over KOFFP (prime) and KOFF
        for (int KOFFP = 1; KOFFP <= NSPL; KOFFP++) {
            if (ILLOC(LTOCS + 4*KOFFP) < 0) continue;
            for (int KOFF = 1; KOFF <= NSPL; KOFF++) {
                if (ILLOC(LTOCS + 4*KOFF) < 0) continue;

                bool SKIPSW = false;

                // Fortran: LAPCO = Z(IAPCO) - QMIN_k + NQ*(KOFF-1 + NSPL*(KOFFP-1))
                int LAPCO = Z_ptr[IAPCO] - QMIN_k + NQ * ((KOFF-1) + NSPL*(KOFFP-1));

                // For diagonal K=KCROS: accumulate cross section
                if (KOFF == KOFFP && K == KCROS) {
                    double T_cs = 2.0;
                    int MX = (ILLOC(LTOCS + 4*KOFF - 3) + ILLOC(LTOCS + 4*KOFF - 2) + 1) / 2;
                    if (MX == 0) T_cs = 1.0;
                    for (int NANG = 1; NANG <= NUMANG; NANG++) {
                        int I = (LF + 2*NUMANG*KOFF) + 2*NANG;
                        ALLOC(LFTF + NANG) = ALLOC(I)*ALLOC(I) + ALLOC(I+1)*ALLOC(I+1);
                        ALLOC(LCROS + NANG) += T_cs * ALLOC(LFTF + NANG);
                    }
                }

                // Print APCO coefficients if debugging
                if (PBUGSW) {
                    std::printf(" APCO%2d%2d", KOFF, KOFFP);
                    for (int Q = QMIN_k; Q <= K; Q++)
                        std::printf("%12.5g", ALLOC(LAPCO + Q));
                    std::printf("\n");
                }

                // Loop over Q, looking for non-zero coefficient
                for (int Q = QMIN_k; Q <= K; Q++) {
                    if (ALLOC(LAPCO + Q) == 0.0) continue;

                    // Compute F*F'* products at all angles if needed
                    if (!SKIPSW) {
                        for (int NANG = 1; NANG <= NUMANG; NANG++) {
                            int I  = (LF + 2*NUMANG*KOFF)  + 2*NANG;
                            int II = (LF + 2*NUMANG*KOFFP) + 2*NANG;
                            if (USEREAL)
                                // Real part: Re(F * F'*) = Re(F)*Re(F') + Im(F)*Im(F')
                                ALLOC(LFTF + NANG) = ALLOC(I)*ALLOC(II) + ALLOC(I+1)*ALLOC(II+1);
                            else
                                // Imag part: Im(F'* * F) = Re(F')*Im(F) - Re(F)*Im(F')
                                ALLOC(LFTF + NANG) = ALLOC(II)*ALLOC(I+1) - ALLOC(I)*ALLOC(II+1);
                        }
                        SKIPSW = true;
                    }

                    // Accumulate analyzing powers at this Q
                    int LAP = LAPOW_base + NUMANG * Q;
                    double T_ap = ALLOC(LAPCO + Q);
                    for (int NANG = 1; NANG <= NUMANG; NANG++)
                        ALLOC(LAP + NANG) += T_ap * ALLOC(LFTF + NANG);
                } // Q loop
            } // KOFF loop
        } // KOFFP loop
    } // K loop

    // Divide analyzing powers by cross section
    {
        int I = Z_ptr[IAPOW];
        for (int II = 1; II <= NUMAP; II++) {
            for (int NANG = 1; NANG <= NUMANG; NANG++) {
                ALLOC(I) /= ALLOC(LCROS + NANG);
                I++;
            }
        }
    }

    // Free temporary arrays
    Z_ptr[IAPCO]  = -Z_ptr[IAPCO];
    Z_ptr[IAPSCR] = -Z_ptr[IAPSCR];
    Z_ptr[IFTF]   = -Z_ptr[IFTF];
    Z_ptr[ICROS]  = -Z_ptr[ICROS];

    if (!PRNTSW) return;

    // Determine block size for blank-line spacing
    double ANGBLK = 300.0;
    if (NUMANG > 40) {
        ANGBLK = 10.0;
        if (ANGSTP < 0.99) {
            ANGBLK = 5.0;
            if (ANGSTP <= 0.201) ANGBLK = 1.0;
        }
    }
    int NUMBLK = (int)(ANGBLK / ANGSTP + 0.5);
    const int LINEMX = 57;

    // Build column label lists
    // Max columns = 12 (from FORMAT '12( 6X, A1, 'T', 2I1 )')
    int LSTK[13] = {}, LSTQ[13] = {};
    char LSTI[13] = {};   // ' ' or 'i'

    int IIS = Z_ptr[IAPOW] - 1;
    int KS  = 1;

    // Outer print loop: one page per group of K values
    while (KS <= JA) {
        int KE = MIN0(JA, 4);
        if (KS == 5) KE = MIN0(JA, 6);
        if (KS  > 5) KE = KS;

        // Build column labels for K = KS..KE
        int Ncols = 0;
        for (int K = KS; K <= KE; K++) {
            int QMIN_k = K % 2;
            for (int Q = QMIN_k; Q <= K; Q++) {
                Ncols++;
                LSTK[Ncols] = K;
                LSTQ[Ncols] = Q;
                LSTI[Ncols] = (QMIN_k != 0) ? 'i' : ' ';
            }
        }

        int IIE  = IIS + (Ncols - 1) * NUMANG;
        int LINECT = 1000;  // Force header on first line

        // Loop through angles
        for (int NANG = 1; NANG <= NUMANG; NANG++) {
            double ANGLE = ANGMIN + ANGSTP * (NANG - 1);
            double ANGLE_print = ANGLE;
            if (LABANG != 0 && ANGLE > KANDM.ABAR)
                ANGLE_print = 2.0 * KANDM.ABAR - ANGLE;

            // Print header if needed
            if (LINECT >= LINEMX) {
                // FORMAT 603: new page header
                // '1' = new page (form feed equivalent)
                std::printf("1%54sP T O L E M Y\n", "");
                std::printf("%10sANALYZING POWERS FOR THE  >>>> %-8s <<<<  CHANNEL\n",
                            "", CHNM);
                // '0' = blank line
                std::printf("0%.45s%7.2f MEV     %.65s\n",
                            &HEDCOM.REACT[1], ELAB, &HEDCOM.HEADER[1]);
                // FORMAT 613: column labels
                // '0' = blank line
                std::printf("0  %-4.4s ", LABCM);
                for (int i = 1; i <= Ncols; i++)
                    std::printf("      %cT%d%d", LSTI[i], LSTK[i], LSTQ[i]);
                std::printf("\n");
                // FORMAT 614
                std::printf("  ANGLE\n\n");
                LINECT = 8;
            }

            // FORMAT 633: 1X, F6.2, 2X, 12F10.5
            std::printf(" %6.2f  ", ANGLE_print);
            for (int II = IIS; II <= IIE; II += NUMANG)
                std::printf("%10.5f", ALLOC(II + NANG));
            std::printf("\n");
            LINECT++;

            // Blank line after ANGBLK-th angle
            double angle_mod = DMOD(ANGLE + ANGSTP + 1.e-7, ANGBLK);
            if (angle_mod <= 1.e-5) {
                std::printf("\n");  // FORMAT 633 with no values (blank line)
                LINECT++;
                if (LINECT + NUMBLK > LINEMX) LINECT = 1000;
            }
        } // angle loop

        IIS = IIE + NUMANG;
        KS  = KE + 1;
    } // page group loop
}
