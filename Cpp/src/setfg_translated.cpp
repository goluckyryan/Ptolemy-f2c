// SETFG — translated from source.f lines 32481-32574
// Computes all F and G Coulomb functions for coupled channels

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include "ptolemy_alloc.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_forward.h"
#include "ptolemy_io.h"
#include <cstdio>

void SETFG(int NCHNDF, int MCHNDF, int* CHNDEF, float* RCHNDF,
           int NBASDF, int MBASDF, int* BASDEF, int MCHNVL, double* CHNVAL,
           int NUMLIS, double* LIS_arr, double* FG_arr,
           int L1, int L2P1, double* FGWORK, double R1, double R2, int& IRET)
{
    // CHNDEF(MCHNDF, NCHNDF) — integer 2D
    // RCHNDF(MCHNDF, NCHNDF) — real*4 2D
    // BASDEF(MBASDF, NBASDF) — integer 2D
    // CHNVAL(MCHNVL, NCHNDF) — real*8 2D
    // LIS(NUMLIS) — integer stored as double (from ALLOC)
    // FG(4, NBASDF, NUMLIS) — real*4 3D
    // FGWORK(L2P1, 6) — real*8 2D

    #define CHNDEF_F(i,j) CHNDEF[((j)-1)*MCHNDF + (i)-1]
    #define RCHNDF_F(i,j) RCHNDF[((j)-1)*MCHNDF + (i)-1]
    #define BASDEF_F(i,j) BASDEF[((j)-1)*MBASDF + (i)-1]
    #define CHNVAL_F(i,j) CHNVAL[((j)-1)*MCHNVL + (i)-1]
    #define FGWORK_F(i,j) FGWORK[((j)-1)*L2P1 + (i)-1]
    // FG is real*4, accessed as float*
    float* FG_f = reinterpret_cast<float*>(FG_arr);
    #define FG_F(i,ib,ili) FG_f[((ili)-1)*NBASDF*4 + ((ib)-1)*4 + (i)-1]

    int L2 = L2P1 - 1;
    int IFLAG = 0;

    // Loop through the channels array
    for (int NC = 1; NC <= NCHNDF; NC++) {
        double AK = CHNVAL_F(1, NC);
        double ETA = CHNVAL_F(2, NC);
        int IB1 = CHNDEF_F(8, NC);
        int IB2 = CHNDEF_F(9, NC) + IB1 - 1;

        double RHO = AK * R1;
        RCWFN(RHO, ETA, L1, L2, &FGWORK_F(1, 1),
              &FGWORK_F(1, 5), &FGWORK_F(1, 2), &FGWORK_F(1, 6),
              1.0e-12, IFLAG);
        if (IFLAG != 0) goto L990;

        RHO = AK * R2;
        RCWFN(RHO, ETA, L1, L2, &FGWORK_F(1, 3),
              &FGWORK_F(1, 5), &FGWORK_F(1, 4), &FGWORK_F(1, 6),
              1.0e-12, IFLAG);
        if (IFLAG != 0) goto L990;

        // Loop through the L-J's for this channel
        for (int IB = IB1; IB <= IB2; IB++) {
            int LMJ = BASDEF_F(3, IB);

            // Loop through the required J's and store F and G
            for (int ILI = 1; ILI <= NUMLIS; ILI++) {
                // Note: J and not 2*J is stored for now
                // LIS is stored in ALLOC as doubles, but contains integers
                int J = 2 * (int)LIS_arr[ILI - 1];
                int LL = LMJ + J;
                if (LL < 0) continue;
                LL = LL >> 1;  // ishft(L, -1)

                for (int II = 1; II <= 4; II++) {
                    FG_F(II, IB, ILI) = (float)FGWORK_F(LL + 1, II);
                }
            }
        }
    }

    IRET = 1;
    return;

    // Error in RCWFN
L990:
    std::printf("\n**** RCWFN ERROR IN SETFG:\n%8d%20.10g%20.10g%8d%8d\n",
                IFLAG, 0.0, 0.0, L1, L2);
    IRET = 0;
    return;

    #undef CHNDEF_F
    #undef RCHNDF_F
    #undef BASDEF_F
    #undef CHNVAL_F
    #undef FGWORK_F
    #undef FG_F
}
