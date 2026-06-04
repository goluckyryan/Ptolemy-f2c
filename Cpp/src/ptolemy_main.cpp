//
// ptolemy_main.cpp — PROGRAM PTOLEMY
//
// Main program for the DWBA code.
// Verbatim transliteration of ptolemy.f (251 lines).
//
// This acts as a switcher for CONTRL.  CONTRL is called and reads input.
// Each time a substantial calculation is to be done, CONTRL sets IGOTO
// and JGOTO and returns here.  IGOTO is then used to switch to the
// appropriate routine.
//
// 1/18/75 — made from CONTRL — S. Pieper
//

#include "ptolemy_types.h"
#include "ptolemy_alloc.h"
#include "ptolemy_commons.h"
#include "ptolemy_forward.h"
#include "ptolemy_io.h"
#include "ptolemy_intrinsics.h"

#include <cstdio>
#include <cstdlib>

int main()
{
    // Local variables
    int IGOTO, JGOTO, IPRM1, IPRM2, IPRM3;
    double PRM1;
    bool ONESW;

    // Initialize with all the defaults
    JGOTO = 3;

    // ---------------------------------------------------------------
    // Main loop: start or continue the input process
    // ---------------------------------------------------------------
label_90:
    CONTRL(IGOTO, JGOTO, IPRM1, IPRM2, IPRM3, PRM1);

    if (IGOTO == 0) {
        // Normal termination
        return 0;   // Fortran STOP
    }

    // What is to be done next
    // IGOTO = 4 available for recycling

    ONESW = true;

    switch (IGOTO) {
    case 1:  goto label_100;
    case 2:  goto label_200;
    case 3:  goto label_300;
    case 4:  goto label_90;
    case 5:  goto label_500;
    case 6:  goto label_600;
    case 7:  goto label_620;
    case 8:  goto label_650;
    case 9:  goto label_670;
    case 10: goto label_700;
    case 11: goto label_720;
    case 12: goto label_740;
    case 13: goto label_800;
    case 14: goto label_850;
    case 15: goto label_870;
    default: goto label_90;
    }

    // ---------------------------------------------------------------
    // IGOTO = 1: after PROJECTILE or TARGET
    // ---------------------------------------------------------------
label_100:
    BOUND(IPRM1, IPRM2, IPRM3);
    goto label_90;

    // ---------------------------------------------------------------
    // IGOTO = 2: after FIT or FITELASTIC — fitting feature removed
    // ---------------------------------------------------------------
label_200:
    printf("\n**** FIT KEYWORD NOT AVAILABLE IN this VERSION (fitting feature removed)\n");
    IPRM1 = 0;
    goto label_90;

    // ---------------------------------------------------------------
    // IGOTO = 3: after CALLTSO
    // ---------------------------------------------------------------
label_300:
    printf("\n**** TSO KEYWORD NOT AVAILABLE IN this VERSION\n");
    goto label_90;

    // ---------------------------------------------------------------
    // IGOTO = 5: after INCOMING or OUTGOING
    // ---------------------------------------------------------------
label_500:
    WAVSET(IPRM1, false, false);
    goto label_90;

    // ---------------------------------------------------------------
    // GRIDSETUP
    // IGOTO = 6: after GRIDSETUP
    // ---------------------------------------------------------------
label_600:
    PRBPRT(IPRM1);
    if (IPRM1 == 0) goto label_90;
    if (SWITCH.PROBLM != 24) goto label_605;

    // Coupled channels calls
    BASLBL();
    // Note: Fortran main passes IPRM1 to BASLBL but subroutine takes no args.
    // IPRM1 retains its value from PRBPRT.
    if (IPRM1 == 0) goto label_90;
    BASCPL(IPRM1);
    if (IPRM1 == 0) goto label_90;

label_605:
    GETSCT(IPRM1);
    if (IPRM1 == 0) goto label_90;
    if (INTRNL.ISTRIP != 0) GRDSET(IPRM1);

    // Following done for all calculations
    if (INTRNL.ISTRIP == 0) INGRST(IPRM1);
    if (IPRM1 == 0) goto label_90;
    ANGSET(IPRM1);
    if (IPRM1 == 0) goto label_90;
    if (INTRNL.ISTRIP == 0) goto label_615;
    if (ONESW) goto label_90;
    goto label_620;

    // Setup Coulomb integrals for inelastic and CC
label_615:
    COULST(IPRM1);

    if (IPRM1 == 0 || ONESW) goto label_90;

    // ---------------------------------------------------------------
    // IGOTO = 7: after RADIALINTEGRALS
    // ---------------------------------------------------------------
label_620:
    if (SWITCH.PROBLM == 20) INELDC(IPRM1);
    if (SWITCH.PROBLM == 22) INRDIN(IPRM1);
    if (SWITCH.PROBLM == 24) COUPLN(IPRM1);
    if (IPRM1 == 0 || ONESW) goto label_90;

    // ---------------------------------------------------------------
    // IGOTO = 8: after LINTERPOLATION
    // ---------------------------------------------------------------
label_650:
    IPRM1 = 1;
    if (SWITCH.PROBLM == 24) FFACST(IPRM1);
    if (IPRM1 == 0) goto label_90;
    LINTRP();
    if (ONESW) goto label_90;

    // ---------------------------------------------------------------
    // IGOTO = 9: after CROSSSECTION
    // ---------------------------------------------------------------
label_670:
    XSECTN(IPRM1);
    goto label_90;

    // ---------------------------------------------------------------
    // IGOTO = 10, 11, 12: after DOTEMP1, DOTEMP2, DOTEMP3
    // ---------------------------------------------------------------
label_700:
    TEMP1();
    goto label_90;

label_720:
    TEMP2();
    goto label_90;

label_740:
    TEMP3();
    goto label_90;

    // ---------------------------------------------------------------
    // IGOTO = 13: after SCATTERING
    // ---------------------------------------------------------------
label_800:
    WAVEF(IPRM1);
    goto label_90;

    // ---------------------------------------------------------------
    // IGOTO = 14: after all channels set up —
    // full transfer or inelastic calculation
    // ---------------------------------------------------------------
label_850:
    ONESW = false;
    goto label_600;

    // ---------------------------------------------------------------
    // IGOTO = 15: after CROSSSECTION if LINTRPOLATION not done —
    // both LINTRP and XSECTN
    // ---------------------------------------------------------------
label_870:
    ONESW = false;
    goto label_650;
}
