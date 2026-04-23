// clrchn_translated.cpp — CLRCHN subroutine
// Translated from source.f lines 9004-9347
// Clears channel-specific parameters between channel calculations

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include <cstdio>

void CLRCHN(int ICHAN)
{
    //
    // CLEAR THE DEFINITIONS OF CHANNEL VARIABLES ( P, T, AND POTS )
    //
    // THIS ROUTINE IS CALLED AFTER EACH BOUNDSTATE AND OPTICAL
    // MODEL SCATTERING STATE.  IT SETS CHANNEL VARIABLES THAT ARE
    // NOT LIKELY TO BE SIGNIFICANT FOR LATER CALCULATIONS TO
    // UNDEFINED STATUS.
    //
    // FOR STAND ALONES, THINGS ARE SET TO UNDEFINED ONLY IF SOMETHING
    // ELSE CAN STILL DEFINE THEM (E.G., R IS SET UNDEFINED IF R0,
    // MP AND MT ARE ALL DEFINED).
    //

    // Local references for convenience
    double& A      = FLOAT_common.A;
    double& AI     = FLOAT_common.AI;
    double& ASO    = FLOAT_common.ASO;
    double& ASOI   = FLOAT_common.ASOI;
    double& ASYMPT = FLOAT_common.ASYMPT;
    double& E      = FLOAT_common.E;
    double& ECM    = FLOAT_common.ECM;
    double& ELAB   = FLOAT_common.ELAB;
    double& AM     = FLOAT_common.AM;
    double& AMP    = FLOAT_common.AMP;
    double& AMT    = FLOAT_common.AMT;
    double& R      = FLOAT_common.R;
    double& R0     = FLOAT_common.R0;
    double& RC     = FLOAT_common.RC;
    double& RC0    = FLOAT_common.RC0;
    double& RCP    = FLOAT_common.RCP;
    double& RCT    = FLOAT_common.RCT;
    double& RC0P   = FLOAT_common.RC0P;
    double& RC0T   = FLOAT_common.RC0T;
    double& RI     = FLOAT_common.RI;
    double& RI0    = FLOAT_common.RI0;
    double& RSO    = FLOAT_common.RSO;
    double& RSO0   = FLOAT_common.RSO0;
    double& RSOI   = FLOAT_common.RSOI;
    double& RSOI0  = FLOAT_common.RSOI0;
    double& RSI    = FLOAT_common.RSI;
    double& RSI0   = FLOAT_common.RSI0;
    double& SPAM   = FLOAT_common.SPAM;
    double& TAU    = FLOAT_common.TAU;
    double& TAUI   = FLOAT_common.TAUI;
    double& V      = FLOAT_common.V;
    double& VI     = FLOAT_common.VI;
    double& VSO    = FLOAT_common.VSO;
    double& VSOI   = FLOAT_common.VSOI;
    double& VSI    = FLOAT_common.VSI;
    double& VE     = FLOAT_common.VE;
    double& VIE    = FLOAT_common.VIE;
    double& VESQ   = FLOAT_common.VESQ;
    double& VIESQ  = FLOAT_common.VIESQ;
    double& ASI    = FLOAT_common.ASI;
    double& R0E    = FLOAT_common.R0E;
    double& RI0E   = FLOAT_common.RI0E;
    double& AE     = FLOAT_common.AE;
    double& AIE    = FLOAT_common.AIE;
    double& R0ESQ  = FLOAT_common.R0ESQ;
    double& RI0ESQ = FLOAT_common.RI0ESQ;
    double& AESQ   = FLOAT_common.AESQ;
    double& AIESQ  = FLOAT_common.AIESQ;
    double& POWRL  = FLOAT_common.POWRL;
    double& POWIM  = FLOAT_common.POWIM;

    int& L       = INTGER.L;
    int& NODES   = INTGER.NODES;
    int& IZP     = INTGER.IZP;
    int& IZT     = INTGER.IZT;
    int& PARITY  = INTGER.PARITY;

    double& J    = JBLOCK.J;
    double& JP   = JBLOCK.JP;
    double& JSP  = JBLOCK.JSP;
    double& JST  = JBLOCK.JST;

    int& IDONE   = INTRNL.IDONE;
    int& ISTRIP  = INTRNL.ISTRIP;
    double UNDEF = INTRNL.UNDEF;

    // EQUIVALENCE arrays (local pointers named with _p suffix to avoid
    // shadowing the inline accessor functions from ptolemy_commons.h)
    double* VTEN_p  = VTEN();   // VTEN(1)..VTEN(6) => VTR..VTPI
    double* RTEN_p  = RTEN();   // RTEN(1)..RTEN(6) => RTR..RTPI
    double* RTEN0_p = RTEN0();  // RTEN0(1)..RTEN0(6) => RTR0..RTPI0
    double* ATEN_p  = ATEN();   // ATEN(1)..ATEN(6) => ATR..ATPI

    if (ICHAN < 5) goto L200;

    //
    // CLEAR ONLY THINGS THAT CAN BE REDEFINED LATER
    //
    if (AMP == UNDEF || AMT == UNDEF) goto L150;
    AM = UNDEF;
    if (R0  != UNDEF) R    = UNDEF;
    if (RI0 != UNDEF) RI   = UNDEF;
    if (RSO0  != UNDEF) RSO  = UNDEF;
    if (RSOI0 != UNDEF) RSOI = UNDEF;
    if (RSI0  != UNDEF) RSI  = UNDEF;
    if (RC0  != UNDEF) RC  = UNDEF;
    if (RC0P != UNDEF) RCP = UNDEF;
    if (RC0T != UNDEF) RCT = UNDEF;
    for (int I = 1; I <= 6; I++) {
        if (RTEN0_p[I] != UNDEF) RTEN_p[I] = UNDEF;
    }
    if (ELAB != UNDEF) E = UNDEF;

L150:
    if (ECM != UNDEF) E = UNDEF;
    goto L400;

    //
    // WHEN DOING DWBA; FREE LOTS
    //
L200:
    E    = UNDEF;
    AM   = UNDEF;
    AMP  = UNDEF;
    AMT  = UNDEF;
    R    = UNDEF;
    RC   = UNDEF;
    RCP  = UNDEF;
    RCT  = UNDEF;
    RI   = UNDEF;
    RSO  = UNDEF;
    RSOI = UNDEF;
    RSI  = UNDEF;
    for (int I = 1; I <= 6; I++) {
        RTEN_p[I] = UNDEF;
    }
    SPAM = UNDEF;
    L     = NOTDEF_INT;
    NODES = NOTDEF_INT;
    IZP   = NOTDEF_INT;
    IZT   = NOTDEF_INT;
    J     = INTRNL.NOTDEF;
    JP    = INTRNL.NOTDEF;
    JSP   = INTRNL.NOTDEF;
    JST   = INTRNL.NOTDEF;
    PARITY = 0;
    for (int I = 1; I <= 2; I++) {
        INTGER.PARIPT[I] = 0;
        FLOAT_common.AMXGPT[I] = UNDEF;
        FLOAT_common.EXSPT[I]  = UNDEF;
        INTRNL.LSPCPT[I] = NOTDEF_INT;
        INTRNL.NODEPT[I] = NOTDEF_INT;
    }

    //
    // DO NOT SET POTENTIALS TO 0 WHEN GOING FROM ONE SCAT. STATE TO
    // THE NEXT.
    //
    if (ICHAN >= 3 && IDONE / 4 != 3) goto L400;
    TAU  = 0;
    TAUI = 0;
    VSO  = 0;
    VSOI = 0;
    VSI  = 0;
    VE   = 0;
    VIE  = 0;
    VESQ = 0;
    VIESQ = 0;
    for (int I = 1; I <= 6; I++) {
        VTEN_p[I] = 0;
    }
    V  = 0;
    VI = 0;

    //
    // CANCEL LINKULE REFERENCES EXCEPT FROM ONE SCAT STATE TO ANOTHER
    //
    // ALSO DO NOT CANCEL PRIOR TO EFFECTIVE EXCITATION POTENTIAL
    //
    if (ISTRIP == 0) goto L300;

    for (int I = 1; I <= LNKBLK.NUMLNK; I++) {
        LNKBLK.LNKADR[I][3] = 0;
    }

    //
    // WE LEAVE THE A'S DEFINED FROM ONE BOUND STATE TO ANOTHER AND
    // FROM ONE SCATTERING TO ANOTHER BUT NOT FROM BOUND TO SCATTER
    //
    if (ICHAN == 0) goto L310;
L300:
    if (IDONE == (1 << (ICHAN - 1))) goto L400;
    if (IDONE == 3 + (1 << (ICHAN - 1))) goto L400;
    if (IDONE == 12 + (1 << (ICHAN - 1))) goto L400;
L310:
    A     = UNDEF;
    AI    = UNDEF;
    ASO   = UNDEF;
    ASOI  = UNDEF;
    ASI   = UNDEF;
    R0    = UNDEF;
    RI0   = UNDEF;
    RSO0  = UNDEF;
    RSOI0 = UNDEF;
    RSI0  = UNDEF;
    R0E   = 0;
    RI0E  = 0;
    R0ESQ = 0;
    RI0ESQ = 0;
    RC0   = UNDEF;
    RC0P  = UNDEF;
    RC0T  = UNDEF;
    AE    = 0;
    AIE   = 0;
    AESQ  = 0;
    AIESQ = 0;
    for (int I = 1; I <= 6; I++) {
        RTEN0_p[I] = UNDEF;
        ATEN_p[I]  = UNDEF;
    }
    ASYMPT = UNDEF;

    //
    // WE ALWAYS UNDEFINE PARAMETERS THAT WERE SET TO OTHER PARAMETERS
    // BY SETPOT
    //
L400:
    if (!INTRNL.WASSET[1]) goto L410;
    RI0   = UNDEF;
    RI0E  = 0;
    RI0ESQ = 0;
L410:
    if (!INTRNL.WASSET[2]) goto L420;
    AI    = UNDEF;
    AIE   = 0;
    AIESQ = 0;
L420:
    if (INTRNL.WASSET[3])  RSO    = UNDEF;
    if (INTRNL.WASSET[4])  ASO    = UNDEF;
    if (INTRNL.WASSET[5])  RSOI   = UNDEF;
    if (INTRNL.WASSET[6])  ASOI   = UNDEF;
    if (INTRNL.WASSET[7])  RSI    = UNDEF;
    if (INTRNL.WASSET[8])  ASI    = UNDEF;
    if (INTRNL.WASSET[9])  RC0    = UNDEF;
    if (INTRNL.WASSET[9])  RC     = UNDEF;
    if (INTRNL.WASSET[10]) RCP    = UNDEF;
    if (INTRNL.WASSET[11]) RCT    = UNDEF;
    if (INTRNL.WASSET[12]) JP     = INTRNL.NOTDEF;
    if (INTRNL.WASSET[13]) J      = INTRNL.NOTDEF;
    if (INTRNL.WASSET[14]) L      = NOTDEF_INT;
    if (INTRNL.WASSET[15]) ASYMPT = UNDEF;
    if (INTRNL.WASSET[16]) ELAB   = UNDEF;
    for (int I = 1; I <= 6; I++) {
        if (INTRNL.WASSET[2*I+15]) RTEN_p[I]  = UNDEF;
        if (INTRNL.WASSET[2*I+16]) ATEN_p[I]  = UNDEF;
    }
    if (INTRNL.WASSET[29]) PARITY = 0;

    //
    // IF THIS IS INELASTIC SCATTERING, SET BACK TO FIRST CHANNEL
    // AFTER BOTH CHANNELS ARE DONE.
    //
    if (ISTRIP != 0 || IDONE != 12) goto L800;

    for (int I = 1; I <= LNKBLK.NUMLNK; I++) {
        for (int II = 1; II <= 6; II++) {
            LNKBLK.LNKADR[I][II] = LOCFIT.linkule.LNKAD2[I][1][II];
        }
    }
    R     = KANDM.RSCTS[1];
    RI    = WAVCOM.RIS[1];
    RC    = KANDM.RCSCTS[1];
    A     = KANDM.ASCTS[1];
    AI    = WAVCOM.AIS[1];
    V     = WAVCOM.V0RS[1];
    VI    = WAVCOM.V0IS[1];
    RSI   = WAVCOM.RSIS[1];
    ASI   = WAVCOM.ASIS[1];
    VSI   = WAVCOM.V0SIS[1];
    POWRL = WAVCOM.POWRLS[1];
    POWIM = WAVCOM.POWIMS[1];

L800:
    return;
}
