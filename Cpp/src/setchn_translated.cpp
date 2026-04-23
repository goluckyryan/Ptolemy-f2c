// setchn_translated.cpp — SETCHN subroutine
// Translated from source.f lines 32002-32479
// Sets up parameters for a given channel

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include "ptolemy_alloc.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_forward.h"
#include "ptolemy_io.h"
#include <cstdio>
#include <cstring>
#include <cmath>

void SETCHN(int ICHAN, int& IRET)
{
    // ELECTRON MASS FOR CONVERTING ATOMIC TO NUCLEAR MASSES
    static const double EMASS = 0.5110034;

    IRET = 1;

    int JPOLD = (int)JBLOCK.JP;
    int JOLD = (int)JBLOCK.J;
    int LOLD = INTGER.L;
    int IPOLD = INTGER.PARITY;

    // Declare all variables at function scope to avoid goto-crosses-init
    double SIGNMX = INTRNL.UNDEF;
    int ISGNZX = NOTDEF_INT;
    int IP = 0, IT = 0;
    double TMP = 0, TMT = 0, CMLAB = 1.0;
    double AMXCP = 0, AMXCT = 0;

    // WASSET tracks things set to other things; CLRCHN uses it
    for (int I = 1; I <= 32; I++)
        INTRNL.WASSET[I] = 0;

    // IBND = 0 means not bound state, > 0 means bound state
    int IBND = 1;
    if (ICHAN == 5) goto L300;
    IBND = 0;
    if (ICHAN == 6) goto L310;

    // Compute mass excess = mass excess (ground state) + E*
    for (int I = 1; I <= 5; I++) {
        if (FLOAT_common.AMXCS[I] != INTRNL.UNDEF) continue;
        if (FLOAT_common.AMXCGS[I] == INTRNL.UNDEF) continue;
        FLOAT_common.AMXCS[I] = FLOAT_common.AMXCGS[I] + FLOAT_common.EXS[I];
    }

    // Find ZX, MX and the signed ZX, MX

    // If charge transfer is 0 or MX = 0 (inelastic), no ambiguity
    if (DABS(FLOAT_common.AMX) < 0.3) SIGNMX = 0;
    if (INTGER.IZS[5] == 0) ISGNZX = 0;

    for (int I = 1; I <= 3; I += 2) {
        if (INTGER.IZS[I] == NOTDEF_INT || INTGER.IZS[I+1] == NOTDEF_INT)
            goto L140;
        ISGNZX = (INTGER.IZS[I] - INTGER.IZS[I+1]) * (2 - I);
        if (INTGER.IZS[5] == NOTDEF_INT) INTGER.IZS[5] = IABS(ISGNZX);
    L140:
        if (AMS()[I] == INTRNL.UNDEF || AMS()[I+1] == INTRNL.UNDEF)
            continue;
        SIGNMX = (AMS()[I] - AMS()[I+1]) * (2 - I);
        if (AMS()[5] == INTRNL.UNDEF) AMS()[5] = DABS(SIGNMX);
    }
    if (DABS(SIGNMX) < 0.3) SIGNMX = 0;

    // Supply missing mass or Z values if possible
    for (int I = 1; I <= 4; I++) {
        int I2 = I - 1 + 2 * MOD(I, 2);
        int KSIGN = IABS(5 - 2*I) - 2;
        if (ISGNZX == NOTDEF_INT) goto L180;
        if (INTGER.IZS[I] != NOTDEF_INT) goto L180;
        if (INTGER.IZS[I2] == NOTDEF_INT) goto L180;
        INTGER.IZS[I] = INTGER.IZS[I2] + KSIGN * ISGNZX;
    L180:
        if (SIGNMX == INTRNL.UNDEF) continue;
        if (AMS()[I] != INTRNL.UNDEF) continue;
        if (AMS()[I2] == INTRNL.UNDEF) continue;
        AMS()[I] = AMS()[I2] + KSIGN * SIGNMX;
    }

    // If requested convert atomic numbers to true masses
    if (SWITCH.MASTYP != 2) goto L220;
    for (int I = 1; I <= 5; I++) {
        if (FLOAT_common.AMXCS[I] == INTRNL.UNDEF) continue;
        int M = (int)AMS()[I];
        if (AMS()[I] != M) continue;
        AMS()[I] = AMS()[I] + FLOAT_common.AMXCS[I] / CNSTNT.AMUMEV;
    }

    // Determine which is P and T
L220:
    INTRNL.ISTRIP = NOTDEF_INT;
    if (ISGNZX != NOTDEF_INT) INTRNL.ISTRIP = ISGNZX;
    if (SIGNMX != INTRNL.UNDEF) INTRNL.ISTRIP = (int)SIGNMX;
    if (INTRNL.ISTRIP == NOTDEF_INT)
        std::printf("\n***** WARNING:  NOT ENOUGH DATA TO DETERMIN IF"
                    " REACTION IS STRIPPING OR PICKUP - STRIPPING ASSUMED.\n");
    if (INTRNL.ISTRIP == NOTDEF_INT) INTRNL.ISTRIP = 1;
    if (INTRNL.ISTRIP != 0) INTRNL.ISTRIP = ISIGN(1, INTRNL.ISTRIP);

    if (ICHAN >= 3) goto L250;

    // Bound state, projectile is always X
    IP = 5;
    if (INTRNL.ISTRIP < 0) goto L240;
    IT = ICHAN + 1;
    IBND = IT + 2*ICHAN - 3;
    goto L280;
L240:
    IBND = ICHAN + 1;
    IT = IBND + 2*ICHAN - 3;
    goto L280;

    // Scattering states - no composite
L250:
    IBND = 0;
    IP = ICHAN - 2;
    IT = IP + 2;

    // Define undefined P and T stuff
L280:
    if (INTGER.IZP == NOTDEF_INT) INTGER.IZP = INTGER.IZS[IP];
    if (INTGER.IZT == NOTDEF_INT) INTGER.IZT = INTGER.IZS[IT];
    if (FLOAT_common.AMP == INTRNL.UNDEF) FLOAT_common.AMP = AMS()[IP];
    if (FLOAT_common.AMT == INTRNL.UNDEF) FLOAT_common.AMT = AMS()[IT];
    if (JBLOCK.JSP == INTRNL.NOTDEF) JBLOCK.JSP = JBLOCK.JS[IP];
    if (JBLOCK.JST == INTRNL.NOTDEF) JBLOCK.JST = JBLOCK.JS[IT];

    // Copy E* for P and T
    if (FLOAT_common.EXSPT[1] == INTRNL.UNDEF) FLOAT_common.EXSPT[1] = FLOAT_common.EXS[IP];
    if (FLOAT_common.EXSPT[2] == INTRNL.UNDEF) FLOAT_common.EXSPT[2] = FLOAT_common.EXS[IT];
    if (FLOAT_common.AMXGPT[1] == INTRNL.UNDEF) FLOAT_common.AMXGPT[1] = FLOAT_common.AMXCGS[IP];
    if (FLOAT_common.AMXGPT[2] == INTRNL.UNDEF) FLOAT_common.AMXGPT[2] = FLOAT_common.AMXCGS[IT];
    if (INTGER.PARIPT[1] == 0) INTGER.PARIPT[1] = INTGER.PARITS[IP];
    if (INTGER.PARIPT[2] == 0) INTGER.PARIPT[2] = INTGER.PARITS[IT];
    INTRNL.LSPCPT[2] = INTRNL.LSPECS[IT];
    INTRNL.NODEPT[2] = INTRNL.NODESP[IT];

    // For scattering allow ranges of L, JP, and J
    if (IBND == 0) goto L310;

    // Following page or so is for bound states only
    if (JBLOCK.J == INTRNL.NOTDEF) JBLOCK.J = JBLOCK.JS[IBND];
    if (INTGER.PARITY == 0) INTGER.PARITY = INTGER.PARITS[IBND];
    if (INTGER.L == NOTDEF_INT) INTGER.L = INTRNL.LSPECS[IBND];
    if (INTGER.NODES == NOTDEF_INT) INTGER.NODES = INTRNL.NODESP[IBND];

    // The following code is also used when not doing a DWBA
L300:
    if (INTGER.L == NOTDEF_INT) INTGER.L = INTRNL.LSPCPT[2];
    if (INTGER.NODES == NOTDEF_INT) INTGER.NODES = INTRNL.NODEPT[2];

    // If some angular momenta are zero then we can define others
    if (JBLOCK.JP == INTRNL.NOTDEF && JBLOCK.JST == 0) JBLOCK.JP = JBLOCK.J;
    if (JBLOCK.JP == INTRNL.NOTDEF && JBLOCK.J == 0) JBLOCK.JP = JBLOCK.JST;
    if (JBLOCK.JP == INTRNL.NOTDEF && INTGER.L == 0) JBLOCK.JP = JBLOCK.JSP;
    if (JBLOCK.JP == INTRNL.NOTDEF && JBLOCK.JSP == 0 &&
        INTGER.L != NOTDEF_INT) JBLOCK.JP = 2*INTGER.L;
    if (JBLOCK.J == INTRNL.NOTDEF && JBLOCK.JST == 0) JBLOCK.J = JBLOCK.JP;
    if (INTGER.L != NOTDEF_INT) goto L305;
    if (JBLOCK.JP != INTRNL.NOTDEF && JBLOCK.JSP == 0) INTGER.L = (int)JBLOCK.JP / 2;
    if (JBLOCK.JSP != INTRNL.NOTDEF && JBLOCK.JP == 0) INTGER.L = (int)JBLOCK.JSP / 2;

    // If only two possible L's, parity can tell us
    if (INTGER.L != NOTDEF_INT) goto L305;
    { int I = INTGER.PARITY * INTGER.PARIPT[1] * INTGER.PARIPT[2];
    if (I == 0) goto L310;
    I = (I + 3) / 2;
    if (JBLOCK.JP == INTRNL.NOTDEF || JBLOCK.JSP == INTRNL.NOTDEF) goto L310;
    int L1 = IABS((int)JBLOCK.JP - (int)JBLOCK.JSP) / 2;
    L1 = L1 + MOD(I + L1, 2);
    int L2 = ((int)JBLOCK.JP + (int)JBLOCK.JSP) / 2;
    if (L2 - L1 > 1) goto L310;
    INTGER.L = L1; }

    // L is defined; we may be able to patch up a parity
L305:
    { int I = +1;
    if (MOD(INTGER.L, 2) != 0) I = -1;
    if (INTGER.PARITY == 0) INTGER.PARITY = INTGER.PARIPT[1] * INTGER.PARIPT[2] * I;
    if (INTGER.PARIPT[1] == 0) INTGER.PARIPT[1] = I * INTGER.PARITY * INTGER.PARIPT[2];
    if (INTGER.PARIPT[2] == 0) INTGER.PARIPT[2] = I * INTGER.PARITY * INTGER.PARIPT[1]; }

    // Setup kinematic projectile and target masses
L310:
    AMXCP = 0;
    if (FLOAT_common.AMXGPT[1] != INTRNL.UNDEF)
        AMXCP = (FLOAT_common.AMXGPT[1] + FLOAT_common.EXSPT[1]) / CNSTNT.AMUMEV;
    AMXCT = 0;
    if (FLOAT_common.AMXGPT[2] != INTRNL.UNDEF)
        AMXCT = (FLOAT_common.AMXGPT[2] + FLOAT_common.EXSPT[2]) / CNSTNT.AMUMEV;

    // MP and MT must be separately defined
    if (FLOAT_common.AMP != INTRNL.UNDEF && FLOAT_common.AMT != INTRNL.UNDEF) goto L330;
    std::printf("\n**** ERROR:  BOTH MP AND MT ARE REQUIRED.\n");
    IRET = 0;
    goto L400;

    // If necessary convert to floating mass
L330:
    if (SWITCH.MASTYP != 2) goto L340;
    { int M_ = (int)FLOAT_common.AMP;
    if (M_ == FLOAT_common.AMP) FLOAT_common.AMP = FLOAT_common.AMP + AMXCP - INTGER.IZP*(EMASS/CNSTNT.AMUMEV); }
    { int M_ = (int)FLOAT_common.AMT;
    if (M_ == FLOAT_common.AMT) FLOAT_common.AMT = FLOAT_common.AMT + AMXCT - INTGER.IZT*(EMASS/CNSTNT.AMUMEV); }

    // Now the kinematic masses
L340:
    TMP = FLOAT_common.AMP;
    TMT = FLOAT_common.AMT;
    if (SWITCH.MASTYP != 0) goto L370;
    { int M_ = (int)TMP;
    if (M_ == TMP) TMP = TMP + AMXCP - INTGER.IZP*(EMASS/CNSTNT.AMUMEV); }
    { int M_ = (int)TMT;
    if (M_ == TMT) TMT = TMT + AMXCT - INTGER.IZT*(EMASS/CNSTNT.AMUMEV); }

    // Define reduced mass if necessary
L370:
    if (FLOAT_common.AM != INTRNL.UNDEF) goto L380;
    FLOAT_common.AM = CNSTNT.AMUMEV * TMP * TMT / (TMP + TMT);

    // Find ELAB to ECM conversion
L380:
    INTRNL.RATMAS = TMP / TMT;
    CMLAB = 1 + INTRNL.RATMAS;

    // Find E if necessary
L400:
    if (ICHAN >= 6) goto L450;
    if (FLOAT_common.E != INTRNL.UNDEF) goto L500;
    switch (ICHAN) {
        case 1: goto L410;
        case 2: goto L410;
        case 3: goto L450;
        case 4: goto L480;
        case 5: goto L430;
        case 6: goto L450;
    }

    // Bound state; first try Q and then mass excess
L410:
    if (FLOAT_common.Q == INTRNL.UNDEF || INTRNL.EBNDS[3-ICHAN] == INTRNL.UNDEF)
        goto L420;
    FLOAT_common.E = (3 - 2*ICHAN) * INTRNL.ISTRIP * FLOAT_common.Q + INTRNL.EBNDS[3-ICHAN];
    goto L500;
L420:
    if (FLOAT_common.AMXCS[IP] == INTRNL.UNDEF || FLOAT_common.AMXCS[IT] == INTRNL.UNDEF
        || FLOAT_common.AMXCS[IBND] == INTRNL.UNDEF) goto L430;
    FLOAT_common.E = FLOAT_common.AMXCS[IBND] - FLOAT_common.AMXCS[IP] - FLOAT_common.AMXCS[IT];
    goto L500;
L430:
    std::printf("\n**** E MUST BE DEFINED.\n");
    IRET = 0;
    goto L500;

    // Incoming scattering state; ECM or ELAB applies
L450:
    if (FLOAT_common.ECM != INTRNL.UNDEF) FLOAT_common.E = FLOAT_common.ECM;
    if (FLOAT_common.ELAB == INTRNL.UNDEF) goto L465;
    FLOAT_common.E = FLOAT_common.ELAB / CMLAB;
    goto L500;
L465:
    if (FLOAT_common.E != INTRNL.UNDEF) goto L500;
    std::printf("\n**** E, ECM OR ELAB MUST BE DEFINED.\n");
    IRET = 0;
    goto L500;

    // Outgoing scattering state; apply Q to ECM
L480:
    if (FLOAT_common.ECM == INTRNL.UNDEF) goto L430;
    if (FLOAT_common.Q != INTRNL.UNDEF) goto L495;
    // Attempt to get Q from total mass excesses
    for (int I = 1; I <= 4; I++) {
        if (FLOAT_common.AMXCS[I] == INTRNL.UNDEF) goto L430;
    }
    FLOAT_common.Q = FLOAT_common.AMXCS[1] - FLOAT_common.AMXCS[2] + FLOAT_common.AMXCS[3] - FLOAT_common.AMXCS[4];
L495:
    FLOAT_common.E = FLOAT_common.ECM + FLOAT_common.Q;
    goto L500;

    // At last E is defined. For incoming wave define ECM, ELAB
L500:
    if (ICHAN != 3 && ICHAN != 6) goto L550;
    FLOAT_common.ECM = FLOAT_common.E;
    INTRNL.WASSET[16] = (FLOAT_common.ELAB == INTRNL.UNDEF && ICHAN == 6) ? 1 : 0;
    FLOAT_common.ELAB = FLOAT_common.E * CMLAB;

    // For outgoing channel, define Q
L550:
    if (ICHAN != 4) goto L600;
    if (FLOAT_common.ECM != INTRNL.UNDEF) FLOAT_common.Q = FLOAT_common.E - FLOAT_common.ECM;

    // Store P and T stuff back into 1-5 stuff
L600:
    if (ICHAN > 4) goto L700;
    AMS()[IP] = FLOAT_common.AMP;
    AMS()[IT] = FLOAT_common.AMT;
    INTGER.IZS[IP] = INTGER.IZP;
    INTGER.IZS[IT] = INTGER.IZT;
    JBLOCK.JS[IP] = JBLOCK.JSP;
    JBLOCK.JS[IT] = JBLOCK.JST;
    INTGER.PARITS[IP] = INTGER.PARIPT[1];
    INTGER.PARITS[IT] = INTGER.PARIPT[2];
    if (IBND != 0) JBLOCK.JS[IBND] = JBLOCK.J;
    if (IBND != 0) INTGER.PARITS[IBND] = INTGER.PARITY;
    if (FLOAT_common.SPAM == INTRNL.UNDEF) goto L700;
    if (ICHAN == 1) FLOAT_common.SPAMP = FLOAT_common.SPAM;
    if (ICHAN == 2) FLOAT_common.SPAMT = FLOAT_common.SPAM;

L700:
    INTRNL.WASSET[12] = (JPOLD != (int)JBLOCK.JP) ? 1 : 0;
    INTRNL.WASSET[13] = (JOLD != (int)JBLOCK.J) ? 1 : 0;
    INTRNL.WASSET[14] = (LOLD != INTGER.L) ? 1 : 0;
    INTRNL.WASSET[29] = (INTGER.PARITY != IPOLD) ? 1 : 0;

    // Set asymptopia if needed
    if (FLOAT_common.ASYMPT != INTRNL.UNDEF) goto L750;
    INTRNL.WASSET[15] = 1;
    FLOAT_common.ASYMPT = FLOAT_common.BNDASY;
    if (IBND == 0) FLOAT_common.ASYMPT = DABS(FLOAT_common.SCTASY);

L750:
    return;
}
