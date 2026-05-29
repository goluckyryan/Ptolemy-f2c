// source_channels.cpp — CC input parsing and setup routines
// BASLBL → baslbl_translated.cpp, SETFG → setfg_translated.cpp
// CHANEL → chanel_translated.cpp, SETCHN → setchn_translated.cpp (when ready)

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include "ptolemy_alloc.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_forward.h"
#include "ptolemy_io.h"
#include <cstdio>
#include <cstring>
#include <cmath>

extern int GIVEAL(int IWRDS);
extern int CHOPIT(int IWRDS);

// ============================================================================
// BASCPL — Constructs basis coupling array
// Translated from source.f lines 2164-3214
// ============================================================================
void BASCPL(int& IRET)
{
    // Local variables
    int MBASCP = 11;
    int LPL = 3;
    int IFACS[5] = {0, 1, 2, -1, -2}; // 1-based

    bool PBUGSW, SYMSW, FLIPSW;

    // Fortran DIMENSION JNUC(2,2) = column-major flat array, 1-based
    // Access: JNUC(i,j) = JNUC_flat[(j-1)*2 + (i-1)]
    int JNUC_flat[5] = {}; // 0-based: positions 0..3 used for (1,1),(2,1),(1,2),(2,2)
    double ENUC_flat[5] = {};
    int KNUC_flat[5] = {};
    #define JNUC(i,j) JNUC_flat[((j)-1)*2 + ((i)-1)]
    #define ENUC(i,j) ENUC_flat[((j)-1)*2 + ((i)-1)]
    #define KNUC(i,j) KNUC_flat[((j)-1)*2 + ((i)-1)]
    int SCHN[3]; // 1-based: SCHN(2)
    #define JP1 JNUC(1,1)
    #define JP2 JNUC(1,2)
    #define JT1 JNUC(2,1)
    #define JT2 JNUC(2,2)
    #define EP1 ENUC(1,1)
    #define EP2 ENUC(1,2)
    #define ET1 ENUC(2,1)
    #define ET2 ENUC(2,2)

    int JPT12_flat[5] = {};
    double EPT12_flat[5] = {};
    int KPT12_flat[5] = {};
    #define JPT12(i,j) JPT12_flat[((j)-1)*2 + ((i)-1)]
    #define EPT12(i,j) EPT12_flat[((j)-1)*2 + ((i)-1)]
    #define KPT12(i,j) KPT12_flat[((j)-1)*2 + ((i)-1)]

    int NAME12[3]; // 1-based: NAME12(2)
    #define NAME1 NAME12[1]
    #define NAME2 NAME12[2]

    double RNUC_local[3]; // 1-based
    double RCOUL[3]; // 1-based

    float BETACS[3]; // 1-based: BETACS(2)
    float BETANS[3]; // 1-based
    int LXX2S[3]; // 1-based
    double xcs[3], xns[3]; // 1-based
    int IBETAC[3], IPTS[3]; // 1-based
    int LXS[3], IBETAN[3]; // 1-based
    double RMENUC[3], RMECOU[3]; // 1-based
    double SCLFAC[3]; // 1-based
    double VALS[9][3]; // 1-based: VALS(8,2)
    double DELTAN[3]; // 1-based

    char8 NAME5("TEMPCUP");
    char8 IPTWRD[3]; // 1-based
    IPTWRD[1] = char8("PROJ");
    IPTWRD[2] = char8("TARG");
    char8 ORDWRD[3];
    ORDWRD[1] = char8("FIRST");
    ORDWRD[2] = char8("SECOND");
    char8 MODNAM;
    char8 MODELS_local[8];
    MODELS_local[1] = char8("ROTATION");
    MODELS_local[2] = char8("VIBRATIO");
    MODELS_local[3] = char8("DEFORMED");
    MODELS_local[4] = char8("HOTFUDGE");
    MODELS_local[5] = char8("MODEL5");
    MODELS_local[6] = char8("MODEL6");
    MODELS_local[7] = char8("MODEL7");

    int NUMFLP = 6;
    int IFLIPS[7] = {0, 4, 6, 8, 11, 15, 17}; // 1-based

    int I, II, III, IC;
    int ICH1, ICH2, ICP;
    int ICH1PT, ICH2PT, ICPPT;
    int IPARI1, IPARI2, IPARTY;
    int IB1ST, IB1END, IB2ST, IB2END;
    int IB1, IB2, IB1PT, IB2PT;
    int IORDER, ITYPE, IMUL, IPT, MODEL;
    int LXTOT, LXTOP, ICFLIP;
    int MINLAM, MAXLAM, LXX2, LX;
    int ICOUNT, ICNT2;
    int IPT1, IRTN, IORD, IUP, ITEMP;
    int IFAC, LDIFX2, LDIF;
    int LCHNDF, NCHNDF, MCHNDF;
    int LCHNCP, NCHNCP, MCHNCP;
    int LBASDF_local, NBASDF, MBASDF;
    int LTRANS, MTRANS;
    int ITMPCP, LBASCP, LTMPCP;
    int IWRDS, IL;
    int IBCPT, IBC, LTEMP, LTMP2;
    int N5;
    double SCALE, FAC, RMESGN, TEMP, BETTMP;
    double XN, XC, RMEGN, RMEGC;
    double X;

    PBUGSW = MOD(INTGER.IPRINT, 10) >= 5;

    // Get the deformation radii
    RCOUL[1] = FLOAT_common.RCP;
    RCOUL[2] = FLOAT_common.RCT;
    I = 0;
    if (FLOAT_common.RCP == INTRNL.UNDEF || FLOAT_common.RCP == 0 ||
        FLOAT_common.RCT == INTRNL.UNDEF) I = 1;
    RNUC_local[1] = FLOAT_common.R0 * std::pow(FLOAT_common.AMA, 1.0/3.0);
    if (I == 1)
        RCOUL[1] = FLOAT_common.RC0 * std::pow(FLOAT_common.AMA, 1.0/3.0);
    RNUC_local[2] = FLOAT_common.R0 * std::pow(FLOAT_common.AMBIGA, 1.0/3.0);
    if (I == 1)
        RCOUL[2] = FLOAT_common.RC0 * std::pow(FLOAT_common.AMBIGA, 1.0/3.0);
    std::printf("\n-ELECTROMAGNETIC MATRIX ELEMENTS ARE RELATED TO BETAS UP TO %.6s ORDER.\n",
        ORDWRD[SWITCH.MEORD+1].data);
    std::printf("\n%16sDEFORMED RADII\n", "");
    std::printf("%15sNUCLEAR   COULOMB\n", "");
    std::printf(" PROJECTILE%1s%10.4f%10.4f\n", "", RNUC_local[1], RCOUL[1]);
    std::printf(" TARGET%5s%10.4f%10.4f\n", "", RNUC_local[2], RCOUL[2]);

    // We need the LOG(FACTORIAL) table
    I = 8 * INELCM.LXMAX;
    SETLOG(I);

    // Identical particle processing
    SYMSW = WAVCOM.ISTATS[2] < 3;
    if (SYMSW) std::printf("\n SYMMETRIZED COUPLINGS ARE BEING GENERATED FOR IDENTICAL NUCLEAR STATES.\n");

    // Start with a temporary array
    N5 = NAMLOC(NAME5.data);
    if (N5 == 0) LOCPTRS.Z[N5] = -LOCPTRS.Z[N5];
    ITMPCP = GIVEAL(IWRDS);
    LBASCP = LOCPTRS.Z[ITMPCP];
    NAMCOM.NAMES[ITMPCP] = NAME5;
    LBASCP = ALLOCS.FACFR4 * (LBASCP - 1) + 1;

    SETLOG(I);
    LCHNDF = LOCPTRS.Z[CCBLK.ICHNDF];
    LCHNDF = ALLOCS.FACFR4 * (LCHNDF - 1) + 1;
    NCHNDF = ILLOC(LCHNDF);
    MCHNDF = ILLOC(LCHNDF + 1);
    LCHNDF = ILLOC(LCHNDF + 2) + LCHNDF;
    LCHNCP = LOCPTRS.Z[CCBLK.ICHNCP];
    LCHNCP = ALLOCS.FACFR4 * (LCHNCP - 1) + 1;
    NCHNCP = ILLOC(LCHNCP);
    MCHNCP = ILLOC(LCHNCP + 1);
    LCHNCP = ILLOC(LCHNCP + 2) + LCHNCP;
    LBASDF_local = LOCPTRS.Z[CCBLK.IBASDF];
    LBASDF_local = ALLOCS.FACFR4 * (LBASDF_local - 1) + 1;
    NBASDF = ILLOC(LBASDF_local);
    MBASDF = ILLOC(LBASDF_local + 1);
    LBASDF_local = LBASDF_local + ILLOC(LBASDF_local + 2);
    LTRANS = 0;
    if (CCBLK.ITRANS == 0) goto L20;
    LTRANS = LOCPTRS.Z[CCBLK.ITRANS];
    LTRANS = ALLOCS.FACFR4 * (LTRANS - 1) + 1;
    MTRANS = ILLOC(LTRANS + 1);

L20:
    CCBLK.IORDX = 1;

    // Begin the guts
    {
        int NBASCP_local = 0;
        std::printf("\n-COUPLINGS USED IN THE CALCULATION:\n");
        std::printf("\n  CHANNEL%4sMODEL ORDER NUC. LX%4sSCALE"
            "%4sBETAN%2sDELTAN%2sR.M.E."
            "%7sBETAC%2sDELTAC%8sR.M.E.%8sBE UP"
            "%9sQ%8sMEANLIFE\n",
            "", "", "", "", "", "", "", "", "", "", "", "");
        std::printf("  FROM   TO%26sFACTOR%10sFM%5sALPHA"
            "%18sFM"
            "%7sE B**LX/2%4sE**2 B**LX%5sE BARN"
            "%8sSEC.\n",
            "", "", "", "", "", "", "", "");

        // Loop over all channels on the left side
        for (ICH1 = 1; ICH1 <= NCHNDF; ICH1++) {
            ICH1PT = (ICH1 - 1) * MCHNDF + LCHNDF;
            NAME1 = ILLOC(ICH1PT);
            JPT12(2,1) = IABS(ILLOC(ICH1PT + 1));
            EPT12(2,1) = ALLOC4(ICH1PT + 2);
            JPT12(1,1) = IABS(ILLOC(ICH1PT + 3));
            IPARI1 = ISIGN(1, ILLOC(ICH1PT + 1)) * ISIGN(1, ILLOC(ICH1PT + 3));
            EPT12(1,1) = ALLOC4(ICH1PT + 4);
            IB1ST = ILLOC(ICH1PT + 7);
            IB1END = IB1ST - 1 + ILLOC(ICH1PT + 8);
            KPT12(2,1) = ILLOC(ICH1PT + 15);
            KPT12(1,1) = ILLOC(ICH1PT + 16);
            std::printf(" \n");

            // Loop over channels on the right side
            for (ICH2 = 1; ICH2 <= NCHNDF; ICH2++) {
                ICH2PT = (ICH2 - 1) * MCHNDF + LCHNDF;
                NAME2 = ILLOC(ICH2PT);
                JPT12(2,2) = IABS(ILLOC(ICH2PT + 1));
                EPT12(2,2) = ALLOC4(ICH2PT + 2);
                JPT12(1,2) = IABS(ILLOC(ICH2PT + 3));
                IPARI2 = ISIGN(1, ILLOC(ICH2PT + 1)) * ISIGN(1, ILLOC(ICH2PT + 3));
                EPT12(1,2) = ALLOC4(ICH2PT + 4);
                IB2ST = ILLOC(ICH2PT + 7);
                IB2END = IB2ST - 1 + ILLOC(ICH2PT + 8);
                KPT12(2,2) = ILLOC(ICH2PT + 15);
                KPT12(1,2) = ILLOC(ICH2PT + 16);

                IPARTY = 0;
                if (IPARI1 * IPARI2 < 0) IPARTY = 1;

                // Loop over all channel couplings
                for (ICP = 1; ICP <= NCHNCP; ICP++) {
                    ICOUNT = 0;
                    ICNT2 = 0;
                    NAME1 = IABS(NAME1);
                    NAME2 = IABS(NAME2);
                    ICPPT = (ICP - 1) * MCHNCP + LCHNCP;
                    ICFLIP = 2;
                    if (IABS(ILLOC(ICPPT)) == NAME2 &&
                        ILLOC(ICPPT + 1) == NAME1) goto L150;
                    if (ILLOC(ICPPT + 3) == 2) goto L949;
                    ICFLIP = 1;
                    if (ILLOC(ICPPT) != NAME1 ||
                        ILLOC(ICPPT + 1) != NAME2) {
                        goto L949;
                    }

                L150:
                    IORDER = ILLOC(ICPPT + 2);
                    ITYPE = ILLOC(ICPPT + 3);
                    SCALE = 1;
                    for (I = 1; I <= IORDER; I++) {
                        LXX2S[I] = ILLOC(ICPPT + 3 + I);
                        BETACS[I] = ALLOC4(ICPPT + 5 + I);
                        BETANS[I] = ALLOC4(ICPPT + 7 + I);
                        IBETAC[I] = ILLOC(ICPPT + 10 + I);
                        IBETAN[I] = ILLOC(ICPPT + 14 + I);
                        SCLFAC[I] = ALLOC4(ICPPT + 16 + I);
                        SCALE = SCALE * SCLFAC[I];
                    }
                    IMUL = ILLOC(ICPPT + 4);
                    if (IORDER >= 2) IMUL = NOTDEF_INT;
                    IPT = ILLOC(ICPPT + 10);
                    ILLOC(ICPPT + 13) = 1;
                    MODEL = ILLOC(ICPPT + 14);
                    MODNAM = MODELS_local[MODEL];
                    LXTOT = ILLOC(ICPPT + 19);
                    LXTOP = ILLOC(ICPPT + 20);

                    // The leftside channel may have been reversed
                    FLIPSW = ILLOC(ICPPT) < 0;
                    if (FLIPSW) NAME12[ICFLIP] = -NAME12[ICFLIP];
                    for (IC = 1; IC <= 2; IC++) {
                        for (I = 1; I <= 2; I++) {
                            II = I;
                            if (FLIPSW && IC == ICFLIP) II = 3 - I;
                            JNUC(II,IC) = JPT12(I,IC);
                            ENUC(II,IC) = EPT12(I,IC);
                            KNUC(II,IC) = KPT12(I,IC);
                        }
                    }

                    // For first order define IPT
                    if (IORDER >= 2) goto L500;
                    if (IPT != NOTDEF_INT) goto L450;

                    // Which nucleus is being coupled in first order?
                    IPT = 1;
                    if (JP1 != JP2 || EP1 != EP2) goto L450;
                    IPT = 2;
                    if (JT1 != JT2 || ET1 != ET2) goto L450;

                    // Diagonal case
                    if (EP1 == 0 || JP1 == 0) goto L450;
                    IPT = 1;
                    if (ET1 == 0 || JT1 == 0) goto L450;
                    std::printf("\n**** CANNOT DETERMIN WHICH NUCLEUS (PROJECTILE OR TARGET) TO COUPLE IN FIRST ORDER.\n");
                L445:
                    std::printf("      NUCLEI SPINS AND ENERGIES ARE:\n");
                    std::printf("      CHANNEL%4d, PROJ. AND TARGET:%4d/2%7.3f  %4d/2%7.3f\n",
                        NAME1, JP1, EP1, JT1, ET1);
                    std::printf("      CHANNEL%4d, PROJ. AND TARGET:%4d/2%7.3f  %4d/2%7.3f\n",
                        NAME2, JP2, EP2, JT2, ET2);
                    IRET = 0;
                    goto L949;

                    // Make sure it is a legal 1st order coupling
                L450:
                    if (JNUC(3-IPT,1) == JNUC(3-IPT,2) &&
                        ENUC(3-IPT,1) == ENUC(3-IPT,2)) goto L470;
                    std::printf("\n*** BOTH NUCLEII CANNOT CHANGE STATE IN A  FIRST ORDER COUPLING.\n");
                    { char buf[5]; std::memcpy(buf, IPTWRD[IPT].data, 4); buf[4] = '\0';
                    std::printf("      OR THE WRONG NUCLEUS (%s) WAS SPECIFIED AS BEING EXCITED.\n", buf); }
                    goto L445;

                L470:
                    ILLOC(ICPPT + 10) = IPT;
                    MINLAM = IABS(JNUC(IPT,1) - JNUC(IPT,2));
                    MAXLAM = JNUC(IPT,1) + JNUC(IPT,2);
                    if (LXTOP != NOTDEF_INT) MAXLAM = MIN0(LXTOP, MAXLAM);

                    MINLAM = MINLAM + MOD(MINLAM + 2 * IPARTY, 4);
                    MAXLAM = MAXLAM - MOD(MAXLAM + 2 * IPARTY, 4);

                    if (MINLAM <= MAXLAM) goto L475;
                    std::printf("\n*** CHANNEL%5d IS TO BE COUPLED TO%5d IN FIRST ORDER BUT SPINS AND PARITIES CANNOT BE.\n",
                        NAME1, NAME2);
                    std::printf("      MINLAM, MAXLAM =%5d/2%5d/2\n", MINLAM, MAXLAM);
                    goto L445;

                L475:
                    if (IMUL != NOTDEF_INT || LXTOT != NOTDEF_INT) goto L480;

                    // LX not explicitly given
                    if (MINLAM == 0 && MAXLAM >= 4) MINLAM = 4;
                    if (BETACS[1] == undef4() && BETANS[1] == undef4()) goto L600;
                    if (MINLAM != MAXLAM)
                        std::printf("\n*** WARNING, CHANNEL%4d CAN BE COUPLED TO CHANNEL%4d"
                            " IN FIRST ORDER WITH MULTIPOLES FROM%4d/2 TO%4d/2.\n"
                            "      YOU HAVE SPECIFIED BETAS BUT NOT A PARTICULAR"
                            " MULTIPOLE, ONLY THE LOWEST WILL BE USED.\n\n",
                            NAME1, NAME2, MINLAM, MAXLAM);
                    MAXLAM = MINLAM;
                    IMUL = MINLAM;
                    goto L600;

                L480:
                    if (IMUL == NOTDEF_INT) IMUL = LXTOT;
                    if (IMUL >= MINLAM && IMUL <= MAXLAM) goto L485;
                    std::printf("\n*** CHANNEL%5d IS TO BE COUPLED TO%5d IN FIRST ORDER WITH MULTIPOLARITY%3d/2\n",
                        NAME1, NAME2, IMUL);
                    std::printf("      BUT ONLY%3d/2 =< LX =<%3d/2 IS POSSIBLE.\n", MINLAM, MAXLAM);
                    goto L949;
                L485:
                    MINLAM = IMUL;
                    MAXLAM = IMUL;
                    goto L600;

                    // Second order coupling
                L500:
                    MINLAM = IABS(LXX2S[1] - LXX2S[2]);
                    MAXLAM = LXX2S[1] + LXX2S[2];
                    if (LXTOP != NOTDEF_INT) MAXLAM = MIN0(LXTOP, MAXLAM);

                    // For mutual excitation
                    if (IPT != 3) goto L550;
                    for (I = 1; I <= 2; I++) {
                        if (LXX2S[I] < IABS(JNUC(I,1) - JNUC(I,2))
                            || LXX2S[I] > JNUC(I,1) + JNUC(I,2)) goto L560;
                    }
                    goto L600;

                L550:
                    MINLAM = MAX0(MINLAM, IABS(JNUC(IPT,1) - JNUC(IPT,2)));
                    MAXLAM = MIN0(MAXLAM, JNUC(IPT,1) + JNUC(IPT,2));
                    if (MINLAM <= MAXLAM) goto L570;
                L560:
                    std::printf("\n*** CHANNEL%5d IS TO BE COUPLED TO%5d IN SECOND ORDER WITH MULTIPOLES%3d/2%3d/2\n",
                        NAME1, NAME2, LXX2S[1], LXX2S[2]);
                    std::printf("      BUT THESE CANNOT CONNECT THE NUCLEAR SPINS.\n");
                    goto L445;

                L570: {
                    int LXX2_local = LXTOT;
                    if (LXX2_local != NOTDEF_INT) goto L572;
                    LXX2_local = (int)(2 * BETACS[1]);
                    if (IBETAC[1] == 3) goto L572;
                    if (IBETAN[1] != 3) goto L580;
                    LXX2_local = (int)(2 * BETANS[1]);
                L572:
                    if (LXX2_local >= MINLAM && LXX2_local <= MAXLAM
                        && MOD(LXX2_local + MINLAM, 4) == 0) goto L575;
                    std::printf("\n*** CHANNEL%5d IS TO BE COUPLED TO%5d IN 2ND ORDER WITH TOTAL LX =%3d/2\n",
                        NAME1, NAME2, LXX2_local);
                    std::printf("     WHICH IS OUTSIDE OF THE ALLOWED RANGE:%3d/2%3d/2\n", MINLAM, MAXLAM);
                    goto L445;
                L575:
                    MINLAM = LXX2_local;
                    MAXLAM = LXX2_local;
                    goto L600;
                }

                L580:
                    if (MINLAM == 0 && MAXLAM >= 4) MINLAM = 4;

                    // Loop over multipolarities
                L600:
                    for (LXX2 = MINLAM; LXX2 <= MAXLAM; LXX2 += 4) {

                        if (IPT == 3 && LXX2 > MINLAM) goto L700;

                        // Get the values of the betas
                        if (IORDER == 1) LXX2S[1] = LXX2;
                        RMESGN = +1;
                        for (I = 1; I <= IORDER; I++) {
                            II = IPT;
                            IORD = IORDER;
                            if (IPT == 3) IORD = 1;
                            if (IPT == 3) II = I;
                            IPTS[I] = II;
                            LX = LXX2S[I];
                            LXS[I] = LX / 2;

                            // Merge from transitions table
                            if (BETACS[I] != undef4() && BETANS[I] != undef4())
                                goto L630;
                            if (LTRANS == 0) goto L630;
                            {
                                int L1 = LTRANS + ILLOC(LTRANS + 5 - II);
                                if (L1 == LTRANS) goto L630;
                                int L2 = L1 + MTRANS * (ILLOC(LTRANS + 7 - II) - 1);
                                for (int I2 = L1; I2 <= L2; I2 += MTRANS) {
                                    if (ILLOC(I2) != LX) continue;
                                    if (BETACS[I] != undef4() ||
                                        ALLOC4(I2 + 1) == undef4()) goto L615;
                                    // Cannot accept BE for 2nd order
                                    if (IORD == 2 && ILLOC(LTRANS + 9 - II) != 0)
                                        goto L615;
                                    IBETAC[I] = ILLOC(LTRANS + 9 - II);
                                    BETACS[I] = ALLOC4(I2 + 1);
                                L615:
                                    if (BETANS[I] != undef4() ||
                                        ALLOC4(I2 + 2) == undef4()) goto L630;
                                    IBETAN[I] = 0;
                                    BETANS[I] = ALLOC4(I2 + 2);
                                    goto L630;
                                }
                            }

                        L630:
                            if (BETANS[I] == undef4()) goto L640;
                            if (IBETAN[I] != 1) goto L635;
                            BETANS[I] = BETANS[I] / (float)RNUC_local[II];
                            IBETAN[I] = 0;
                        L635:
                            if (BETACS[I] != undef4()) goto L645;
                            BETACS[I] = (float)(RNUC_local[II] / RCOUL[II]) * BETANS[I];
                            IBETAC[I] = IBETAN[I];
                            goto L645;
                        L640:
                            IBETAN[I] = -99;
                            if (BETACS[I] != undef4()) goto L645;
                            if (IORD == 2 && I == 1) goto L829;
                            IBETAC[I] = -99;
                            goto L650;
                        L645:
                            if (IBETAC[1] != 1) goto L647;
                            BETACS[I] = BETACS[I] / (float)RCOUL[II];
                            IBETAC[I] = 0;
                        L647:
                            if (BETANS[I] != undef4() || IBETAC[I] != 0) goto L650;
                            BETANS[I] = (float)(RCOUL[II] / RNUC_local[II]) * BETACS[I];

                        L650:
                            if (IORD == 2 && I == 1) goto L659_cont;
                            IUP = 1;
                            if (ENUC(II,2) > ENUC(II,1)) IUP = 2;
                            if (IUP == 2) {
                                int pw = (JNUC(II,1) - JNUC(II,2)) / 2;
                                RMESGN = RMESGN * ((pw % 2 == 0) ? 1.0 : -1.0);
                            }
                            ITEMP = LXX2;
                            if (IORD == 1) ITEMP = LXX2S[I];
                            TEMP = SCALE;
                            if (IORD == 1) TEMP = SCLFAC[I];

                            // For non-mutual we cannot allow Coulomb for LXTOT=0
                            BETTMP = BETACS[I];
                            if (MODEL == 4) BETACS[I] = 0;
                            if (IPT == 3) goto L655;
                            if (LXX2 != 0) goto L655;
                            BETACS[I] = 0;
                            IBETAC[I] = 2;

                        L655:
                            for (III = 1; III <= 2; III++) {
                                xcs[III] = BETACS[III];
                                xns[III] = BETANS[III];
                            }
                            {
                                int IZPT_local[3]; // 1-based
                                IZPT_local[1] = INTGER.IZS[2*II-1];
                                IZPT_local[2] = INTGER.IZS[2*II];
                                // Note: IZPT(II) in Fortran is just IZS(2*II-1) for IPT dimension
                                // Actually IZPT uses EQUIVALENCE (IZPT(1), IZS(2))
                                // So IZPT(1) = IZP, IZPT(2) = IZT
                                // But we need single IZ value for BEBETA: IZPT(II)
                                int IZ_val = (&INTGER.IZS[2])[II-1]; // IZPT(II)
                                BEBETA(MODEL, IORD, I, IBETAC[I], SWITCH.MEORD, &xcs[0],
                                    ITEMP, &LXX2S[0], IZ_val,
                                    AMS(2*II-1), RCOUL[II], ENUC(II,IUP) - ENUC(II,3-IUP),
                                    JNUC(II,IUP), JNUC(II,3-IUP), KNUC(II,IUP),
                                    KNUC(II,3-IUP),
                                    IBETAN[I], &xns[0], RNUC_local[II], TEMP, II,
                                    RMENUC[I], RMECOU[I], IRTN, &VALS[0][I]);
                            }
                            for (III = 1; III <= 2; III++) {
                                BETACS[III] = (float)xcs[III];
                                BETANS[III] = (float)xns[III];
                            }
                            if (MODEL == 4) BETACS[I] = (float)BETTMP;
                            if (PBUGSW)
                                std::printf(" BEBETA%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%3d%13.4g%13.4g%13.4g%13.4g%13.4g%13.4g%2d\n",
                                    NAME1, NAME2, I, II, IORD, IORDER,
                                    ITEMP, IUP, LXX2, LXX2S[1], LXX2S[2], IBETAC[I], IBETAN[I],
                                    (&INTGER.IZS[2])[II-1],
                                    (double)BETACS[1], (double)BETACS[2], (double)BETANS[1], (double)BETANS[2],
                                    RMENUC[I], RMECOU[I], IRTN);
                            if (IRTN != 0) goto L445;

                            // If both RME's are still undefined, must give up
                            if (DABS(RMENUC[I]) + DABS(RMECOU[I]) == 0) goto L829;

                            if (IORD == 2) goto L659_cont;
                            IBETAN[I] = 0;
                            IBETAC[I] = 0;

                            DELTAN[I] = BETANS[I] * (float)RNUC_local[II];

                        L659_cont:;
                        } // end for I = 1, IORDER

                        // Print coupling lines
                        if (IORD == 2) goto L665;
                        for (I = 1; I <= IORDER; I++) {
                            char iptwrd_buf[5]; std::memcpy(iptwrd_buf, IPTWRD[IPTS[I]].data, 4); iptwrd_buf[4] = '\0';
                            if (I == 1) {
                                std::printf(" %5d%5d %.8s%3d", NAME2, NAME1, MODNAM.data, IORDER);
                            }
                            std::printf("%23s%s%3d%9.4f%9.4f%8.4f%11.4g%9.4f%8.4f%12.4g%12.4g%9.4f%13.4g\n",
                                "", iptwrd_buf, LXS[I], SCLFAC[I],
                                (double)BETANS[I], (double)DELTAN[I], RMENUC[I],
                                (double)BETACS[I], VALS[1][I], VALS[2][I], VALS[3][I], VALS[4][I], VALS[5][I]);
                        }
                        goto L700;
                    L665:
                        if (LXX2 > MINLAM) goto L670;
                        { char iptwrd_buf[5]; std::memcpy(iptwrd_buf, IPTWRD[IPTS[1]].data, 4); iptwrd_buf[4] = '\0';
                        std::printf(" %5d%5d %.8s%3d%23s%s%9s%9.4f%3s LX1 =%2d%28sLX2 =%2d\n",
                            NAME2, NAME1, MODNAM.data, IORDER, "", iptwrd_buf, "", SCALE, "", LXS[1], "", LXS[2]); }
                        if (IBETAN[1] == 0)
                            std::printf("+%53sBETAN =%7.4f%30sBETAN =%7.4f\n",
                                "", (double)BETANS[1], "", (double)BETANS[2]);
                        if (IBETAC[1] == 0)
                            std::printf("+%69sBETAC =%7.4f%30sBETAC =%7.4f\n",
                                "", (double)BETACS[1], "", (double)BETACS[2]);
                    L670:
                        LX = LXX2 / 2;
                        std::printf("%31s%3d  (TOTAL)%29s%12.4g%16s%12.4g%12.4g%9.4f%13.4g\n",
                            "", LX, "", RMENUC[2], "", VALS[2][2], VALS[3][2], VALS[4][2], VALS[5][2]);

                    L700:
                        CUPAB(IORDER, IPT, &LXX2S[0], LXX2, &RMENUC[0], &RMECOU[0],
                            RNUC_local, &INTGER.IZS[2], RMEGN, RMEGC, IPT1);

                        // Setup the INGRST indicator
                        IPT1 = IORDER;
                        if (MODEL <= 2) goto L710;
                        if (IPT == 1) IPT1 = 3 + 50 * LXX2S[1];
                        if (IPT == 2) IPT1 = 4 + 5000 * LXX2S[1];
                        if (IPT == 3) IPT1 = 5 + 50 * LXX2S[1] + 5000 * LXX2S[2];
                        if (MODEL == 4) IPT1 = 6 + (int)(100 * BETTMP);
                    L710:
                        if (PBUGSW)
                            std::printf(" CUPAB%15.4g%15.4g%15.4g%15.4g%15.4g%15.4g%15.4g%5d\n",
                                RMENUC[1], RMENUC[2], RMECOU[1], RMECOU[2], RMEGN, RMEGC,
                                RMESGN, IPT1);

                        RMEGN = RMESGN * RMEGN;
                        RMEGC = RMESGN * RMEGC;

                        // Loop over basis states for channel 1
                        for (IB1 = IB1ST; IB1 <= IB1END; IB1++) {
                            IB1PT = LBASDF_local + (IB1 - 1) * MBASDF;
                            SCHN[1] = ILLOC(IB1PT + 1);
                            int L1 = ILLOC(IB1PT + 2);

                            // Loop over all basis states coupled to
                            for (IB2 = IB2ST; IB2 <= IB2END; IB2++) {
                                IB2PT = LBASDF_local + (IB2 - 1) * MBASDF;
                                SCHN[2] = ILLOC(IB2PT + 1);
                                int L2 = ILLOC(IB2PT + 2);

                                LDIFX2 = L2 - L1;
                                if (IABS(LDIFX2) > LXX2 ||
                                    IABS(SCHN[1] - SCHN[2]) > LXX2) continue;
                                if (SCHN[1] + SCHN[2] < LXX2) continue;

                                // Setup the symmetry indicator
                                FAC = 1;
                                if (!SYMSW) goto L740;
                                if (ILLOC(IB1PT + 6) == 2) goto L740;
                                if (ILLOC(IB2PT + 6) != 2) goto L730;
                                FAC = 2;
                                if (IPT != 3) goto L740;
                                I = LXX2S[1] + LXX2S[2] + SCHN[2];
                                if (BTEST(I, 1)) continue;
                                goto L740;

                            L730:
                                if (IPT != 3 || (IPT == 3 && (LXX2S[1] != LXX2S[2])))
                                    FAC = 2;

                            L740:
                                if (!FLIPSW) goto L750;
                                I = SCHN[ICFLIP];
                                if (BTEST(I, 1)) FAC = -FAC;

                            L750:
                                LDIF = LDIFX2 / 2;
                                IFAC = MOD(LDIF + 4000, 4) + 1;
                                IFAC = IFACS[IFAC];

                                CUPSPN(IRTN, LXX2, IPT, SCHN,
                                    LXX2S, JNUC_flat, RMEGN, RMEGC, XN, XC);

                                if (PBUGSW)
                                    std::printf(" CUPSPN%3d%3d%3d%3d%3d%3d%3d/2%3d/2%3d/2%3d/2%3d/2%3d/2%3d/2%3d/2%2d%2d%3d%2d%3d/2%3d/2%3d/2%3d%3d%15.5g%15.5g%3.0f\n",
                                        NAME1, NAME2, IB1, IB2, L1, L2,
                                        JT1, JP1, JT2, JP2, SCHN[1], SCHN[2],
                                        IPARI1, IPARI2, ICP, IORDER,
                                        LXX2S[1], LXX2S[2], LXX2, IRTN, IPT1, XN, XC, FAC);
                                if (IRTN != 0) goto L770;
                                IRET = 0;
                                continue; // goto L789

                            L770:
                                ICNT2 = ICNT2 + 1;
                                { int IZP_val = INTGER.IZS[2]; // IZPT(1)
                                  int IZT_val = INTGER.IZS[3]; // IZPT(2)
                                  if (IZP_val * IZT_val == 0) XC = 0; }
                                if (DABS(XN) + DABS(XC) == 0) continue;
                                ICOUNT = ICOUNT + 1;
                                ILLOC(LBASCP) = IB1;
                                ILLOC(LBASCP + 1) = IB2;
                                ILLOC(LBASCP + 2) = IORDER;
                                ILLOC(LBASCP + 3) = LXX2;
                                ALLOC4(LBASCP + 4) = (float)(FAC * XN);
                                ALLOC4(LBASCP + 5) = (float)(FAC * XC);
                                ILLOC(LBASCP + 6) = IFAC;
                                ILLOC(LBASCP + 7) = IPT1;
                                ILLOC(LBASCP + 8) = NOTDEF_INT;
                                ILLOC(LBASCP + 9) = NOTDEF_INT;
                                ILLOC(LBASCP + 10) = NOTDEF_INT;

                                // Accumulate maximum order for INGRST
                                if (IORDER == 2) CCBLK.IORDX = 2;

                                LBASCP = LBASCP + MBASCP;
                                NBASCP_local = NBASCP_local + 1;

                                // For C.C. we want LXTOP = MAX(LOUT-LIN, MULTIPOLE)
                                LXTOP = MAX0(LXTOP, LXX2 / 2);

                            } // end for IB2
                        } // end for IB1

                        if (IORDER != 1) goto L829;

                        // For first order, betas must be found again for each LX
                        BETANS[1] = undef4();
                        BETACS[1] = undef4();

                    L829:;
                    } // end for LXX2

                    if (ICNT2 != 0) goto L840;
                    std::printf("\n**** CHANNEL%5d IS TO BE COUPLED TO%5d BY ORDER =%2d AND MULTIPOLE =%3d/2\n",
                        NAME1, NAME2, IORDER, IMUL);
                    std::printf("      (****/2 MEANS MULTIPOLE WAS NOT SPECIFIED; ALL POSSIBILITIES WERE TRIED.)\n");
                    std::printf("      BUT THERE IS NO DEFINITION OF THE REQUIRED BETAS, OR THE SPECIFIED MULTIPOLE IS IMPOSSIBLE TO USE.\n");
                    std::printf("      (CHECK TRANSISTIONS SPECIFICATIONS.)\n");
                    std::printf("      (OR PERHAPS PARITY IS VIOLATED)\n");
                    std::printf("      YOU NEVER KNOW FOR SURE.\n");
                    IRET = 0;
                    goto L949;
                L840:
                    if (ICOUNT != 0) goto L949;
                    std::printf("\n**** CHANNEL%5d IS TO BE COUPLED TO%5d BY ORDER =%2d AND MULTIPOLE =%3d/2\n",
                        NAME1, NAME2, IORDER, IMUL);
                    std::printf("      (****/2 MEANS MULTIPOLE WAS NOT SPECIFIED; ALL POSSIBILITIES WERE TRIED.)\n");
                    std::printf("      BUT ALL COUPLING CONSTANTS WERE ZERO (SUSPISCOUS)\n");

                L949:
                } // end for ICP

            } // end for ICH2
        } // end for ICH1

        X = DFLOAT(NBASCP_local) / NBASDF;
        std::printf("\n THERE ARE A TOTAL OF%6d COUPLINGS OF BASIS STATES FOR AN AVERAGE OF%6.1f COUPLINGS PER BASIS STATE.\n",
            NBASCP_local, X);

        IWRDS = (MBASCP * NBASCP_local) / ALLOCS.FACFR4 + 1;
        IL = CHOPIT(IWRDS);

        // Was every coupling specification used
        for (ICP = 1; ICP <= NCHNCP; ICP++) {
            ICPPT = (ICP - 1) * MCHNCP + LCHNCP;
            if (ILLOC(ICPPT + 13) != 0) continue;
            std::printf("\n*** CHANNEL%5d IS TO BE COUPLED TO%5d BUT ONE OF THESE CHANNELS IS NOT DEFINED.\n",
                ILLOC(ICPPT), ILLOC(ICPPT + 1));
            IRET = 0;
        }

        // Prepare to sort the basis coupling array
        if (IRET == 0) return;
        int ITEMP_sort = NALLOC("SORTINDX", NBASCP_local);
        CCBLK.IBASCP = NALLOC("BASISCUP", (NBASCP_local * MBASCP + LPL + 1) / ALLOCS.FACFR4);
        LBASCP = ALLOCS.FACFR4 * (LOCPTRS.Z[CCBLK.IBASCP] - 1) + 1;
        ILLOC(LBASCP) = NBASCP_local;
        ILLOC(LBASCP + 1) = MBASCP;
        ILLOC(LBASCP + 2) = LPL;
        LBASCP = LBASCP + LPL;
        LTMPCP = ALLOCS.FACFR4 * (LOCPTRS.Z[ITMPCP] - 1) + 1;
        LTEMP = ALLOCS.FACFR4 * (LOCPTRS.Z[ITEMP_sort] - 1);
        LTMP2 = (LBASCP + 1) / ALLOCS.FACFR4;

        // Construct the sort parameter
        for (IBC = 1; IBC <= NBASCP_local; IBC++) {
            IBCPT = LTMPCP + (IBC - 1) * MBASCP;
            ALLOC(LTMP2 + IBC) = ILLOC(IBCPT + 3)
                + 100 * (ILLOC(IBCPT + 2) + 10 * (ILLOC(IBCPT + 1)
                + 5000 * ILLOC(IBCPT)));
        }

        IXSORT(&ILLOC(LTEMP), &ALLOC(LTMP2), NBASCP_local);

        // Straighten it out
        LCHNDF = LOCPTRS.Z[CCBLK.ICHNDF];
        LCHNDF = ALLOCS.FACFR4 * (LCHNDF - 1) + 1;
        LCHNDF = ILLOC(LCHNDF + 2) + LCHNDF;
        LBASDF_local = LOCPTRS.Z[CCBLK.IBASDF];
        LBASDF_local = ALLOCS.FACFR4 * (LBASDF_local - 1) + 1;
        LBASDF_local = LBASDF_local + ILLOC(LBASDF_local + 2);
        IB1 = 1;
        ILLOC(LBASDF_local + 3) = 1;
        for (IBC = 1; IBC <= NBASCP_local; IBC++) {
            IBCPT = LTMPCP + (ILLOC(LTEMP + IBC) - 1) * MBASCP;
            for (I = 1; I <= MBASCP; I++) {
                ILLOC(LBASCP + I - 1) = ILLOC(IBCPT + I - 1);
            }
            LBASCP = LBASCP + MBASCP;

            // When we change basis states, update pointers
            if (IB1 == ILLOC(IBCPT)) continue;
            ILLOC(LBASDF_local + 4) = IBC - ILLOC(LBASDF_local + 3);
            if (ILLOC(LBASDF_local + 4) > 0) goto L1350;
            std::printf("\n*** BASIS STATE%5d OF CHANNEL%5d IS NOT COUPLED TO ANYTHING.\n",
                IB1, ILLOC(LCHNDF + (ILLOC(LBASDF_local) - 1) * MCHNDF));
            IRET = 0;
        L1350:
            IB1 = IB1 + 1;
            LBASDF_local = LBASDF_local + MBASDF;
            ILLOC(LBASDF_local + 3) = IBC;
            if (IB1 >= ILLOC(IBCPT)) continue;
            std::printf("\n*** BASIS STATE%5d OF CHANNEL%5d IS NOT COUPLED TO ANYTHING.\n",
                IB1, ILLOC(LCHNDF + (ILLOC(LBASDF_local) - 1) * MCHNDF));
            IRET = 0;
            goto L1350;
        }
        ILLOC(LBASDF_local + 4) = NBASCP_local - ILLOC(LBASDF_local + 3) + 1;
        if (IB1 == NBASDF) goto L1400;
        IB1 = IB1 + 1;
        std::printf("\n*** BASIS STATES%5d THROUGH%5d ARE NOT COUPLED TO ANYTHING.\n",
            IB1, NBASDF);
        IRET = 0;

    L1400:
        if (MOD(INTGER.IPRINT / 100, 10) > 0) return;
        LOCPTRS.Z[ITEMP_sort] = -LOCPTRS.Z[ITEMP_sort];
        LOCPTRS.Z[ITMPCP] = -LOCPTRS.Z[ITMPCP];
    } // end of NBASCP_local scope

    #undef JP1
    #undef JP2
    #undef JT1
    #undef JT2
    #undef EP1
    #undef EP2
    #undef ET1
    #undef ET2
    #undef NAME1
    #undef NAME2
    return;
}

// ============================================================================
// CCHAN — Channel specification input parser
// Translated from source.f lines 6543-6759
// ============================================================================
void CCHAN(int& IRET)
{
    extern int GIVEAL(int IWRDS);
    extern int CHOPIT(int IWRDS);

    const char PLUS  = '+';
    const char MINUS = '-';

    int NUMKEY = 7;
    char8 KEYS[8]; // 1-based
    KEYS[1] = char8("NAME");
    KEYS[2] = char8("JT");
    KEYS[3] = char8("ET");
    KEYS[4] = char8("JP");
    KEYS[5] = char8("EP");
    KEYS[6] = char8("KT");
    KEYS[7] = char8("KP");

    int LENDAT = 17;
    int LEKEY = 3;
    int LAL = 5;
    char8 NAME8("CHANNELS");
    char8 AKEY;

    int INAME, IJT, IET, IJP, IEP, IKP, IKT, ICH;
    int IORDER, I, INT_val, IWRDS, IL, LABEL, LDATA, LDATA4;
    double TEMP;
    char IPM;
    int r;

    INTRNL.ICD3I = 0;
    IRET = 0;
    INAME = 0;
    IJT = 0;
    IET = 0;
    IJP = 0;
    IEP = 0;
    IKP = 0;
    IKT = 0;
    ICH = 0;
    I = NAMLOC(NAME8.data);
    if (I != 0) LOCPTRS.Z[I] = -LOCPTRS.Z[I];
    CCBLK.ICHNDF = GIVEAL(IWRDS);
    LDATA = LOCPTRS.Z[CCBLK.ICHNDF];
    NAMCOM.NAMES[CCBLK.ICHNDF] = NAME8;
    LDATA4 = ALLOCS.FACFR4 * (LDATA - 1) + 1;
    ILLOC(LDATA4 + 1) = LENDAT;
    ILLOC(LDATA4 + 2) = LAL;
    ILLOC(LDATA4 + 3) = 0;
    ILLOC(LDATA4 + 4) = 0;
    LABEL = LDATA4 + LAL;
    r = FITKEY(AKEY.data);
    if (r == 1) goto L90;
    if (r == 2 || r == 3 || r == 4) goto L900;
    goto L900;

    // ( ENCOUNTERED LETS GO
L90:
    IORDER = 0;
L100:
    r = NXWORD(AKEY.data);
    if (r == 1) goto L699;
    if (r == 2 || r == 3) goto L900;
    for (I = 1; I <= NUMKEY; I++) {
        if (KEYS[I] != AKEY) continue;
        switch (I) {
            case 1: goto L200;
            case 2: goto L250;
            case 3: goto L300;
            case 4: goto L350;
            case 5: goto L400;
            case 6: goto L410;
            case 7: goto L420;
        }
    }
    goto L900;

L200:
    IORDER = IORDER + 1;
    INAME = IORDER;
    if (INAME != 1) goto L930;
    goto L100;
L250:
    IORDER = IORDER + 1;
    IJT = IORDER;
    goto L100;
L300:
    IORDER = IORDER + 1;
    IET = IORDER;
    goto L100;
L350:
    IORDER = IORDER + 1;
    IJP = IORDER;
    goto L100;
L400:
    IORDER = IORDER + 1;
    IEP = IORDER;
    goto L100;
L410:
    IORDER = IORDER + 1;
    IKT = IORDER;
    goto L100;
L420:
    IORDER = IORDER + 1;
    IKP = IORDER;
    goto L100;

L699:
    if (IORDER < LEKEY) goto L940;
    LABEL = LABEL - LENDAT;

    // Advance to the next channel
L700:
    LABEL = LABEL + LENDAT;
    ILLOC(LABEL + 1) = 0;
    ILLOC(LABEL + 3) = 0;
    ALLOC4(LABEL + 2) = 0.0f;
    ALLOC4(LABEL + 4) = 0.0f;
    ALLOC4(LABEL + 5) = undef4();
    ALLOC4(LABEL + 6) = undef4();
    for (I = 8; I <= LENDAT; I++) {
        ILLOC(LABEL + I - 1) = NOTDEF_INT;
    }
    // KT and KP are zero by default
    ILLOC(LABEL + 15) = 0;
    ILLOC(LABEL + 16) = 0;
    for (I = 1; I <= IORDER; I++) {
        if (I == INAME) goto L720;
        if (I == IJT) goto L730;
        if (I == IJP) goto L740;
        if (I == IKT) goto L750;
        if (I == IKP) goto L760;
        r = NXVALF(TEMP);
        if (r != 0) goto L910;
        if (I == IET) ALLOC4(LABEL + 2) = (float)TEMP;
        if (I == IEP) ALLOC4(LABEL + 4) = (float)TEMP;
        goto L799;
    L720:
        { int iret2;
          NXINT(INT_val, iret2);
          if (iret2 == 1 || iret2 == 2) goto L910;
          if (iret2 == 3) goto L800;
        }
        ILLOC(LABEL) = INT_val;
        goto L799;
    L730:
        r = NXHINT(INT_val);
        if (r != 0) goto L910;
        IPM = INPBUF.IBUF[INPBUF.INCH - 1];
        if (IPM == PLUS) ILLOC(LABEL + 1) = INT_val;
        if (IPM == MINUS) ILLOC(LABEL + 1) = -INT_val;
        goto L799;
    L740:
        r = NXHINT(INT_val);
        if (r != 0) goto L910;
        IPM = INPBUF.IBUF[INPBUF.INCH - 1];
        if (IPM == PLUS) ILLOC(LABEL + 3) = INT_val;
        if (IPM == MINUS) ILLOC(LABEL + 3) = -INT_val;
        goto L799;
    L750:
        r = NXHINT(ILLOC(LABEL + 15));
        if (r != 0) goto L910;
        goto L799;
    L760:
        r = NXHINT(ILLOC(LABEL + 16));
        if (r != 0) goto L910;
    L799:;
    }
    ICH = ICH + 1;
    goto L700;

L800:
    r = FITKEY(AKEY.data);
    if (r == 1 || r == 2 || r == 4) goto L900;
    if (r == 3) goto L810;
    goto L900;
L810:
    ILLOC(LDATA4) = ICH;
    IWRDS = (LABEL - LDATA4) / ALLOCS.FACFR4 + 1;
    IL = CHOPIT(IWRDS);
    IRET = 1;
    return;

L900:
    std::printf("\n****  CHANNEL DEFINITIONS NOT ENCLOSED IN PARENTHESES OR A KEYWORD IS MISPELLED; CHANNELS KEYWORDS ARE:\n\n");
    for (I = 1; I <= NUMKEY; I++) std::printf(" %.8s ", KEYS[I].data);
    std::printf("\n");
    return;
L910:
    std::printf("\n****  ERROR IN READING%5d NUMERICAL DATA IN SET.\n", I);
    return;
L930:
    std::printf("\n****  CHANNEL SPECIFICATION MUST BE FIRST KEYWORD\n");
    return;
L940:
    std::printf("\n*** AT LEAST THREE KEYWORDS ARE REQUIRED.\n");
    return;
}

// ============================================================================
// CCOUP — Coupling scheme input parser
// Translated from source.f lines 7300-7951
// ============================================================================
void CCOUP(int& IRET)
{
    extern int GIVEAL(int IWRDS);
    extern int CHOPIT(int IWRDS);

    int NUMKEY = 16;
    char8 KEYS[17]; // 1-based
    KEYS[1] = char8("COUPLED");
    KEYS[2] = char8("ORDER");
    KEYS[3] = char8("ORDERMAX");
    KEYS[4] = char8("TYPE");
    KEYS[5] = char8("MULTIPOL");
    KEYS[6] = char8("LX");
    KEYS[7] = char8("NUCLEUS");
    KEYS[8] = char8("MODEL");
    KEYS[9] = char8("SCALE");
    KEYS[10] = char8("SCALEFAC");
    KEYS[11] = char8("TOTALLX");
    KEYS[12] = char8("TOTALMUL");
    KEYS[13] = char8("MAXLX");
    KEYS[14] = char8("MAXIMUML");
    KEYS[15] = char8("MAXMULTI");
    KEYS[16] = char8("MAXIMUMM");

    int NMBCKY = 21;
    char8 BCKEYS[22]; // 1-based
    BCKEYS[1] = char8("BETAC"); BCKEYS[2] = char8("BETACOUL");
    BCKEYS[3] = char8("DELTAC"); BCKEYS[4] = char8("DELTACOU");
    BCKEYS[5] = char8("BE"); BCKEYS[6] = char8("BELX"); BCKEYS[7] = char8("BEUP");
    BCKEYS[8] = char8("Q"); BCKEYS[9] = char8("QMOMENT"); BCKEYS[10] = char8("QUADRUPO");
    BCKEYS[11] = char8("HALFLIFE"); BCKEYS[12] = char8("MEANLIFE");
    BCKEYS[13] = char8("BEDOWN"); BCKEYS[14] = char8("RMEC"); BCKEYS[15] = char8("RMECOULO");
    BCKEYS[16] = char8("BEUPWU"); BCKEYS[17] = char8("BEDOWNWU");
    BCKEYS[18] = char8("RMECALPH"); BCKEYS[19] = char8("LX");
    BCKEYS[20] = char8("TOTALLX"); BCKEYS[21] = char8("WIDTH");
    int IBCKYS[22] = {0, 0, 0, 1, 1, 10, 10, 10, 11, 11, 11, 4, 5, 6, 7, 7, 8, 9, 2, 3, 3, 12}; // 1-based

    int NMBNKY = 9;
    char8 BNKEYS[10]; // 1-based
    BNKEYS[1] = char8("BETAN"); BNKEYS[2] = char8("BETANUCL");
    BNKEYS[3] = char8("DELTAN"); BNKEYS[4] = char8("DELTANUC");
    BNKEYS[5] = char8("RMEN"); BNKEYS[6] = char8("RMENUCLE");
    BNKEYS[7] = char8("RMENALPH"); BNKEYS[8] = char8("LX");
    BNKEYS[9] = char8("TOTALLX");
    int IBNKYS[10] = {0, 0, 0, 1, 1, 2, 2, 2, 3, 3}; // 1-based

    int LENDAT = 21;
    int LAL = 3;
    char8 NAME8("CHANCOUP");

    int LISTN = 4;
    char8 KEY1[5]; // 1-based
    KEY1[1] = char8("UPDOWN"); KEY1[2] = char8("UPONLY");
    KEY1[3] = char8("DIAGONAL"); KEY1[4] = char8("BORN");
    int KEY1VL[5] = {0, 1, 2, 1, 2}; // 1-based

    int LISTN1 = 10;
    char8 KEY2[11]; // 1-based
    KEY2[1] = char8("P"); KEY2[2] = char8("T"); KEY2[3] = char8("PT");
    KEY2[4] = char8("TP"); KEY2[5] = char8("PP"); KEY2[6] = char8("TT");
    KEY2[7] = char8("PROJ"); KEY2[8] = char8("PROJECTI");
    KEY2[9] = char8("TARG"); KEY2[10] = char8("TARGET");
    int KEY2VL[11] = {0, 1, 2, 3, 4, 1, 2, 1, 1, 2, 2}; // 1-based

    int NUMMOD = 7;
    char8 MODELS_local[8]; // 1-based
    MODELS_local[1] = char8("ROTATION"); MODELS_local[2] = char8("VIBRATIO");
    MODELS_local[3] = char8("DEFORMED"); MODELS_local[4] = char8("HOTFUDGE");
    MODELS_local[5] = char8("MODEL5"); MODELS_local[6] = char8("MODEL6");
    MODELS_local[7] = char8("MODEL7");

    int NUMFLP = 6;
    int IFLIPS[7] = {0, 4, 6, 8, 11, 15, 17}; // 1-based

    float DFLTS[31]; // 1-based
    int* IDFLTS = reinterpret_cast<int*>(DFLTS); // EQUIVALENCE
    char8 BLANK(" ");
    char8 AKEY, AKEY2, AKEY3;
    bool EQULSW, SCNDSW;

    int ICOP, IOR, ITY, IMUL, IBC, IBN, INC, ICH, IMD, ISC, ITOTMU, IMAXMU;
    int IORDER, I, I1, I2, II, J, INT_val, IWRDS, IL, LABEL, LDATA, LDATA4, JORDER;
    int IPT, ITEMP, L;
    double TEMP, VAL;
    int r;

    IRET = 0;
    ICOP = 0; IOR = 0; ITY = 0; IMUL = 0;
    IBC = 0; IBN = 0; INC = 0; ICH = 0;
    IMD = 0; ISC = 0; ITOTMU = 0; IMAXMU = 0;

    // Initialize the default values
    for (I = 1; I <= 30; I++) {
        IDFLTS[I] = NOTDEF_INT;
    }
    // Defined defaults
    IDFLTS[2] = 1; // ORDER=1
    IDFLTS[3] = 1; // TYPE=UPDOWN
    IDFLTS[13] = 0; // USE=0
    IDFLTS[14] = 1; // MODEL=ROTATIONAL
    DFLTS[17] = 1.0f; // SCALE=1
    DFLTS[18] = 1.0f;

    I = NAMLOC(NAME8.data);
    if (I != 0) LOCPTRS.Z[I] = -LOCPTRS.Z[I];
    CCBLK.ICHNCP = GIVEAL(IWRDS);
    LDATA = LOCPTRS.Z[CCBLK.ICHNCP];
    NAMCOM.NAMES[CCBLK.ICHNCP] = NAME8;
    LDATA4 = ALLOCS.FACFR4 * (LDATA - 1) + 1;
    ILLOC(LDATA4 + 2) = LAL;
    LABEL = LDATA4 + LAL;
    r = FITKEY(AKEY.data);
    if (r == 1) goto L90;
    if (r == 2 || r == 3 || r == 4) goto L900;
    goto L900;

    // ( ENCOUNTERED LETS GO
L90:
    IORDER = 0;
    AKEY2 = BLANK;

L100:
    AKEY = AKEY2;
    if (AKEY2 == BLANK) {
        r = NXWORD(AKEY.data);
        if (r == 1) goto L400;
        if (r == 2 || r == 3) goto L900;
    }

    // Look for possible = sign after this keyword
    EQULSW = false;
    r = FITKEY(AKEY2.data);
    if (r == 1) goto L900;
    if (r == 2) goto L120;
    if (r == 3) goto L900;
    if (r == 4) goto L130;
    goto L150;

L120:
    EQULSW = true;
L130:
    AKEY2 = BLANK;

    // Identify the keyword
L150:
    for (I = 1; I <= NUMKEY; I++) {
        if (KEYS[I] != AKEY) continue;
        switch (I) {
            case 1: goto L200;
            case 2: case 3: goto L210;
            case 4: goto L220;
            case 5: case 6: goto L230;
            case 7: goto L250;
            case 8: goto L260;
            case 9: case 10: goto L270;
            case 11: case 12: goto L280;
            case 13: case 14: case 15: case 16: goto L290;
        }
    }

    // Is it a default coulomb beta specification
    for (I = 1; I <= NMBCKY; I++) {
        if (AKEY != BCKEYS[I]) continue;
        IDFLTS[11] = IBCKYS[I];
        IDFLTS[12] = IBCKYS[I];
        IORDER = IORDER + 1;
        IBC = IORDER;
        if (!EQULSW) goto L100;
        { double tmp6 = DFLTS[6]; r = NXVAL(tmp6); DFLTS[6] = (float)tmp6; }
        if (r != 0) goto L945;
        DFLTS[7] = DFLTS[6];
        goto L100;
    }

    // Is it a default nuclear beta
    for (I = 1; I <= NMBNKY; I++) {
        if (AKEY != BNKEYS[I]) continue;
        IDFLTS[15] = IBNKYS[I];
        IDFLTS[16] = IBNKYS[I];
        IORDER = IORDER + 1;
        IBN = IORDER;
        if (!EQULSW) goto L100;
        { double tmp8 = DFLTS[8]; r = NXVAL(tmp8); DFLTS[8] = (float)tmp8; }
        if (r != 0) goto L945;
        DFLTS[9] = DFLTS[8];
        goto L100;
    }

    goto L900;

    // COUPLED
L200:
    IORDER = IORDER + 1;
    ICOP = IORDER;
    if (ICOP != 1) goto L930;
    goto L100;

    // ORDER
L210:
    IORDER = IORDER + 1;
    IOR = IORDER;
    if (IOR != 2) goto L940;
    if (EQULSW) {
        int iret2;
        NXINT(IDFLTS[2], iret2);
        if (iret2 != 0) goto L945;
    }
    goto L100;

    // TYPE
L220:
    IORDER = IORDER + 1;
    ITY = IORDER;
    if (!EQULSW) goto L100;
    r = NXWORD(AKEY3.data);
    if (r != 0) goto L945;
    AKEY = AKEY3;
    goto L472;

    // MULTIPOLARITY
L230:
    IORDER = IORDER + 1;
    IMUL = IORDER;
    if (!EQULSW) goto L100;
    r = NXHINT(INT_val);
    if (r != 0) goto L945;
    IDFLTS[4] = INT_val;
    IDFLTS[5] = INT_val;
    goto L100;

    // NUCLEUS
L250:
    IORDER = IORDER + 1;
    INC = IORDER;
    if (!EQULSW) goto L100;
    r = NXWORD(AKEY3.data);
    if (r != 0) goto L945;
    AKEY = AKEY3;
    goto L555;

    // MODEL
L260:
    IORDER = IORDER + 1;
    IMD = IORDER;
    if (!EQULSW) goto L100;
    r = NXWORD(AKEY3.data);
    if (r != 0) goto L945;
    AKEY = AKEY3;
    goto L572;

    // SCALE FACTOR
L270:
    IORDER = IORDER + 1;
    ISC = IORDER;
    if (!EQULSW) goto L100;
    { double dval; r = NXVAL(dval); DFLTS[17] = (float)dval; }
    if (r != 0) goto L945;
    DFLTS[18] = DFLTS[17];
    goto L100;

    // TOTALLX
L280:
    IORDER = IORDER + 1;
    ITOTMU = IORDER;
    if (!EQULSW) goto L100;
    r = NXHINT(INT_val);
    if (r != 0) goto L945;
    IDFLTS[19] = INT_val;
    goto L100;

    // MAXIMUMLX
L290:
    IORDER = IORDER + 1;
    IMAXMU = IORDER;
    if (!EQULSW) goto L100;
    r = NXHINT(INT_val);
    if (r != 0) goto L945;
    IDFLTS[20] = INT_val;
    goto L100;

    // Done with keywords, go on to data
L400:
    ILLOC(LDATA4 + 1) = LENDAT;
    LABEL = LABEL - LENDAT;
    EQULSW = false;

    // Go to next batch of data
L410:
    LABEL = LABEL + LENDAT;
    II = 1;
L420_loop:
    { int iret2;
      NXINT(INT_val, iret2);
      if (iret2 == 1) goto L910;
      if (iret2 == 2) goto L420_loop;
      if (iret2 == 3) goto L800;
    }
    ILLOC(LABEL) = INT_val;
    I1 = 2;
    II = 2;
    { int iret2;
      NXINT(INT_val, iret2);
      if (iret2 != 0) goto L910;
    }
    ILLOC(LABEL + 1) = -INT_val;
    if (INT_val > 0) goto L915;
    ICH = ICH + 1;

    // Move in the defaults
    for (I = 2; I <= 5; I++) {
        ILLOC(LABEL + I) = IDFLTS[I];
    }
    for (I = 6; I <= 9; I++) {
        ALLOC4(LABEL + I) = DFLTS[I];
    }
    L = LENDAT - 1;
    for (I = 10; I <= L; I++) {
        ILLOC(LABEL + I) = IDFLTS[I];
    }
    JORDER = ILLOC(LABEL + 2);
    SCNDSW = JORDER == 2;

    // Process the coupling order if specified
    if (IOR == 0) goto L440;
    II = 3;
    I1 = 3;
    { int iret2;
      NXINT(JORDER, iret2);
      if (iret2 == 1) goto L910;
      if (iret2 == 2) goto L410;
      if (iret2 == 3) goto L810;
    }
    ILLOC(LABEL + 2) = JORDER;
L440:
    if (JORDER != 1 && JORDER != 2) goto L925;
    SCNDSW = JORDER == 2;
    if (SCNDSW && (IMUL == 0)) goto L950;
    if (SCNDSW && (INC == 0)) goto L980;

    // Process the keywords after order
    if (I1 > IORDER) goto L410;
L450:
    II = II + 1;
    if (I1 == IMUL) goto L490;
    if (I1 == IBC) goto L500;
    if (I1 == IBN) goto L530;
    if (I1 == INC) goto L550;
    if (I1 == IMD) goto L570;
    if (I1 == ITY) goto L470;
    if (I1 == ISC) goto L580;
    if (I1 == ITOTMU) goto L590;
    if (I1 == IMAXMU) goto L600;

    std::printf("\n*** THE%3d INPUT FIELD CANNOT BE CLASSIFIED.  DID YOU SPECIFY THE SAME TYPE OF DATA TWICE IN THE KEYWORD LIST?\n", I1);
    goto L908;

    // COUPLING TYPE
L470:
    r = NXWORD(AKEY.data);
    if (r == 1) goto L960;
    if (r == 2) goto L700;
    if (r == 3) goto L810;
L472:
    for (J = 1; J <= LISTN; J++) {
        if (KEY1[J] != AKEY) continue;
        INT_val = KEY1VL[J];
        if (EQULSW) { IDFLTS[3] = INT_val; goto L100; } // goto L225
        ILLOC(LABEL + 3) = INT_val;
        goto L689;
    }
    goto L960;

    // MULTIPOLARITY
L490:
    r = NXHINT(INT_val);
    if (r == 1) goto L910;
    if (r == 2) goto L700;
    if (r == 3) goto L810;
    ILLOC(LABEL + 4) = INT_val;
    if (INT_val < 0) goto L935;
    if (!SCNDSW) goto L689;
    II = II + 1;
    r = NXHINT(INT_val);
    if (r == 1) goto L910;
    if (r == 2 || r == 3) goto L950;
    ILLOC(LABEL + 5) = INT_val;
    if (INT_val < 0) goto L935;
    goto L689;

    // BETA COULOMB
L500:
    for (J = 1; J <= JORDER; J++) {
    L505:
        r = NXVALF(TEMP);
        if (r == 0) goto L525;
        if (r == 1) goto L510;
        if (r == 2) goto L595;
        if (r == 3) goto L805;
        if (r == 4) goto L805;
        goto L525;
    L510:
        r = NXWORD(AKEY.data);
        if (r == 1) goto L910;
        if (r == 2) goto L595;
        if (r == 3) goto L805;
        for (I2 = 1; I2 <= NMBCKY; I2++) {
            if (AKEY != BCKEYS[I2]) continue;
            ILLOC(LABEL + 10 + J) = IBCKYS[I2];
            goto L505;
        }
        std::printf("\n*** INVALID KEYWORD FOR SPECIFYING COULOMB DEFORMATION: %.8s,  VALID KEYWORDS ARE:\n", AKEY.data);
        for (I2 = 1; I2 <= NMBCKY; I2++) std::printf("      %.8s ", BCKEYS[I2].data);
        std::printf("\n");
        return;
    L525:
        ALLOC4(LABEL + 5 + J) = (float)TEMP;
        II = II + 1;
    }
    goto L689;

    // NUCLEAR DEFORMATION OR BETA
L530:
    for (J = 1; J <= JORDER; J++) {
    L535:
        r = NXVALF(TEMP);
        if (r == 0) goto L545;
        if (r == 1) goto L537;
        if (r == 2) goto L595;
        if (r == 3) goto L805;
        if (r == 4) goto L805;
        goto L545;
    L537:
        r = NXWORD(AKEY.data);
        if (r == 1) goto L910;
        if (r == 2) goto L595;
        if (r == 3) goto L805;
        for (I2 = 1; I2 <= NMBNKY; I2++) {
            if (AKEY != BNKEYS[I2]) continue;
            ILLOC(LABEL + 14 + J) = IBNKYS[I2];
            goto L535;
        }
        std::printf("\n*** INVALID KEYWORD FOR SPECIFYING NUCLEAR DEFORMATION: %.8s,  VALID KEYWORDS ARE:\n", AKEY.data);
        for (I2 = 1; I2 <= NMBNKY; I2++) std::printf("      %.8s ", BNKEYS[I2].data);
        std::printf("\n");
        return;
    L545:
        ALLOC4(LABEL + 7 + J) = (float)TEMP;
        II = II + 1;
    }
    goto L689;

    // NUCLEUS
L550:
    r = NXWORD(AKEY.data);
    if (r == 1) goto L970;
    if (r == 2) goto L700;
    if (r == 3) goto L810;
L555:
    for (J = 1; J <= LISTN1; J++) {
        if (KEY2[J] != AKEY) continue;
        INT_val = KEY2VL[J];
        if (EQULSW) { IDFLTS[10] = INT_val; goto L100; } // goto L255
        ILLOC(LABEL + 10) = INT_val;
        goto L689;
    }
    goto L970;

    // MODEL
L570:
    r = NXWORD(AKEY.data);
    if (r == 1) goto L975;
    if (r == 2) goto L700;
    if (r == 3) goto L810;
L572:
    for (J = 1; J <= NUMMOD; J++) {
        if (AKEY != MODELS_local[J]) continue;
        INT_val = J;
        if (EQULSW) { IDFLTS[14] = INT_val; goto L100; } // goto L265
        ILLOC(LABEL + 14) = INT_val;
        goto L689;
    }
    goto L975;

    // SCALE FACTOR
L580:
    { double dval; r = NXVALF(dval); VAL = dval; }
    if (r == 1) goto L910;
    if (r == 2) goto L700;
    if (r == 3) goto L810;
    ALLOC4(LABEL + 17) = (float)VAL;
    if (!SCNDSW) goto L689;
    II = II + 1;
    { double dval; r = NXVALF(dval); VAL = dval; }
    if (r == 1) goto L910;
    if (r == 2) goto L700;
    if (r == 3) goto L810;
    ALLOC4(LABEL + 18) = (float)VAL;
    goto L689;

    // TOTALLX
L590:
    r = NXHINT(INT_val);
    if (r == 1) goto L910;
    if (r == 2) goto L700;
    if (r == 3) goto L810;
    ILLOC(LABEL + 19) = INT_val;
    if (INT_val < 0) goto L935;
    goto L689;

    // MAXIMUMLX
L600:
    r = NXHINT(INT_val);
    if (r == 1) goto L910;
    if (r == 2) goto L700;
    if (r == 3) goto L810;
    ILLOC(LABEL + 20) = INT_val;
    if (INT_val < 0) goto L935;
    goto L689;

L689:
    I1 = I1 + 1;
    if (I1 <= IORDER) goto L450;
    goto L700;

    // All data for one coupling has been read — check consistency
L595:
    if (J != 1) goto L950;
L700:
    IPT = ILLOC(LABEL + 10);
    if (!SCNDSW) goto L750;
    if (IPT == NOTDEF_INT) goto L980;
    if (ILLOC(LABEL + 4) == NOTDEF_INT) goto L950;

    // Convert TP to PT input
    if (IPT != 4) goto L790;
    for (II = 1; II <= NUMFLP; II++) {
        I = IFLIPS[II];
        ITEMP = ILLOC(LABEL + I);
        ILLOC(LABEL + I) = ILLOC(LABEL + I + 1);
        ILLOC(LABEL + I + 1) = ITEMP;
    }
    ILLOC(LABEL + 10) = 3;
    goto L790;

    // First order cannot have PT
L750:
    if (IPT <= 2 || IPT == NOTDEF_INT) goto L790;
    std::printf("\n*** NUCLEUS SPECIFICATION OF PT OR TP IS INVALID FOR FIRST ODER.\n");
    return;

L790:
    if (I1 < 0) goto L840;
    goto L410;

    // Check for ) to stop the input
L800:
    r = FITKEY(AKEY.data);
    if (r == 1 || r == 2 || r == 4) goto L920;
    if (r == 3) goto L850;
    goto L920;
L805:
    if (J != 1) goto L950;
L810:
    I1 = -1;
    r = FITKEY(AKEY.data);
    if (r == 1 || r == 2 || r == 4) goto L920;
    if (r == 3) goto L700;
    goto L920;

L840:
    LABEL = LABEL + LENDAT;
L850:
    ILLOC(LDATA4) = ICH;
    IWRDS = (LABEL - LDATA4) / ALLOCS.FACFR4 + 1;
    IL = CHOPIT(IWRDS);
    IRET = 1;
    return;

    // Error messages
L900:
    std::printf(" *** DATA NOT ENCLOSED IN PARENTHESIS OR KEYWORD MISPELLED.   DATA KEYWORDS ARE:\n\n");
    for (I = 1; I <= NUMKEY; I++) std::printf(" %.8s ", KEYS[I].data);
    std::printf("\n");
    for (I = 1; I <= NMBCKY; I++) std::printf(" %.8s ", BCKEYS[I].data);
    std::printf("\n");
    for (I = 1; I <= NMBNKY; I++) std::printf(" %.8s ", BNKEYS[I].data);
    std::printf("\n");
L905:
    std::printf("\n*** THIS OR LAST KEYWORD READ IS: %.8s\n", AKEY.data);
L908:
    std::printf("\n*** INPUT LINE UP TO ERROR IS:  ");
    for (I = 1; I <= INPBUF.INCH; I++) std::printf("%c", INPBUF.IBUF[I]);
    std::printf("\n");
    return;
L910:
    std::printf("\n*** ERROR IN READING%5d DATA VALUE IN SET\n", II);
    goto L908;
L915:
    std::printf("\n*** CHANNELS TO BE COUPLED MUST BE SPECIFIED AS CH1-CH2.  YOU HAVE LEFT OUT THE \"-\".\n");
    goto L908;
L920:
    std::printf("\n*** INVALID OR MISPLACED PUNCTUATION.\n");
    goto L908;
L925:
    std::printf("\n*** INVALID ORDER OF COUPLING:%12d\n", JORDER);
    goto L908;
L930:
    std::printf("\n*** COUPLED MUST BE FIRST KEYWORD\n");
    goto L905;
L935:
    std::printf("\n*** ILLEGAL MULTIPOLARITY:%12d\n", INT_val);
    goto L908;
L940:
    std::printf("\n*** IF USING ORDER MUST BE SECOND KEYWORD\n");
    goto L905;
L945:
    std::printf("\n*** INCORRECTLY SPECIFIED DEFAULT VALUE.\n");
    goto L905;
L950:
    std::printf("\n*** FOR SECOND ORDER COUPLING MUST DEFINE TWO MULTIPOLES.\n");
    goto L908;
L960:
    std::printf("\n*** COUPLING TYPE MISPELLED OR MISPLACED.  VALID KEYWORDS ARE:\n\n");
    for (I = 1; I <= LISTN; I++) std::printf(" %.8s ", KEY1[I].data);
    std::printf("\n");
    goto L905;
L970:
    std::printf("\n*** NUCLEUS KEYWORD MISPELLED OR MISPLACED.  VALID KEYWORDS ARE:\n\n");
    for (I = 1; I <= LISTN1; I++) std::printf(" %.8s ", KEY2[I].data);
    std::printf("\n");
    goto L905;
L975:
    std::printf("\n*** MODEL SPECIFICATION MISSPELLED OR MISPLACED.  VALID MODELS ARE\n\n");
    for (I = 1; I <= NUMMOD; I++) std::printf(" %.8s ", MODELS_local[I].data);
    std::printf("\n");
    goto L905;
L980:
    std::printf("\n*** FOR SECOND ORDER COUPLING MUST USE NUCLEUS KEY\n");
    goto L905;
}

// CHNSCN — Scans an input line for channel specification
// Translated from source.f lines 8267-8523
void CHNSCN(char8* GUY, double* ESTARS, int* NODVLS, int* LVALS, int* JVALS,
            int* IPARYS, int& IRTN)
{
    // Character constants
    const char SLASH = '/';
    const char PLUS  = '+';
    const char MINUS = '-';
    const char EQUAL = '=';
    const char BLANK = ' ';
    const char LEFTP = '(';

    int I, II;
    int INCHST;
    bool PREDSW, EREDSW;
    int ITOKEN, IVALU, IMESS;
    char STPCHR;
    char8 CVALU;
    double VALU, DUMVAL;
    char dummy[4] = {};   // unused MESS parameter for MSCAN

    IRTN    = 0;
    INCHST  = INPBUF.INCH;
    PREDSW  = false;
    EREDSW  = false;

    // Initialise all 3 slots (1-based)
    char8 BLANK8;
    std::memset(BLANK8.data, ' ', 8);
    for (I = 1; I <= 3; I++) {
        GUY[I]    = BLANK8;
        NODVLS[I] = NOTDEF_INT;
        LVALS[I]  = NOTDEF_INT;
        IPARYS[I] = 0;
        JVALS[I]  = NOTDEF_INT;
        ESTARS[I] = INTRNL.UNDEF;
    }

    I = 1;

// -----------------------------------------------------------------------
//  LOOP ON REACTION PARTICIPANTS
// -----------------------------------------------------------------------
L100:
    MSCAN(-20, ITOKEN, GUY[I], DUMVAL, IVALU, dummy, IMESS, STPCHR);
    if (IMESS > 5)  goto L9004;
    if (IMESS == 0) goto L9005;
    INPBUF.INCH = INPBUF.INCH + 1;   // advance past stop char

    // If stopped by '(' -> process excited-state / J specification
    if (STPCHR != LEFTP) goto L400;

// -----------------------------------------------------------------------
//  PROCESS EXCITED STATE INFO  (E*, J, parity, level descriptor)
// -----------------------------------------------------------------------
L200:
    MSCAN(18, ITOKEN, CVALU, VALU, IVALU, dummy, IMESS, STPCHR);
    if (ITOKEN == -9) goto L9900;
    if (ITOKEN == -8) goto L9007;
    if (ITOKEN != -2) goto L230;

    // It is E*
    if (ESTARS[I] != INTRNL.UNDEF) {
        std::printf("0**** ATTEMPTING TO ENTER 2 VALUES FOR E*\n");
        goto L9900;
    }
    ESTARS[I] = VALU;
    goto L200;

L230:
    if (ITOKEN != -5) goto L270;

    // It is J (half-integer token)
    if (JVALS[I] != NOTDEF_INT) {
        std::printf("0**** ATTEMPTING TO ENTER 2 VALUES FOR J\n");
        goto L9900;
    }
    JVALS[I] = 2 * IVALU;
    if (STPCHR != SLASH) goto L250;

    // Bypass /2: advance past '/', then verify the '2'
    INPBUF.INCH = INPBUF.INCH + 1;
    JVALS[I] = IVALU;
    MSCAN(18, ITOKEN, CVALU, VALU, IVALU, dummy, IMESS, STPCHR);
    if (ITOKEN == -9) goto L9900;
    if (ITOKEN == -8) goto L9007;
    if (!(ITOKEN == -5 && IVALU == 2)) {
        std::printf("0**** INVALID J FIELD\n");
        goto L9900;
    }

L250:
    // Parity may be given by terminal + or -
    if (STPCHR == PLUS)  { IPARYS[I] = +1; INPBUF.INCH++; goto L200; }
    if (STPCHR == MINUS) { IPARYS[I] = -1; INPBUF.INCH++; goto L200; }
    goto L200;

L270:
    if (ITOKEN != -3) goto L280;

    // It is a nodes-and-L level descriptor (e.g. 2P3/2-)
    if (LVALS[I] != NOTDEF_INT) {
        std::printf("0**** ATTEMPTING TO ENTER TWO LEVEL DESCRIPTORS\n");
        goto L9900;
    }
    LVALS[I]  = IVALU;
    NODVLS[I] = (int)VALU;
    goto L200;

L280:
    // Is it the closing ')'?
    if (ITOKEN != 11) {
        std::printf("0**** UNRECOGNIZABLE EXCITED STATE DESCRIPTOR (TOKEN =%4d)\n", ITOKEN);
        goto L9900;
    }
    // Done: bypass the ')' already consumed by MSCAN; set STPCHR to char after it
    STPCHR         = INPBUF.IBUF[INPBUF.INCH];
    INPBUF.INCH    = INPBUF.INCH + 1;

// -----------------------------------------------------------------------
//  READY TO PROCESS THE NEXT FIELD (or end of channel spec)
// -----------------------------------------------------------------------
L400:
    INPBUF.INCH = INPBUF.INCH - 1;
    if (PREDSW && EREDSW) goto L800;

    // Skip blanks
L420:
    INPBUF.INCH = INPBUF.INCH + 1;
    if (STPCHR != BLANK) goto L450;
    if (INPBUF.INCH >= INPBUF.NOCH) goto L460;
    STPCHR = INPBUF.IBUF[INPBUF.INCH];
    goto L420;

    // Only '+' and '=' are valid delimiters inside the specification
L450:
    if (STPCHR == EQUAL) goto L500;
    if (STPCHR == PLUS)  goto L550;

    // End of channel spec: STPCHR will be reconsidered by CONTRL
    INPBUF.INCH = INPBUF.INCH - 1;
L460:
    if (PREDSW) goto L800;
    std::printf("0**** CHANNEL SPECIFICATION REQUIRES TWO"
                " NUCLEAR SYMBOLS SEPARATED BY A \"+\".\n");
    goto L9900;

// -----------------------------------------------------------------------
//  '=' SIGN -- bound-state result
// -----------------------------------------------------------------------
L500:
    if (EREDSW) goto L9020;
    EREDSW = true;
    // Always store result of '=' in 3rd slot
    if (I == 2) {
        I = 3;
    } else {
        GUY[3]    = GUY[1];
        JVALS[3]  = JVALS[1];
        ESTARS[3] = ESTARS[1];
        // I stays the same (1), re-scan GUY[1]
    }
    goto L100;

// -----------------------------------------------------------------------
//  '+' SIGN
// -----------------------------------------------------------------------
L550:
    if (PREDSW) goto L9022;
    PREDSW = true;
    I = I + 1;
    goto L100;

// -----------------------------------------------------------------------
//  SUCCESS
// -----------------------------------------------------------------------
L800:
    IRTN = 2;
    if (EREDSW) IRTN = 3;
    return;

// -----------------------------------------------------------------------
//  ERROR MESSAGES
// -----------------------------------------------------------------------
L9004:
    std::printf("0**** TOO MANY CHARACTERS IN SYMBOL AND ATOMIC MASS.\n");
    goto L9900;
L9005:
    std::printf("0**** UNEXPECTED DELIMITOR OR END OF LINE ENCOUNTERED"
                " WHEN LOOKING FOR ELEMENT SYMBOL.\n");
    goto L9900;
L9007:
    std::printf("0**** COMPLETE SPECIFICATION MUST BE ON SAME LINE AS"
                " CHANNEL KEYWORD.\n");
    goto L9900;
L9020:
    std::printf("0**** TWO \"=\" SIGNS ARE NOT ALLOWED\n");
    goto L9900;
L9022:
    std::printf("0**** TWO \"+\" SIGNS ARE NOT ALLOWED.\n");
    goto L9900;

L9900:
    {
        int INCHM1 = INPBUF.INCH - 1;
        std::printf("      WAS PROCESSING FIELD%2d OF THE SPECIFICATION\n"
                    "      SPECIFICATION UP TO THE ERROR IS: ", I);
        for (II = INCHST; II <= INCHM1 && II <= IBUFSIZE; II++)
            std::printf("%c", INPBUF.IBUF[II]);
        std::printf("\n");
    }
    IRTN = -1;
    return;
}

// CLRCHN -> clrchn_translated.cpp

// ============================================================================
// CTRAN — Transition parameters input parser
// Translated from source.f lines 12805-12992
// ============================================================================
void CTRAN(int& IRET)
{
    extern int GIVEAL(int IWRDS);
    extern int CHOPIT(int IWRDS);

    int NUMKEY = 8;
    char8 KEYS[9]; // 1-based
    KEYS[1] = char8("TARGET");
    KEYS[2] = char8("PROJECTI");
    KEYS[3] = char8("MULTIPOL");
    KEYS[4] = char8("BELX");
    KEYS[5] = char8("BEUP");
    KEYS[6] = char8("BE");
    KEYS[7] = char8("BETAC");
    KEYS[8] = char8("BETAN");

    int LENDAT = 3;
    int LAL = 9;
    char8 NAME8("TRANSITI");
    char8 AKEY;

    int ICG, ICH1, ICH, IBN, IBC, IMT;
    int IORDER, IPT, I, INT_val, IWRDS, IL, LABEL, LTRANS;
    double TEMP;
    int r;

    INTRNL.ICD3I = 1;
    IRET = 0;
    ICG = 1;
    ICH1 = 0;

L10:
    ICH = 0;
    IBN = 0;
    IBC = 0;
    IMT = 0;
    if (ICG == 2) goto L90;
    I = NAMLOC(NAME8.data);
    if (I != 0) LOCPTRS.Z[I] = -LOCPTRS.Z[I];
    CCBLK.ITRANS = GIVEAL(IWRDS);
    { int LDATA = LOCPTRS.Z[CCBLK.ITRANS];
    NAMCOM.NAMES[CCBLK.ITRANS] = NAME8;
    LTRANS = ALLOCS.FACFR4 * (LDATA - 1) + 1; }
    LABEL = LTRANS + LAL;
    ILLOC(LTRANS + 3) = 0;
    ILLOC(LTRANS + 4) = 0;
    ILLOC(LTRANS + 5) = 0;
    ILLOC(LTRANS + 6) = 0;
    ILLOC(LTRANS + 7) = NOTDEF_INT;
    ILLOC(LTRANS + 8) = NOTDEF_INT;

L20:
    r = FITKEY(AKEY.data);
    if (r == 1) goto L90;
    if (r == 2 || r == 3 || r == 4) goto L900;
    goto L900;

    // ( ENCOUNTERED LETS GO
L90:
    IORDER = 0;
    IPT = 0;
L100_ct:
    r = NXWORD(AKEY.data);
    if (r == 1) goto L699;
    if (r == 2 || r == 3) goto L900;
    for (I = 1; I <= NUMKEY; I++) {
        if (KEYS[I] != AKEY) continue;
        if (I >= 3 && IPT == 0) goto L930;
        switch (I) {
            case 1: case 2: goto L200_ct;
            case 3: goto L300_ct;
            case 4: case 5: case 6: goto L350_ct;
            case 7: goto L400_ct;
            case 8: goto L450_ct;
        }
    }
    goto L900;

L200_ct:
    ILLOC(LTRANS + 2 + I) = LABEL - LTRANS;
    if (IORDER != 0) goto L930;
    IPT = I;
    goto L100_ct;

L300_ct:
    IORDER = IORDER + 1;
    IMT = IORDER;
    goto L100_ct;

    // BE
L350_ct:
    IORDER = IORDER + 1;
    IBC = IORDER;
    ILLOC(LTRANS + 6 + IPT) = -1;
    goto L100_ct;

    // BETAC
L400_ct:
    IORDER = IORDER + 1;
    IBC = IORDER;
    ILLOC(LTRANS + 6 + IPT) = 0;
    goto L100_ct;

    // BETAN
L450_ct:
    IORDER = IORDER + 1;
    IBN = IORDER;
    goto L100_ct;

    // End of keywords, go on to data
L699:
    LABEL = LABEL - LENDAT;

    // Read the next set of data
L700_ct:
    LABEL = LABEL + LENDAT;
    ILLOC(LABEL) = NOTDEF_INT;
    ALLOC4(LABEL + 1) = undef4();
    ALLOC4(LABEL + 2) = undef4();
    for (I = 1; I <= IORDER; I++) {
        if (I == IMT) goto L720;
        if (I == IBC) goto L730;
        r = NXVALF(TEMP);
        if (r == 1) goto L910;
        if (r == 2) goto L800_ct;
        if (r == 3) goto L910;
        if (r == 4) goto L890;
        ALLOC4(LABEL + 2) = (float)TEMP;
        goto L750;
    L720:
        r = NXHINT(INT_val);
        if (r == 1) goto L910;
        if (r == 2) goto L800_ct;
        if (r == 3) goto L885;
        ILLOC(LABEL) = INT_val;
        goto L750;
    L730:
        r = NXVALF(TEMP);
        if (r == 1) goto L910;
        if (r == 2) goto L800_ct;
        if (r == 3) goto L910;
        if (r == 4) goto L890;
        ALLOC4(LABEL + 1) = (float)TEMP;
    L750:;
    }
    ICH = ICH + 1;
    goto L700_ct;

    // Semi-colon found, go to next nucleus
L800_ct:
    ILLOC(LTRANS + 4 + IPT) = ICH;
    ICH1 = ICH;
    ICG = ICG + 1;
    goto L10;

L885:
    r = FITKEY(AKEY.data);
    if (r == 1 || r == 2 || r == 4) goto L900;
    if (r == 3) goto L890;
    goto L900;

L890:
    ILLOC(LTRANS + 4 + IPT) = ICH;
    ILLOC(LTRANS) = ICH1 + ICH;
    ILLOC(LTRANS + 1) = LENDAT;
    ILLOC(LTRANS + 2) = LAL;
    IWRDS = (LABEL - LTRANS) / ALLOCS.FACFR4 + 1;
    IL = CHOPIT(IWRDS);
    IRET = 1;
    return;

L900:
    std::printf("\n**** DATA NOT ENCLOSED IN PARENTHESIS OR KEYWORD MISSPELLED, DATA KEYWORDS ARE:\n\n");
    for (I = 1; I <= NUMKEY; I++) std::printf(" %.8s ", KEYS[I].data);
    std::printf("\n");
    return;
L910:
    std::printf("\n**** ERROR IN READING THE%5d NUMERICAL DATA VALUE IN A SET\n", I);
    return;
L930:
    std::printf("\n**** MUST USE PROJECTILE OR TARGET KEYWORD FIRST\n");
    return;
}

// SETCHN -> setchn_translated.cpp
