// xsectn_translated.cpp — XSECTN
// Translated from source.f lines 38258-39158
// Computes transfer/inelastic cross sections from S-matrix elements

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include "ptolemy_alloc.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_forward.h"
#include "ptolemy_io.h"
#include <cstdio>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <string>

extern double second();

// Fortran Gw.d format emulation.
// F-mode (0.1 <= |x| < 10^d): F(w-4).(d-n) right-justified + 4 trailing spaces.
//   n = floor(log10(|x|))+1 for |x|>=1, 0 for 0.1<=|x|<1, 1 for x=0 (gfortran behavior).
// E-mode otherwise: [sign]0.dddddE±ee right-justified in field w.
static std::string fG(double x, int w, int d) {
    double absval = std::fabs(x);
    double limit = std::pow(10.0, (double)d);
    bool f_mode;
    int n;
    if (x == 0.0) {
        f_mode = true; n = 1;
    } else if (absval >= 0.1 && absval < limit) {
        f_mode = true;
        n = (int)std::floor(std::log10(absval)) + 1;
        if (n < 0) n = 0;
    } else {
        f_mode = false; n = 0;
    }
    char buf[128];
    if (f_mode) {
        int fw = w - 4, fd = d - n;
        if (fd < 0) fd = 0;
        if (fd == 0) {
            // Fortran Fw.0 prints a trailing decimal point; C's %.0f does not
            char num[64];
            std::snprintf(num, sizeof(num), "%.0f.", x);
            int numlen = (int)std::strlen(num);
            int pad = fw - numlen; if (pad < 0) pad = 0;
            return std::string(pad, ' ') + num + "    ";
        }
        std::snprintf(buf, sizeof(buf), "%*.*f", fw, fd, x);
        return std::string(buf) + "    ";
    } else {
        // E-mode: Fortran notation 0.ddddd x 10^e_f
        int e_f = (absval > 0.0) ? (int)std::floor(std::log10(absval)) + 1 : 0;
        double m_f = (absval > 0.0) ? absval / std::pow(10.0, (double)e_f) : 0.0;
        if (m_f >= 1.0)      { e_f++; m_f /= 10.0; }
        if (m_f > 0.0 && m_f < 0.1) { e_f--; m_f *= 10.0; }
        char mant[64];
        std::snprintf(mant, sizeof(mant), "%.*f", d, m_f);
        if (mant[0] == '1') { e_f++; m_f /= 10.0; std::snprintf(mant, sizeof(mant), "%.*f", d, m_f); }
        char sign_c = (x < 0.0) ? '-' : ' ';
        char exp_s  = (e_f >= 0)  ? '+' : '-';
        char content[64];
        std::snprintf(content, sizeof(content), "%c%sE%c%02d",
                      sign_c, mant, exp_s, std::abs(e_f));
        int clen = (int)std::strlen(content);
        int pad = w - clen; if (pad < 0) pad = 0;
        return std::string(pad, ' ') + content;
    }
}

void XSECTN(int& IRTN)
{
    // COMMON block references
    auto* Z      = LOCPTRS.Z;
    auto& FACFR4 = ALLOCS.FACFR4;
    auto& PI     = CNSTNT.PI;
    auto& RADIAN = CNSTNT.RADIAN;
    auto& DEGREE = CNSTNT.DEGREE;
    auto& SMLNUM = CNSTNT.SMLNUM;
    auto& LMIN   = INTGER.LMIN;
    auto& LMAX   = INTGER.LMAX;
    auto& IPRINT = INTGER.IPRINT;
    auto& LEBACK = INTGER.LEBACK;
    auto& ANGMIN = FLOAT_common.ANGMIN;
    auto& ANGMAX = FLOAT_common.ANGMAX;
    auto& ANGSTP = FLOAT_common.ANGSTP;
    auto& ELAB   = FLOAT_common.ELAB;
    auto& ECM    = FLOAT_common.ECM;
    auto& AMA    = FLOAT_common.AMA;
    auto& AMB    = FLOAT_common.AMB;
    auto& AMBIGA = FLOAT_common.AMBIGA;
    auto& AMBIGB = FLOAT_common.AMBIGB;
    auto& Q      = FLOAT_common.Q;
    auto& TAU_f  = FLOAT_common.TAU;
    auto& PROBLM = SWITCH.PROBLM;
    auto& LABANG = SWITCH.LABANG;
    auto& ISAVMX = SWITCH.ISAVMX;
    auto& ISAVB  = SWITCH.ISAVB;
    auto& IEXTYP = SWITCH.IEXTYP;
    auto& IBRNSB = SWITCH.IBRNSB;
    auto& ISTRIP = INTRNL.ISTRIP;
    auto& LIMOST = INTRNL.LIMOST;
    auto& AKI    = KANDM.AKI;
    auto& AKO    = KANDM.AKO;
    auto& LOMOST = KANDM.LOMOST;
    auto& LXMIN  = INELCM.LXMIN;
    auto& LXMAX  = INELCM.LXMAX;
    auto& LXSTEP = INELCM.LXSTEP;
    auto& NSPLI  = INELCM.NSPLI;
    auto& ITOCS  = INELCM.ITOCS;
    auto& IRDINT = INELCM.IRDINT;
    auto& IDENSW = INELCM.IDENSW;
    auto& FORMF_ref = FORMF;

    // CHARACTER*8 data
    char8 BETNAM("BETAS   ");
    char8 RUTNAM[3]  = {char8(), char8("TORUTHIN"), char8("TORUTHOT")};
    char8 FENAME[3]  = {char8(), char8("FIN     "), char8("FOUT    ")};
    char8 RUTHNM[3]  = {char8(), char8("RUTHIN  "), char8("RUTHOUT ")};
    char8 SIGNAM[3]  = {char8(), char8("SIGMAIN "), char8("SIGMAOUT")};
    char8 INHEAD[3]  = {char8(), char8("INCOMING"), char8(" COUPLED")};
    char8 OUTHED[3]  = {char8(), char8("OUTGOING"), char8(" OPTICAL")};
    char8 LABWRD("LAB.    ");
    char8 CMWRD("C.M.    ");
    char8 BLANK("        ");
    char8 RUTWRD("/RUTH   ");
    char SIGNS[3] = {0, '+', '-'};

    // Local variables
    double TSTART, T1, T2, T3, T4, T5, T6, T7, T8, T9;
    double TBET, TCROSS, CRSFLP;
    bool TRANSW, CCSW, PBUGSW, PCONSW, PHEDSW, FKEPSW;
    int IPRNT, JMOST, NUMANG, LINECT, LINEMX = 57;
    int NCHNDF, MCHN, LTOCOF, LSOFF, LFOFF;
    double RUTHFC, TAU_v;
    auto& ABAR = KANDM.ABAR;  // Fortran: ABAR is in KANDM common block
    int IMAG, IPHASE, ISIGIN, ISIGOT;
    int ITORT1, ITORT2, IFIN, IFOUT, IRUTH1, IRUTH2;
    double REACIN, REACOT, SIGTOT;
    int ICROSS, IFEPSL, IF_var;
    int IDNPAR;
    char8 LABCM;
    int LXMN, LXMX, LDELMX, JPMX, LXRNGE, NMLX, NTEMPS;
    int ITEMPS, ILXTOT, IMXTOT;
    int LTOCS, LMAG, LPHASE, LRDINT_l, LTOCS2;
    int IPLM, ITMPCN, ILXCRO, IMXCRO, IMXCON, ILXCNT, IFTEMP, IFLO, IFHI, IFERR;
    int ISPACE;
    int LLXCNT, LPLM, LTMPCN, LLXCRO, LMXCRO, LFTEMP, LFEND;
    int LFLO, LFHI, LFERR, LF, LMXCON;
    int LTORT1_l, LTORT2_l, LSIGIN, LRUTH1, LCROSS, LLXTOT, LMXTOT;
    double ANGBLK;
    int NUMBLK, LINADD, NMLXH;

    TSTART = (float)second();
    IRTN = 0;
    CRSFLP = 0;
    TRANSW = (ISTRIP != 0);
    CCSW = (PROBLM == 24);
    JMOST = LIMOST;
    if (CCSW) JMOST = LOMOST - LXMAX;

    IPRNT = IPRINT % 10;
    PBUGSW = (IPRNT >= 4);
    PCONSW = (IPRNT >= 3);

    // Print phase shifts (allocates SMAG/SPHASE as side-effect for XSECTN)
    PHSPRT(IPRNT >= 1);
    T1 = (float)second();
    TBET = 0;
    TCROSS = 0;

    // Print header
    std::printf("1%54sP T O L E M Y\n"
                "%10sCOMPUTATION OF CROSS SECTIONS"
                "%38sSOPHISTICATED PROGRAMS CONTAIN SUBTLE ERRORS\n"
                "0%.45s%7.2f MEV     %.65s\n\n",
                "", "", "",
                &HEDCOM.REACT[1], ELAB, &HEDCOM.HEADER[1]);
    LINECT = 4;

    PHEDSW = true;

    LABCM = CMWRD;
    if (LABANG != 0) LABCM = LABWRD;

    // Locate magnitudes and phases (allocated by PHSPRT)
    IMAG = NAMLOC("SMAG    ");
    IPHASE = NAMLOC("SPHASE  ");

    if (CCSW) Z[IRDINT] = -Z[IRDINT];

    // Lab angle range limits
    ABAR = 1.0e+5;
    if (LABANG == 0) goto L100;

    {
        double X = AMBIGA / AMA;
        double QMIN_l = Q;
        if (CCSW) QMIN_l = KANDM.QMIN;
        else QMIN_l = Q;
        double TEMP = sqrt(X * AMBIGB * (AMA + AMBIGA) / (AMB * (AMB + AMBIGB))
                         * (ECM + QMIN_l) / ECM);
        if (PBUGSW) std::printf("\nLIMITS ON SIN(THETALAB):%15.5G%15.5G\n", X, TEMP);
        X = std::min(X, TEMP);
        if (X > 1) goto L100;
        double ANGBAR = RADIAN * asin(X);
        ABAR = ANGSTP * ((int)((ANGBAR - ANGMIN) / ANGSTP + 1.0e-10) + 0.5) + ANGMIN;
        double TABAR = 2 * ABAR;
        ANGMAX = std::min(ANGMAX, TABAR);
    }

L100:
    NUMANG = (int)((ANGMAX - ANGMIN) / ANGSTP + 1.5);

    // Compute elastic cross sections
    FKEPSW = (ISAVB == 1) || WAVCOM.SOSWS[1];
    ELDCS(AKI, KANDM.ETAS[1], ANGMIN, ANGMAX, ANGSTP,
          LMIN, JMOST, WAVCOM.LSKIPS[1], WAVCOM.ISTATS[1], WAVCOM.JSPS[1],
          KANDM.ISMATS[1], WAVCOM.ITOCE[1], WAVCOM.NUMJS[1], KANDM.ISIGS[1],
          ELAB, KANDM.RATMSS[1], 0 /*false*/,
          ITORT1, RUTNAM[1].data, FKEPSW ? 1 : 0, IFIN, FENAME[1].data,
          REACIN, IRUTH1, RUTHNM[1].data);
    ISIGIN = NAMLOC("CROSSSEC");
    NAMCOM.NAMES[ISIGIN] = SIGNAM[1];

    // 2nd elastic channel
    int IPROB = 2;
    if (PROBLM >= 23) IPROB = 1;
    FKEPSW = (ISAVB == 1) || WAVCOM.SOSWS[IPROB];
    double AKS_p = (IPROB == 1) ? AKI : AKO;
    ELDCS(AKS_p, KANDM.ETAS[IPROB], ANGMIN, ANGMAX, ANGSTP,
          LMIN, JMOST, WAVCOM.LSKIPS[IPROB], WAVCOM.ISTATS[IPROB],
          WAVCOM.JSPS[IPROB], KANDM.ISMATS[2], WAVCOM.ITOCE[IPROB],
          WAVCOM.NUMJS[IPROB], KANDM.ISIGS[IPROB], ELAB, KANDM.RATMSS[IPROB],
          0 /*false*/,
          ITORT2, RUTNAM[2].data, FKEPSW ? 1 : 0, IFOUT, FENAME[2].data,
          REACOT, IRUTH2, RUTHNM[2].data);
    ISIGOT = NAMLOC("CROSSSEC");
    NAMCOM.NAMES[ISIGOT] = SIGNAM[2];

    T2 = (float)second();

    // Setup for calculations
    NCHNDF = 2;
    LTOCOF = 0;
    LSOFF = 0;
    MCHN = 2;
    IF_var = 0;
    LFOFF = 0;

    // Rutherford conversion factor
    RUTHFC = 1;
    if (KANDM.ETAS[1] != 0)
        RUTHFC = KANDM.ETAS[2] * AKI / (KANDM.ETAS[1] * AKO);

    if (CCSW) {
        // CC setup (stub for now)
        goto L200;
    }

L200:
    ICROSS = NALLOC("CROSSSEC", NUMANG * (NCHNDF - 1));
    IFEPSL = NALLOC("FEPSILON", 2 * LEBACK + 1);

    IDNPAR = 0;
    if (IDENSW) IDNPAR = 1;

    // Main channel loop
    for (int NCHN = 2; NCHN <= NCHNDF; NCHN++) {
        if (CCSW) {
            // CC channel setup (stub)
        }

        // Scan TOCS for extrema
        LXMN = 1000;
        LXMX = -1000;
        LDELMX = -1000;
        JPMX = -1000;
        LTOCS = Z[ITOCS] * FACFR4 - FACFR4 + LTOCOF;
        for (int KOFFS = 1; KOFFS <= NSPLI; KOFFS++) {
            if (ILLOC(LTOCS + 4 * KOFFS) < 0) continue;
            int LX = ILLOC(LTOCS + 4 * KOFFS - 2);
            LXMN = std::min(LXMN, LX);
            LXMX = std::max(LXMX, LX);
            int LDEL = ILLOC(LTOCS + 4 * KOFFS - 3);
            LDELMX = std::max(LDELMX, LDEL);
            JPMX = std::max(JPMX, ILLOC(LTOCS + 4 * KOFFS - 1));
        }
        LXRNGE = LXMX - LXMN + 1;
        NMLX = (LXMX - LXMN + LXSTEP) / LXSTEP;
        if (PBUGSW) std::printf(" LXMN, LXMX, LDELMX, JPMX, NMLX:%8d%8d%8d%8d%8d\n",
                                LXMN, LXMX, LDELMX, JPMX, NMLX);

        // Allocate scratch
        NTEMPS = (LXMX + 1) * LXRNGE * (LDELMX + 1);
        ITEMPS = NALLOC("TEMPS   ", NTEMPS);
        ILXTOT = NALLOC("LXTOTS  ", LXMX + 1);
        IMXTOT = NALLOC("MXTOTS  ", NSPLI);

        // Compute BETA's
        T3 = (float)second();
        SETLOG(2 * (LOMOST + LXMAX));

        LMAG = FACFR4 * (Z[IMAG] - 1) + 1 + LSOFF;
        LPHASE = FACFR4 * (Z[IPHASE] - 1) + 1 + LSOFF;
        LRDINT_l = FACFR4 * (Z[IRDINT] - 1) + 1;
        LTOCS = Z[ITOCS] * FACFR4 - FACFR4 + LTOCOF;

        if (CCSW) DSGMAL(KANDM.ETAS[2], LIMOST + LXMAX,
                         &ALLOC(Z[KANDM.ISIGS[2]]));

        {
            int ntot_sp = NSPLI * (LIMOST + 1);
            for (int qq = 0; qq < ntot_sp; qq++) {
                if (std::isnan((double)ALLOC4(LMAG + qq))) ALLOC4(LMAG + qq) = 0.0f;
                if (std::isnan((double)ALLOC4(LPHASE + qq))) ALLOC4(LPHASE + qq) = 0.0f;
            }
        }
        double DUMMY_d = 0;
        double DUMMY8[6] = {};  // ALOWFC(2,2)=4 elems, FCOUL(2,3)=6 elems
        BETCAL(0 /*false*/, AKI, WAVCOM.JSPS[1], NSPLI, LMIN, LIMOST,
               INELCM.LSKIP,
               LXMN, LXMX, LDELMX, NTEMPS, WAVCOM.ISTATS[1], PBUGSW ? 1 : 0,
               SIGTOT, DUMMY_d,
               &ILLOC(LTOCS + 1), DUMMY8, &ALLOC4(LMAG), &ALLOC4(LPHASE),
               &ALLOC(Z[KANDM.ISIGS[1]]), &ALLOC(Z[KANDM.ISIGS[2]]),
               &ALLOC4(LRDINT_l), &ALLOC(Z[ILXTOT]), &ALLOC(Z[IMXTOT]),
               DUMMY8, &ALLOC(Z[ITEMPS]));

        Z[ITEMPS] = -Z[ITEMPS];

        T4 = (float)second();
        TBET = TBET + T4 - T3;

        // Compute TAU for lab/CM conversion
        TAU_v = sqrt(AMA * AMB * (AMB + AMBIGB) * ECM /
                    (AMBIGA * AMBIGB * (AMA + AMBIGA) * (ECM + Q)));

        // Allocate space for PLM
        ISPACE = LIMOST + 1 + ((2 * LIMOST + 1 - LXMX) * LXMX) / 2;
        IPLM = NALLOC("PLM     ", ISPACE);
        ITMPCN = NALLOC("TMPCONTR", 2 * (LIMOST + 1));

        // Space for cross sections
        ILXCRO = NALLOC("LXCROSSS", NUMANG * LXRNGE);
        IMXCRO = 0;
        if (ISAVMX == 1) IMXCRO = NALLOC("MXCROSSS", NUMANG * NSPLI);
        IMXCON = NALLOC("MXCONS  ", NSPLI);
        if (!CCSW) {
            FKEPSW = (ISAVB == 1) || (std::min({(int)JBLOCK.JS[1], JPMX, LXMX}) > 0);
            if (FKEPSW) IF_var = NALLOC("F       ", 2 * NUMANG * NSPLI);
        }
        ILXCNT = NALLOC("LXCONTRB", LXMX + 1);
        IFTEMP = NALLOC("FTEMP   ", 2 * NSPLI);
        IFLO = NALLOC("FLO     ", 2 * NSPLI);
        IFHI = NALLOC("FHI     ", 2 * NSPLI);
        IFERR = NALLOC("FERROR  ", NSPLI);
        INIT8(IFLO, 0.0);

        SETLOG(2 * (LOMOST + LXMAX));

        // Get locations
        LLXCNT = Z[ILXCNT];
        LRDINT_l = FACFR4 * (Z[IRDINT] - 1) + 1;
        LPLM = Z[IPLM];
        LTMPCN = Z[ITMPCN];
        LLXCRO = Z[ILXCRO] - LXMN * NUMANG - 1;
        LMXCRO = (IMXCRO != 0) ? Z[IMXCRO] - NUMANG - 1 : 0;
        LFTEMP = Z[IFTEMP] - 2;
        LFEND = LFTEMP + 2 * NSPLI - 1;
        LFLO = Z[IFLO] - 2;
        LFHI = Z[IFHI] - 2;
        LFERR = Z[IFERR] - 1;
        LF = (IF_var != 0) ? Z[IF_var] - 2 * NUMANG - 2 + LFOFF : 0;
        LMXCON = Z[IMXCON] - 1;
        LTORT1_l = Z[ITORT1] - 1;
        LTORT2_l = Z[ITORT2] - 1;
        if (INTGER.IZP * INTGER.IZT == 0) LTORT2_l = Z[ISIGOT] - 1;
        LSIGIN = Z[ISIGIN] - 1;
        LRUTH1 = Z[IRUTH1] - 1;
        LCROSS = Z[ICROSS] - 1 + NUMANG * (NCHN - 2);
        LLXTOT = Z[ILXTOT];
        LMXTOT = Z[IMXTOT];
        LTOCS = Z[ITOCS] * FACFR4 - FACFR4 + LTOCOF;

        // Setup print spacing
        ANGBLK = 300;
        if (NUMANG > 40) {
            ANGBLK = 10;
            if (ANGSTP < 0.99) {
                ANGBLK = 5;
                if (ANGSTP <= 0.201) ANGBLK = 1;
            }
        }

        // Compute I_v and NMLXH before header
        int I_v = 3;
        if (NCHN > 2) I_v = 5;
        NMLXH = std::min(NMLX, I_v);
        // I2 = NMLXH+1, except I2=1 when NMLX==1
        int I2_v = (NMLX == 1) ? 1 : NMLXH + 1;
        int N_v = CCSW ? 2 : 1;
        // TEMPwd for FORMAT 536
        const char* TEMPwd = (INTGER.IZP * INTGER.IZT == 0) ? "     " : "/RUTH";

        // LINADD: number of subheader lines per page (depends on NMLX and I_v)
        LINADD = (NMLX + I_v - 1) / I_v;

        // Loop over angles
        NUMBLK = (int)(ANGBLK / ANGSTP + 0.5);
        NUMBLK = NUMBLK * LINADD;

        for (int NANG = 1; NANG <= NUMANG; NANG++) {
            // Page break check (Fortran: IF (LINECT .LT. LINEMX) GO TO 520)
            if (LINECT >= LINEMX) {
                // FORMAT 3: new page header
                std::printf("1%54sP T O L E M Y\n"
                            "%10sCOMPUTATION OF CROSS SECTIONS"
                            "%38sSOPHISTICATED PROGRAMS CONTAIN SUBTLE ERRORS\n"
                            "0%.45s%7.2f MEV     %.65s\n\n",
                            "", "", "",
                            &HEDCOM.REACT[1], ELAB, &HEDCOM.HEADER[1]);
                LINECT = 4;
                PHEDSW = true;
                if (CCSW) {
                    // FORMAT 203 for CC channel (not yet implemented)
                    LINECT += 2;
                }
            }
            // Subheader (Fortran label 520: IF (.NOT. PHEDSW) GO TO 550)
            if (PHEDSW) {
                // FORMAT 523: first column header line
                std::printf("0  %.4s  REACTION     REACTION   LOW L  HIGH L   %% FROM\n", LABCM.data);
                // FORMAT 526 (NCHN==2) / FORMAT 527 (NCHN>2)
                if (NCHN == 2) {
                    std::printf("+%58s%-8.8s ELASTIC%6s%-8.8s ",
                                "", INHEAD[N_v].data, "", OUTHED[N_v].data);
                    for (int i = 0; i < NMLXH; i++) {
                        if (i < NMLXH - 1)
                            std::printf("   REACTION   ");
                        else
                            std::printf("   REACTION");
                    }
                    std::printf("\n");
                } else if (NCHN > 2) {
                    std::printf("+ %53s", "");
                    for (int i = 0; i < NMLXH; i++) {
                        if (i < NMLXH - 1)
                            std::printf("   REACTION   ");
                        else
                            std::printf("   REACTION");
                    }
                    std::printf("\n");
                }
                // FORMAT 533: second column header line
                std::printf("  ANGLE  %.4s MB        /RUTH     %%/L   %% ERROR  L>LMAX\n", LABCM.data);
                if (NCHN <= 2) {
                    // FORMAT 536
                    std::printf("+%58s%.4s MB      /RUTH      %.5s\n",
                                "", LABCM.data, TEMPwd);
                    // FORMAT 538 (NMLX>1)
                    if (NMLX > 1) {
                        bool fl = true;
                        int cnt = 0;
                        for (int LX = LXMN; LX <= LXMX; LX += LXSTEP) {
                            if (cnt % 3 == 0) {
                                if (fl) { std::printf("+%89s", ""); fl = false; }
                                else      std::printf("\n %89s", "");
                            }
                            std::printf("    LX =%2d", LX);
                            if (cnt % 3 < 2 && (LX + LXSTEP) <= LXMX)
                                std::printf("    ");
                            cnt++;
                        }
                        std::printf("\n");
                    }
                } else {
                    // FORMAT 543 (NCHN>2, NMLX>1)
                    if (NMLX > 1) {
                        bool fl = true;
                        int cnt = 0;
                        for (int LX = LXMN; LX <= LXMX; LX += LXSTEP) {
                            if (cnt % 5 == 0) {
                                if (fl) { std::printf("+%55s", ""); fl = false; }
                                else      std::printf("\n %55s", "");
                            }
                            std::printf("    LX =%2d", LX);
                            if (cnt % 5 < 4 && (LX + LXSTEP) <= LXMX)
                                std::printf("    ");
                            cnt++;
                        }
                        std::printf("\n");
                    }
                }
                // FORMAT 548: blank line
                std::printf(" \n");
                LINECT += LINADD + 1;
                PHEDSW = false;
            }
            // label 550: process angle
            double ANGLE = ANGMIN + (NANG - 1) * ANGSTP;
            double ANGCM = ANGLE;
            double AJACOB = 1.0;

            if (LABANG != 0) {
                // Lab to CM conversion
                double ANG = ANGLE;
                if (ANGLE > ABAR) ANGLE = 2 * ABAR - ANGLE;
                double AN = DEGREE * ANGLE;
                double TEMP = 1 - (TAU_v * sin(AN)) * (TAU_v * sin(AN));
                if (TEMP < 0) TEMP = 0;
                TEMP = cos(AN) * sqrt(TEMP);
                if (ANG > ABAR) TEMP = -TEMP;
                ANGCM = acos(TEMP - TAU_v * sin(AN) * sin(AN));
                TEMP = TAU_v * (TAU_v + 2 * cos(ANGCM)) + 1;
                if (TEMP < 0) TEMP = 0;
                AJACOB = TEMP * sqrt(TEMP) / fabs(TAU_v * cos(ANGCM) + 1);
                ANGCM = RADIAN * ANGCM;
            }

            // Calculate transition amplitudes
            int LHI = LMAX;
            int NFLPAM;
            AMPCAL(ANGCM, 0 /*false*/, NSPLI, LMIN, LIMOST, INELCM.LSKIP,
                   LXMX, IDNPAR, 0.0, 0.0, 0.0, 1 /*true*/, LHI,
                   &ILLOC(LTOCS + 1), &ALLOC4(LRDINT_l), DUMMY8,
                   &ALLOC(LFTEMP + 2), &ALLOC(LFLO + 2), &ALLOC(LFHI + 2),
                   &ALLOC(LFERR + 1), DUMMY8,
                   &ALLOC(LPLM), &ALLOC(LTMPCN), &ALLOC(LTMPCN + LIMOST + 1),
                   LEBACK, &ALLOC(Z[IFEPSL]), NFLPAM);
            CRSFLP += NFLPAM;

            double SIGMA = 0;
            double SUMLOW = 0, SUMHI = 0, SUMERR = 0;

            for (int LX = 0; LX <= LXMX; LX++)
                ALLOC(LLXCNT + LX) = 0;

            // Sum over KOFFS = JT, JP, LX, MX
            for (int KOFFS = 1; KOFFS <= NSPLI; KOFFS++) {
                if (ILLOC(LTOCS + 4 * KOFFS) < 0) continue;
                int LX = ILLOC(LTOCS + 4 * KOFFS - 2);
                int MX = (ILLOC(LTOCS + 4 * KOFFS - 3) + LX + 1) / 2;

                double FMNEG = 1.0;
                if (MX != 0) FMNEG = 2.0;

                double FR = ALLOC(LFTEMP + 2 * KOFFS);
                double FI = ALLOC(LFTEMP + 2 * KOFFS + 1);

                // Save F if requested
                if (FKEPSW && IF_var != 0) {
                    ALLOC(LF + 2 * NANG + 2 * NUMANG * KOFFS) = FR;
                    ALLOC(LF + 2 * NANG + 2 * NUMANG * KOFFS + 1) = FI;
                }

                double CONTRI = (10 * AJACOB) * (FR * FR + FI * FI);
                ALLOC(LMXCON + KOFFS) = CONTRI;
                if (ISAVMX == 1)
                    ALLOC(LMXCRO + NANG + NUMANG * KOFFS) = CONTRI;
                CONTRI = CONTRI * FMNEG;
                ALLOC(LLXCNT + LX) = CONTRI + ALLOC(LLXCNT + LX);
                SIGMA += CONTRI;

                SUMLOW += (10 * FMNEG * AJACOB) *
                    (ALLOC(LFLO + 2 * KOFFS) * ALLOC(LFLO + 2 * KOFFS)
                   + ALLOC(LFLO + 2 * KOFFS + 1) * ALLOC(LFLO + 2 * KOFFS + 1));
                SUMHI += (10 * FMNEG * AJACOB) *
                    (ALLOC(LFHI + 2 * KOFFS) * ALLOC(LFHI + 2 * KOFFS)
                   + ALLOC(LFHI + 2 * KOFFS + 1) * ALLOC(LFHI + 2 * KOFFS + 1));
                SUMERR += ALLOC(LFERR + KOFFS) * CONTRI;

                ALLOC(LLXCRO + NANG + NUMANG * LX) = ALLOC(LLXCNT + LX);
            }

            ALLOC(LCROSS + NANG) = SIGMA;

            // Compute errors
            double SAFE = std::max(SIGMA, SMLNUM);
            double ERRLOW = 100 * (SIGMA - SUMLOW) / (SAFE + SUMLOW);
            double HICONT = 100 * (SIGMA - SUMHI) / SAFE;
            SUMERR = 200 * SUMERR / SAFE;

            // Print cross sections
            double RUTH = ALLOC(LRUTH1 + NANG);
            double ELAS = ALLOC(LSIGIN + NANG);
            double SIGMAR = 0;
            if (RUTH != 0) SIGMAR = SIGMA / (RUTHFC * AJACOB * RUTH);

            if (NCHN <= 2) {
                // FORMAT 853: F7.2, G13.5, F10.6, 2F8.2, F8.1, G14.5, F10.6, F12.6, (T91,3G14.5)
                std::printf("%7.2f%s%10.6f%8.2f%8.2f%8.1f%s%10.6f%12.6f",
                            ANGLE, fG(SIGMA, 13, 5).c_str(),
                            SIGMAR, ERRLOW, SUMERR, HICONT,
                            fG(ELAS, 14, 5).c_str(),
                            ALLOC(LTORT1_l + NANG), ALLOC(LTORT2_l + NANG));
                if (NMLX > 1) {
                    // (T91, 3G14.5): first group at col 91 (already there), overflow → new line
                    int lx_count = 0;
                    for (int LX = LXMN; LX <= LXMX; LX += LXSTEP) {
                        if (lx_count > 0 && lx_count % 3 == 0)
                            std::printf("\n %89s", "");  // ' ' CC + 89 spaces → T91
                        std::printf("%s", fG(ALLOC(LLXCNT + LX), 14, 5).c_str());
                        lx_count++;
                    }
                }
                std::printf("\n");
            } else {
                // FORMAT 863: F7.2, G13.5, F10.6, 2F8.2, F8.1, (T56, 5G14.5)
                std::printf("%7.2f%s%10.6f%8.2f%8.2f%8.1f",
                            ANGLE, fG(SIGMA, 13, 5).c_str(),
                            SIGMAR, ERRLOW, SUMERR, HICONT);
                if (NMLX > 1) {
                    std::printf(" ");  // T56: 1 space from col 55 to col 56
                    int lx_count = 0;
                    for (int LX = LXMN; LX <= LXMX; LX += LXSTEP) {
                        if (lx_count > 0 && lx_count % 5 == 0)
                            std::printf("\n %54s", "");  // ' ' CC + 54 spaces → T56
                        std::printf("%s", fG(ALLOC(LLXCNT + LX), 14, 5).c_str());
                        lx_count++;
                    }
                }
                std::printf("\n");
            }

            LINECT += LINADD;
            if (PCONSW) LINECT += NSPLI;
            // FORMAT 548: blank separator at ANGBLK block boundaries (Fortran label 880)
            double origAngle = ANGMIN + (NANG - 1) * ANGSTP;
            if (fmod(origAngle + ANGSTP + 1e-7, ANGBLK) <= 1e-5) {
                if (LINECT < LINEMX) std::printf(" \n");
                LINECT++;
                if (LINECT + NUMBLK > LINEMX) LINECT = 1000;
            }
        } // end angle loop

        T5 = (float)second();
        TCROSS += T5 - T4;
        CRSFLP += 26 * NUMANG;

        // FORMAT 923: / '0TOTAL:', G15.5
        std::printf("\n0TOTAL:%s\n", fG(SIGTOT, 15, 5).c_str());
        double SAFE = std::max(SIGTOT, SMLNUM);
        if (NCHN <= 2) {
            // FORMAT 927: '+', T55, G14.5, G24.5
            std::printf("+%53s%s%s\n", "", fG(REACIN, 14, 5).c_str(), fG(REACOT, 24, 5).c_str());
            if (NMLX > 1) {
                // FORMAT 933: '+', (T91, 3G14.5)
                int count = 0;
                for (int LX = LXMN; LX <= LXMX; LX += LXSTEP) {
                    if (count % 3 == 0) {
                        if (count == 0) std::printf("+%89s", "");  // '+' + T91
                        else            std::printf("\n %89s", ""); // ' ' CC + T91
                    }
                    std::printf("%s", fG(ALLOC(LLXTOT + LX), 14, 5).c_str());
                    count++;
                }
                std::printf("\n");
            }
        } else if (NMLX > 1) {
            // FORMAT 938: '+', (T55, 5G14.5)
            int count = 0;
            for (int LX = LXMN; LX <= LXMX; LX += LXSTEP) {
                if (count % 5 == 0) {
                    if (count == 0) std::printf("+%53s", "");   // '+' + T55
                    else            std::printf("\n %53s", ""); // ' ' CC + T55
                }
                std::printf("%s", fG(ALLOC(LLXTOT + LX), 14, 5).c_str());
                count++;
            }
            std::printf("\n");
        }

        // FORMAT 950/955/963/968: convert LX totals to percent and print
        if (NMLX > 1) {
            for (int LX = LXMN; LX <= LXMX; LX += LXSTEP)
                ALLOC(LLXTOT + LX) = 100.0 * ALLOC(LLXTOT + LX) / SAFE;
            if (NCHN == 2) {
                // FORMAT 963: (T87, 3(F13.2, '%'))  — ' ' CC + 85 spaces = T87
                int count = 0;
                for (int LX = LXMN; LX <= LXMX; LX += LXSTEP) {
                    if (count % 3 == 0) {
                        if (count > 0) std::printf("\n");
                        std::printf(" %85s", "");
                    }
                    std::printf("%13.2f%%", ALLOC(LLXTOT + LX));
                    count++;
                }
                std::printf("\n");
            } else if (NCHN > 2) {
                // FORMAT 968: (T51, 5(F13.2, '%'))  — ' ' CC + 49 spaces = T51
                int count = 0;
                for (int LX = LXMN; LX <= LXMX; LX += LXSTEP) {
                    if (count % 5 == 0) {
                        if (count > 0) std::printf("\n");
                        std::printf(" %49s", "");
                    }
                    std::printf("%13.2f%%", ALLOC(LLXTOT + LX));
                    count++;
                }
                std::printf("\n");
            }
        }

        // FORMAT 983/993: JP/JT/LX/MX breakdown
        if (NSPLI > 1) {
            // FORMAT 983: / '0 JP  JT  LX MX      TOTAL      PERCENT' 5X '(...)' /
            std::printf("\n0 JP  JT  LX MX      TOTAL      PERCENT"
                        "%5s(VALUES FOR M > 0 NOT DOUBLED.)\n\n", "");
            for (int KOFFS = 1; KOFFS <= NSPLI; KOFFS++) {
                int JTT = ILLOC(LTOCS + 4 * KOFFS);
                if (JTT < 0) continue;
                int LX = ILLOC(LTOCS + 4 * KOFFS - 2);
                int MX = (ILLOC(LTOCS + 4 * KOFFS - 3) + LX + 1) / 2;
                double GR = 100 * ALLOC(LMXTOT + KOFFS - 1) / SAFE;
                // FORMAT 993: 2(I3,'/2'), 2I3, G15.5, F8.2
                std::printf("%3d/2%3d/2%3d%3d%s%8.2f\n",
                            ILLOC(LTOCS + 4 * KOFFS - 1), JTT,
                            LX, MX, fG(ALLOC(LMXTOT + KOFFS - 1), 15, 5).c_str(), GR);
            }
        }

        LFOFF += 2 * NSPLI * NUMANG;
    } // end NCHN loop

    // Rename S array to BETAS
    NAMCOM.NAMES[IRDINT] = BETNAM;
    Z[IPLM] = -Z[IPLM];
    Z[IMXCON] = -Z[IMXCON];
    Z[ILXCNT] = -Z[ILXCNT];
    Z[IFTEMP] = -Z[IFTEMP];
    Z[IFLO] = -Z[IFLO];
    Z[IFHI] = -Z[IFHI];

    // Analyzing powers
    T6 = (float)second();
    T7 = T6;
    T8 = T7;

    if (JBLOCK.JS[1] == 0) goto L2150;

    // Reaction analyzing powers
    if (LXMX == 0 || JPMX == 0) {
        std::printf("\n---- ALL REACTION ANALYZING POWERS ARE ZERO ----\n");
        goto L2100;
    }
    {
        int IPARIT = INTGER.PARITS[1] * INTGER.PARITS[2] * INTGER.PARITS[3] * INTGER.PARITS[4];
        if (IPARIT == 0 && TRANSW)
            IPARIT = 1 - 2 * ((FORMF_ref.LBP + FORMF_ref.LBT) % 2);
        if (IPARIT == 0) {
            std::printf("\n**** REACTION ANALYZING POWERS CANNOT BE CALCULATED"
                        " BECAUSE THE PARITY CHANGE IS UNKNOWN ****\n");
            goto L2100;
        }
        int IANPOW;
        ANAPOW(ANGMIN, ANGMAX, ANGSTP, (int)JBLOCK.JS[1], (int)JBLOCK.JS[2],
               (int)JBLOCK.JS[3], (int)JBLOCK.JS[4], NSPLI, IPARIT,
               1 /*true*/, PBUGSW ? 1 : 0, "ANPOW   ", ELAB, "REACTION",
               IF_var, ITOCS, IANPOW);
    }
    T7 = (float)second();

L2100:
    // Incoming elastic analyzing powers
    if (!WAVCOM.SOSWS[1]) {
        std::printf("\n---- ALL INCOMING ELASTIC ANALYZING POWERS ARE ZERO ----\n");
        goto L2150;
    }
    {
        int IAPOWA;
        ANAPOW(ANGMIN, ANGMAX, ANGSTP, (int)JBLOCK.JS[1], (int)JBLOCK.JS[1],
               (int)JBLOCK.JS[3], (int)JBLOCK.JS[3], WAVCOM.NUMJS[1], 1,
               1 /*true*/, PBUGSW ? 1 : 0, "ANPOWA  ", ELAB, "INCOMING",
               IFIN, WAVCOM.ITOCE[1], IAPOWA);
    }

L2150:
    // Outgoing elastic analyzing powers
    if (JBLOCK.JS[2] != 0) {
        if (WAVCOM.SOSWS[2]) {
            double EELAB = (1 + KANDM.RATMSS[2]) * (WAVCOM.ES[2] + Q);
            int IAPOWB;
            ANAPOW(ANGMIN, ANGMAX, ANGSTP, (int)JBLOCK.JS[2], (int)JBLOCK.JS[2],
                   (int)JBLOCK.JS[4], (int)JBLOCK.JS[4], WAVCOM.NUMJS[2], 1,
                   1 /*true*/, PBUGSW ? 1 : 0, "ANPOWB  ", EELAB, "OUTGOING",
                   IFOUT, WAVCOM.ITOCE[2], IAPOWB);
        } else {
            std::printf("\n---- ALL OUTGOING ELASTIC ANALYZING POWERS ARE ZERO ----\n");
        }
    }
    T8 = (float)second();

    // Print timing
    {
        double TIME = second() - TSTART;
        double TPHS = T1 - TSTART;
        double T1_v = T2 - T1;
        T6 = T7 - T6;
        T7 = T8 - T7;
        double T9 = TIME - TPHS - T1_v - TBET - TCROSS - T6 - T7;
        CRSFLP = 1.0e-6 * CRSFLP;
        double CRSFPS = (TCROSS > 0) ? CRSFLP / TCROSS : 0;
        // FORMAT 2583: '0'/'0'/'0CPU TIMES...'/'0PRINTING...'/.../F17.3//`0`
        // T30 = F10.3/F17.3 starts at col 30 (29 chars before it)
        std::printf("0\n0\n0CPU TIMES FOR CROSS SECTIONS (SECONDS):\n"
                    "0PRINTING PHASE SHIFTS%7s%10.3f\n"
                    " ELASTIC CROSS SECTIONS%6s%10.3f\n"
                    " BETA'S (SUM ON LI)%10s%10.3f\n"
                    " REACTION CROSS SECTIONS%5s%10.3f%10s(%7.2f MFLOPS =%7.3f MFLOP/SEC )\n"
                    " REACTION ANALYZING POWERS%3s%10.3f\n"
                    " ELASTIC ANALYZING POWERS%4s%10.3f\n"
                    " ALL OTHER TIME%14s%10.3f\n"
                    " TOTAL TIME%18s%17.3f\n\n0\n",
                    "", TPHS, "", T1_v, "", TBET,
                    "", TCROSS, "", CRSFLP, CRSFPS,
                    "", T6, "", T7, "", T9, "", TIME);
    }

    IRTN = 1;
    return;
}
