// phiffer.cpp -- translated from phiffer.f
//
// PHIFFER - linkule for phiffer-computed bound states
// PHIFER  - generates single-particle radial function
// STOP_F  - print error message and stop (renamed to avoid clash with C stop)
//
// 2/13/07 - new linkule based on av18 somewhat
//
// Translated from Fortran to C++ preserving all logic and labels.

#include "ptolemy_types.h"
#include "ptolemy_commons.h"
#include "ptolemy_forward.h"
#include "ptolemy_intrinsics.h"
#include "ptolemy_io.h"
#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <algorithm>

// Forward declarations for routines in this file
void PHIFER(int itype, int ievopt, int nodes, int l, int jtwo,
            double h2o2mu, double* params, bool prntsw, int npts,
            double rmax, double* rgrid, double* phi_by_rl,
            int nvgrid, double* vongrid, double& vrgrid);
void STOP_F(const char* msg, int icode);

// ============================================================================
// SUBROUTINE PHIFFER
// ============================================================================

void phiffer(char8 alias, real4* savint, int ipotyp, int ireque, int& iretur,
             int iunit, int& numout, int l, double& jp, double rstart,
             double stepsz, int numpts, double* array1, double* array2,
             double* flt, double* tempvs, double* param, int* intger,
             double* jblock, int* iswtch, double* intrnl, double* fintrn,
             double* cnstnt, double* wavbl, double* fkandm,
             char* id, double* alloc, int* illoc, int facfr4,
             integer4* loc, int* length,
             int (*nalloc)(const char*, int), int (*namloc)(const char*),
             int* iparam)
{
    // implicit real*8 (a-h, o-z)

    // Local persistent variables (SAVE)
    static double params_loc[15];  // 1-based: params(14) -> [15]
    static int nodes, izp, izt;
    static double v, rv, rc, av, e, am;
    static double rhowb, alphawb;
    static double vso, rso, aso;

    // Local
    double undef, rt, rend, rt4pi, hbarc, hb2o2m, afine, alpha_coul;
    int notdef;
    bool wfsw;
    int itype, ievopt;
    double x;

    iretur = 0;
    numout = 0;

    undef = fintrn[1];
    notdef = (int)intrnl[3];

    if (ireque - 2 < 0) goto L100;
    if (ireque - 2 == 0) goto L200;
    goto L300;

    // =========================================================================
    // Setup (ireque = 1): check for errors
    // =========================================================================
L100:
    if (ipotyp != 6) {
        printf(" **** bad ipotyp in %.8s %d\n", alias.data, ipotyp);
        numout = 1;
        iretur = -1;
        return;
    }

    // set jp, nodes, v, r, a, e, j, jsp, jst, spam
    if (l == notdef || (int)jp == notdef) {
        printf(" L or JP not defined: %d %d\n", l, (int)jp);
        numout = 1;
        iretur = -1;
        return;
    }

    nodes = intger[12];
    v = flt[62];
    if (v == undef) v = 60.0;
    rv = flt[43];
    rc = flt[45];
    av = flt[1];
    e = flt[12];
    am = flt[24];
    if (e == undef || am == undef) {
        printf(" E or red. mass undefined: %g %g\n", e, am);
        numout = 1;
        iretur = -1;
        return;
    }
    if (rv == undef || av == undef) {
        printf(" R or A undefined: %g %g\n", rv, av);
        numout = 1;
        iretur = -1;
        return;
    }
    if (rc == undef) rc = rv;
    rhowb = param[1];
    alphawb = param[2];
    if (alphawb == undef) alphawb = 0.0;
    if (alphawb == 0.0) rhowb = 1.0;
    if (rhowb == undef) {
        printf(" rho(wb) undefined, rho, alpha: %g %g\n", param[1], param[2]);
        numout = 1;
        iretur = -1;
        return;
    }

    izp = intger[24];
    izt = intger[25];
    if (izp == notdef || izt == notdef) {
        printf(" zp or zt not defined %d %d\n", izp, izt);
        numout = 1;
        iretur = -1;
        return;
    }

    vso = flt[64];
    rso = flt[49];
    aso = flt[7];
    if (vso == undef) vso = 0.0;
    if (vso != 0.0) {
        if (rso == undef || aso == undef) {
            printf(" rso or vso undef\" vso, rso, aso = %g %g %g\n", vso, rso, aso);
            numout = 1;
            iretur = -1;
            return;
        }
    } else {
        if (rso == undef) rso = 1.0;
        if (aso == undef) aso = 1.0;
    }

    return;

    // =========================================================================
    // Printing (ireque = 2)
    // =========================================================================
L200:
    if (ipotyp == 6) {
        // Printout for wave function calculation.
        printf(" Phiffer calculation of wave function\n");
        printf("L, nodes, jp =%4d%4d%4d/2\n", l, nodes, (int)jp);
        printf("E, mu =%10.3f%10.3f\n", e, am);
        printf("V(guess), V(convrg) =%10.3f%10.3f\n", v, -params_loc[10]);
        printf("R, A =%10.3f%10.3f\n", rv, av);
        printf("w.b. rho, alpha =%10.3f%10.3f\n", rhowb, alphawb);
        printf("Vso, Rso, Aso =%10.3f%10.3f%10.3f\n", vso, rso, aso);
        printf("Zp, Zt, Rc =%4d%4d%10.3f\n", izp, izt, rc);
        numout = 7;
    } else {
        printf(" ******* we should not be here****\n");
        numout = 1;
    }

    return;

    // =========================================================================
    // Calculation (ireque = 3)
    // =========================================================================
L300:
    am = flt[24];
    e = flt[12];
    rt4pi = cnstnt[2];
    hbarc = cnstnt[6];

    hb2o2m = hbarc * hbarc / (2.0 * am);
    afine = cnstnt[8];
    alpha_coul = izp * izt * hbarc / afine;

    wfsw = (ipotyp == 6);

    rt = rstart;
    rend = rstart + (numpts - 1) * stepsz;

    if (rt != 0.0) {
        printf(" rstart not 0: %g\n", rt);
        iretur = -1;
        return;
    }

    for (int i = 1; i <= 14; i++) params_loc[i] = 0.0;
    params_loc[1] = rv;
    params_loc[2] = av;
    params_loc[3] = e;
    params_loc[4] = rhowb;
    params_loc[5] = alphawb;
    params_loc[6] = rso;
    params_loc[7] = aso;
    params_loc[8] = vso;
    params_loc[10] = v;
    params_loc[11] = alpha_coul;
    params_loc[12] = rc;

    // allocate rgrid
    double* rgrid = new double[numpts + 2];  // 1-based

    itype = 11;
    ievopt = 1;

    PHIFER(itype, ievopt, nodes, l, (int)jp, hb2o2m,
           params_loc, true, numpts, rend, rgrid,
           array1, numpts, array2, x);

    // phifer returns phi/r^L. We want u = r^(L+1) * phi/r^L
    for (int i = 1; i <= numpts; i++) {
        array1[i] = std::pow(rgrid[i], l + 1) * array1[i];
    }

    flt[62] = -params_loc[10];

    delete[] rgrid;
    return;
}


// ============================================================================
// SUBROUTINE PHIFER
//
// Generates single-particle radial function.
// Returns phi_by_rl = phi / r**L
// Integral r^2 phi^2 dr = 1
// ============================================================================

void PHIFER(int itype, int ievopt, int nodes, int l, int jtwo,
            double h2o2mu, double* params, bool prntsw, int npts,
            double rmax, double* rgrid, double* phi_by_rl,
            int nvgrid, double* vongrid, double& vrgrid)
{
    // implicit real*8 (a-h, o-z)

    bool infwellsw, varye, scatsw, eminsw, emaxsw, gridsw;
    int nodval[3];    // 1-based
    double xval[3];   // 1-based

    // Allocatable arrays (1-based)
    double* ws     = nullptr;
    double* psir   = nullptr;
    double* centrg = nullptr;
    double* splb   = nullptr;
    double* splc   = nullptr;
    double* spld   = nullptr;

    double wsr, wsa, wse, rho, alpha, rso_l, aso_l, vso_l;
    double alpha_fermi, wsv, alpha_coul, rc, rb, xa, xb;
    double sml, dr, dx, dl, consso, deriv;
    double estrt, vstrt, vme;
    double ak2, akap, gg1, gg2, ff1, ff2, dldif;
    double dldif1, dldif2, wse1, wse2, wsv1, wsv2;
    double efac, efac0, emin, emax;
    double sum1, sum2, anorm, rms, fac;
    double r, vc, so, x_loc, f_loc, pi;
    int iu, it, ih, nptsm, nptsm2;
    int match, matchp, nodlop, nodelast, nodhave, nod;
    int maxgrid, loop;

    if (prntsw) {
        printf("phifer input:%5d%3d%3d%3d%4d/2%6d%6.1f%10.4f\n",
               itype, ievopt, nodes, l, jtwo, npts, rmax, h2o2mu);
        printf(" %10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n",
               params[1], params[2], params[3], params[4],
               params[5], params[6], params[7], params[8]);
        printf(" %10.4f%10.4f%10.4f%10.4f%10.4f\n",
               params[9], params[10], params[11], params[12], params[13]);
    }

    maxgrid = npts;
    // gridsw not yet set, initialize
    iu = itype % 10;
    gridsw = (iu == 9);
    if (gridsw) maxgrid = std::max(maxgrid, nvgrid);

    // Allocate work arrays (1-based)
    ws     = new double[maxgrid + 1];
    psir   = new double[maxgrid + 1];
    centrg = new double[maxgrid + 1];
    splb   = new double[maxgrid + 1];
    splc   = new double[maxgrid + 1];
    spld   = new double[maxgrid + 1];

    it = (itype / 10) % 10;
    ih = (itype / 100) % 10;
    infwellsw = (ih == 1);
    scatsw = (ih == 3);
    varye = (ievopt == 0);

    sml = 1.e-10;
    dr = rmax / (npts - 1);
    dx = dr * dr / (12.0 * h2o2mu);
    nptsm = npts - 1;
    nptsm2 = npts - 2;
    for (int i = 1; i <= npts; i++) {
        rgrid[i] = dr * (i - 1);
    }

    wsr = params[1];
    wsa = params[2];
    wse = params[3];
    rho = params[4];
    alpha = params[5];
    rso_l = params[6];
    aso_l = params[7];
    vso_l = params[8];
    alpha_fermi = params[9];
    wsv = params[10];
    alpha_coul = params[11];
    rc = params[12];
    rb = params[13];
    xa = rc / std::sqrt(12.0);
    xb = rb / std::sqrt(12.0);

    if (iu == 8) goto L800;

    if (scatsw) deriv = params[14];
    if (rho == 0.0) rho = 1.0;
    dl = l * (l + 1) * dr * dr / 12.0;
    consso = 0.0;
    if (jtwo != -1) {
        consso = -2.0 * vso_l * (jtwo * (jtwo + 2) / 4.0 - l * (l + 1) - 0.75) / aso_l;
    }

    if (gridsw) {
        SPLNCB(nvgrid, &vrgrid, vongrid, splb + 1, splc + 1, spld + 1);
        INTRPC(nvgrid, &vrgrid, vongrid, splb + 1, splc + 1, spld + 1,
               npts, rgrid + 1, ws + 1);
    }

    // Build potential and centrifugal arrays
    for (int i = 1; i <= npts; i++) {
        r = dr * (i - 1);

        switch (iu) {
            case 1:
                ws[i] = 1.0 / (1.0 + std::exp((r - wsr) / wsa));
                if (wsr > 0.0) ws[i] = ws[i] * (1.0 + (alpha_fermi / (wsr * wsr)) * r * r);
                ws[i] = ws[i] - alpha * std::exp(-(r / rho) * (r / rho));
                break;
            case 2:
                ws[i] = std::exp(-((r - wsr) / wsa) * ((r - wsr) / wsa));
                break;
            case 9:
                // potential already in ws from spline interpolation
                break;
            default:
                STOP_F("PHIFER: invalid itype", itype);
                break;
        }

        ws[i] = dx * ws[i];
        centrg[i] = dl / ((r + 1.e-20) * (r + 1.e-20));
        if (consso != 0.0) {
            x_loc = std::exp((r - rso_l) / aso_l);
            so = consso * x_loc / ((r + 1.e-20) * (1.0 + x_loc) * (1.0 + x_loc));
            centrg[i] = centrg[i] + dx * so;
        }

        switch (it) {
            case 0:
                vc = 0.0;
                break;
            case 1:
                if (r > rc) {
                    vc = alpha_coul / r;
                } else {
                    vc = (alpha_coul / (2.0 * rc)) * (3.0 - (r / rc) * (r / rc));
                }
                break;
            case 2:
                if (r > 1.e-3) {
                    vc = alpha_coul * (1.0
                         - (0.5 / ((1.0 - (xa / xb) * (xa / xb)) * (1.0 - (xa / xb) * (xa / xb))))
                           * std::exp(-r / xb) * (2.0 + (r / xb) + 4.0 / (1.0 - (xb / xa) * (xb / xa)))
                         - (0.5 / ((1.0 - (xb / xa) * (xb / xa)) * (1.0 - (xb / xa) * (xb / xa))))
                           * std::exp(-r / xa) * (2.0 + (r / xa) + 4.0 / (1.0 - (xa / xb) * (xa / xb)))
                         ) / r;
                } else {
                    vc = (alpha_coul / (2.0 * std::pow(xb * xb - xa * xa, 3)))
                       * (xb * xb * xb * (xb * xb - 5.0 * xa * xa)
                        + xa * xa * xa * (5.0 * xb * xb - xa * xa));
                }
                break;
            default:
                STOP_F("PHIFER: invalid itype", itype);
                break;
        }
        centrg[i] = centrg[i] + dx * vc;
    }

    if (ievopt == 9) {
        for (int i = 1; i <= npts; i++) {
            phi_by_rl[i] = (wsv * ws[i] + centrg[i]) / dx;
        }
        goto L900;
    }

    nodlop = 1;
    efac0 = 1.3;
    efac = efac0;
    nodelast = 9999;
    nodhave = 0;
    eminsw = false;
    emaxsw = false;
    emin = 0.0;
    emax = 0.0;
    if (!infwellsw) {
        if (!gridsw) wsv = -std::fabs(wsv);
        wse = std::fabs(wse);
    }

    // ===== Loop back (label 15) =====
L15:
    estrt = wse;
    vstrt = wsv;
    vme = (wsv * ws[1] / dx + wse) / h2o2mu;
    match = (int)(std::sqrt(-l * (l + 1) / vme) / dr + 0.5);
    match = std::max(match, npts / 20);
    matchp = match + 1;

    // Converge on the eigenvalue
    for (loop = 1; loop <= 100; loop++) {

        // Integrate from outside in
        ak2 = -dx * wse;
        akap = std::sqrt(+wse / h2o2mu);
        if (infwellsw) {
            psir[npts] = 0.0;
            psir[nptsm] = 0.01;
        } else if (scatsw) {
            ak2 = -ak2;
            if (deriv == 0.0) {
                psir[npts] = 0.0;
                psir[nptsm] = 0.01;
            } else {
                psir[npts] = 1.0;
                double d1 = deriv + 1.0 / rmax;
                double d2 = -akap * akap + l * (l + 1) / (rmax * rmax);
                psir[nptsm] = 1.0 - d1 * dr + 0.5 * d2 * dr * dr;
            }
        } else {
            psir[npts] = std::exp(-akap * rgrid[npts]);
            psir[nptsm] = std::exp(-akap * rgrid[nptsm]);
        }

        double xa_l = psir[npts];
        double xb_l = psir[nptsm];
        double xd = wsv * ws[npts] + centrg[npts] - ak2;
        double xe = wsv * ws[nptsm] + centrg[nptsm] - ak2;
        for (int j = nptsm2; j >= match; j--) {
            double xf = wsv * ws[j] + centrg[j] - ak2;
            psir[j] = ((2.0 + 10.0 * xe) * xb_l - (1.0 - xd) * xa_l) / (1.0 - xf);
            xa_l = xb_l;
            xb_l = psir[j];
            xd = xe;
            xe = xf;
        }
        gg1 = xb_l;
        gg2 = xa_l;

        // Integrate from inside out
        psir[1] = 0.0;
        psir[2] = (1.0 + 12.0 * (wsv * ws[2] - ak2) / (4 * l + 6)) * std::pow(dr, l + 1);
        xb_l = psir[2];
        xe = wsv * ws[2] + centrg[2] - ak2;
        double xf = wsv * ws[3] + centrg[3] - ak2;
        psir[3] = (2.0 + 10.0 * xe) * xb_l / (1.0 - xf);
        if (l == 1) psir[3] = psir[3] + 2.0 * h2o2mu * dx / (1.0 - xf);
        xa_l = xb_l;
        xb_l = psir[3];
        xd = xe;
        xe = xf;
        for (int j = 4; j <= matchp; j++) {
            xf = wsv * ws[j] + centrg[j] - ak2;
            psir[j] = ((2.0 + 10.0 * xe) * xb_l - (1.0 - xd) * xa_l) / (1.0 - xf);
            xa_l = xb_l;
            xb_l = psir[j];
            xd = xe;
            xe = xf;
        }
        ff1 = xa_l;
        ff2 = xb_l;
        dldif = (ff2 * gg1 - ff1 * gg2) / (std::fabs(ff1 * gg1) + std::fabs(ff2 * gg2));

        if (loop == 1) {
            dldif2 = dldif;
            wse2 = wse;
            wsv2 = wsv;
            if (varye) {
                wse = wse2 * 1.2;
            } else {
                wsv = wsv2 * 1.2;
            }
            continue;  // go to 80
        } else if (loop == 2) {
            dldif1 = dldif;
            wse1 = wse;
            wsv1 = wsv;
        } else {
            if (std::fabs(dldif) <= sml) goto L90;
            // Arithmetic IF replacements (three-way branch)
            if (dldif * dldif1 <= 0.0) {
                // L54
                if (dldif * dldif2 <= 0.0) {
                    // L58
                    if (std::fabs(dldif1) - std::fabs(dldif2) <= 0.0) {
                        // L60
                        goto L60;
                    } else {
                        // L65
                        goto L65;
                    }
                } else {
                    // L60
                    goto L60;
                }
            } else {
                // L56
                if (dldif * dldif2 <= 0.0) {
                    // L65
                    goto L65;
                } else {
                    // L58
                    if (std::fabs(dldif1) - std::fabs(dldif2) <= 0.0) {
                        // L60
                        goto L60;
                    } else {
                        // L65
                        goto L65;
                    }
                }
            }

        L60:
            dldif2 = dldif1;
            wse2 = wse1;
            wsv2 = wsv1;
            dldif1 = dldif;
            wse1 = wse;
            wsv1 = wsv;
            goto L70;

        L65:
            dldif1 = dldif;
            wse1 = wse;
            wsv1 = wsv;
        }

    L70:
        if (varye) {
            wse = (dldif2 * wse1 - dldif1 * wse2) / (dldif2 - dldif1);
            if (wse < 0.0 && !infwellsw) {
                wse = 0.5 * std::min(wse1, wse2);
                if (wse < 0.001) {
                    if (prntsw)
                        printf(" wse tried for 0, start= %g\n", estrt);
                    wse = 2.0 * estrt;
                    goto L15;
                }
            }
            if (std::min(std::fabs(wse - wse1), std::fabs(wse - wse2))
                / std::fabs(wse1 - wse2) < 0.1)
                wse = 0.5 * (wse1 + wse2);
            if (eminsw) {
                if (wse < emin) {
                    wse = 0.5 * (emin + std::min(wse1, wse2));
                }
            }
            if (emaxsw) {
                if (wse < emax) {
                    wse = 0.5 * (emax + std::max(wse1, wse2));
                }
            }
        } else {
            wsv = (dldif2 * wsv1 - dldif1 * wsv2) / (dldif2 - dldif1);
            // Don't let wsv change sign, but it could be either sign
            if (wsv * vstrt <= 0.0) {
                wsv = std::copysign(0.5 * std::min(std::fabs(wsv1), std::fabs(wsv2)), vstrt);
            }
            if (std::min(std::fabs(wsv - wsv1), std::fabs(wsv - wsv2))
                / std::fabs(wsv1 - wsv2) < 0.1)
                wsv = 0.5 * (wsv1 + wsv2);
        }
    }  // end loop 80

    // Did not converge
    if (std::fabs(dldif) > 10000.0 * sml) {
        printf("\n*********** did not converge on wavefunction ********************\n");
        printf("%25.15g%25.15g\n", wse, wsv);
        printf("%25.15g%25.15g%25.15g\n", wse1, wsv1, dldif1);
        printf("%25.15g%25.15g%25.15g\n", wse2, wsv2, dldif2);
        printf("%25.15g%25.15g%25.15g%25.15g\n", ff1, gg1, ff2, gg2);
        printf("******** phiffer could not converge ********\n");
    } else {
        if (prntsw)
            printf("phifer accepts poor converg%15.5g%15.5g\n", dldif, sml);
    }

    // Normalize inside r(matchp)
L90:
    fac = gg2 / psir[matchp];
    for (int i = 2; i <= matchp; i++) {
        psir[i] = fac * psir[i];
    }

    // Normalize and count nodes
    sum1 = 0.0;
    sum2 = 0.0;
    f_loc = psir[4];
    nod = 0;
    for (int i = 2; i <= npts; i++) {
        sum1 = sum1 + psir[i] * psir[i];
        sum2 = sum2 + psir[i] * psir[i] * rgrid[i] * rgrid[i];
        phi_by_rl[i] = psir[i] / std::pow(rgrid[i], l + 1);
        if (psir[i] * f_loc <= 0.0 && i > 4 && i < npts) {
            nod = nod + 1;
            f_loc = psir[i + 1];
        }
    }
    if (nod == nodes) goto L130;

    // Wrong number of nodes; try again
    if (prntsw) {
        printf(" got%3d nodes:  e =%15.5g%15.5g efac=%15.5g     v =%10.5f%10.5f%5d matching loops\n",
               nod, estrt, wse, efac, vstrt, wsv, loop);
    }
    nodlop = nodlop + 1;
    if (nodlop > 40) {
        printf(" *** node loop limit reached: %d %g %g %g %g %g\n",
               nod, estrt, wse, efac, vstrt, wsv);
        goto L130;
    }

    if (varye) {
        // Varying e
        if (nod == nodelast) {
            efac = 1.1 * efac;
        } else {
            efac = efac0;
        }
        nodelast = nod;

        // Set barriers on e based on what we have found
        double xdel = 0.1 * std::fabs(wse);
        if (nod > nodes) {
            if (!eminsw || (eminsw && wse - xdel > emin)) {
                eminsw = true;
                emin = wse - xdel;
            }
        } else {
            if (!emaxsw || (emaxsw && wse + xdel < emax)) {
                emaxsw = true;
                emax = wse + xdel;
            }
        }

        // Accumulate eigenvalues and number of nodes
        if (nodhave == 0) {
            nodhave = 1;
            nodval[1] = nod;
            xval[1] = wse;
            goto L120;
        } else if (nodhave == 1) {
            if (nod == nodval[1]) goto L120;
            nodhave = 2;
            nodval[2] = nod;
            xval[2] = wse;
        } else {
            if (nod == nodval[1] || nod == nodval[2]) goto L120;
            // Replace the worse number of nodes
            int idx;
            if (std::abs(nodval[1] - nodes) > std::abs(nodval[2] - nodes)) {
                idx = 1;
            } else {
                idx = 2;
            }
            if (std::abs(nodval[idx] - nodes) > std::abs(nod - nodes)) {
                nodval[idx] = nod;
                xval[idx] = wse;
            } else {
                goto L120;
            }
        }

        // Predict the new eigenval based on number of nodes found
        x_loc = xval[1] + ((xval[2] - xval[1]) / (double)(nodval[2] - nodval[1]))
              * (nodes - nodval[1]);
        wse = x_loc;
        goto L125;

    L120:
        if (nod < nodes) {
            wse = std::min(wse, estrt);
            if (wse > 0.0) {
                wse = wse / efac;
            } else {
                wse = wse * efac;
            }
        } else {
            wse = std::max(wse, estrt);
            if (wse > 0.0) {
                wse = wse * efac;
            } else {
                wse = wse / efac;
            }
        }
    } else {
        // Varying v
        if (nod < nodes) {
            wsv = std::copysign(1.3 * std::max(std::fabs(wsv), std::fabs(vstrt)), vstrt);
        } else {
            wsv = std::copysign(std::min(std::fabs(wsv), std::fabs(vstrt)) / 1.2, vstrt);
        }
    }
L125:
    goto L15;

    // Have correct number of nodes
L130:
    phi_by_rl[1] = (rgrid[3] * rgrid[3] * phi_by_rl[2]
                   - rgrid[2] * rgrid[2] * phi_by_rl[3])
                 / (rgrid[3] * rgrid[3] - rgrid[2] * rgrid[2]);
    anorm = std::sqrt(sum1 * dr);
    for (int i = 1; i <= npts; i++) {
        phi_by_rl[i] = phi_by_rl[i] / anorm;
    }
    rms = std::sqrt(sum2 / sum1);
    if (prntsw) {
        printf(" woods-saxon strength v =%15.10f     e(sep) =%15.10f;   rms radius =%8.5f%5d iters;\n",
               wsv, wse, rms, loop);
        printf(" convrg =%12.3g     phi/r**l:%18.9g%18.9g%18.9g%18.9g%18.9g\n",
               dldif, phi_by_rl[1], phi_by_rl[2], phi_by_rl[3], phi_by_rl[4], phi_by_rl[5]);
    }

    // If we completely failed to converge, stop
    if (std::fabs(dldif) > 1000.0 * sml) {
        printf(" could not converge: %g %g\n", dldif, 1000.0 * sml);
    }
    params[3] = wse;
    params[10] = wsv;

    // PTOLEMY: return potential shape
    for (int i = 1; i <= nvgrid; i++) {
        vongrid[i] = ws[i] / dx;
    }

    goto L900;

    // Return a h.o. solution
L800:
    if (nodes != 0) {
        std::exit(9876);
    }
    pi = std::acos(-1.0);
    anorm = std::pow(2.0 / wsa, 2 * l + 3) * std::sqrt(2.0 / pi);
    for (int ll = 1; ll <= 2 * l + 1; ll += 2) {
        anorm = anorm / ll;
    }
    anorm = std::sqrt(anorm);
    printf(" H.O.; L= %d  A= %g  Norm= %g\n", l, wsa, anorm);
    for (int i = 1; i <= npts; i++) {
        x_loc = rgrid[i] / wsa;
        // We leave out the r**l because phifer returns phi/r**l
        phi_by_rl[i] = anorm * std::exp(-(x_loc * x_loc));
    }
    delete[] ws; delete[] psir; delete[] centrg;
    delete[] splb; delete[] splc; delete[] spld;
    return;

L900:
    delete[] ws; delete[] psir; delete[] centrg;
    delete[] splb; delete[] splc; delete[] spld;
    return;
}


// ============================================================================
// SUBROUTINE STOP (renamed STOP_F to avoid C/C++ conflict)
//
// Print an error message and stop.
// ============================================================================

void STOP_F(const char* msg, int icode)
{
    printf("\n %s\n", std::string(78, '*').c_str());
    printf(" *%20sStopping because of an error  !!%24s*\n", "", "");
    printf(" * %s : %10d%*s*\n", msg, icode, (int)(74 - strlen(msg) - 13), "");
    printf(" %s\n\n", std::string(78, '*').c_str());

    fprintf(stderr, "\n");
    fprintf(stderr, " %s\n", std::string(78, '*').c_str());
    fprintf(stderr, " *%20sStopping because of an error  !!%24s*\n", "", "");
    fprintf(stderr, " * %s : %10d%*s*\n", msg, icode, (int)(74 - strlen(msg) - 13), "");
    fprintf(stderr, " %s\n\n", std::string(78, '*').c_str());

    std::exit(9999);
}
