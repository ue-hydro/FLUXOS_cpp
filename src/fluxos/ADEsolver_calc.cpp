
// Copyright 1992: Cornel Beffa and Roland Faeh
// Copyright 2013: Kashif Shaad and Diogo Costa
// Copyright 2019, Diogo Costa

// This program, FLUXOS, is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) aNCOLS later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <iostream>
#include <armadillo>
#include <memory>

#include "GlobVar.h"
#include "ADEsolver_calc.h"

// ADE solver - optimized version
void adesolver_calc(
    GlobVar& ds,
    int it,
    int ichem)
{
    // Cache frequently used values
    const double dx = ds.dxy;
    const double dy = ds.dxy;
    const double dtfl = ds.dtfl;
    const double hdry = ds.hdry;
    const double D_coef = ds.D_coef;
    const double area = ds.arbase;
    const unsigned int NROWS = ds.NROWS;
    const unsigned int NCOLS = ds.NCOLS;
    const double nt = 1.0;  // eddy viscosity (m2/s)
    const double sigc = 0.5;

    // Get raw pointers for faster access (avoid repeated unique_ptr dereference)
    arma::Mat<double>& conc = (*ds.conc_SW)[ichem];
    const arma::Mat<double>& h_mat = *ds.h;
    const arma::Mat<double>& h0_mat = *ds.h0;
    const arma::Mat<float>& ldry_mat = *ds.ldry;
    const arma::Mat<float>& ldry_prev_mat = *ds.ldry_prev;
    const arma::Mat<double>& fe_1_mat = *ds.fe_1;
    const arma::Mat<double>& fn_1_mat = *ds.fn_1;

    // Local matrices for bounds calculation
    arma::mat cmaxr(ds.MROWS, ds.MCOLS);
    arma::mat cminr(ds.MROWS, ds.MCOLS);

    if(it > 1) {
        // ADJUST CONCENTRATION TO NEW DEPTH - parallelized
        #pragma omp parallel for collapse(2) schedule(static)
        for(unsigned int iy = 1; iy <= NCOLS; iy++) {
            for(unsigned int ix = 1; ix <= NROWS; ix++) {
                // Calculate bounds from neighbors
                cmaxr(ix, iy) = std::fmax(
                    conc(ix-1, iy),
                    std::fmax(
                        conc(ix+1, iy),
                        std::fmax(conc(ix, iy-1), conc(ix, iy+1))));

                cminr(ix, iy) = std::fmin(
                    conc(ix-1, iy),
                    std::fmin(conc(ix+1, iy),
                        std::fmin(conc(ix, iy-1), conc(ix, iy+1))));

                double hnew = h_mat(ix, iy);

                if(ldry_mat(ix, iy) == 0 && ldry_prev_mat(ix, iy) == 0) {
                    conc(ix, iy) = conc(ix, iy) * h0_mat(ix, iy) / hnew;
                } else if(ldry_mat(ix, iy) == 1) {
                    conc(ix, iy) = 0.0;
                }
            }
        }
    } else {
        return;
    }

    // Row-wise flux storage for better cache locality
    // Process column by column (Y direction) to enable some parallelization
    #pragma omp parallel
    {
        // Thread-local storage for flux propagation
        arma::vec qfcds_local(NROWS + 2, arma::fill::zeros);

        #pragma omp for schedule(static)
        for(unsigned int iy = 1; iy <= NCOLS; iy++) {
            double pfe = 0.0;
            qfcds_local.zeros();

            for(unsigned int ix = 1; ix <= NROWS; ix++) {
                double pfw, qfs, qfn, pfce, pfde, fe, fp, hne, fem;
                double he, hp, hn, qxl, qyl, fw, fee, fs, fn, fnn;
                double fnm, qfcn, qfdn, cvolrate, cf, cbilan, dc, cvolpot, cvolrat, con;
                double hnue, hnn;

                const unsigned int is = iy - 1;
                const unsigned int in = iy + 1;
                const unsigned int inn = std::min(iy + 2, NCOLS + 1);
                const unsigned int iw = ix - 1;
                const unsigned int ie = ix + 1;
                const unsigned int iee = std::min(ix + 2, NROWS + 1);

                // BC at start of row
                if(ix == 1) {
                    pfce = conc(0, iy) * fe_1_mat(0, iy) * dy;
                    hp = std::fmax(h_mat(1, iy), hdry);
                    he = std::fmax(h_mat(2, iy), hdry);
                    fp = conc(0, iy);
                    fe = conc(1, iy);
                    hne = std::sqrt(hp * nt * he * nt) / sigc / std::fabs(dx) * dy * D_coef;
                    pfde = 0.0;
                    pfe = pfce;
                }

                // CHECK IF THE DOMAIN IS DRY
                if(ldry_mat(ix, iy) == 1) {
                    pfe = 0.0;
                    qfcds_local(ix) = 0.0;
                    conc(ix, iy) = 0.0;
                    continue;
                }

                // INITIALIZATION
                hp = h_mat(ix, iy);
                he = h_mat(ie, iy);
                hn = h_mat(ix, in);
                qxl = fe_1_mat(ix, iy);
                qyl = fn_1_mat(ix, iy);
                fw = conc(iw, iy);
                fp = conc(ix, iy);
                fe = conc(ie, iy);
                fee = conc(iee, iy);
                fs = conc(ix, is);
                fn = conc(ix, in);
                fnn = conc(ix, inn);

                // FLUXES OVER WEST AND SOUTH FACES
                pfw = pfe;
                qfs = qfcds_local(ix);

                // X-DIRECTION
                if(ldry_mat(ie, iy) == 0) {
                    hnue = std::fmax(hp * nt * he * nt, 0.0001);
                    hne = std::sqrt(hnue) / sigc / dx * dy * D_coef;
                    pfde = -hne * (fe - fp);

                    if(qxl > 0.0) {
                        if(ldry_mat(iw, iy) == 0) {
                            fem = -0.125 * fw + 0.75 * fp + 0.375 * fe;
                        } else {
                            fem = 0.5 * fp + 0.5 * fe;
                        }
                    } else {
                        if(ldry_mat(iee, iy) == 0) {
                            fem = 0.375 * fp + 0.75 * fe - 0.125 * fee;
                        } else {
                            fem = 0.5 * fp + 0.5 * fe;
                        }
                    }
                } else {
                    fem = 0.0;
                    pfde = 0.0;
                }

                fem = std::fmax(0.0, fem);

                if(ix == NROWS) {
                    fem = conc(NROWS + 1, iy);
                }

                pfce = qxl * fem * dy;
                pfe = pfce + pfde;

                if(pfe < 0) {
                    if(ldry_mat(ie, iy) == 0) {
                        cvolrate = -(fe * he) * area / dtfl;
                        pfe = std::fmax(pfe, cvolrate);
                    } else {
                        pfe = 0.0;
                    }
                }

                // Y-DIRECTION
                if(ldry_mat(ix, in) == 0) {
                    hnue = std::fmax(0.0001, hp * nt * hn * nt);
                    hnn = std::sqrt(hnue) / sigc / dy * dx * D_coef;
                    qfdn = -hnn * (fn - fp);

                    if(qyl > 0.0) {
                        if(ldry_mat(ix, is) == 0) {
                            fnm = -0.125 * fs + 0.75 * fp + 0.375 * fn;
                        } else {
                            fnm = 0.5 * fp + 0.5 * fn;
                        }
                    } else {
                        if(ldry_mat(ix, inn) == 0) {
                            fnm = 0.375 * fp + 0.75 * fn - 0.125 * fnn;
                        } else {
                            fnm = 0.5 * fp + 0.5 * fn;
                        }
                    }
                } else {
                    fnm = 0.0;
                    qfdn = 0.0;
                }

                fnm = std::fmax(0.0, fnm);

                if(iy == NCOLS) {
                    fnm = conc(ix, NCOLS + 1);
                }

                qfcn = qyl * fnm * dx;
                qfn = qfcn + qfdn;

                if(qfn < 0) {
                    if(ldry_mat(ix, in) == 0) {
                        cvolrate = -(fn * hn) * area / dtfl;
                        qfn = std::fmax(qfn, cvolrate);
                    } else {
                        qfn = 0.0;
                    }
                }

                // CHECK AVAILABLE MATERIAL
                cvolpot = (fp * hp) * area;
                cvolrat = cvolpot / dtfl + pfw + qfs;

                if(cvolrat > 0.0) {
                    if(pfe > 0.0 && qfn > 0.0) {
                        if(pfe + qfn > cvolrat) {
                            cf = qfn / (pfe + qfn);
                            pfe = (1.0 - cf) * cvolrat;
                            qfn = cf * cvolrat;
                        }
                    } else if(pfe > 0.0) {
                        pfe = std::fmin(pfe, (cvolrat - qfn));
                    } else if(qfn > 0.0) {
                        qfn = std::fmin(qfn, (cvolrat - pfe));
                    }
                } else {
                    if(pfe >= 0.0 && qfn < 0.0) {
                        cbilan = cvolrat - qfn;
                        if(cbilan > 0.0) {
                            pfe = std::fmin(pfe, cbilan);
                        } else {
                            pfe = 0.0;
                        }
                    } else if(pfe < 0.0 && qfn >= 0.0) {
                        cbilan = cvolrat - pfe;
                        if(cbilan > 0.0) {
                            qfn = std::fmin(qfn, cbilan);
                        } else {
                            qfn = 0.0;
                        }
                    } else if(pfe >= 0.0 && qfn >= 0.0) {
                        pfe = 0.0;
                        qfn = 0.0;
                    }
                }

                // CALCULATE NEW CONCENTRATION
                dc = (pfw - pfe + qfs - qfn) * dtfl / area;
                con = conc(ix, iy) + dc / hp;
                con = std::fmin(cmaxr(ix, iy), con);
                con = std::fmax(cminr(ix, iy), con);
                conc(ix, iy) = con;

                qfcds_local(ix) = qfn;
            }
        }
    }
}
