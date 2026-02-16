

// Copyright 1992: Cornel Beffa and Roland Faeh
// Copyright 2013: Kashif Shaad and Diogo Costa
// Copyright 2019: Diogo Costa

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
#include <cmath>

#include "GlobVar.h"
#include "solver_drydomain.h"

void solver_dry(
    GlobVar& ds,
    unsigned int irow,
    unsigned int icol) {

    // Cache global constants
    const unsigned int NROWSl = ds.NROWS;
    const unsigned int NCOLSl = ds.NCOLS;
    const double gaccl = ds.gacc;
    const double hdryl = ds.hdry;
    const double dtl = ds.dtfl;
    const double dxy = ds.dxy;

    // Pre-compute neighbor indices
    const unsigned int is = icol - 1;
    const unsigned int in = icol + 1;
    const unsigned int iw = irow - 1;
    const unsigned int ie = irow + 1;

    // Get raw references for faster access
    const arma::Mat<float>& ldry_ref = *ds.ldry;
    const arma::Mat<double>& zb_ref = *ds.zb;
    const arma::Mat<double>& z_ref = *ds.z;
    arma::Mat<double>& qx_ref = *ds.qx;
    arma::Mat<double>& qy_ref = *ds.qy;
    arma::Mat<double>& dh_ref = *ds.dh;
    arma::Mat<double>& dqx_ref = *ds.dqx;
    arma::Mat<double>& dqy_ref = *ds.dqy;
    arma::Mat<double>& qxf_ref = *ds.qxf;
    arma::Mat<double>& qyf_ref = *ds.qyf;
    arma::Mat<double>& fn_1_ref = *ds.fn_1;
    arma::Mat<double>& fn_2_ref = *ds.fn_2;
    arma::Mat<double>& fn_3_ref = *ds.fn_3;
    arma::Mat<double>& fe_1_ref = *ds.fe_1;
    arma::Mat<double>& fe_2_ref = *ds.fe_2;
    arma::Mat<double>& fe_3_ref = *ds.fe_3;

    const float ldw = ldry_ref(iw, icol);
    const float ldp = ldry_ref(irow, icol);
    const float lde = ldry_ref(ie, icol);
    const float lds = ldry_ref(irow, is);
    const float ldn = ldry_ref(irow, in);

    // CHECK IF ALL NEIGHBOUR CELLS ARE DRY
    if(ldw == 1 && ldp == 1 && lde == 1 && lds == 1 && ldn == 1)
    {
        dh_ref(irow, icol) = 0.0;
        dqx_ref(irow, icol) = 0.0;
        dqy_ref(irow, icol) = 0.0;
        qxf_ref(irow, icol) = 0.0;
        qyf_ref(irow, icol) = 0.0;
        qx_ref(irow, icol) = 0.0;
        qy_ref(irow, icol) = 0.0;
        fn_1_ref(irow, icol) = 0.0;
        fn_2_ref(irow, icol) = 0.0;
        fn_3_ref(irow, icol) = 0.0;
        fe_1_ref(irow, icol) = 0.0;
        fe_2_ref(irow, icol) = 0.0;
        fe_3_ref(irow, icol) = 0.0;
        return;
    }

    // CELL CENTER VALUES
    const double zbp = zb_ref(irow, icol);
    const double zbe = zb_ref(ie, icol);
    const double zbn = zb_ref(irow, in);
    const double zp = z_ref(irow, icol);
    const double ze = z_ref(ie, icol);
    const double zn = z_ref(irow, in);
    const double hp = std::fmax(0.0, zp - zbp);
    const double he = std::fmax(0.0, ze - zbe);
    const double hn = std::fmax(0.0, zn - zbn);
    const double qp = qx_ref(irow, icol);
    const double qe = qx_ref(ie, icol);
    const double rp = qy_ref(irow, icol);
    const double rn = qy_ref(irow, in);

    // CELL FACE VALUES — hydrostatic reconstruction (Audusse et al. 2004)
    const double zbpe = std::fmax(zbe, zbp);
    const double zbpn = std::fmax(zbn, zbp);
    double hme = 0.5 * (hp + he);
    double qme = 0.5 * (qp + qe);
    double hmn = 0.5 * (hp + hn);
    double rmn = 0.5 * (rp + rn);
    double dze = ze - zp;
    double dzn = zn - zp;

    // INITIATION
    double fe1 = 0.0, fe2 = 0.0, fe3 = 0.0;
    double fn1 = 0.0, fn2 = 0.0, fn3 = 0.0;
    double volrat = 0.0;

    // CELLS WITH SOME DRY NEIGHBOURS
    if(lde == 0.0f)
    {
        hme = std::fmax(0.0, ze - zbpe);

        if(hme > hdryl)
        {
            if(ze <= zp)
            {
                double fe2p = 0.5 * gaccl * hme * hme;
                fe1 = 0.0;
                fe2 = fe2p;
                fe3 = 0.0;
            }
            else
            {
                dze = std::fmin(std::fabs(dze), hme);
                qme = 0.296 * dze * sqrt(gaccl * dze);      // Ritter solution
                hme = 0.444 * dze;                          // at cell side
                volrat = dxy * he / dtl;                    // available volume rate per m of cell E
                qme = -std::fmin(qme, volrat);
                double ume = qme / hme;                     // from cell center
                fe1 = qme;
                fe2 = qme * ume + 0.5 * gaccl * hme * hme;
                fe3 = 0.0;
            }
        }
        else
        {
            fe1 = 0.0;
            fe2 = 0.0;
            fe3 = 0.0;
        }
    }

    if(ldn == 0.0f)
    {
        hmn = std::fmax(0.0, zn - zbpn);
        if(hmn > hdryl)
        {
            if(zn <= zp)
            {
                double fn3p = 0.5 * gaccl * hmn * hmn;
                fn1 = 0.0;
                fn2 = 0.0;
                fn3 = fn3p;
            }
            else
            {
                dzn = std::fmin(std::fabs(dzn), hmn);
                rmn = 0.296 * dzn * sqrt(gaccl * dzn);        // Ritter solution
                hmn = 0.444 * dzn;                            // at cell side
                volrat = dxy * hn / dtl;                      // available volume rate per m of cell N
                rmn = -std::fmin(rmn, volrat);
                double vmn = rmn / hmn;
                fn1 = rmn;
                fn2 = 0.0;
                fn3 = rmn * vmn + 0.5 * gaccl * hmn * hmn;
            }
        }
        else
        {
            fn1 = 0.0;
            fn2 = 0.0;
            fn3 = 0.0;
        }
    }

    // BOUNDARY CONDITIONS (WEIR DISCHARGE RATE)
    if(icol == 1 || icol == NCOLSl)
    {
        fn1 = std::fmin(volrat, sqrt(gaccl) * pow(std::fmax(hp, 0.0), 1.5));
    }
    if(irow == 1 || irow == NROWSl)
    {
        fe1 = std::fmin(volrat, sqrt(gaccl) * pow(std::fmax(hp, 0.0), 1.5));
    }

    // SAVE MASS AND MOMENTUM FLUXES
    fn_1_ref(irow, icol) = fn1;
    fn_2_ref(irow, icol) = fn2;
    fn_3_ref(irow, icol) = fn3;
    fe_1_ref(irow, icol) = fe1;
    fe_2_ref(irow, icol) = fe2;
    fe_3_ref(irow, icol) = fe3;
}
