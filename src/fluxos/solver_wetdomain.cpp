

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
#include "solver_wetdomain.h"

void solver_wet(
    GlobVar& ds,
    unsigned int irow,
    unsigned int icol) {

    // Cache global constants
    const unsigned int NROWSl = ds.NROWS;
    const unsigned int NCOLSl = ds.NCOLS;
    const double hdryl = ds.hdry;
    const double gaccl = ds.gacc;
    const double nueml = ds.nuem;
    const double cvdefl = ds.cvdef;
    const double dtl = ds.dtfl;
    const unsigned int dx = ds.dxy;
    const unsigned int dy = ds.dxy;
    const double arbase = ds.arbase;

    // Pre-compute neighbor indices
    const unsigned int is = icol - 1;
    const unsigned int in = icol + 1;
    const unsigned int inn = std::min(icol + 2, NCOLSl + 1);
    const unsigned int iw = irow - 1;
    const unsigned int ie = irow + 1;

    // Get raw references for faster access
    const arma::Mat<double>& zb_ref = *ds.zb;
    const arma::Mat<double>& z_ref = *ds.z;
    const arma::Mat<double>& qx_ref = *ds.qx;
    const arma::Mat<double>& qy_ref = *ds.qy;
    const arma::Mat<double>& h_ref = *ds.h;
    const arma::Mat<float>& ldry_ref = *ds.ldry;
    const arma::Mat<double>& us_ref = *ds.us;
    const arma::Mat<double>& ks_ref = *ds.ks;
    arma::Mat<double>& fn_1_ref = *ds.fn_1;
    arma::Mat<double>& fn_2_ref = *ds.fn_2;
    arma::Mat<double>& fn_3_ref = *ds.fn_3;
    arma::Mat<double>& fe_1_ref = *ds.fe_1;
    arma::Mat<double>& fe_2_ref = *ds.fe_2;
    arma::Mat<double>& fe_3_ref = *ds.fe_3;

    const double kspl = ks_ref(irow, icol);
    const float ldp = ldry_ref(irow, icol);
    const float lde = ldry_ref(ie, icol);
    const float ldn = ldry_ref(irow, in);

    // CELL CENTER VALUES - batch read for cache efficiency
    const double zbw = zb_ref(iw, icol);
    const double zbp = zb_ref(irow, icol);
    const double zbe = zb_ref(ie, icol);
    const double zbs = zb_ref(irow, is);
    const double zbn = zb_ref(irow, in);
    const double zbnn = zb_ref(irow, inn);
    const double zw = z_ref(iw, icol);
    const double zp = z_ref(irow, icol);
    const double ze = z_ref(ie, icol);
    const double zs = z_ref(irow, is);
    double zn = z_ref(irow, in);
    const double znn = z_ref(irow, inn);
    double qp = qx_ref(irow, icol);
    double qe = qx_ref(ie, icol);
    double qn = qx_ref(irow, in);
    const double rw = qy_ref(iw, icol);
    double rp = qy_ref(irow, icol);
    double re = qy_ref(ie, icol);
    double rn = qy_ref(irow, in);

    // Hydrostatic reconstruction: use fmax of bed elevations at face
    // (Audusse et al. 2004 — preserves lake-at-rest / C-property)
    const double zbpe = std::fmax(zbe, zbp);
    const double zbpn = std::fmax(zbn, zbp);

    // Cell-center water depths (original, before reconstruction)
    double hp = std::fmax(0.0, zp - zbp);
    const double hp0 = std::fmax(std::fmax(hdryl, hp), kspl);
    const double hw = std::fmax(0.0, zw - zbw);
    double he = std::fmax(0.0, ze - zbe);
    const double hs = std::fmax(0.0, zs - zbs);
    double hn = std::fmax(0.0, zn - zbn);
    double hnn_orig = std::fmax(0.0, znn - zbnn);

    // Reconstructed face depths (water depth as seen from each side of the face)
    const double hp_e = std::fmax(0.0, zp - zbpe);   // P side of east face
    const double he_e = std::fmax(0.0, ze - zbpe);   // E side of east face
    const double hp_n = std::fmax(0.0, zp - zbpn);   // P side of north face
    const double hn_n = std::fmax(0.0, zn - zbpn);   // N side of north face

    // CELL FACE VALUES (use reconstructed depths)
    double hme = 0.5 * (hp_e + he_e);
    double qme = 0.5 * (qp + qe);
    double rme = 0.5 * (rp + re);
    double hmn = 0.5 * (hp_n + hn_n);
    double qmn = 0.5 * (qp + qn);
    double rmn = 0.5 * (rp + rn);

    // Hydrostatic pressure source correction for cell P
    // Mismatch between original cell depth and reconstructed face depth
    // S = 0.5 * g * (h_orig^2 - h_recon^2)
    const double src_e = 0.5 * gaccl * (hp * hp - hp_e * hp_e);
    const double src_n = 0.5 * gaccl * (hp * hp - hp_n * hp_n);

    double dze = ze - zp;
    double dqe = qe - qp;
    double dre = re - rp;
    double dzn = zn - zp;
    double drn = rn - rp;
    double dqn = qn - qp;

    bool lroe = true;
    double volrat;

    // CELLS WITH SOME DRY NEIGHBOURS
    if(lde == 1.0f) {
        hme = std::fmax(0.0, zp - zbpe);
        he = 0.0;
        qe = 0.0;
        re = 0.0;
        if(hme > hdryl) {
            if(ze >= zp) {
                qme = 0.0;
                rme = 0.0;
            } else {
                dze = fmin(fabs(dze), hme);
                qme = 0.296 * dze * sqrt(gaccl * dze);  // Ritter solution
                volrat = 0.5 * dx * hp / dtl;           // available volume rate per m of cell P
                qme = fmin(qme, volrat);
                rme = 0.0;
            }
        } else {
            qme = 0.0;
            rme = 0.0;
            hme = 0.0;
        }
        lroe = false;
    }

    if(ldn == 1.0f) {
        hmn = std::fmax(0.0, zp - zbpn);
        hn = 0.0;
        qn = 0.0;
        rn = 0.0;
        if(hmn > hdryl) {
            if(zn >= zp) {
                qmn = 0.0;
                rmn = 0.0;
            } else {
                dzn = fmin(fabs(dzn), hmn);
                rmn = 0.296 * dzn * sqrt(gaccl * dzn);  // Ritter solution
                qmn = 0.0;
                volrat = 0.5 * dy * hp / dtl;           // available volume rate per m of cell P
                rmn = fmin(rmn, volrat);
            }
        } else {
            qmn = 0.0;
            rmn = 0.0;
            hmn = 0.0;
        }
        lroe = false;
    }

    // CALC TURBULENT STRESS
    const double cme = sqrt(gaccl * hme);
    const double cmn = sqrt(gaccl * hmn);
    const double ume = qme / std::fmax(hme, hdryl);
    const double vme = rme / std::fmax(hme, hdryl);
    const double umn = qmn / std::fmax(hmn, hdryl);
    const double vmn = rmn / std::fmax(hmn, hdryl);

    const double cnp = cvdefl * us_ref(irow, icol) * hp + nueml;
    const double cne = cvdefl * us_ref(ie, icol) * he + nueml;
    const double cnn = cvdefl * us_ref(irow, in) * hn + nueml;
    const double hne = 0.5 * (cnp + cne) * sqrt(hp * he);
    double hnn = 0.5 * (cnp + cnn) * sqrt(hp * hn);

    const double up = qp / hp0;
    const double un = qn / std::fmax(std::fmax(hn, hdryl), ks_ref(irow, in));
    const double us0 = qx_ref(irow, is) / std::fmax(std::fmax(hs, hdryl), ks_ref(irow, is));
    const double use0 = qx_ref(ie, is) / std::fmax(std::fmax(h_ref(ie, is), hdryl), ks_ref(ie, is));
    const double une = qx_ref(ie, in) / std::fmax(std::fmax(h_ref(ie, in), hdryl), ks_ref(ie, in));
    const double vp = rp / hp0;
    const double ve = re / std::fmax(std::fmax(he, hdryl), ks_ref(ie, icol));
    const double vw = rw / std::fmax(std::fmax(hw, hdryl), ks_ref(iw, icol));
    const double vwn = qy_ref(iw, in) / std::fmax(std::fmax(h_ref(iw, in), hdryl), ks_ref(iw, in));
    const double ven = qy_ref(ie, in) / std::fmax(std::fmax(h_ref(ie, in), hdryl), ks_ref(ie, in));

    const double txye = hne * ((ve - vp) / fabs(dx) + 0.25 * (un + une - us0 - use0) / dy);
    const double txyn = hnn * ((un - up) / fabs(dy) + 0.25 * (ve + ven - vw - vwn) / dx);

    // CALC OF CONVECTION FLUXES
    double fe1c = qme;
    double fe2c = qme * ume;
    double fe3c = qme * vme - txye;
    double fn1c = rmn;
    double fn2c = rmn * umn - txyn;
    double fn3c = rmn * vmn;

    // ROE's DISSIPATION
    double fe1r = 0.0, fe2r = 0.0, fe3r = 0.0;
    double fn1r = 0.0, fn2r = 0.0, fn3r = 0.0;

    if(lroe) {
        if(hme > hdryl) {
            double dzea = fabs(dze);
            if(dzea > 0.5 * hme) {
                double dhea = fabs(he - hp);
                dzea = fmin(dzea, dhea);
                dze = std::copysign(dzea, dze);
            }
            double cc = 0.25 / cme;
            double c1 = ume;
            double c2 = ume + cme;
            double c3 = ume - cme;
            double c1a = fabs(c1);
            double c2a = fabs(c2);
            double c3a = fabs(c3);
            double a11 = c2 * c3a - c2a * c3;
            double a12 = c2a - c3a;
            double a21 = c2 * c3 * (c3a - c2a);
            double a22 = c2a * c2 - c3a * c3;
            double a31 = vme * (c2 * c3a - 2.0 * cme * c1a - c2a * c3);
            double a32 = vme * (c2a - c3a);
            double a33 = 2.0 * cme * c1a;

            fe1r = -(a11 * dze + a12 * dqe) * cc;
            fe2r = -(a21 * dze + a22 * dqe) * cc;
            fe3r = -(a31 * dze + a32 * dqe + a33 * dre) * cc;
        }

        if(ldp == 0.0f && hmn > hdryl) {
            double dzna = fabs(dzn);
            double dhna = fabs(hn - hp);
            dzna = fmin(dzna, dhna);
            dzn = std::copysign(dzna, dzn);
            double cc = 0.25 / cmn;
            double c1 = vmn;
            double c2 = vmn + cmn;
            double c3 = vmn - cmn;
            double c1a = std::fabs(c1);
            double c2a = std::fabs(c2);
            double c3a = std::fabs(c3);
            double a11 = c2 * c3a - c2a * c3;
            double a13 = c2a - c3a;
            double a21 = umn * (c2 * c3a - 2.0 * cmn * c1a - c2a * c3);
            double a22 = 2.0 * cmn * c1a;
            double a23 = umn * (c2a - c3a);
            double a31 = c2 * c3 * (c3a - c2a);
            double a33 = c2a * c2 - c3a * c3;

            fn1r = -(a11 * dzn + a13 * drn) * cc;
            fn2r = -(a21 * dzn + a22 * dqn + a23 * drn) * cc;
            fn3r = -(a31 * dzn + a33 * drn) * cc;
        }
    }

    // PRESSURE AT CELL SIDE
    const double fe2p = 0.5 * gaccl * (hme * hme);
    const double fn3p = 0.5 * gaccl * (hmn * hmn);

    // SUM OF ALL FLUXES (include hydrostatic pressure source correction)
    // The src_e/src_n terms absorb the bed-slope source from hydrostatic
    // reconstruction into the momentum flux (Audusse et al. 2004)
    double fe1 = fe1c + fe1r;
    double fe2 = fe2c + fe2r + fe2p - src_e;
    double fe3 = fe3c + fe3r;
    double fn1 = fn1c + fn1r;
    double fn2 = fn2c + fn2r;
    double fn3 = fn3c + fn3r + fn3p - src_n;

    // BOUNDARY CONDITIONS (WEIR DISCHARGE RATE)
    if(icol == 1 || icol == NCOLSl) {
        fn1 = std::fmin(volrat, sqrt(gaccl) * pow(std::fmax(hp, 0.0), 1.5));
    }
    if(irow == 1 || irow == NROWSl) {
        fe1 = std::fmin(volrat, sqrt(gaccl) * pow(std::fmax(hp, 0.0), 1.5));
    }

    // CHECK MASS BALANCE (restrict outflow flux to available water)
    double volpot = arbase * hp;    // volume in cell P [m3]
    volrat = volpot / dtl;          // max flux rate

    if(volrat > 0.0) {  // cell has water
        if(fe1 > 0.0 && fn1 > 0.0) {
            if(fe1 * dy + fn1 * dx > volrat) {
                double cf = fn1 * dx / (fe1 * dy + fn1 * dx);
                fe1 = (1.0 - cf) * volrat / dy;
                fn1 = cf * volrat / dx;
            }
        } else if(fe1 > 0.0) {
            fe1 = fmin(fe1 * dy, volrat) / dy;
        } else if(fn1 > 0.0) {
            fn1 = fmin(fn1 * dx, volrat) / dx;
        }
    } else {  // cell has no water
        fe1 = 0.0;
        fn1 = 0.0;
    }

    // SAVE MASS AND MOMENTUM FLUXES
    fn_1_ref(irow, icol) = fn1;
    fn_2_ref(irow, icol) = fn2;
    fn_3_ref(irow, icol) = fn3;
    fe_1_ref(irow, icol) = fe1;
    fe_2_ref(irow, icol) = fe2;
    fe_3_ref(irow, icol) = fe3;
}
