// Copyright 2019, Diogo Costa
// CUDA device port of solver_wet (solver_wetdomain.cpp)
// Column-major indexing: element(irow, icol) = ptr[irow + icol * MROWS]

#ifndef CUDA_SOLVER_WET_CUH
#define CUDA_SOLVER_WET_CUH

#ifdef USE_CUDA

#include <cuda_runtime.h>

__device__ inline double d_fmax(double a, double b) { return fmax(a, b); }
__device__ inline double d_fmin(double a, double b) { return fmin(a, b); }

__device__ void solver_wet_device(
    const double* __restrict__ d_zb,
    const double* __restrict__ d_z,
    const double* __restrict__ d_qx,
    const double* __restrict__ d_qy,
    const double* __restrict__ d_h,
    const float*  __restrict__ d_ldry,
    const double* __restrict__ d_us,
    const double* __restrict__ d_ks,
    double* __restrict__ d_fe_1,
    double* __restrict__ d_fe_2,
    double* __restrict__ d_fe_3,
    double* __restrict__ d_fn_1,
    double* __restrict__ d_fn_2,
    double* __restrict__ d_fn_3,
    unsigned int irow, unsigned int icol,
    unsigned int MROWS, unsigned int NROWS, unsigned int NCOLS,
    double hdry, double gacc, double nuem, double cvdef,
    double dtfl, unsigned int dxy, double arbase)
{
    // Column-major index macro
    #define IDX(r, c) ((r) + (c) * MROWS)

    const unsigned int dx = dxy;
    const unsigned int dy = dxy;

    // Neighbor indices
    const unsigned int is = icol - 1;
    const unsigned int in = icol + 1;
    const unsigned int inn = (icol + 2 <= NCOLS + 1) ? icol + 2 : NCOLS + 1;
    const unsigned int iw = irow - 1;
    const unsigned int ie = irow + 1;

    const double kspl = d_ks[IDX(irow, icol)];
    const float ldp = d_ldry[IDX(irow, icol)];
    const float lde = d_ldry[IDX(ie, icol)];
    const float ldn = d_ldry[IDX(irow, in)];

    // Cell center values
    const double zbw = d_zb[IDX(iw, icol)];
    const double zbp = d_zb[IDX(irow, icol)];
    const double zbe = d_zb[IDX(ie, icol)];
    const double zbs = d_zb[IDX(irow, is)];
    const double zbn = d_zb[IDX(irow, in)];
    const double zbnn = d_zb[IDX(irow, inn)];
    const double zw = d_z[IDX(iw, icol)];
    const double zp = d_z[IDX(irow, icol)];
    const double ze = d_z[IDX(ie, icol)];
    const double zs = d_z[IDX(irow, is)];
    double zn = d_z[IDX(irow, in)];
    const double znn = d_z[IDX(irow, inn)];
    double qp = d_qx[IDX(irow, icol)];
    double qe = d_qx[IDX(ie, icol)];
    double qn = d_qx[IDX(irow, in)];
    const double rw = d_qy[IDX(iw, icol)];
    double rp = d_qy[IDX(irow, icol)];
    double re = d_qy[IDX(ie, icol)];
    double rn = d_qy[IDX(irow, in)];

    const double zbpe = 0.5 * (zbe + zbp);
    const double zbpn = 0.5 * (zbn + zbp);
    double hp = d_fmax(0.0, zp - zbp);
    const double hp0 = d_fmax(d_fmax(hdry, hp), kspl);
    const double hw = d_fmax(0.0, zw - zbw);
    double he = d_fmax(0.0, ze - zbe);
    const double hs = d_fmax(0.0, zs - zbs);
    double hn = d_fmax(0.0, zn - zbn);
    double hnn_orig = d_fmax(0.0, znn - zbnn);

    // Cell face values
    double hme = 0.5 * (hp + he);
    double qme = 0.5 * (qp + qe);
    double rme = 0.5 * (rp + re);
    double hmn = 0.5 * (hp + hn);
    double qmn = 0.5 * (qp + qn);
    double rmn = 0.5 * (rp + rn);

    double dze = ze - zp;
    double dqe = qe - qp;
    double dre = re - rp;
    double dzn = zn - zp;
    double drn = rn - rp;
    double dqn = qn - qp;

    bool lroe = true;
    double volrat;

    // Cells with some dry neighbours (east)
    if (lde == 1.0f) {
        hme = d_fmax(0.0, zp - zbpe);
        he = 0.0;
        qe = 0.0;
        re = 0.0;
        if (hme > hdry) {
            if (ze >= zp) {
                qme = 0.0;
                rme = 0.0;
            } else {
                dze = d_fmin(fabs(dze), hme);
                qme = 0.296 * dze * sqrt(gacc * dze);
                volrat = 0.5 * dx * hp / dtfl;
                qme = d_fmin(qme, volrat);
                rme = 0.0;
            }
        } else {
            qme = 0.0;
            rme = 0.0;
            hme = 0.0;
        }
        lroe = false;
    }

    // Cells with some dry neighbours (north)
    if (ldn == 1.0f) {
        hmn = d_fmax(0.0, zp - zbpn);
        hn = 0.0;
        qn = 0.0;
        rn = 0.0;
        if (hmn > hdry) {
            if (zn >= zp) {
                qmn = 0.0;
                rmn = 0.0;
            } else {
                dzn = d_fmin(fabs(dzn), hmn);
                rmn = 0.296 * dzn * sqrt(gacc * dzn);
                qmn = 0.0;
                volrat = 0.5 * dy * hp / dtfl;
                rmn = d_fmin(rmn, volrat);
            }
        } else {
            qmn = 0.0;
            rmn = 0.0;
            hmn = 0.0;
        }
        lroe = false;
    }

    // Turbulent stress
    const double cme = sqrt(gacc * hme);
    const double cmn = sqrt(gacc * hmn);
    const double ume = qme / d_fmax(hme, hdry);
    const double vme = rme / d_fmax(hme, hdry);
    const double umn = qmn / d_fmax(hmn, hdry);
    const double vmn = rmn / d_fmax(hmn, hdry);

    const double cnp = cvdef * d_us[IDX(irow, icol)] * hp + nuem;
    const double cne = cvdef * d_us[IDX(ie, icol)] * he + nuem;
    const double cnn = cvdef * d_us[IDX(irow, in)] * hn + nuem;
    const double hne = 0.5 * (cnp + cne) * sqrt(hp * he);
    double hnn = 0.5 * (cnp + cnn) * sqrt(hp * hn);

    const double up = qp / hp0;
    const double un = qn / d_fmax(d_fmax(hn, hdry), d_ks[IDX(irow, in)]);
    const double us0 = d_qx[IDX(irow, is)] / d_fmax(d_fmax(hs, hdry), d_ks[IDX(irow, is)]);
    const double use0 = d_qx[IDX(ie, is)] / d_fmax(d_fmax(d_h[IDX(ie, is)], hdry), d_ks[IDX(ie, is)]);
    const double une = d_qx[IDX(ie, in)] / d_fmax(d_fmax(d_h[IDX(ie, in)], hdry), d_ks[IDX(ie, in)]);
    const double vp = rp / hp0;
    const double ve = re / d_fmax(d_fmax(he, hdry), d_ks[IDX(ie, icol)]);
    const double vw = rw / d_fmax(d_fmax(hw, hdry), d_ks[IDX(iw, icol)]);
    const double vwn = d_qy[IDX(iw, in)] / d_fmax(d_fmax(d_h[IDX(iw, in)], hdry), d_ks[IDX(iw, in)]);
    const double ven = d_qy[IDX(ie, in)] / d_fmax(d_fmax(d_h[IDX(ie, in)], hdry), d_ks[IDX(ie, in)]);

    const double txye = hne * ((ve - vp) / fabs((double)dx) + 0.25 * (un + une - us0 - use0) / dy);
    const double txyn = hnn * ((un - up) / fabs((double)dy) + 0.25 * (ve + ven - vw - vwn) / dx);

    // Convection fluxes
    double fe1c = qme;
    double fe2c = qme * ume;
    double fe3c = qme * vme - txye;
    double fn1c = rmn;
    double fn2c = rmn * umn - txyn;
    double fn3c = rmn * vmn;

    // Roe's dissipation
    double fe1r = 0.0, fe2r = 0.0, fe3r = 0.0;
    double fn1r = 0.0, fn2r = 0.0, fn3r = 0.0;

    if (lroe) {
        if (hme > hdry) {
            double dzea = fabs(dze);
            if (dzea > 0.5 * hme) {
                double dhea = fabs(he - hp);
                dzea = d_fmin(dzea, dhea);
                dze = copysign(dzea, dze);
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

        if (ldp == 0.0f && hmn > hdry) {
            double dzna = fabs(dzn);
            double dhna = fabs(hn - hp);
            dzna = d_fmin(dzna, dhna);
            dzn = copysign(dzna, dzn);
            double cc = 0.25 / cmn;
            double c1 = vmn;
            double c2 = vmn + cmn;
            double c3 = vmn - cmn;
            double c1a = fabs(c1);
            double c2a = fabs(c2);
            double c3a = fabs(c3);
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

    // Pressure at cell side
    const double fe2p = 0.5 * gacc * (hme * hme);
    const double fn3p = 0.5 * gacc * (hmn * hmn);

    // Sum of all fluxes
    double fe1 = fe1c + fe1r;
    double fe2 = fe2c + fe2r + fe2p;
    double fe3 = fe3c + fe3r;
    double fn1 = fn1c + fn1r;
    double fn2 = fn2c + fn2r;
    double fn3 = fn3c + fn3r + fn3p;

    // Boundary conditions (weir discharge rate)
    if (icol == 1 || icol == NCOLS) {
        fn1 = d_fmin(volrat, sqrt(gacc) * pow(d_fmax(hp, 0.0), 1.5));
    }
    if (irow == 1 || irow == NROWS) {
        fe1 = d_fmin(volrat, sqrt(gacc) * pow(d_fmax(hp, 0.0), 1.5));
    }

    // Mass balance check
    double volpot = arbase * hp;
    volrat = volpot / dtfl;

    if (volrat > 0.0) {
        if (fe1 > 0.0 && fn1 > 0.0) {
            if (fe1 * dy + fn1 * dx > volrat) {
                double cf = fn1 * dx / (fe1 * dy + fn1 * dx);
                fe1 = (1.0 - cf) * volrat / dy;
                fn1 = cf * volrat / dx;
            }
        } else if (fe1 > 0.0) {
            fe1 = d_fmin(fe1 * dy, volrat) / dy;
        } else if (fn1 > 0.0) {
            fn1 = d_fmin(fn1 * dx, volrat) / dx;
        }
    } else {
        fe1 = 0.0;
        fn1 = 0.0;
    }

    // Save fluxes
    d_fn_1[IDX(irow, icol)] = fn1;
    d_fn_2[IDX(irow, icol)] = fn2;
    d_fn_3[IDX(irow, icol)] = fn3;
    d_fe_1[IDX(irow, icol)] = fe1;
    d_fe_2[IDX(irow, icol)] = fe2;
    d_fe_3[IDX(irow, icol)] = fe3;

    #undef IDX
}

#endif // USE_CUDA
#endif // CUDA_SOLVER_WET_CUH
