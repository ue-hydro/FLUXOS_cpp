// Copyright 2019, Diogo Costa
// CUDA device port of solver_dry (solver_drydomain.cpp)
// Column-major indexing: element(irow, icol) = ptr[irow + icol * MROWS]

#ifndef CUDA_SOLVER_DRY_CUH
#define CUDA_SOLVER_DRY_CUH

#ifdef USE_CUDA

#include <cuda_runtime.h>

__device__ void solver_dry_device(
    const double* __restrict__ d_zb,
    const double* __restrict__ d_z,
    double* __restrict__ d_qx,
    double* __restrict__ d_qy,
    const double* __restrict__ d_h,
    const float*  __restrict__ d_ldry,
    double* __restrict__ d_dh,
    double* __restrict__ d_dqx,
    double* __restrict__ d_dqy,
    double* __restrict__ d_qxf,
    double* __restrict__ d_qyf,
    double* __restrict__ d_fe_1,
    double* __restrict__ d_fe_2,
    double* __restrict__ d_fe_3,
    double* __restrict__ d_fn_1,
    double* __restrict__ d_fn_2,
    double* __restrict__ d_fn_3,
    unsigned int irow, unsigned int icol,
    unsigned int MROWS, unsigned int NROWS, unsigned int NCOLS,
    double gacc, double hdry, double dtfl, unsigned int dxy)
{
    #define IDX(r, c) ((r) + (c) * MROWS)

    // Neighbor indices
    const unsigned int is = icol - 1;
    const unsigned int in = icol + 1;
    const unsigned int iw = irow - 1;
    const unsigned int ie = irow + 1;

    const float ldw = d_ldry[IDX(iw, icol)];
    const float ldp = d_ldry[IDX(irow, icol)];
    const float lde = d_ldry[IDX(ie, icol)];
    const float lds = d_ldry[IDX(irow, is)];
    const float ldn = d_ldry[IDX(irow, in)];

    // All neighbours dry -> zero everything
    if (ldw == 1 && ldp == 1 && lde == 1 && lds == 1 && ldn == 1) {
        d_dh[IDX(irow, icol)] = 0.0;
        d_dqx[IDX(irow, icol)] = 0.0;
        d_dqy[IDX(irow, icol)] = 0.0;
        d_qxf[IDX(irow, icol)] = 0.0;
        d_qyf[IDX(irow, icol)] = 0.0;
        d_qx[IDX(irow, icol)] = 0.0;
        d_qy[IDX(irow, icol)] = 0.0;
        d_fn_1[IDX(irow, icol)] = 0.0;
        d_fn_2[IDX(irow, icol)] = 0.0;
        d_fn_3[IDX(irow, icol)] = 0.0;
        d_fe_1[IDX(irow, icol)] = 0.0;
        d_fe_2[IDX(irow, icol)] = 0.0;
        d_fe_3[IDX(irow, icol)] = 0.0;
        return;
    }

    // Cell center values
    const double zbp = d_zb[IDX(irow, icol)];
    const double zbe = d_zb[IDX(ie, icol)];
    const double zbn = d_zb[IDX(irow, in)];
    const double zp = d_z[IDX(irow, icol)];
    const double ze = d_z[IDX(ie, icol)];
    const double zn = d_z[IDX(irow, in)];
    const double hp = fmax(0.0, d_z[IDX(irow, icol)] - d_z[IDX(irow, icol)]);
    const double he = fmax(0.0, ze - zbe);
    const double hn = fmax(0.0, zn - zbn);
    const double qp = d_qx[IDX(irow, icol)];
    const double qe = d_qx[IDX(ie, icol)];
    const double rp = d_qy[IDX(irow, icol)];
    const double rn = d_qy[IDX(irow, in)];

    // Cell face values
    const double zbpe = 0.5 * (zbe + zbp);
    const double zbpn = 0.5 * (zbn + zbp);
    double hme = 0.5 * (hp + he);
    double qme = 0.5 * (qp + qe);
    double hmn = 0.5 * (hp + hn);
    double rmn = 0.5 * (rp + rn);
    double dze = ze - zp;
    double dzn = zn - zp;

    double fe1 = 0.0, fe2 = 0.0, fe3 = 0.0;
    double fn1 = 0.0, fn2 = 0.0, fn3 = 0.0;
    double volrat = 0.0;

    // East neighbor is wet
    if (lde == 0.0f) {
        hme = fmax(0.0, ze - zbpe);

        if (hme > hdry) {
            if (ze <= zp) {
                double fe2p = 0.5 * gacc * hme * hme;
                fe1 = 0.0;
                fe2 = fe2p;
                fe3 = 0.0;
            } else {
                dze = fmin(fabs(dze), hme);
                qme = 0.296 * dze * sqrt(gacc * dze);
                hme = 0.444 * dze;
                volrat = dxy * he / dtfl;
                qme = -fmin(qme, volrat);
                double ume = qme / hme;
                fe1 = qme;
                fe2 = qme * ume + 0.5 * gacc * hme * hme;
                fe3 = 0.0;
            }
        } else {
            fe1 = 0.0;
            fe2 = 0.0;
            fe3 = 0.0;
        }
    }

    // North neighbor is wet
    if (ldn == 0.0f) {
        hmn = fmax(0.0, zn - zbpn);
        if (hmn > hdry) {
            if (zn <= zp) {
                double fn3p = 0.5 * gacc * hmn * hmn;
                fn1 = 0.0;
                fn2 = 0.0;
                fn3 = fn3p;
            } else {
                dzn = fmin(fabs(dzn), hmn);
                rmn = 0.296 * dzn * sqrt(gacc * dzn);
                hmn = 0.444 * dzn;
                volrat = dxy * hn / dtfl;
                rmn = -fmin(rmn, volrat);
                double vmn = rmn / hmn;
                fn1 = rmn;
                fn2 = 0.0;
                fn3 = rmn * vmn + 0.5 * gacc * hmn * hmn;
            }
        } else {
            fn1 = 0.0;
            fn2 = 0.0;
            fn3 = 0.0;
        }
    }

    // Boundary conditions (weir discharge rate)
    if (icol == 1 || icol == NCOLS) {
        fn1 = fmin(volrat, sqrt(gacc) * pow(fmax(hp, 0.0), 1.5));
    }
    if (irow == 1 || irow == NROWS) {
        fe1 = fmin(volrat, sqrt(gacc) * pow(fmax(hp, 0.0), 1.5));
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
#endif // CUDA_SOLVER_DRY_CUH
