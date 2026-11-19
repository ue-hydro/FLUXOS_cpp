// Copyright 2026, Diogo Costa, diogo.costa@uevora.pt
// This file is part of the FLUXOS model.

// This program, FLUXOS, is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
#include <cmath>

#include "GlobVar.h"
#include "hydrodynamics_calc.h"
#include "solver_drydomain.h"
#include "solver_wetdomain.h"
#include "mpi_domain.h"

void hydrodynamics_calc(
    GlobVar& ds)
{
//-----------------------------------------------------------------------
// Solves shallow water equation for one time-step - Ritter solver
// Pressure term excluded from numerical flux
// Discretized as central difference.
// Adjustment of source term
// MUSCL-approach with limiter in roe dissipation
// fw1= mass flux per unit width
// fw2= momentum flux per unit width in x-direction
// fw3= momentum flux per unit width in y-direction
//
// MPI+OpenMP hybrid parallelization:
// - Domain is decomposed across MPI processes
// - Each process uses OpenMP for local parallelism
// - Ghost cells are exchanged between MPI neighbors
//-----------------------------------------------------------------------

    // Cache frequently used values for faster access
    const unsigned int NCOLS = ds.NCOLS;
    const unsigned int NROWS = ds.NROWS;
    const double hdry = ds.hdry;
    const double dxy = ds.dxy;
    const double dtl = ds.dtfl;

    // Get raw references to avoid repeated unique_ptr dereference
    arma::Mat<double>& z_ref = *ds.z;
    arma::Mat<double>& zb_ref = *ds.zb;
    arma::Mat<double>& h_ref = *ds.h;
    arma::Mat<double>& qx_ref = *ds.qx;
    arma::Mat<double>& qy_ref = *ds.qy;
    arma::Mat<double>& us_ref = *ds.us;
    arma::Mat<float>& ldry_ref = *ds.ldry;
    arma::Mat<double>& twetimetracer_ref = *ds.twetimetracer;
    arma::Mat<float>& innerNeumannBCWeir_ref = *ds.innerNeumannBCWeir;
    arma::Mat<double>& dh_ref = *ds.dh;
    arma::Mat<double>& dqx_ref = *ds.dqx;
    arma::Mat<double>& dqy_ref = *ds.dqy;
    arma::Mat<double>& qxf_ref = *ds.qxf;
    arma::Mat<double>& qyf_ref = *ds.qyf;
    arma::Mat<double>& fe_1_ref = *ds.fe_1;
    arma::Mat<double>& fe_2_ref = *ds.fe_2;
    arma::Mat<double>& fe_3_ref = *ds.fe_3;
    arma::Mat<double>& fn_1_ref = *ds.fn_1;
    arma::Mat<double>& fn_2_ref = *ds.fn_2;
    arma::Mat<double>& fn_3_ref = *ds.fn_3;

    // Pre-compute constants
    const double dtl_div_3600 = dtl / 3600.0;
    const double inv_dxy = 1.0 / dxy;
    const double manning_n = ds.manning_n;
    const double gn2 = ds.gacc * manning_n * manning_n;
    const double sqrt_g = std::sqrt(ds.gacc);
    const double qxf_scale = (manning_n > 0.0) ? 1.0 : dtl;

#ifdef USE_MPI
    // Exchange ghost cells for input fields before computation
    // This ensures each process has up-to-date boundary data from neighbors
    mpi_domain.exchange_ghost_cells(z_ref);
    mpi_domain.exchange_ghost_cells(zb_ref);
    mpi_domain.exchange_ghost_cells(qx_ref);
    mpi_domain.exchange_ghost_cells(qy_ref);
    mpi_domain.exchange_ghost_cells(h_ref);
    mpi_domain.exchange_ghost_cells(ldry_ref);
#endif

    // Begin OpenMP parallel region
    #pragma omp parallel
    {
        unsigned int irow, icol;
        double hp;
        float cell_neumann;

        // GET hp AND CHECK IF DRY OR WET
        #pragma omp for schedule(static) collapse(2)
        for(icol = 1; icol <= NCOLS; icol++)
        {
            for(irow = 1; irow <= NROWS; irow++)
            {
                hp = std::fmax(0.0, z_ref(irow, icol) - zb_ref(irow, icol));
                h_ref(irow, icol) = hp;

                if(hp <= hdry)
                {
                    qx_ref(irow, icol) = 0.0;
                    qy_ref(irow, icol) = 0.0;
                    us_ref(irow, icol) = 0.0;
                    ldry_ref(irow, icol) = 1.0f;
                }
                else
                {
                    ldry_ref(irow, icol) = 0.0f;
                    twetimetracer_ref(irow, icol) += dtl_div_3600;
                }
            }
        }

#ifdef USE_MPI
        // Single-threaded ghost exchange after ldry update
        #pragma omp single
        {
            mpi_domain.exchange_ghost_cells(ldry_ref);
            mpi_domain.exchange_ghost_cells(h_ref);
        }
#endif

        // CALL FLOW SOLVERS (compute mass and momentum fluxes)
        #pragma omp for schedule(static) collapse(2)
        for(icol = 1; icol <= NCOLS; icol++)
        {
            for(irow = 1; irow <= NROWS; irow++)
            {
                cell_neumann = innerNeumannBCWeir_ref(irow, icol);
                if(ldry_ref(irow, icol) == 1.0f && cell_neumann == 0.0f)
                {
                    solver_dry(ds, irow, icol);
                }
                else if(cell_neumann == 0.0f)
                {
                    solver_wet(ds, irow, icol);
                }
            }
        }

#ifdef USE_MPI
        // Exchange flux ghost cells before derivative computation
        #pragma omp single
        {
            mpi_domain.exchange_ghost_cells(fe_1_ref);
            mpi_domain.exchange_ghost_cells(fe_2_ref);
            mpi_domain.exchange_ghost_cells(fe_3_ref);
            mpi_domain.exchange_ghost_cells(fn_1_ref);
            mpi_domain.exchange_ghost_cells(fn_2_ref);
            mpi_domain.exchange_ghost_cells(fn_3_ref);
        }
#endif

        // CALCULATE TOTAL MASS AND MOMENTUM DERIVATIVE
        #pragma omp for schedule(static) collapse(2)
        for(icol = 1; icol <= NCOLS; icol++)
        {
            for(irow = 1; irow <= NROWS; irow++)
            {
                cell_neumann = innerNeumannBCWeir_ref(irow, icol);
                if(cell_neumann == 0.0f)
                {
                    // Use cached references and pre-computed inverse
                    const double fe1_w = fe_1_ref(irow - 1, icol);
                    const double fe1_p = fe_1_ref(irow, icol);
                    const double fn1_s = fn_1_ref(irow, icol - 1);
                    const double fn1_p = fn_1_ref(irow, icol);

                    dh_ref(irow, icol) = ((fe1_w - fe1_p) * inv_dxy + (fn1_s - fn1_p) * inv_dxy) * dtl;

                    const double fe2_w = fe_2_ref(irow - 1, icol);
                    const double fe2_p = fe_2_ref(irow, icol);
                    const double fn2_s = fn_2_ref(irow, icol - 1);
                    const double fn2_p = fn_2_ref(irow, icol);

                    dqx_ref(irow, icol) = ((fe2_w - fe2_p) * inv_dxy + (fn2_s - fn2_p) * inv_dxy) * dtl;

                    const double fe3_w = fe_3_ref(irow - 1, icol);
                    const double fe3_p = fe_3_ref(irow, icol);
                    const double fn3_s = fn_3_ref(irow, icol - 1);
                    const double fn3_p = fn_3_ref(irow, icol);

                    dqy_ref(irow, icol) = ((fe3_w - fe3_p) * inv_dxy + (fn3_s - fn3_p) * inv_dxy) * dtl;

                    // Hydrostatic-reconstruction bed-slope source, applied per
                    // cell (Audusse et al. 2004): 0.5 g (h_e*^2 - h_w*^2) with
                    // h_f* = max(0, z_p - max(zb_p, zb_f)). Keeps lake-at-rest
                    // balanced against steps/walls on every side.
                    {
                        const double zp  = z_ref(irow, icol);
                        const double zbp = zb_ref(irow, icol);
                        const double he_s = std::fmax(0.0, zp - std::fmax(zbp, zb_ref(irow + 1, icol)));
                        const double hw_s = std::fmax(0.0, zp - std::fmax(zbp, zb_ref(irow - 1, icol)));
                        const double hn_s = std::fmax(0.0, zp - std::fmax(zbp, zb_ref(irow, icol + 1)));
                        const double hs_s = std::fmax(0.0, zp - std::fmax(zbp, zb_ref(irow, icol - 1)));
                        dqx_ref(irow, icol) += 0.5 * ds.gacc * (he_s * he_s - hw_s * hw_s) * inv_dxy * dtl;
                        dqy_ref(irow, icol) += 0.5 * ds.gacc * (hn_s * hn_s - hs_s * hs_s) * inv_dxy * dtl;
                    }

                    // Face discharges for the momentum smoothing below. The legacy
                    // form stores fe1*dt (a volume, m2) and blends it into q (m2/s),
                    // which scales q by ~(0.8 + 0.2 dt) per step: damping for
                    // dt < 1 s, amplifying for dt > 1 s. Bed friction slows the flow
                    // so the CFL step grows past 1 s, so with MANNING_N the
                    // unit-consistent face discharge fe1 is used. Without it the
                    // legacy form is kept (results unchanged).
                    qxf_ref(irow, icol) = fe1_p * qxf_scale;
                    qyf_ref(irow, icol) = fn1_p * qxf_scale;
                }
                else
                {
                    dh_ref(irow, icol) = 0.0;
                    dqx_ref(irow, icol) = 0.0;
                    dqy_ref(irow, icol) = 0.0;
                    qxf_ref(irow, icol) = 0.0;
                    qyf_ref(irow, icol) = 0.0;
                }
            }
        }

#ifdef USE_MPI
        // Exchange qxf/qyf ghost cells for smoothing step
        #pragma omp single
        {
            mpi_domain.exchange_ghost_cells(qxf_ref);
            mpi_domain.exchange_ghost_cells(qyf_ref);
        }
#endif

        // CALCULATE NEW VALUES
        #pragma omp for schedule(static) collapse(2)
        for(icol = 1; icol <= NCOLS; icol++)
        {
            for(irow = 1; irow <= NROWS; irow++)
            {
                z_ref(irow, icol) = z_ref(irow, icol) + dh_ref(irow, icol);
                hp = std::fmax(0.0, z_ref(irow, icol) - zb_ref(irow, icol));
                h_ref(irow, icol) = hp;
                cell_neumann = innerNeumannBCWeir_ref(irow, icol);

                if(hp < hdry || cell_neumann == 1.0f)
                {
                    qx_ref(irow, icol) = 0.0;
                    qy_ref(irow, icol) = 0.0;
                    us_ref(irow, icol) = 0.0;
                    ldry_ref(irow, icol) = 1.0f;
                }
                else
                {
                    // Update momentum at cell center
                    qx_ref(irow, icol) = qx_ref(irow, icol) + dqx_ref(irow, icol);
                    qy_ref(irow, icol) = qy_ref(irow, icol) + dqy_ref(irow, icol);

                    // Smoothing with neighboring fluxes
                    const double qxf_w = qxf_ref(irow - 1, icol);
                    const double qxf_p = qxf_ref(irow, icol);
                    qx_ref(irow, icol) = 0.1 * qxf_w + 0.8 * qx_ref(irow, icol) + 0.1 * qxf_p;

                    const double qyf_s = qyf_ref(irow, icol - 1);
                    const double qyf_p = qyf_ref(irow, icol);
                    qy_ref(irow, icol) = 0.1 * qyf_s + 0.8 * qy_ref(irow, icol) + 0.1 * qyf_p;

                    // Manning bed friction (MANNING_N > 0 only), point-implicit:
                    //   dq/dt = -g n^2 |q| q / h^(7/3)   (= -g h S_f, S_f = n^2|u|u/h^(4/3))
                    //   q_new = q* / (1 + dt g n^2 |q*| / h^(7/3))
                    // Unconditionally stable, never reverses the flow. Also sets the
                    // shear velocity u* = sqrt(g) n |u| / h^(1/6), which feeds the
                    // existing eddy-viscosity term (cvdef u* h) in solver_wet.
                    if(manning_n > 0.0)
                    {
                        const double qxp = qx_ref(irow, icol);
                        const double qyp = qy_ref(irow, icol);
                        const double qmag = std::sqrt(qxp * qxp + qyp * qyp);
                        const double hf = std::fmax(hp, hdry);
                        const double cfric = gn2 * qmag / std::pow(hf, 7.0 / 3.0);
                        const double ffac = 1.0 / (1.0 + dtl * cfric);
                        qx_ref(irow, icol) = qxp * ffac;
                        qy_ref(irow, icol) = qyp * ffac;
                        us_ref(irow, icol) = sqrt_g * manning_n * (qmag * ffac / hf) / std::pow(hf, 1.0 / 6.0);
                    }

                    ldry_ref(irow, icol) = 0.0f;
                }
            }
        }

    } // end of OpenMP parallel region

#ifdef USE_MPI
    // Final ghost exchange for updated state variables
    mpi_domain.exchange_ghost_cells(z_ref);
    mpi_domain.exchange_ghost_cells(h_ref);
    mpi_domain.exchange_ghost_cells(qx_ref);
    mpi_domain.exchange_ghost_cells(qy_ref);
    mpi_domain.exchange_ghost_cells(ldry_ref);
#endif
}
