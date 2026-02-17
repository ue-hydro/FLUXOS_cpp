
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
#include "mpi_domain.h"

// ADE solver — explicit, fully-parallel, mass-conservative formulation.
//
// Two-pass approach (no sequential sweep, no race conditions):
//
//   Pass 1:  Compute face fluxes at ALL east and north faces.
//            Each face flux depends only on the OLD concentration field
//            (read-only), so this is fully parallel.
//
//   Pass 2:  For each cell, gather the 4 face fluxes, compute net mass
//            change, convert to concentration, and apply monotonicity clamp.
//
// Advection:  first-order upwind  (unconditionally TVD)
// Diffusion:  Fickian,  flux = −D · h_face · edge_len · (c_R − c_L) / dx
void adesolver_calc(
    GlobVar& ds,
    int it,
    int ichem)
{
    const double dx    = ds.dxy;
    const double dy    = ds.dxy;
    const double dtfl  = ds.dtfl;
    const double D_coef = ds.D_coef;
    const double area  = ds.arbase;          // dx * dy
    const unsigned int NROWS = ds.NROWS;
    const unsigned int NCOLS = ds.NCOLS;

    // Minimum depth for ADE transport.  Cells thinner than this
    // are treated as "too shallow for meaningful solute transport"
    // and have their concentration set to zero.  This prevents the
    // c = M/h blowup at the wetting front where h → 0.
    const double h_min = 0.001;              // 1 mm

    arma::Mat<double>& conc     = (*ds.conc_SW)[ichem];
    const arma::Mat<double>& h  = *ds.h;     // depth AFTER hydro update
    const arma::Mat<double>& h0 = *ds.h0;    // depth BEFORE hydro update
    const arma::Mat<float>&  dry = *ds.ldry;
    const arma::Mat<double>& qx = *ds.fe_1;  // discharge per unit width x (m²/s)
    const arma::Mat<double>& qy = *ds.fn_1;  // discharge per unit width y (m²/s)

#ifdef USE_MPI
    mpi_domain.exchange_ghost_cells(conc);
#endif

    if (it <= 1) return;

    // Clamp negative concentrations and zero out thin/dry cells
    #pragma omp parallel for collapse(2) schedule(static)
    for (unsigned int iy = 1; iy <= NCOLS; iy++) {
        for (unsigned int ix = 1; ix <= NROWS; ix++) {
            if (dry(ix, iy) == 1 || h0(ix, iy) < h_min) {
                conc(ix, iy) = 0.0;
            } else {
                conc(ix, iy) = std::fmax(conc(ix, iy), 0.0);
            }
        }
    }

    // ────────────────────────────────────────────────────────
    // PASS 1: Compute face fluxes (east-face and north-face)
    //
    //   flux_e(ix,iy) = mass flux through the EAST  face of cell (ix,iy)
    //   flux_n(ix,iy) = mass flux through the NORTH face of cell (ix,iy)
    //
    //   Positive = flow in +x / +y direction.
    //   Units: [g/s]  (mass flow rate,  = c [g/m³] × q [m²/s] × dy [m])
    //
    //   All reads are from the OLD concentration field — fully parallel.
    // ────────────────────────────────────────────────────────
    arma::mat flux_e(ds.MROWS, ds.MCOLS, arma::fill::zeros);
    arma::mat flux_n(ds.MROWS, ds.MCOLS, arma::fill::zeros);

    #pragma omp parallel for collapse(2) schedule(static)
    for (unsigned int iy = 1; iy <= NCOLS; iy++) {
        for (unsigned int ix = 1; ix <= NROWS; ix++) {

            double cp = conc(ix, iy);
            double hp = h0(ix, iy);

            // ═══ EAST FACE ═══
            // Shared between cell (ix,iy) and cell (ix+1,iy)
            {
                const unsigned int ie = ix + 1;
                double qxl = qx(ix, iy);           // hydro discharge/width at east face

                // Boundary: outflow at domain edge
                if (ix == NROWS) {
                    if (qxl > 0.0 && hp >= h_min) {
                        flux_e(ix, iy) = qxl * cp * dy;
                    }
                }
                // Both cells wet and deep enough
                else if (dry(ie, iy) == 0 && hp >= h_min && h0(ie, iy) >= h_min) {
                    double ce = conc(ie, iy);
                    double he = h0(ie, iy);

                    // First-order upwind advection
                    double c_up = (qxl >= 0.0) ? cp : ce;
                    double adv = qxl * c_up * dy;

                    // Fickian diffusion
                    double h_face = 0.5 * (hp + he);
                    double diff = -D_coef * h_face * dy * (ce - cp) / dx;

                    flux_e(ix, iy) = adv + diff;
                }
                // Current cell wet, east cell dry: advective outflow only
                else if (dry(ie, iy) == 1 && hp >= h_min && qxl > 0.0) {
                    flux_e(ix, iy) = qxl * cp * dy;
                }
                // Current cell dry/thin, east cell wet: advective inflow only
                else if (dry(ix, iy) == 0 && hp < h_min && h0(ie, iy) >= h_min && qxl < 0.0) {
                    double ce = conc(ie, iy);
                    flux_e(ix, iy) = qxl * ce * dy;  // negative flux = inflow from east
                }
            }

            // ═══ NORTH FACE ═══
            // Shared between cell (ix,iy) and cell (ix,iy+1)
            {
                const unsigned int in = iy + 1;
                double qyl = qy(ix, iy);           // hydro discharge/width at north face

                // Boundary: outflow at domain edge
                if (iy == NCOLS) {
                    if (qyl > 0.0 && hp >= h_min) {
                        flux_n(ix, iy) = qyl * cp * dx;
                    }
                }
                // Both cells wet and deep enough
                else if (dry(ix, in) == 0 && hp >= h_min && h0(ix, in) >= h_min) {
                    double cn = conc(ix, in);
                    double hn = h0(ix, in);

                    double c_up = (qyl >= 0.0) ? cp : cn;
                    double adv = qyl * c_up * dx;

                    double h_face = 0.5 * (hp + hn);
                    double diff = -D_coef * h_face * dx * (cn - cp) / dy;

                    flux_n(ix, iy) = adv + diff;
                }
                // Current cell wet, north cell dry: outflow only
                else if (dry(ix, in) == 1 && hp >= h_min && qyl > 0.0) {
                    flux_n(ix, iy) = qyl * cp * dx;
                }
                // Current cell dry/thin, north cell wet: inflow only
                else if (dry(ix, iy) == 0 && hp < h_min && h0(ix, in) >= h_min && qyl < 0.0) {
                    double cn = conc(ix, in);
                    flux_n(ix, iy) = qyl * cn * dx;
                }
            }
        }
    }

    // ────────────────────────────────────────────────────────
    // PASS 2: Update mass and convert to concentration
    //
    //   For each cell, gather the 4 face fluxes:
    //     West inflow  = +flux_e(ix-1, iy)    [positive east flux from west neighbor]
    //     East outflow = +flux_e(ix,   iy)    [our east face]
    //     South inflow = +flux_n(ix, iy-1)    [positive north flux from south neighbor]
    //     North outflow= +flux_n(ix,   iy)    [our north face]
    //
    //   Net flux in = (west - east) + (south - north)
    //   M_new = max(0, M_old + net_flux * dt / area)
    //   c_new = M_new / h_new   (clamped by monotonicity)
    // ────────────────────────────────────────────────────────
    #pragma omp parallel for collapse(2) schedule(static)
    for (unsigned int iy = 1; iy <= NCOLS; iy++) {
        for (unsigned int ix = 1; ix <= NROWS; ix++) {
            if (dry(ix, iy) == 1) {
                conc(ix, iy) = 0.0;
                continue;
            }

            double h_new = h(ix, iy);

            // If the cell is too shallow after the hydro step,
            // zero out concentration to prevent blowup
            if (h_new < h_min) {
                conc(ix, iy) = 0.0;
                continue;
            }

            double hp = std::fmax(h0(ix, iy), h_min);
            double cp = conc(ix, iy);
            double M_old = cp * hp;               // mass per unit area [g/m²]

            // Gather face fluxes [g/s]
            double F_west  = flux_e(ix - 1, iy);  // flux INTO cell from west face
            double F_east  = flux_e(ix,     iy);  // flux OUT of cell through east face
            double F_south = flux_n(ix, iy - 1);  // flux INTO cell from south face
            double F_north = flux_n(ix,     iy);  // flux OUT of cell through north face

            double net_flux = (F_west - F_east) + (F_south - F_north);

            double M_new = M_old + net_flux * dtfl / area;

            // Clamp mass to non-negative
            M_new = std::fmax(M_new, 0.0);

            // Convert to concentration
            double c_new = M_new / h_new;

            // ═══ MONOTONICITY CLAMP ═══
            // The new concentration must not exceed the maximum of
            // the cell's own old value and its face-neighbor old values.
            // This prevents artificial extrema (especially at shallow
            // wetting fronts where M/h_new can explode).
            double c_max = cp;
            c_max = std::fmax(c_max, conc(ix - 1, iy));      // west
            c_max = std::fmax(c_max, conc(ix + 1, iy));      // east
            c_max = std::fmax(c_max, conc(ix, iy - 1));      // south
            c_max = std::fmax(c_max, conc(ix, iy + 1));      // north

            conc(ix, iy) = std::fmin(c_new, c_max);
        }
    }

#ifdef USE_MPI
    mpi_domain.exchange_ghost_cells(conc);
#endif
}
