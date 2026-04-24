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

// Copyright 2019, Diogo Costa
// Forcing routines for triangular mesh

#include "tri_forcing.h"
#include <cmath>
#include <iostream>

// ============================================================================
// Add meteo forcing (rainfall/snowmelt) to all cells
// Same logic as regular mesh add_meteo but loops over cells
// ============================================================================
bool tri_add_meteo(
    GlobVar& ds,
    const TriMesh& mesh,
    TriSolution& sol,
    int nchem)
{
    bool errflag = false;

    // Return if no METEO_FILE provided
    if (ds.meteo_file.empty())
        return errflag;

    try {
        unsigned int meteo_rowi = 0;
        for (unsigned int a = 0; a <= (*ds.meteo).col(0).n_elem; a++) {
            meteo_rowi = a;
            if ((*ds.meteo).at(a, 0) > ds.tim) {
                break;
            }
        }

        // Convert mm/day to m/s then multiply by dt
        double meteoi = (*ds.meteo).at(meteo_rowi, 1) / (1000.0 * 3600.0 * 24.0) * ds.dtfl;

        // Get chem data
        std::vector<double> meteo_conci(nchem);
        for (int ichem = 0; ichem < nchem; ichem++) {
            meteo_conci[ichem] = (*ds.meteo).at(meteo_rowi, ichem + 2);
        }

        #pragma omp parallel for schedule(static)
        for (int ci = 0; ci < mesh.num_cells; ci++) {
            if (sol.innerNeumannBCWeir[ci] != 1.0f) {
                double hp = std::fmax(sol.z[ci] - sol.zb[ci], 0.0);

                sol.z[ci] += meteoi;
                sol.h[ci] = std::fmax(sol.z[ci] - sol.zb[ci], 0.0);

                if (sol.h[ci] <= ds.hdry) {
                    sol.ldry[ci] = 1.0f;
                } else {
                    sol.ldry[ci] = 0.0f;
                }
                sol.h0[ci] = sol.h[ci];

                // Concentration adjustment for mass balance
                if (ds.ade_solver && hp > 0.0) {
                    for (int ichem = 0; ichem < nchem; ichem++) {
                        sol.conc_SW[ichem][ci] =
                            (sol.conc_SW[ichem][ci] * hp + meteoi * meteo_conci[ichem])
                            / sol.h[ci];
                    }
                }
            }
        }
    } catch (...) {
        std::cout << "problem in 'tri_add_meteo' module" << std::endl;
        errflag = true;
    }

    return errflag;
}

// ============================================================================
// Add inflow at discharge point
// Uses inflow_nrow as cell index (set during initialization)
// ============================================================================
bool tri_add_inflow(
    GlobVar& ds,
    const TriMesh& mesh,
    TriSolution& sol,
    int nchem)
{
    bool errflag = false;

    // Return if no INFLOW_FILE provided
    if (ds.inflow_file.empty())
        return errflag;

    try {
        unsigned int inflow_rowi = 0;
        for (unsigned int a = 0; a <= (*ds.inflow).col(0).n_elem; a++) {
            inflow_rowi = a;
            if ((*ds.inflow).at(a, 0) > ds.tim) {
                break;
            }
        }

        // Distribute Q*dt across the injection PATCH built in
        // tri_initiate (a disk of cells around the inflow point, total
        // area ≈ π·dxy²). This decouples the per-cell depth rise from
        // the local mesh refinement, so a 5 m triangle and a 30 m
        // triangle get the same dh — matching the regular-mesh single-
        // cell injection at a coarse dxy scale.
        //
        // Fallback: if the patch is empty (should only happen when
        // tri_initiate could not find the inflow point at all), we
        // simply skip injection this step.
        if (ds.inflow_patch_cells.empty() ||
            ds.inflow_patch_total_area <= 0.0) {
            return errflag;
        }

        const double Q_dt = (*ds.inflow).at(inflow_rowi, 1) * ds.dtfl;
        const double dh   = Q_dt / ds.inflow_patch_total_area;

        // Get chem data
        std::vector<double> inflow_conci(nchem);
        for (int ichem = 0; ichem < nchem; ichem++) {
            inflow_conci[ichem] = (*ds.inflow).at(inflow_rowi, ichem + 2);
        }

        for (int ci : ds.inflow_patch_cells) {
            if (ci < 0 || ci >= mesh.num_cells) continue;
            if (sol.innerNeumannBCWeir[ci] == 1.0f) continue; // NODATA

            if (sol.z[ci] <= sol.zb[ci]) {
                sol.z[ci] = sol.zb[ci];
            }

            double hp = std::fmax(sol.z[ci] - sol.zb[ci], 0.0);
            sol.z[ci] += dh;
            sol.h[ci] = std::fmax(sol.z[ci] - sol.zb[ci], 0.0);

            sol.ldry[ci] = (sol.h[ci] <= ds.hdry) ? 1.0f : 0.0f;
            sol.h0[ci]   = sol.h[ci];

            // Concentration: local mixing of existing water (depth hp)
            // with injected water (depth dh) per cell.
            if (ds.ade_solver && sol.h[ci] > 0.0) {
                for (int ichem = 0; ichem < nchem; ichem++) {
                    sol.conc_SW[ichem][ci] =
                        (sol.conc_SW[ichem][ci] * hp + dh * inflow_conci[ichem])
                        / sol.h[ci];
                }
            }
        }
    } catch (...) {
        std::cout << "problem in 'tri_add_inflow' module" << std::endl;
        errflag = true;
    }

    return errflag;
}
