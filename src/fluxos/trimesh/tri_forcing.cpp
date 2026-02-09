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
                if (ds.ade_solver && hp > 0.0 && ds.openwq == false) {
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

        int ci = ds.inflow_nrow;  // Cell index for triangular mesh

        if (ci < 0 || ci >= mesh.num_cells) {
            return errflag;  // No valid inflow cell
        }

        // Inflow as m3/s, convert to depth change: dh = Q * dt / cell_area
        double inflowi = (*ds.inflow).at(inflow_rowi, 1) * ds.dtfl / mesh.cells[ci].area;

        // Get chem data
        std::vector<double> inflow_conci(nchem);
        for (int ichem = 0; ichem < nchem; ichem++) {
            inflow_conci[ichem] = (*ds.inflow).at(inflow_rowi, ichem + 2);
        }

        if (sol.innerNeumannBCWeir[ci] != 1.0f) {
            if (sol.z[ci] <= sol.zb[ci]) {
                sol.z[ci] = sol.zb[ci];
            }

            double hp = std::fmax(sol.z[ci] - sol.zb[ci], 0.0);
            sol.z[ci] += inflowi;
            sol.h[ci] = std::fmax(sol.z[ci] - sol.zb[ci], 0.0);

            if (sol.h[ci] <= ds.hdry)
                sol.ldry[ci] = 1.0f;
            else
                sol.ldry[ci] = 0.0f;

            sol.h0[ci] = sol.h[ci];

            // Concentration adjustment
            if (ds.ade_solver && hp > 0.0 && ds.openwq == false) {
                for (int ichem = 0; ichem < nchem; ichem++) {
                    sol.conc_SW[ichem][ci] =
                        (sol.conc_SW[ichem][ci] * hp + inflowi * inflow_conci[ichem])
                        / sol.h[ci];
                }
            }
        } else {
            std::cout << "inflow cell is in NODATA_VALUE region" << std::endl;
            errflag = true;
        }
    } catch (...) {
        std::cout << "problem in 'tri_add_inflow' module" << std::endl;
        errflag = true;
    }

    return errflag;
}
