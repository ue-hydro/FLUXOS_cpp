// Copyright 2019, Diogo Costa
// WINTRA soil release solver for triangular mesh

#include "tri_wintra_solver.h"
#include <cmath>

// ============================================================================
// WINTRA solver: cell-based soil mass release
// Same physics as regular mesh, but uses per-cell area instead of constant arbase
// ============================================================================
void tri_wintrasolver_calc(
    GlobVar& ds,
    const TriMesh& mesh,
    TriSolution& sol,
    int ichem)
{
    const int ncells = mesh.num_cells;

    #pragma omp parallel for schedule(static)
    for (int ci = 0; ci < ncells; ci++) {
        double hp = sol.h[ci];
        double zbp = sol.zb[ci];

        if (hp > ds.hdry && sol.innerNeumannBCWeir[ci] != 1.0f) {
            // Mass release: deltam = soil_mass * rate/3600 * dt
            double deltam = sol.soil_mass[ci] * ds.soil_release_rate / 3600.0 * ds.dtfl;

            sol.soil_mass[ci] -= deltam;

            // Concentration increase: deltam / (h * cell_area)
            sol.conc_SW[ichem][ci] += deltam / (hp * mesh.cells[ci].area);
        }
    }
}
