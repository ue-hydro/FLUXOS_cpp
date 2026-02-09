// Copyright 2019, Diogo Costa
// Edge-based advection-dispersion equation solver for triangular mesh

#include "tri_ade_solver.h"
#include <cmath>
#include <algorithm>
#include <vector>

// ============================================================================
// ADE solver for triangular mesh
// Edge-based formulation replaces the Cartesian x/y sweep
// ============================================================================
void tri_adesolver_calc(
    GlobVar& ds,
    const TriMesh& mesh,
    TriSolution& sol,
    int it,
    int ichem)
{
    if (it <= 1) return;

    const int ncells = mesh.num_cells;
    const int nedges = mesh.num_edges;
    const double dtfl = ds.dtfl;
    const double hdry = ds.hdry;
    const double D_coef = ds.D_coef;
    const double nt = 1.0;    // eddy viscosity (m2/s)
    const double sigc = 0.5;  // turbulent Schmidt number

    auto& conc = sol.conc_SW[ichem];

    // ---- Phase 1: Concentration bounds for limiting ----
    std::vector<double> cmaxr(ncells);
    std::vector<double> cminr(ncells);

    #pragma omp parallel for schedule(static)
    for (int ci = 0; ci < ncells; ci++) {
        cmaxr[ci] = conc[ci];
        cminr[ci] = conc[ci];

        for (int e = 0; e < 3; e++) {
            int nj = mesh.cells[ci].neighbors[e];
            if (nj >= 0) {
                cmaxr[ci] = std::fmax(cmaxr[ci], conc[nj]);
                cminr[ci] = std::fmin(cminr[ci], conc[nj]);
            }
        }
    }

    // ---- Phase 2: Adjust concentration to new depth ----
    #pragma omp parallel for schedule(static)
    for (int ci = 0; ci < ncells; ci++) {
        double hnew = sol.h[ci];
        if (sol.ldry[ci] == 0.0f && sol.ldry_prev[ci] == 0.0f) {
            conc[ci] = conc[ci] * sol.h0[ci] / std::fmax(hnew, hdry);
        } else if (sol.ldry[ci] == 1.0f) {
            conc[ci] = 0.0;
        }
    }

    // ---- Phase 3: Compute edge fluxes for concentration ----
    std::vector<double> conc_flux(nedges, 0.0);

    #pragma omp parallel for schedule(static)
    for (int ei = 0; ei < nedges; ei++) {
        const auto& edge = mesh.edges[ei];
        int lc = edge.left_cell;
        int rc = edge.right_cell;

        if (rc < 0) continue;  // Skip boundary edges for now

        bool left_dry = (sol.ldry[lc] == 1.0f);
        bool right_dry = (sol.ldry[rc] == 1.0f);

        if (left_dry && right_dry) {
            conc_flux[ei] = 0.0;
            continue;
        }

        double h_L = sol.h[lc];
        double h_R = sol.h[rc];
        double c_L = conc[lc];
        double c_R = conc[rc];

        // Advective flux: upwind based on mass flux direction
        double mass_flux = sol.flux_mass[ei];  // from hydrodynamics
        double c_face;

        if (mass_flux >= 0.0) {
            // Flow from left to right: upwind from left
            if (!left_dry) {
                c_face = c_L;
            } else {
                c_face = 0.0;
            }
        } else {
            // Flow from right to left: upwind from right
            if (!right_dry) {
                c_face = c_R;
            } else {
                c_face = 0.0;
            }
        }

        double adv_flux = mass_flux * c_face;

        // Dispersive flux: -D_eff * (c_R - c_L) / dist_LR * edge_length
        double disp_flux = 0.0;
        if (!left_dry && !right_dry && edge.dist_LR > 1e-15) {
            double h_face = std::sqrt(std::fmax(h_L * nt * h_R * nt, 1e-8));
            double D_eff = h_face / sigc * D_coef;
            disp_flux = -D_eff * (c_R - c_L) / edge.dist_LR * edge.length;
        }

        conc_flux[ei] = adv_flux + disp_flux;

        // Limit outgoing flux
        if (conc_flux[ei] > 0.0 && !left_dry) {
            double max_flux = c_L * h_L * mesh.cells[lc].area / dtfl;
            conc_flux[ei] = std::fmin(conc_flux[ei], max_flux);
        } else if (conc_flux[ei] < 0.0 && !right_dry) {
            double max_flux = c_R * h_R * mesh.cells[rc].area / dtfl;
            conc_flux[ei] = std::fmax(conc_flux[ei], -max_flux);
        }
    }

    // ---- Phase 4: Accumulate fluxes and update concentration ----
    // Sequential to avoid race conditions
    std::vector<double> dc(ncells, 0.0);

    for (int ei = 0; ei < nedges; ei++) {
        const auto& edge = mesh.edges[ei];
        int lc = edge.left_cell;
        int rc = edge.right_cell;

        if (rc < 0) continue;

        // Left cell loses flux, right cell gains flux
        dc[lc] -= conc_flux[ei] * dtfl / mesh.cells[lc].area;
        dc[rc] += conc_flux[ei] * dtfl / mesh.cells[rc].area;
    }

    // Apply changes with bounding
    #pragma omp parallel for schedule(static)
    for (int ci = 0; ci < ncells; ci++) {
        if (sol.ldry[ci] == 1.0f) {
            conc[ci] = 0.0;
            continue;
        }

        double hp = std::fmax(sol.h[ci], hdry);
        double con = conc[ci] + dc[ci] / hp;

        // Bound by neighbor min/max
        con = std::fmin(cmaxr[ci], con);
        con = std::fmax(cminr[ci], con);
        con = std::fmax(0.0, con);

        conc[ci] = con;
    }
}
