// Copyright 2019, Diogo Costa
// Multi-phase hydrodynamics driver for triangular mesh

#include "tri_hydrodynamics.h"
#include "tri_gradient.h"
#include "tri_solver_wet.h"
#include "tri_solver_dry.h"
#include <cmath>
#include <algorithm>
#include <limits>

// ============================================================================
// Courant condition for triangular mesh
// CFL = cfl * inradius / (|u| + c) for each cell
// ============================================================================
void tri_courant_condition(
    GlobVar& ds,
    const TriMesh& mesh,
    TriSolution& sol,
    double& dtfl_out,
    double& hpall_out)
{
    const int ncells = mesh.num_cells;
    const double hdry = ds.hdry;
    const double gacc = ds.gacc;
    const double cfl = ds.cfl;

    double dtfl_local = 9.0e10;
    double hpall_local = 0.0;

    #pragma omp parallel for schedule(static) reduction(min:dtfl_local) reduction(max:hpall_local)
    for (int ci = 0; ci < ncells; ci++) {
        double hp = sol.h[ci];

        // Save h0 for ADE solver
        sol.h0[ci] = hp;
        sol.ldry_prev[ci] = sol.ldry[ci];

        if (hp > hdry) {
            sol.ldry[ci] = 0.0f;
            hp = std::fmax(hp, hdry);
            hpall_local = std::fmax(hpall_local, hp);

            double c0 = std::sqrt(gacc * hp);
            double u0 = std::fmax(1e-6, std::fabs(sol.qx[ci] / hp));
            double v0 = std::fmax(1e-6, std::fabs(sol.qy[ci] / hp));
            double speed = std::sqrt(u0 * u0 + v0 * v0) + c0;

            // CFL with inradius (characteristic length for triangles)
            double dt_candidate = cfl * mesh.cells[ci].inradius / speed;
            dtfl_local = std::fmin(dt_candidate, dtfl_local);
        } else {
            sol.ldry[ci] = 1.0f;
        }
    }

    dtfl_out = dtfl_local;
    hpall_out = hpall_local;
}

// ============================================================================
// Main hydrodynamics timestep driver
// ============================================================================
void tri_hydrodynamics_calc(
    GlobVar& ds,
    const TriMesh& mesh,
    TriSolution& sol)
{
    const int ncells = mesh.num_cells;
    const int nedges = mesh.num_edges;
    const double hdry = ds.hdry;
    const double gacc = ds.gacc;
    const double cvdef = ds.cvdef;
    const double nuem = ds.nuem;
    const double dtfl = ds.dtfl;

    // ================================================================
    // Phase 1: Wet/dry check and compute water depth
    // ================================================================
    #pragma omp parallel for schedule(static)
    for (int ci = 0; ci < ncells; ci++) {
        double hp = std::fmax(0.0, sol.z[ci] - sol.zb[ci]);
        sol.h[ci] = hp;

        if (hp <= hdry) {
            sol.qx[ci] = 0.0;
            sol.qy[ci] = 0.0;
            sol.us[ci] = 0.0;
            sol.ldry[ci] = 1.0f;
        } else {
            sol.ldry[ci] = 0.0f;
            sol.twetimetracer[ci] += dtfl / 3600.0;
        }
    }

    // ================================================================
    // Phase 2: Compute gradients (least-squares)
    // ================================================================
    tri_compute_gradients(mesh, sol);

    // ================================================================
    // Phase 3: Apply Barth-Jespersen limiter
    // ================================================================
    tri_apply_limiter(mesh, sol);

    // ================================================================
    // Phase 4: Compute edge fluxes
    // For each edge, choose wet or dry solver based on neighbor states
    // ================================================================
    #pragma omp parallel for schedule(static)
    for (int ei = 0; ei < nedges; ei++) {
        const Edge& edge = mesh.edges[ei];
        int lc = edge.left_cell;
        int rc = edge.right_cell;

        bool left_dry = (sol.ldry[lc] == 1.0f);
        bool right_dry = (rc < 0) ? true : (sol.ldry[rc] == 1.0f);
        bool left_neumann = (sol.innerNeumannBCWeir[lc] == 1.0f);
        bool right_neumann = (rc >= 0) ? (sol.innerNeumannBCWeir[rc] == 1.0f) : false;

        // Skip Neumann BC cells
        if (left_neumann && (rc < 0 || right_neumann)) {
            sol.flux_mass[ei] = 0.0;
            sol.flux_momx[ei] = 0.0;
            sol.flux_momy[ei] = 0.0;
            continue;
        }

        if (left_dry && right_dry) {
            // Both dry: zero flux
            sol.flux_mass[ei] = 0.0;
            sol.flux_momx[ei] = 0.0;
            sol.flux_momy[ei] = 0.0;
        } else if (left_dry || right_dry) {
            // One side dry: use Ritter dry-front solver
            tri_edge_flux_dry(mesh, sol, ei, gacc, hdry, dtfl);
        } else {
            // Both wet: use full Roe solver
            tri_edge_flux_wet(mesh, sol, ei, gacc, hdry, cvdef, nuem, dtfl);
        }
    }

    // ================================================================
    // Phase 5: Accumulate fluxes to cells and update state
    // dU_i = -(dt / A_i) * sum_edges(F_edge * sign)
    // ================================================================

    // Zero accumulators
    #pragma omp parallel for schedule(static)
    for (int ci = 0; ci < ncells; ci++) {
        sol.dh[ci] = 0.0;
        sol.dqx[ci] = 0.0;
        sol.dqy[ci] = 0.0;
    }

    // Accumulate edge fluxes to cells
    // Note: This loop is NOT parallelizable with OpenMP due to race conditions
    // (multiple edges write to same cell). Use atomic operations or coloring.
    for (int ei = 0; ei < nedges; ei++) {
        const Edge& edge = mesh.edges[ei];
        int lc = edge.left_cell;
        int rc = edge.right_cell;

        // Left cell: flux goes OUT (subtract)
        double inv_area_L = 1.0 / mesh.cells[lc].area;
        sol.dh[lc] -= sol.flux_mass[ei] * dtfl * inv_area_L;
        sol.dqx[lc] -= sol.flux_momx[ei] * dtfl * inv_area_L;
        sol.dqy[lc] -= sol.flux_momy[ei] * dtfl * inv_area_L;

        // Right cell: flux comes IN (add)
        if (rc >= 0) {
            double inv_area_R = 1.0 / mesh.cells[rc].area;
            sol.dh[rc] += sol.flux_mass[ei] * dtfl * inv_area_R;
            sol.dqx[rc] += sol.flux_momx[ei] * dtfl * inv_area_R;
            sol.dqy[rc] += sol.flux_momy[ei] * dtfl * inv_area_R;
        }
    }

    // ================================================================
    // Phase 6: State update
    // ================================================================
    #pragma omp parallel for schedule(static)
    for (int ci = 0; ci < ncells; ci++) {
        // Update water surface elevation
        sol.z[ci] += sol.dh[ci];

        // Recompute depth
        double hp = std::fmax(0.0, sol.z[ci] - sol.zb[ci]);
        sol.h[ci] = hp;

        if (hp < hdry || sol.innerNeumannBCWeir[ci] == 1.0f) {
            sol.qx[ci] = 0.0;
            sol.qy[ci] = 0.0;
            sol.us[ci] = 0.0;
            sol.ux[ci] = 0.0;
            sol.uy[ci] = 0.0;
            sol.ldry[ci] = 1.0f;
        } else {
            // Update momentum
            sol.qx[ci] += sol.dqx[ci];
            sol.qy[ci] += sol.dqy[ci];

            // Compute velocities (desingularized)
            double h2 = hp * hp;
            double eps2 = hdry * hdry;
            double denom = h2 + std::fmax(h2, eps2);
            sol.ux[ci] = 2.0 * hp * sol.qx[ci] / denom;
            sol.uy[ci] = 2.0 * hp * sol.qy[ci] / denom;

            // Shear stress velocity: us = sqrt(g * h * Sf)
            double speed = std::sqrt(sol.ux[ci] * sol.ux[ci] + sol.uy[ci] * sol.uy[ci]);
            double cf = gacc * sol.ks[ci] * sol.ks[ci] / std::pow(std::fmax(hp, sol.ks[ci]), 4.0 / 3.0);
            sol.us[ci] = std::sqrt(cf) * speed;

            sol.ldry[ci] = 0.0f;
        }
    }
}
