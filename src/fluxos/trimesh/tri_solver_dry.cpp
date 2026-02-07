// Copyright 2019, Diogo Costa
// Dry-front Ritter solution in edge-normal frame for triangular mesh

#include "tri_solver_dry.h"
#include <cmath>
#include <algorithm>

// ============================================================================
// Compute Ritter dry-front flux for an edge where one side is dry
// Ritter solution: analytical dam-break onto dry bed
// q = 0.296 * dz * sqrt(g * dz) in edge-normal direction
// ============================================================================
void tri_edge_flux_dry(
    const TriMesh& mesh,
    TriSolution& sol,
    int edge_id,
    double gacc,
    double hdry,
    double dtfl)
{
    const Edge& edge = mesh.edges[edge_id];
    const int lc = edge.left_cell;
    const int rc = edge.right_cell;

    const double nx = edge.nx;
    const double ny = edge.ny;

    // Determine which side is wet and which is dry
    double h_wet, z_wet, zb_wet;
    int wet_cell;
    double sign;  // +1 if wet cell is on left, -1 if on right

    if (rc < 0) {
        // Boundary edge: left cell is the only cell
        h_wet = sol.h[lc];
        z_wet = sol.z[lc];
        zb_wet = sol.zb[lc];
        wet_cell = lc;
        sign = 1.0;
    } else if (sol.ldry[lc] == 0.0f && sol.ldry[rc] == 1.0f) {
        // Left is wet, right is dry
        h_wet = sol.h[lc];
        z_wet = sol.z[lc];
        zb_wet = sol.zb[lc];
        wet_cell = lc;
        sign = 1.0;
    } else if (sol.ldry[lc] == 1.0f && sol.ldry[rc] == 0.0f) {
        // Left is dry, right is wet
        h_wet = sol.h[rc];
        z_wet = sol.z[rc];
        zb_wet = sol.zb[rc];
        wet_cell = rc;
        sign = -1.0;
    } else {
        // Both dry or both wet: shouldn't be called, zero flux
        sol.flux_mass[edge_id] = 0.0;
        sol.flux_momx[edge_id] = 0.0;
        sol.flux_momy[edge_id] = 0.0;
        return;
    }

    if (h_wet <= hdry) {
        sol.flux_mass[edge_id] = 0.0;
        sol.flux_momx[edge_id] = 0.0;
        sol.flux_momy[edge_id] = 0.0;
        return;
    }

    // Bed elevation at face
    double zb_face;
    if (rc >= 0) {
        zb_face = std::fmax(sol.zb[lc], sol.zb[rc]);
    } else {
        zb_face = zb_wet;
    }

    // Water depth at face from wet side
    double h_face = std::fmax(0.0, z_wet - zb_face);

    if (h_face <= hdry) {
        sol.flux_mass[edge_id] = 0.0;
        sol.flux_momx[edge_id] = 0.0;
        sol.flux_momy[edge_id] = 0.0;
        return;
    }

    // Ritter solution: outgoing flux in normal direction
    // dz = depth difference driving flow
    double dz = std::fmin(h_face, h_wet);

    // Ritter discharge per unit width
    double qn = 0.296 * dz * std::sqrt(gacc * dz);

    // Per-edge volume rate limiter: each edge may drain at most 1/3 of cell volume
    double vol_rate = h_wet * mesh.cells[wet_cell].area / (3.0 * dtfl);
    double edge_flux_rate = qn * edge.length;
    if (edge_flux_rate > vol_rate) {
        qn = vol_rate / edge.length;
    }

    // Ritter depth at front
    double h_ritter = 0.444 * dz;
    double un = (h_ritter > 1e-15) ? (qn / h_ritter) : 0.0;

    // Pressure flux
    double pressure = 0.5 * gacc * h_ritter * h_ritter;

    // Mass flux (positive = from wet to dry)
    double flux_mass_n = sign * qn * edge.length;

    // Momentum flux in normal frame
    double fn2 = sign * (qn * un + pressure);
    double fn3 = 0.0;  // No tangential momentum in dry front

    // Rotate to global frame
    double fx_global = fn2 * nx;
    double fy_global = fn2 * ny;

    sol.flux_mass[edge_id] = flux_mass_n;
    sol.flux_momx[edge_id] = fx_global * edge.length;
    sol.flux_momy[edge_id] = fy_global * edge.length;
}
