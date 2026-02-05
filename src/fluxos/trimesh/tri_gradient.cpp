// Copyright 2019, Diogo Costa
// Least-squares gradient reconstruction + Barth-Jespersen limiter

#include "tri_gradient.h"
#include <cmath>
#include <algorithm>

// ============================================================================
// Compute least-squares gradients for z, qx, qy
// grad(phi)_i = LSQ_inv * sum_j[ w_j * (phi_j - phi_i) * [dx_j, dy_j]^T ]
// ============================================================================
void tri_compute_gradients(
    const TriMesh& mesh,
    TriSolution& sol)
{
    const int ncells = mesh.num_cells;

    #pragma omp parallel for schedule(static)
    for (int ci = 0; ci < ncells; ci++) {
        const auto& cell = mesh.cells[ci];

        // Right-hand side vectors for LSQ
        double rhs_z_x = 0.0, rhs_z_y = 0.0;
        double rhs_qx_x = 0.0, rhs_qx_y = 0.0;
        double rhs_qy_x = 0.0, rhs_qy_y = 0.0;

        for (int e = 0; e < 3; e++) {
            double dx, dy, w;
            double dz, dqx, dqy;

            int nj = cell.neighbors[e];

            if (nj >= 0) {
                // Internal neighbor
                dx = mesh.cells[nj].cx - cell.cx;
                dy = mesh.cells[nj].cy - cell.cy;
                double dist = std::sqrt(dx * dx + dy * dy);
                w = (dist > 1e-15) ? (1.0 / dist) : 0.0;

                dz = sol.z[nj] - sol.z[ci];
                dqx = sol.qx[nj] - sol.qx[ci];
                dqy = sol.qy[nj] - sol.qy[ci];
            } else {
                // Boundary: use edge midpoint with zero-gradient (wall)
                int eidx = cell.edges[e];
                dx = mesh.edges[eidx].mx - cell.cx;
                dy = mesh.edges[eidx].my - cell.cy;
                double dist = std::sqrt(dx * dx + dy * dy);
                w = (dist > 1e-15) ? (1.0 / dist) : 0.0;

                // Zero-gradient BC: ghost value = cell value
                dz = 0.0;
                dqx = 0.0;
                dqy = 0.0;
            }

            rhs_z_x += w * dz * dx;
            rhs_z_y += w * dz * dy;
            rhs_qx_x += w * dqx * dx;
            rhs_qx_y += w * dqx * dy;
            rhs_qy_x += w * dqy * dx;
            rhs_qy_y += w * dqy * dy;
        }

        // Apply precomputed LSQ inverse: grad = LSQ_inv * rhs
        const double* inv = &cell.lsq_inv[0][0];

        sol.grad_z_x[ci] = inv[0] * rhs_z_x + inv[1] * rhs_z_y;
        sol.grad_z_y[ci] = inv[2] * rhs_z_x + inv[3] * rhs_z_y;

        sol.grad_qx_x[ci] = inv[0] * rhs_qx_x + inv[1] * rhs_qx_y;
        sol.grad_qx_y[ci] = inv[2] * rhs_qx_x + inv[3] * rhs_qx_y;

        sol.grad_qy_x[ci] = inv[0] * rhs_qy_x + inv[1] * rhs_qy_y;
        sol.grad_qy_y[ci] = inv[2] * rhs_qy_x + inv[3] * rhs_qy_y;
    }
}

// ============================================================================
// Helper: compute Barth-Jespersen limiter for a single variable
// phi = min over edges of: min(1, alpha_edge)
// where alpha_edge bounds the reconstructed value at edge midpoint
// ============================================================================
static inline double barth_jespersen_phi(
    double val_i,        // cell center value
    double grad_x,       // gradient x
    double grad_y,       // gradient y
    double val_min,      // min of neighbors
    double val_max,      // max of neighbors
    const TriMesh& mesh,
    const TriCell& cell)
{
    double phi = 1.0;
    const double eps = 1e-12;

    for (int e = 0; e < 3; e++) {
        int eidx = cell.edges[e];
        // Vector from centroid to edge midpoint
        double dx = mesh.edges[eidx].mx - cell.cx;
        double dy = mesh.edges[eidx].my - cell.cy;

        // Reconstructed increment
        double delta = grad_x * dx + grad_y * dy;

        double alpha;
        if (delta > eps) {
            alpha = std::min(1.0, (val_max - val_i) / delta);
        } else if (delta < -eps) {
            alpha = std::min(1.0, (val_min - val_i) / delta);
        } else {
            alpha = 1.0;  // No reconstruction -> no limiting needed
        }

        phi = std::min(phi, alpha);
    }

    return std::max(0.0, phi);
}

// ============================================================================
// Apply Barth-Jespersen limiter to gradients
// ============================================================================
void tri_apply_limiter(
    const TriMesh& mesh,
    TriSolution& sol)
{
    const int ncells = mesh.num_cells;

    #pragma omp parallel for schedule(static)
    for (int ci = 0; ci < ncells; ci++) {
        const auto& cell = mesh.cells[ci];

        // Find min/max of z, qx, qy among cell and its neighbors
        double z_min = sol.z[ci], z_max = sol.z[ci];
        double qx_min = sol.qx[ci], qx_max = sol.qx[ci];
        double qy_min = sol.qy[ci], qy_max = sol.qy[ci];

        for (int e = 0; e < 3; e++) {
            int nj = cell.neighbors[e];
            if (nj >= 0) {
                z_min = std::min(z_min, sol.z[nj]);
                z_max = std::max(z_max, sol.z[nj]);
                qx_min = std::min(qx_min, sol.qx[nj]);
                qx_max = std::max(qx_max, sol.qx[nj]);
                qy_min = std::min(qy_min, sol.qy[nj]);
                qy_max = std::max(qy_max, sol.qy[nj]);
            }
        }

        // Compute limiter values
        sol.phi_z[ci] = barth_jespersen_phi(
            sol.z[ci], sol.grad_z_x[ci], sol.grad_z_y[ci],
            z_min, z_max, mesh, cell);

        sol.phi_qx[ci] = barth_jespersen_phi(
            sol.qx[ci], sol.grad_qx_x[ci], sol.grad_qx_y[ci],
            qx_min, qx_max, mesh, cell);

        sol.phi_qy[ci] = barth_jespersen_phi(
            sol.qy[ci], sol.grad_qy_x[ci], sol.grad_qy_y[ci],
            qy_min, qy_max, mesh, cell);
    }
}
