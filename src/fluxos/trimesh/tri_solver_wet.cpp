// Copyright 2019, Diogo Costa
// Rotated Roe flux solver with MUSCL + hydrostatic reconstruction
// Edge-based computation for unstructured triangular mesh

#include "tri_solver_wet.h"
#include <cmath>
#include <algorithm>

// ============================================================================
// Velocity desingularization for wetting/drying
// u = 2*h*q / (h^2 + max(h^2, eps^2))
// ============================================================================
static inline double desing_vel(double h, double q, double eps)
{
    double h2 = h * h;
    double eps2 = eps * eps;
    return 2.0 * h * q / (h2 + std::fmax(h2, eps2));
}

// ============================================================================
// Compute Roe flux for a single edge
// Rotated Roe solver: rotate into edge-normal frame, solve 1D Riemann,
// rotate back
// ============================================================================
void tri_edge_flux_wet(
    const TriMesh& mesh,
    TriSolution& sol,
    int edge_id,
    double gacc,
    double hdry,
    double cvdef,
    double nuem,
    double dtfl)
{
    const Edge& edge = mesh.edges[edge_id];
    const int lc = edge.left_cell;
    const int rc = edge.right_cell;

    // Edge normal and tangent
    const double nx = edge.nx;
    const double ny = edge.ny;
    const double tx = -ny;  // tangent = 90 degrees CW from normal
    const double ty = nx;

    // ---- MUSCL reconstruction with limiter ----
    // Left cell reconstruction to edge midpoint
    double dx_L = edge.mx - mesh.cells[lc].cx;
    double dy_L = edge.my - mesh.cells[lc].cy;

    double z_L = sol.z[lc] + sol.phi_z[lc] * (sol.grad_z_x[lc] * dx_L + sol.grad_z_y[lc] * dy_L);
    double qx_L = sol.qx[lc] + sol.phi_qx[lc] * (sol.grad_qx_x[lc] * dx_L + sol.grad_qx_y[lc] * dy_L);
    double qy_L = sol.qy[lc] + sol.phi_qy[lc] * (sol.grad_qy_x[lc] * dx_L + sol.grad_qy_y[lc] * dy_L);
    double zb_L = sol.zb[lc];

    double z_R, qx_R, qy_R, zb_R;

    if (rc >= 0) {
        // Internal edge: reconstruct from right cell
        double dx_R = edge.mx - mesh.cells[rc].cx;
        double dy_R = edge.my - mesh.cells[rc].cy;

        z_R = sol.z[rc] + sol.phi_z[rc] * (sol.grad_z_x[rc] * dx_R + sol.grad_z_y[rc] * dy_R);
        qx_R = sol.qx[rc] + sol.phi_qx[rc] * (sol.grad_qx_x[rc] * dx_R + sol.grad_qx_y[rc] * dy_R);
        qy_R = sol.qy[rc] + sol.phi_qy[rc] * (sol.grad_qy_x[rc] * dx_R + sol.grad_qy_y[rc] * dy_R);
        zb_R = sol.zb[rc];
    } else {
        // Boundary edge: apply boundary condition
        if (edge.boundary_tag == 2) {
            // Outflow: zero-gradient (transmissive)
            z_R = z_L;
            qx_R = qx_L;
            qy_R = qy_L;
            zb_R = zb_L;
        } else {
            // Wall (default): reflect normal velocity, keep tangential
            z_R = z_L;
            double qn_L = qx_L * nx + qy_L * ny;
            double qt_L = qx_L * tx + qy_L * ty;
            double qn_R = -qn_L;  // reflect normal component
            qx_R = qn_R * nx + qt_L * tx;
            qy_R = qn_R * ny + qt_L * ty;
            zb_R = zb_L;
        }
    }

    // ---- Hydrostatic reconstruction (preserves C-property / lake-at-rest) ----
    double zb_face = std::fmax(zb_L, zb_R);
    double h_L = std::fmax(0.0, z_L - zb_face);
    double h_R = std::fmax(0.0, z_R - zb_face);

    // ---- Rotate velocities into edge-normal frame ----
    double qn_L = qx_L * nx + qy_L * ny;   // normal momentum
    double qt_L = qx_L * tx + qy_L * ty;   // tangential momentum
    double qn_R = qx_R * nx + qy_R * ny;
    double qt_R = qx_R * tx + qy_R * ty;

    // Desingularized velocities
    const double eps = hdry;
    double un_L = desing_vel(h_L, qn_L, eps);
    double ut_L = desing_vel(h_L, qt_L, eps);
    double un_R = desing_vel(h_R, qn_R, eps);
    double ut_R = desing_vel(h_R, qt_R, eps);

    // ---- Both sides dry: zero flux ----
    if (h_L <= hdry && h_R <= hdry) {
        sol.flux_mass[edge_id] = 0.0;
        sol.flux_momx[edge_id] = 0.0;
        sol.flux_momy[edge_id] = 0.0;
        return;
    }

    // ---- Roe averages in normal frame ----
    double sqrt_hL = std::sqrt(std::fmax(h_L, 0.0));
    double sqrt_hR = std::sqrt(std::fmax(h_R, 0.0));
    double denom = sqrt_hL + sqrt_hR;

    double u_roe, v_roe, c_roe, h_roe;
    if (denom > 1e-15) {
        u_roe = (sqrt_hL * un_L + sqrt_hR * un_R) / denom;
        v_roe = (sqrt_hL * ut_L + sqrt_hR * ut_R) / denom;
        h_roe = 0.5 * (h_L + h_R);
        c_roe = std::sqrt(gacc * h_roe);
    } else {
        u_roe = 0.0;
        v_roe = 0.0;
        h_roe = 0.0;
        c_roe = 0.0;
    }

    // ---- Wave speeds ----
    double c_L = std::sqrt(gacc * std::fmax(h_L, 0.0));
    double c_R = std::sqrt(gacc * std::fmax(h_R, 0.0));

    // HLL wave speed estimates
    double s_L = std::fmin(un_L - c_L, u_roe - c_roe);
    double s_R = std::fmax(un_R + c_R, u_roe + c_roe);

    // ---- Compute 1D fluxes in normal frame ----
    // Left flux
    double f1_L = h_L * un_L;
    double f2_L = h_L * un_L * un_L + 0.5 * gacc * h_L * h_L;
    double f3_L = h_L * un_L * ut_L;

    // Right flux
    double f1_R = h_R * un_R;
    double f2_R = h_R * un_R * un_R + 0.5 * gacc * h_R * h_R;
    double f3_R = h_R * un_R * ut_R;

    // State vectors
    double U1_L = h_L;
    double U2_L = h_L * un_L;
    double U3_L = h_L * ut_L;
    double U1_R = h_R;
    double U2_R = h_R * un_R;
    double U3_R = h_R * ut_R;

    // ---- HLL flux ----
    double fn1, fn2, fn3;  // flux in normal frame

    if (s_L >= 0.0) {
        // Supersonic from left
        fn1 = f1_L;
        fn2 = f2_L;
        fn3 = f3_L;
    } else if (s_R <= 0.0) {
        // Supersonic from right
        fn1 = f1_R;
        fn2 = f2_R;
        fn3 = f3_R;
    } else {
        // Subsonic: HLL average
        double inv_ds = 1.0 / (s_R - s_L);
        fn1 = (s_R * f1_L - s_L * f1_R + s_L * s_R * (U1_R - U1_L)) * inv_ds;
        fn2 = (s_R * f2_L - s_L * f2_R + s_L * s_R * (U2_R - U2_L)) * inv_ds;
        fn3 = (s_R * f3_L - s_L * f3_R + s_L * s_R * (U3_R - U3_L)) * inv_ds;
    }

    // ---- Add Roe dissipation (improves resolution of contact discontinuities) ----
    // Eigenvalues and eigenvector correction
    if (c_roe > 1e-15 && h_roe > hdry) {
        double dz = (z_R - zb_face) - (z_L - zb_face);  // = h_R - h_L
        double dqn = qn_R - qn_L;
        double dqt = qt_R - qt_L;

        // Limit jumps
        double dza = std::fabs(dz);
        if (dza > 0.5 * h_roe) {
            double dha = std::fabs(h_R - h_L);
            dza = std::fmin(dza, dha);
            dz = std::copysign(dza, dz);
        }

        double cc = 0.25 / c_roe;
        double c2 = u_roe + c_roe;
        double c3 = u_roe - c_roe;
        double c1a = std::fabs(u_roe);
        double c2a = std::fabs(c2);
        double c3a = std::fabs(c3);

        // Roe dissipation matrix applied to jumps
        double r1 = -(c2 * c3a - c2a * c3) * dz * cc
                     -(c2a - c3a) * dqn * cc;
        double r2 = -(c2 * c3 * (c3a - c2a)) * dz * cc
                     -(c2a * c2 - c3a * c3) * dqn * cc;
        double r3 = -(v_roe * (c2 * c3a - 2.0 * c_roe * c1a - c2a * c3)) * dz * cc
                     -(v_roe * (c2a - c3a)) * dqn * cc
                     -(2.0 * c_roe * c1a) * dqt * cc;

        // Blend HLL with Roe dissipation (50% blend for stability)
        fn1 = 0.5 * (f1_L + f1_R) + 0.5 * r1;
        fn2 = 0.5 * (f2_L + f2_R) + 0.5 * r2;
        fn3 = 0.5 * (f3_L + f3_R) + 0.5 * r3;
    }

    // ---- Hydrostatic pressure correction for well-balancedness ----
    // Source term due to bed step: 0.5*g*(h_L^2 - h_Lstar^2) on left side
    double h_L_orig = std::fmax(0.0, z_L - zb_L);
    double h_R_orig = (rc >= 0) ? std::fmax(0.0, z_R - zb_R) : h_L_orig;
    double pressure_corr = 0.5 * gacc * (h_L_orig * h_L_orig - h_L * h_L
                                        - h_R_orig * h_R_orig + h_R * h_R);
    // This is applied as a source term rather than modifying the flux directly

    // ---- Turbulent stress contribution ----
    double turb_fn2 = 0.0, turb_fn3 = 0.0;
    if (h_roe > hdry && rc >= 0) {
        double h_face = h_roe;
        double nu_L = cvdef * sol.us[lc] * std::fmax(sol.h[lc], 0.0) + nuem;
        double nu_R = cvdef * sol.us[rc] * std::fmax(sol.h[rc], 0.0) + nuem;
        double nu_face = 0.5 * (nu_L + nu_R) * std::sqrt(std::fmax(sol.h[lc], 0.0) * std::fmax(sol.h[rc], 0.0));

        // Gradient-based turbulent stress at edge
        // Use average gradients from left and right cells
        double dudx = 0.5 * (sol.grad_qx_x[lc] + sol.grad_qx_x[rc]);
        double dudy = 0.5 * (sol.grad_qx_y[lc] + sol.grad_qx_y[rc]);
        double dvdx = 0.5 * (sol.grad_qy_x[lc] + sol.grad_qy_x[rc]);
        double dvdy = 0.5 * (sol.grad_qy_y[lc] + sol.grad_qy_y[rc]);

        // Stress tensor: tau_ij = nu_t * (du_i/dx_j + du_j/dx_i)
        // Normal flux of x-momentum: tau_xx * nx + tau_xy * ny
        // Normal flux of y-momentum: tau_xy * nx + tau_yy * ny
        double tau_xx = nu_face * 2.0 * dudx;
        double tau_yy = nu_face * 2.0 * dvdy;
        double tau_xy = nu_face * (dudy + dvdx);

        turb_fn2 = -(tau_xx * nx + tau_xy * ny);
        turb_fn3 = -(tau_xy * nx + tau_yy * ny);
    }

    // ---- Rotate flux back to global frame ----
    // F_global = fn1 (mass unchanged), fn2*n + fn3*t for momentum
    double fx_global = (fn2 + turb_fn2) * nx + (fn3 + turb_fn3) * tx;
    double fy_global = (fn2 + turb_fn2) * ny + (fn3 + turb_fn3) * ty;

    // ---- Mass balance check ----
    // Limit outgoing mass flux to available water
    double mass_flux = fn1 * edge.length;
    if (mass_flux > 0.0) {
        // Outgoing from left cell
        double vol_avail = sol.h[lc] * mesh.cells[lc].area;
        double max_flux = vol_avail / dtfl;
        if (mass_flux > max_flux) {
            double scale = max_flux / mass_flux;
            fn1 *= scale;
            fx_global *= scale;
            fy_global *= scale;
        }
    } else if (mass_flux < 0.0 && rc >= 0) {
        // Outgoing from right cell
        double vol_avail = sol.h[rc] * mesh.cells[rc].area;
        double max_flux = vol_avail / dtfl;
        if (-mass_flux > max_flux) {
            double scale = max_flux / (-mass_flux);
            fn1 *= scale;
            fx_global *= scale;
            fy_global *= scale;
        }
    }

    // ---- Store edge fluxes (multiplied by edge length) ----
    sol.flux_mass[edge_id] = fn1 * edge.length;
    sol.flux_momx[edge_id] = fx_global * edge.length;
    sol.flux_momy[edge_id] = fy_global * edge.length;
}

// ============================================================================
// Compute all edge fluxes
// ============================================================================
void tri_compute_edge_fluxes(
    const TriMesh& mesh,
    TriSolution& sol,
    double gacc,
    double hdry,
    double cvdef,
    double nuem,
    double dtfl)
{
    const int nedges = mesh.num_edges;

    #pragma omp parallel for schedule(static)
    for (int ei = 0; ei < nedges; ei++) {
        tri_edge_flux_wet(mesh, sol, ei, gacc, hdry, cvdef, nuem, dtfl);
    }
}
