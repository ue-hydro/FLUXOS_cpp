// Copyright 2019, Diogo Costa
// CUDA kernels for triangular mesh hydrodynamics
// 1D thread indexing: 1 thread per cell or 1 thread per edge

#ifdef USE_CUDA

#include "tri_cuda_kernels.h"
#include <cuda_runtime.h>
#include <cfloat>
#include <cstdio>

#define BLOCK_SIZE 256

// ============================================================================
// Kernel: Wet/dry check (1 thread per cell)
// ============================================================================
__global__ void kernel_tri_wetdry(
    double* __restrict__ z,
    const double* __restrict__ zb,
    double* __restrict__ h,
    double* __restrict__ qx,
    double* __restrict__ qy,
    double* __restrict__ us,
    float* __restrict__ ldry,
    double* __restrict__ twetimetracer,
    int num_cells, double hdry, double dtfl)
{
    int ci = blockIdx.x * blockDim.x + threadIdx.x;
    if (ci >= num_cells) return;

    double hp = fmax(0.0, z[ci] - zb[ci]);
    h[ci] = hp;

    if (hp <= hdry) {
        qx[ci] = 0.0;
        qy[ci] = 0.0;
        us[ci] = 0.0;
        ldry[ci] = 1.0f;
    } else {
        ldry[ci] = 0.0f;
        twetimetracer[ci] += dtfl / 3600.0;
    }
}

// ============================================================================
// Kernel: LSQ gradient computation (1 thread per cell)
// ============================================================================
__global__ void kernel_tri_gradient(
    const double* __restrict__ z,
    const double* __restrict__ qx,
    const double* __restrict__ qy,
    const double* __restrict__ cell_cx,
    const double* __restrict__ cell_cy,
    const int* __restrict__ cell_neighbors,
    const int* __restrict__ cell_edges,
    const double* __restrict__ cell_lsq_inv,
    const double* __restrict__ edge_mx,
    const double* __restrict__ edge_my,
    double* __restrict__ grad_z_x, double* __restrict__ grad_z_y,
    double* __restrict__ grad_qx_x, double* __restrict__ grad_qx_y,
    double* __restrict__ grad_qy_x, double* __restrict__ grad_qy_y,
    int num_cells)
{
    int ci = blockIdx.x * blockDim.x + threadIdx.x;
    if (ci >= num_cells) return;

    double cx_i = cell_cx[ci];
    double cy_i = cell_cy[ci];

    double rhs_zx = 0.0, rhs_zy = 0.0;
    double rhs_qxx = 0.0, rhs_qxy = 0.0;
    double rhs_qyx = 0.0, rhs_qyy = 0.0;

    for (int e = 0; e < 3; e++) {
        int nj = cell_neighbors[ci * 3 + e];
        double dx, dy, w, dz, dqx_val, dqy_val;

        if (nj >= 0) {
            dx = cell_cx[nj] - cx_i;
            dy = cell_cy[nj] - cy_i;
            double dist = sqrt(dx * dx + dy * dy);
            w = (dist > 1e-15) ? (1.0 / dist) : 0.0;
            dz = z[nj] - z[ci];
            dqx_val = qx[nj] - qx[ci];
            dqy_val = qy[nj] - qy[ci];
        } else {
            int eidx = cell_edges[ci * 3 + e];
            dx = edge_mx[eidx] - cx_i;
            dy = edge_my[eidx] - cy_i;
            double dist = sqrt(dx * dx + dy * dy);
            w = (dist > 1e-15) ? (1.0 / dist) : 0.0;
            dz = 0.0;
            dqx_val = 0.0;
            dqy_val = 0.0;
        }

        rhs_zx += w * dz * dx;
        rhs_zy += w * dz * dy;
        rhs_qxx += w * dqx_val * dx;
        rhs_qxy += w * dqx_val * dy;
        rhs_qyx += w * dqy_val * dx;
        rhs_qyy += w * dqy_val * dy;
    }

    double inv0 = cell_lsq_inv[ci * 4 + 0];
    double inv1 = cell_lsq_inv[ci * 4 + 1];
    double inv2 = cell_lsq_inv[ci * 4 + 2];
    double inv3 = cell_lsq_inv[ci * 4 + 3];

    grad_z_x[ci] = inv0 * rhs_zx + inv1 * rhs_zy;
    grad_z_y[ci] = inv2 * rhs_zx + inv3 * rhs_zy;
    grad_qx_x[ci] = inv0 * rhs_qxx + inv1 * rhs_qxy;
    grad_qx_y[ci] = inv2 * rhs_qxx + inv3 * rhs_qxy;
    grad_qy_x[ci] = inv0 * rhs_qyx + inv1 * rhs_qyy;
    grad_qy_y[ci] = inv2 * rhs_qyx + inv3 * rhs_qyy;
}

// ============================================================================
// Kernel: Barth-Jespersen limiter (1 thread per cell)
// ============================================================================
__device__ double bj_limiter(double val_i, double grad_x, double grad_y,
                              double val_min, double val_max,
                              const double* edge_mx, const double* edge_my,
                              const int* cell_edges, int ci, double cx, double cy)
{
    double phi = 1.0;
    const double eps = 1e-12;

    for (int e = 0; e < 3; e++) {
        int eidx = cell_edges[ci * 3 + e];
        double dx = edge_mx[eidx] - cx;
        double dy = edge_my[eidx] - cy;
        double delta = grad_x * dx + grad_y * dy;

        double alpha;
        if (delta > eps) {
            alpha = fmin(1.0, (val_max - val_i) / delta);
        } else if (delta < -eps) {
            alpha = fmin(1.0, (val_min - val_i) / delta);
        } else {
            alpha = 1.0;
        }
        phi = fmin(phi, alpha);
    }
    return fmax(0.0, phi);
}

__global__ void kernel_tri_limiter(
    const double* __restrict__ z,
    const double* __restrict__ qx,
    const double* __restrict__ qy,
    const double* __restrict__ grad_z_x, const double* __restrict__ grad_z_y,
    const double* __restrict__ grad_qx_x, const double* __restrict__ grad_qx_y,
    const double* __restrict__ grad_qy_x, const double* __restrict__ grad_qy_y,
    const double* __restrict__ cell_cx, const double* __restrict__ cell_cy,
    const int* __restrict__ cell_neighbors,
    const int* __restrict__ cell_edges,
    const double* __restrict__ edge_mx, const double* __restrict__ edge_my,
    double* __restrict__ phi_z,
    double* __restrict__ phi_qx,
    double* __restrict__ phi_qy,
    int num_cells)
{
    int ci = blockIdx.x * blockDim.x + threadIdx.x;
    if (ci >= num_cells) return;

    double z_min = z[ci], z_max = z[ci];
    double qx_min = qx[ci], qx_max = qx[ci];
    double qy_min = qy[ci], qy_max = qy[ci];

    for (int e = 0; e < 3; e++) {
        int nj = cell_neighbors[ci * 3 + e];
        if (nj >= 0) {
            z_min = fmin(z_min, z[nj]); z_max = fmax(z_max, z[nj]);
            qx_min = fmin(qx_min, qx[nj]); qx_max = fmax(qx_max, qx[nj]);
            qy_min = fmin(qy_min, qy[nj]); qy_max = fmax(qy_max, qy[nj]);
        }
    }

    double cx = cell_cx[ci], cy = cell_cy[ci];
    phi_z[ci] = bj_limiter(z[ci], grad_z_x[ci], grad_z_y[ci], z_min, z_max,
                             edge_mx, edge_my, cell_edges, ci, cx, cy);
    phi_qx[ci] = bj_limiter(qx[ci], grad_qx_x[ci], grad_qx_y[ci], qx_min, qx_max,
                              edge_mx, edge_my, cell_edges, ci, cx, cy);
    phi_qy[ci] = bj_limiter(qy[ci], grad_qy_x[ci], grad_qy_y[ci], qy_min, qy_max,
                              edge_mx, edge_my, cell_edges, ci, cx, cy);
}

// ============================================================================
// Device helper: desingularized velocity
// ============================================================================
__device__ double d_desing_vel(double h, double q, double eps)
{
    double h2 = h * h;
    double eps2 = eps * eps;
    return 2.0 * h * q / (h2 + fmax(h2, eps2));
}

// ============================================================================
// Kernel: Edge flux computation (1 thread per edge)
// Rotated Roe/HLL solver with MUSCL + hydrostatic reconstruction
// ============================================================================
__global__ void kernel_tri_edge_flux(
    const double* __restrict__ z, const double* __restrict__ zb,
    const double* __restrict__ h,
    const double* __restrict__ qx_arr, const double* __restrict__ qy_arr,
    const double* __restrict__ us_arr,
    const float* __restrict__ ldry,
    const float* __restrict__ innerNeumann,
    const double* __restrict__ grad_z_x, const double* __restrict__ grad_z_y,
    const double* __restrict__ grad_qx_x, const double* __restrict__ grad_qx_y,
    const double* __restrict__ grad_qy_x, const double* __restrict__ grad_qy_y,
    const double* __restrict__ phi_z_arr,
    const double* __restrict__ phi_qx_arr,
    const double* __restrict__ phi_qy_arr,
    const double* __restrict__ cell_cx, const double* __restrict__ cell_cy,
    const double* __restrict__ cell_area,
    const int* __restrict__ edge_left, const int* __restrict__ edge_right,
    const double* __restrict__ edge_nx_arr, const double* __restrict__ edge_ny_arr,
    const double* __restrict__ edge_length, const double* __restrict__ edge_mx_arr,
    const double* __restrict__ edge_my_arr,
    const int* __restrict__ edge_btag,
    double* __restrict__ flux_mass,
    double* __restrict__ flux_momx,
    double* __restrict__ flux_momy,
    int num_edges, double gacc, double hdry, double cvdef, double nuem, double dtfl)
{
    int ei = blockIdx.x * blockDim.x + threadIdx.x;
    if (ei >= num_edges) return;

    int lc = edge_left[ei];
    int rc = edge_right[ei];

    // Check for Neumann cells
    bool lc_neumann = (innerNeumann[lc] == 1.0f);
    bool rc_neumann = (rc >= 0) ? (innerNeumann[rc] == 1.0f) : false;
    if (lc_neumann && (rc < 0 || rc_neumann)) {
        flux_mass[ei] = 0.0; flux_momx[ei] = 0.0; flux_momy[ei] = 0.0;
        return;
    }

    bool left_dry = (ldry[lc] == 1.0f);
    bool right_dry = (rc < 0) ? true : (ldry[rc] == 1.0f);

    if (left_dry && right_dry) {
        flux_mass[ei] = 0.0; flux_momx[ei] = 0.0; flux_momy[ei] = 0.0;
        return;
    }

    double nx = edge_nx_arr[ei];
    double ny = edge_ny_arr[ei];
    double tx = -ny, ty = nx;
    double emx = edge_mx_arr[ei];
    double emy = edge_my_arr[ei];

    // MUSCL reconstruction
    double dx_L = emx - cell_cx[lc];
    double dy_L = emy - cell_cy[lc];
    double z_L = z[lc] + phi_z_arr[lc] * (grad_z_x[lc] * dx_L + grad_z_y[lc] * dy_L);
    double qx_L = qx_arr[lc] + phi_qx_arr[lc] * (grad_qx_x[lc] * dx_L + grad_qx_y[lc] * dy_L);
    double qy_L = qy_arr[lc] + phi_qy_arr[lc] * (grad_qy_x[lc] * dx_L + grad_qy_y[lc] * dy_L);
    double zb_L = zb[lc];

    double z_R, qx_R, qy_R, zb_R;
    if (rc >= 0) {
        double dx_R = emx - cell_cx[rc];
        double dy_R = emy - cell_cy[rc];
        z_R = z[rc] + phi_z_arr[rc] * (grad_z_x[rc] * dx_R + grad_z_y[rc] * dy_R);
        qx_R = qx_arr[rc] + phi_qx_arr[rc] * (grad_qx_x[rc] * dx_R + grad_qx_y[rc] * dy_R);
        qy_R = qy_arr[rc] + phi_qy_arr[rc] * (grad_qy_x[rc] * dx_R + grad_qy_y[rc] * dy_R);
        zb_R = zb[rc];
    } else {
        if (edge_btag[ei] == 2) { // outflow
            z_R = z_L; qx_R = qx_L; qy_R = qy_L; zb_R = zb_L;
        } else { // wall: reflect
            z_R = z_L; zb_R = zb_L;
            double qn = qx_L * nx + qy_L * ny;
            double qt = qx_L * tx + qy_L * ty;
            qx_R = -qn * nx + qt * tx;
            qy_R = -qn * ny + qt * ty;
        }
    }

    // Hydrostatic reconstruction
    double zb_face = fmax(zb_L, zb_R);
    double h_L = fmax(0.0, z_L - zb_face);
    double h_R = fmax(0.0, z_R - zb_face);

    // Rotate to edge-normal frame
    double qn_L = qx_L * nx + qy_L * ny;
    double qt_L = qx_L * tx + qy_L * ty;
    double qn_R = qx_R * nx + qy_R * ny;
    double qt_R = qx_R * tx + qy_R * ty;

    double un_L = d_desing_vel(h_L, qn_L, hdry);
    double ut_L = d_desing_vel(h_L, qt_L, hdry);
    double un_R = d_desing_vel(h_R, qn_R, hdry);
    double ut_R = d_desing_vel(h_R, qt_R, hdry);

    if (h_L <= hdry && h_R <= hdry) {
        flux_mass[ei] = 0.0; flux_momx[ei] = 0.0; flux_momy[ei] = 0.0;
        return;
    }

    // HLL wave speeds
    double c_L = sqrt(gacc * fmax(h_L, 0.0));
    double c_R = sqrt(gacc * fmax(h_R, 0.0));
    double h_roe = 0.5 * (h_L + h_R);
    double sqrt_hL = sqrt(fmax(h_L, 0.0));
    double sqrt_hR = sqrt(fmax(h_R, 0.0));
    double denom = sqrt_hL + sqrt_hR;
    double u_roe = (denom > 1e-15) ? (sqrt_hL * un_L + sqrt_hR * un_R) / denom : 0.0;
    double c_roe = sqrt(gacc * h_roe);

    double s_L = fmin(un_L - c_L, u_roe - c_roe);
    double s_R = fmax(un_R + c_R, u_roe + c_roe);

    // Fluxes
    double f1_L = h_L * un_L;
    double f2_L = h_L * un_L * un_L + 0.5 * gacc * h_L * h_L;
    double f3_L = h_L * un_L * ut_L;
    double f1_R = h_R * un_R;
    double f2_R = h_R * un_R * un_R + 0.5 * gacc * h_R * h_R;
    double f3_R = h_R * un_R * ut_R;

    double fn1, fn2, fn3;
    if (s_L >= 0.0) {
        fn1 = f1_L; fn2 = f2_L; fn3 = f3_L;
    } else if (s_R <= 0.0) {
        fn1 = f1_R; fn2 = f2_R; fn3 = f3_R;
    } else {
        double inv_ds = 1.0 / (s_R - s_L);
        fn1 = (s_R * f1_L - s_L * f1_R + s_L * s_R * (h_R - h_L)) * inv_ds;
        fn2 = (s_R * f2_L - s_L * f2_R + s_L * s_R * (h_R * un_R - h_L * un_L)) * inv_ds;
        fn3 = (s_R * f3_L - s_L * f3_R + s_L * s_R * (h_R * ut_R - h_L * ut_L)) * inv_ds;
    }

    // Rotate back to global frame
    double elen = edge_length[ei];
    double fx = fn2 * nx + fn3 * tx;
    double fy = fn2 * ny + fn3 * ty;

    // Mass balance limiter
    double mflux = fn1 * elen;
    if (mflux > 0.0) {
        double vavail = h[lc] * cell_area[lc] / dtfl;
        if (mflux > vavail) {
            double s = vavail / mflux;
            fn1 *= s; fx *= s; fy *= s;
        }
    } else if (mflux < 0.0 && rc >= 0) {
        double vavail = h[rc] * cell_area[rc] / dtfl;
        if (-mflux > vavail) {
            double s = vavail / (-mflux);
            fn1 *= s; fx *= s; fy *= s;
        }
    }

    flux_mass[ei] = fn1 * elen;
    flux_momx[ei] = fx * elen;
    flux_momy[ei] = fy * elen;
}

// ============================================================================
// Kernel: Accumulate fluxes to cells (1 thread per edge)
// Uses atomicAdd to handle race conditions
// ============================================================================
__global__ void kernel_tri_accumulate(
    const int* __restrict__ edge_left, const int* __restrict__ edge_right,
    const double* __restrict__ flux_mass,
    const double* __restrict__ flux_momx,
    const double* __restrict__ flux_momy,
    const double* __restrict__ cell_area,
    double* __restrict__ dh, double* __restrict__ dqx, double* __restrict__ dqy,
    int num_edges, double dtfl)
{
    int ei = blockIdx.x * blockDim.x + threadIdx.x;
    if (ei >= num_edges) return;

    int lc = edge_left[ei];
    int rc = edge_right[ei];

    double fm = flux_mass[ei];
    double fmx = flux_momx[ei];
    double fmy = flux_momy[ei];

    double inv_area_L = dtfl / cell_area[lc];
    atomicAdd(&dh[lc], -fm * inv_area_L);
    atomicAdd(&dqx[lc], -fmx * inv_area_L);
    atomicAdd(&dqy[lc], -fmy * inv_area_L);

    if (rc >= 0) {
        double inv_area_R = dtfl / cell_area[rc];
        atomicAdd(&dh[rc], fm * inv_area_R);
        atomicAdd(&dqx[rc], fmx * inv_area_R);
        atomicAdd(&dqy[rc], fmy * inv_area_R);
    }
}

// ============================================================================
// Kernel: State update (1 thread per cell)
// ============================================================================
__global__ void kernel_tri_update(
    double* __restrict__ z, const double* __restrict__ zb,
    double* __restrict__ h,
    double* __restrict__ qx, double* __restrict__ qy,
    double* __restrict__ ux, double* __restrict__ uy,
    double* __restrict__ us_arr, const double* __restrict__ ks,
    float* __restrict__ ldry, const float* __restrict__ innerNeumann,
    const double* __restrict__ dh, const double* __restrict__ dqx, const double* __restrict__ dqy,
    int num_cells, double hdry, double gacc)
{
    int ci = blockIdx.x * blockDim.x + threadIdx.x;
    if (ci >= num_cells) return;

    z[ci] += dh[ci];
    double hp = fmax(0.0, z[ci] - zb[ci]);
    h[ci] = hp;

    if (hp < hdry || innerNeumann[ci] == 1.0f) {
        qx[ci] = 0.0; qy[ci] = 0.0;
        ux[ci] = 0.0; uy[ci] = 0.0;
        us_arr[ci] = 0.0;
        ldry[ci] = 1.0f;
    } else {
        qx[ci] += dqx[ci];
        qy[ci] += dqy[ci];

        double h2 = hp * hp;
        double eps2 = hdry * hdry;
        double den = h2 + fmax(h2, eps2);
        ux[ci] = 2.0 * hp * qx[ci] / den;
        uy[ci] = 2.0 * hp * qy[ci] / den;

        double speed = sqrt(ux[ci] * ux[ci] + uy[ci] * uy[ci]);
        double cf = gacc * ks[ci] * ks[ci] / pow(fmax(hp, ks[ci]), 4.0 / 3.0);
        us_arr[ci] = sqrt(cf) * speed;

        ldry[ci] = 0.0f;
    }
}

// ============================================================================
// Kernel: Courant condition with block reduction
// ============================================================================
__global__ void kernel_tri_courant(
    const double* __restrict__ h,
    const double* __restrict__ qx, const double* __restrict__ qy,
    double* __restrict__ h0, float* __restrict__ ldry, float* __restrict__ ldry_prev,
    const double* __restrict__ cell_inradius,
    double* __restrict__ block_reduce,
    int num_cells, double hdry, double gacc, double cfl)
{
    __shared__ double s_dtfl[BLOCK_SIZE];
    __shared__ double s_hpall[BLOCK_SIZE];

    int ci = blockIdx.x * blockDim.x + threadIdx.x;
    int tid = threadIdx.x;

    double local_dtfl = 9e10;
    double local_hpall = 0.0;

    if (ci < num_cells) {
        double hp = h[ci];
        h0[ci] = hp;
        ldry_prev[ci] = ldry[ci];

        if (hp > hdry) {
            ldry[ci] = 0.0f;
            hp = fmax(hp, hdry);
            local_hpall = hp;

            double c0 = sqrt(gacc * hp);
            double u0 = fmax(1e-6, fabs(qx[ci] / hp));
            double v0 = fmax(1e-6, fabs(qy[ci] / hp));
            double speed = sqrt(u0 * u0 + v0 * v0) + c0;
            local_dtfl = cfl * cell_inradius[ci] / speed;
        } else {
            ldry[ci] = 1.0f;
        }
    }

    s_dtfl[tid] = local_dtfl;
    s_hpall[tid] = local_hpall;
    __syncthreads();

    // Block reduction
    for (int s = BLOCK_SIZE / 2; s > 0; s >>= 1) {
        if (tid < s) {
            s_dtfl[tid] = fmin(s_dtfl[tid], s_dtfl[tid + s]);
            s_hpall[tid] = fmax(s_hpall[tid], s_hpall[tid + s]);
        }
        __syncthreads();
    }

    if (tid == 0) {
        block_reduce[blockIdx.x * 2] = s_dtfl[0];
        block_reduce[blockIdx.x * 2 + 1] = s_hpall[0];
    }
}

// ============================================================================
// Host function: Courant condition
// ============================================================================
void tri_cuda_courant_condition(TriCudaMemoryManager& cmem,
                                 double& dtfl_out,
                                 double& hpall_out)
{
    auto& d = cmem.data;
    int nc = d.num_cells;
    int num_blocks = (nc + BLOCK_SIZE - 1) / BLOCK_SIZE;

    kernel_tri_courant<<<num_blocks, BLOCK_SIZE>>>(
        d.d_h, d.d_qx, d.d_qy, d.d_h0, d.d_ldry, d.d_ldry_prev,
        d.d_cell_inradius, d.d_block_reduce,
        nc, d.hdry, d.gacc, d.cfl);

    // Read back block results
    std::vector<double> h_reduce(num_blocks * 2);
    cudaMemcpy(h_reduce.data(), d.d_block_reduce,
               num_blocks * 2 * sizeof(double), cudaMemcpyDeviceToHost);

    dtfl_out = 9e10;
    hpall_out = 0.0;
    for (int b = 0; b < num_blocks; b++) {
        dtfl_out = fmin(dtfl_out, h_reduce[b * 2]);
        hpall_out = fmax(hpall_out, h_reduce[b * 2 + 1]);
    }
}

// ============================================================================
// Host function: Full hydrodynamics timestep
// ============================================================================
void tri_cuda_hydrodynamics_calc(TriCudaMemoryManager& cmem)
{
    auto& d = cmem.data;
    int nc = d.num_cells;
    int ne = d.num_edges;
    int cell_blocks = (nc + BLOCK_SIZE - 1) / BLOCK_SIZE;
    int edge_blocks = (ne + BLOCK_SIZE - 1) / BLOCK_SIZE;

    // Phase 1: Wet/dry
    kernel_tri_wetdry<<<cell_blocks, BLOCK_SIZE>>>(
        d.d_z, d.d_zb, d.d_h, d.d_qx, d.d_qy, d.d_us,
        d.d_ldry, d.d_twetimetracer,
        nc, d.hdry, d.dtfl);

    // Phase 2: Gradients
    kernel_tri_gradient<<<cell_blocks, BLOCK_SIZE>>>(
        d.d_z, d.d_qx, d.d_qy,
        d.d_cell_cx, d.d_cell_cy,
        d.d_cell_neighbors, d.d_cell_edges, d.d_cell_lsq_inv,
        d.d_edge_mx, d.d_edge_my,
        d.d_grad_z_x, d.d_grad_z_y,
        d.d_grad_qx_x, d.d_grad_qx_y,
        d.d_grad_qy_x, d.d_grad_qy_y,
        nc);

    // Phase 3: Limiter
    kernel_tri_limiter<<<cell_blocks, BLOCK_SIZE>>>(
        d.d_z, d.d_qx, d.d_qy,
        d.d_grad_z_x, d.d_grad_z_y,
        d.d_grad_qx_x, d.d_grad_qx_y,
        d.d_grad_qy_x, d.d_grad_qy_y,
        d.d_cell_cx, d.d_cell_cy,
        d.d_cell_neighbors, d.d_cell_edges,
        d.d_edge_mx, d.d_edge_my,
        d.d_phi_z, d.d_phi_qx, d.d_phi_qy,
        nc);

    // Phase 4: Edge fluxes
    kernel_tri_edge_flux<<<edge_blocks, BLOCK_SIZE>>>(
        d.d_z, d.d_zb, d.d_h, d.d_qx, d.d_qy, d.d_us,
        d.d_ldry, d.d_innerNeumannBCWeir,
        d.d_grad_z_x, d.d_grad_z_y,
        d.d_grad_qx_x, d.d_grad_qx_y,
        d.d_grad_qy_x, d.d_grad_qy_y,
        d.d_phi_z, d.d_phi_qx, d.d_phi_qy,
        d.d_cell_cx, d.d_cell_cy, d.d_cell_area,
        d.d_edge_left, d.d_edge_right,
        d.d_edge_nx, d.d_edge_ny, d.d_edge_length,
        d.d_edge_mx, d.d_edge_my,
        d.d_edge_boundary_tag,
        d.d_flux_mass, d.d_flux_momx, d.d_flux_momy,
        ne, d.gacc, d.hdry, d.cvdef, d.nuem, d.dtfl);

    // Zero accumulators
    cudaMemset(d.d_dh, 0, nc * sizeof(double));
    cudaMemset(d.d_dqx, 0, nc * sizeof(double));
    cudaMemset(d.d_dqy, 0, nc * sizeof(double));

    // Phase 5: Accumulate fluxes
    kernel_tri_accumulate<<<edge_blocks, BLOCK_SIZE>>>(
        d.d_edge_left, d.d_edge_right,
        d.d_flux_mass, d.d_flux_momx, d.d_flux_momy,
        d.d_cell_area,
        d.d_dh, d.d_dqx, d.d_dqy,
        ne, d.dtfl);

    // Phase 6: State update
    kernel_tri_update<<<cell_blocks, BLOCK_SIZE>>>(
        d.d_z, d.d_zb, d.d_h, d.d_qx, d.d_qy,
        d.d_ux, d.d_uy, d.d_us, d.d_ks,
        d.d_ldry, d.d_innerNeumannBCWeir,
        d.d_dh, d.d_dqx, d.d_dqy,
        nc, d.hdry, d.gacc);
}

// ============================================================================
// Host function: ADE adjust (concentration adjustment for depth change)
// ============================================================================
void tri_cuda_ade_adjust(TriCudaMemoryManager& cmem)
{
    // The full ADE solver with edge-based dispersion still runs on CPU
    // This only handles the concentration adjustment for depth changes
    // which is embarrassingly parallel
    auto& d = cmem.data;
    // For simplicity, the ADE adjust is performed on CPU after copying data
    // A dedicated kernel can be added later if needed
}

#endif // USE_CUDA
