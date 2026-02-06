// Copyright 2019, Diogo Costa
// CUDA GPU acceleration for FLUXOS - Kernel implementations

#ifdef USE_CUDA

#include <cuda_runtime.h>
#include <cstdio>
#include <cfloat>
#include "cuda_memory.h"
#include "cuda_kernels.h"
#include "cuda_solver_wet.cuh"
#include "cuda_solver_dry.cuh"

#define CUDA_CHECK(call)                                                       \
    do {                                                                        \
        cudaError_t err = call;                                                 \
        if (err != cudaSuccess) {                                               \
            fprintf(stderr, "CUDA error at %s:%d: %s\n", __FILE__, __LINE__,   \
                    cudaGetErrorString(err));                                    \
            exit(EXIT_FAILURE);                                                 \
        }                                                                       \
    } while (0)

// Column-major index
#define IDX(r, c, MROWS) ((r) + (c) * (MROWS))

// Block dimensions for all 2D kernels
static constexpr int BLOCK_X = 16;
static constexpr int BLOCK_Y = 16;
static constexpr int BLOCK_SIZE = BLOCK_X * BLOCK_Y;

// ============================================================================
// Helper: ceil division
// ============================================================================
static inline unsigned int ceil_div(unsigned int a, unsigned int b) {
    return (a + b - 1) / b;
}

// ============================================================================
// KERNEL: Courant condition (CFL timestep calculation)
// ============================================================================
__global__ void kernel_courant(
    const double* __restrict__ d_h,
    double* __restrict__ d_h0,
    const float* __restrict__ d_ldry,
    float* __restrict__ d_ldry_prev,
    const double* __restrict__ d_qx,
    const double* __restrict__ d_qy,
    double* __restrict__ d_block_min,
    double* __restrict__ d_block_max,
    unsigned int MROWS, unsigned int NROWS, unsigned int NCOLS,
    double hdry, double gacc, double cfl, unsigned int dxy)
{
    unsigned int irow = blockIdx.x * blockDim.x + threadIdx.x + 1;
    unsigned int icol = blockIdx.y * blockDim.y + threadIdx.y + 1;
    unsigned int tid = threadIdx.x + threadIdx.y * blockDim.x;

    __shared__ double s_min[BLOCK_SIZE];
    __shared__ double s_max[BLOCK_SIZE];

    double local_dtfl = 9.0e10;
    double local_hpall = 0.0;

    if (irow <= NROWS && icol <= NCOLS) {
        unsigned int idx = IDX(irow, icol, MROWS);
        double hp = d_h[idx];

        // Save for ADE solver
        d_h0[idx] = hp;
        d_ldry_prev[idx] = d_ldry[idx];

        if (hp > hdry) {
            hp = fmax(hp, hdry);
            local_hpall = d_h[idx];
            double c0 = sqrt(gacc * d_h[idx]);
            double u0 = fmax(0.000001, fabs(d_qx[idx] / hp));
            double v0 = fmax(0.000001, fabs(d_qy[idx] / hp));
            double dt_u = cfl * dxy / (u0 + c0);
            double dt_v = cfl * dxy / (v0 + c0);
            local_dtfl = fmin(dt_u, dt_v);
        }
    }

    s_min[tid] = local_dtfl;
    s_max[tid] = local_hpall;
    __syncthreads();

    // Block-level reduction
    for (int stride = BLOCK_SIZE / 2; stride > 0; stride >>= 1) {
        if (tid < stride) {
            s_min[tid] = fmin(s_min[tid], s_min[tid + stride]);
            s_max[tid] = fmax(s_max[tid], s_max[tid + stride]);
        }
        __syncthreads();
    }

    if (tid == 0) {
        unsigned int bid = blockIdx.x + blockIdx.y * gridDim.x;
        d_block_min[bid] = s_min[0];
        d_block_max[bid] = s_max[0];
    }
}

// Second-stage reduction kernel
__global__ void kernel_reduce_final(
    const double* __restrict__ d_block_min,
    const double* __restrict__ d_block_max,
    double* __restrict__ d_result, // [0] = min, [1] = max
    unsigned int num_blocks)
{
    extern __shared__ double s_data[];
    double* s_min = s_data;
    double* s_max = s_data + blockDim.x;

    unsigned int tid = threadIdx.x;
    double local_min = 9.0e10;
    double local_max = 0.0;

    // Grid-stride loop over all blocks
    for (unsigned int i = tid; i < num_blocks; i += blockDim.x) {
        local_min = fmin(local_min, d_block_min[i]);
        local_max = fmax(local_max, d_block_max[i]);
    }

    s_min[tid] = local_min;
    s_max[tid] = local_max;
    __syncthreads();

    for (int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
        if (tid < stride) {
            s_min[tid] = fmin(s_min[tid], s_min[tid + stride]);
            s_max[tid] = fmax(s_max[tid], s_max[tid + stride]);
        }
        __syncthreads();
    }

    if (tid == 0) {
        d_result[0] = s_min[0];
        d_result[1] = s_max[0];
    }
}

void cuda_courant_condition(CudaMemoryManager& cmem,
                            double& dtfl_out,
                            double& hpall_out)
{
    CudaGridData& g = cmem.grid;
    dim3 block(BLOCK_X, BLOCK_Y);
    dim3 grid(ceil_div(g.NROWS, BLOCK_X), ceil_div(g.NCOLS, BLOCK_Y));
    unsigned int num_blocks = grid.x * grid.y;

    // Phase 1: per-block reduction
    kernel_courant<<<grid, block>>>(
        g.d_h, g.d_h0, g.d_ldry, g.d_ldry_prev,
        g.d_qx, g.d_qy,
        g.d_block_reduce, g.d_block_reduce + num_blocks,
        g.MROWS, g.NROWS, g.NCOLS,
        g.hdry, g.gacc, g.cfl, g.dxy);

    // Phase 2: final reduction (single block, 256 threads)
    unsigned int reduce_threads = 256;
    if (reduce_threads > num_blocks) reduce_threads = num_blocks;
    // Round up to power of 2
    unsigned int p = 1;
    while (p < reduce_threads) p <<= 1;
    reduce_threads = p;
    if (reduce_threads < 1) reduce_threads = 1;

    // Use d_block_reduce[0..1] for final result
    double* d_result;
    CUDA_CHECK(cudaMalloc(&d_result, 2 * sizeof(double)));

    kernel_reduce_final<<<1, reduce_threads, 2 * reduce_threads * sizeof(double)>>>(
        g.d_block_reduce, g.d_block_reduce + num_blocks,
        d_result, num_blocks);

    double result[2];
    CUDA_CHECK(cudaMemcpy(result, d_result, 2 * sizeof(double), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaFree(d_result));

    dtfl_out = result[0];
    hpall_out = result[1];
}

// ============================================================================
// KERNEL: Hydrodynamics Phase 1 - Wet/Dry check
// ============================================================================
__global__ void kernel_hydro_phase1(
    const double* __restrict__ d_z,
    const double* __restrict__ d_zb,
    double* __restrict__ d_h,
    double* __restrict__ d_qx,
    double* __restrict__ d_qy,
    double* __restrict__ d_us,
    float* __restrict__ d_ldry,
    double* __restrict__ d_twetimetracer,
    unsigned int MROWS, unsigned int NROWS, unsigned int NCOLS,
    double hdry, double dtfl)
{
    unsigned int irow = blockIdx.x * blockDim.x + threadIdx.x + 1;
    unsigned int icol = blockIdx.y * blockDim.y + threadIdx.y + 1;

    if (irow > NROWS || icol > NCOLS) return;

    unsigned int idx = IDX(irow, icol, MROWS);
    double hp = fmax(0.0, d_z[idx] - d_zb[idx]);
    d_h[idx] = hp;

    if (hp <= hdry) {
        d_qx[idx] = 0.0;
        d_qy[idx] = 0.0;
        d_us[idx] = 0.0;
        d_ldry[idx] = 1.0f;
    } else {
        d_ldry[idx] = 0.0f;
        d_twetimetracer[idx] += dtfl / 3600.0;
    }
}

// ============================================================================
// KERNEL: Hydrodynamics Phase 2 - Flux computation (calls solver_wet/dry)
// ============================================================================
__global__ void kernel_hydro_phase2(
    const double* __restrict__ d_zb,
    const double* __restrict__ d_z,
    double* __restrict__ d_qx,
    double* __restrict__ d_qy,
    const double* __restrict__ d_h,
    const float* __restrict__ d_ldry,
    const double* __restrict__ d_us,
    const double* __restrict__ d_ks,
    const float* __restrict__ d_innerNeumannBCWeir,
    double* __restrict__ d_dh,
    double* __restrict__ d_dqx,
    double* __restrict__ d_dqy,
    double* __restrict__ d_qxf,
    double* __restrict__ d_qyf,
    double* __restrict__ d_fe_1,
    double* __restrict__ d_fe_2,
    double* __restrict__ d_fe_3,
    double* __restrict__ d_fn_1,
    double* __restrict__ d_fn_2,
    double* __restrict__ d_fn_3,
    unsigned int MROWS, unsigned int NROWS, unsigned int NCOLS,
    double hdry, double gacc, double nuem, double cvdef,
    double dtfl, unsigned int dxy, double arbase)
{
    unsigned int irow = blockIdx.x * blockDim.x + threadIdx.x + 1;
    unsigned int icol = blockIdx.y * blockDim.y + threadIdx.y + 1;

    if (irow > NROWS || icol > NCOLS) return;

    unsigned int idx = IDX(irow, icol, MROWS);
    float cell_neumann = d_innerNeumannBCWeir[idx];

    if (d_ldry[idx] == 1.0f && cell_neumann == 0.0f) {
        solver_dry_device(
            d_zb, d_z, d_qx, d_qy, d_h, d_ldry,
            d_dh, d_dqx, d_dqy, d_qxf, d_qyf,
            d_fe_1, d_fe_2, d_fe_3,
            d_fn_1, d_fn_2, d_fn_3,
            irow, icol, MROWS, NROWS, NCOLS,
            gacc, hdry, dtfl, dxy);
    } else if (cell_neumann == 0.0f) {
        solver_wet_device(
            d_zb, d_z, d_qx, d_qy, d_h, d_ldry, d_us, d_ks,
            d_fe_1, d_fe_2, d_fe_3,
            d_fn_1, d_fn_2, d_fn_3,
            irow, icol, MROWS, NROWS, NCOLS,
            hdry, gacc, nuem, cvdef, dtfl, dxy, arbase);
    }
}

// ============================================================================
// KERNEL: Hydrodynamics Phase 3 - Flux derivatives
// ============================================================================
__global__ void kernel_hydro_phase3(
    const double* __restrict__ d_fe_1,
    const double* __restrict__ d_fe_2,
    const double* __restrict__ d_fe_3,
    const double* __restrict__ d_fn_1,
    const double* __restrict__ d_fn_2,
    const double* __restrict__ d_fn_3,
    const float* __restrict__ d_innerNeumannBCWeir,
    double* __restrict__ d_dh,
    double* __restrict__ d_dqx,
    double* __restrict__ d_dqy,
    double* __restrict__ d_qxf,
    double* __restrict__ d_qyf,
    unsigned int MROWS, unsigned int NROWS, unsigned int NCOLS,
    double dtfl, double inv_dxy)
{
    unsigned int irow = blockIdx.x * blockDim.x + threadIdx.x + 1;
    unsigned int icol = blockIdx.y * blockDim.y + threadIdx.y + 1;

    if (irow > NROWS || icol > NCOLS) return;

    unsigned int idx = IDX(irow, icol, MROWS);
    float cell_neumann = d_innerNeumannBCWeir[idx];

    if (cell_neumann == 0.0f) {
        double fe1_w = d_fe_1[IDX(irow - 1, icol, MROWS)];
        double fe1_p = d_fe_1[idx];
        double fn1_s = d_fn_1[IDX(irow, icol - 1, MROWS)];
        double fn1_p = d_fn_1[idx];

        d_dh[idx] = ((fe1_w - fe1_p) * inv_dxy + (fn1_s - fn1_p) * inv_dxy) * dtfl;

        double fe2_w = d_fe_2[IDX(irow - 1, icol, MROWS)];
        double fe2_p = d_fe_2[idx];
        double fn2_s = d_fn_2[IDX(irow, icol - 1, MROWS)];
        double fn2_p = d_fn_2[idx];

        d_dqx[idx] = ((fe2_w - fe2_p) * inv_dxy + (fn2_s - fn2_p) * inv_dxy) * dtfl;

        double fe3_w = d_fe_3[IDX(irow - 1, icol, MROWS)];
        double fe3_p = d_fe_3[idx];
        double fn3_s = d_fn_3[IDX(irow, icol - 1, MROWS)];
        double fn3_p = d_fn_3[idx];

        d_dqy[idx] = ((fe3_w - fe3_p) * inv_dxy + (fn3_s - fn3_p) * inv_dxy) * dtfl;

        d_qxf[idx] = fe1_p * dtfl;
        d_qyf[idx] = fn1_p * dtfl;
    } else {
        d_dh[idx] = 0.0;
        d_dqx[idx] = 0.0;
        d_dqy[idx] = 0.0;
        d_qxf[idx] = 0.0;
        d_qyf[idx] = 0.0;
    }
}

// ============================================================================
// KERNEL: Hydrodynamics Phase 4 - State update + smoothing
// ============================================================================
__global__ void kernel_hydro_phase4(
    double* __restrict__ d_z,
    const double* __restrict__ d_zb,
    double* __restrict__ d_h,
    double* __restrict__ d_qx,
    double* __restrict__ d_qy,
    double* __restrict__ d_us,
    float* __restrict__ d_ldry,
    const double* __restrict__ d_dh,
    const double* __restrict__ d_dqx,
    const double* __restrict__ d_dqy,
    const double* __restrict__ d_qxf,
    const double* __restrict__ d_qyf,
    const float* __restrict__ d_innerNeumannBCWeir,
    unsigned int MROWS, unsigned int NROWS, unsigned int NCOLS,
    double hdry)
{
    unsigned int irow = blockIdx.x * blockDim.x + threadIdx.x + 1;
    unsigned int icol = blockIdx.y * blockDim.y + threadIdx.y + 1;

    if (irow > NROWS || icol > NCOLS) return;

    unsigned int idx = IDX(irow, icol, MROWS);

    d_z[idx] = d_z[idx] + d_dh[idx];
    double hp = fmax(0.0, d_z[idx] - d_zb[idx]);
    d_h[idx] = hp;
    float cell_neumann = d_innerNeumannBCWeir[idx];

    if (hp < hdry || cell_neumann == 1.0f) {
        d_qx[idx] = 0.0;
        d_qy[idx] = 0.0;
        d_us[idx] = 0.0;
        d_ldry[idx] = 1.0f;
    } else {
        // Update momentum
        d_qx[idx] = d_qx[idx] + d_dqx[idx];
        d_qy[idx] = d_qy[idx] + d_dqy[idx];

        // Smoothing with neighboring fluxes
        double qxf_w = d_qxf[IDX(irow - 1, icol, MROWS)];
        double qxf_p = d_qxf[idx];
        d_qx[idx] = 0.1 * qxf_w + 0.8 * d_qx[idx] + 0.1 * qxf_p;

        double qyf_s = d_qyf[IDX(irow, icol - 1, MROWS)];
        double qyf_p = d_qyf[idx];
        d_qy[idx] = 0.1 * qyf_s + 0.8 * d_qy[idx] + 0.1 * qyf_p;

        d_ldry[idx] = 0.0f;
    }
}

// ============================================================================
// HOST: Launch hydrodynamics (4 phases)
// ============================================================================
void cuda_hydrodynamics_calc(CudaMemoryManager& cmem) {
    CudaGridData& g = cmem.grid;
    dim3 block(BLOCK_X, BLOCK_Y);
    dim3 grid_dim(ceil_div(g.NROWS, BLOCK_X), ceil_div(g.NCOLS, BLOCK_Y));

    double inv_dxy = 1.0 / (double)g.dxy;

    // Phase 1: Wet/dry check
    kernel_hydro_phase1<<<grid_dim, block>>>(
        g.d_z, g.d_zb, g.d_h, g.d_qx, g.d_qy, g.d_us,
        g.d_ldry, g.d_twetimetracer,
        g.MROWS, g.NROWS, g.NCOLS, g.hdry, g.dtfl);

    // Phase 2: Flux computation
    kernel_hydro_phase2<<<grid_dim, block>>>(
        g.d_zb, g.d_z, g.d_qx, g.d_qy, g.d_h, g.d_ldry,
        g.d_us, g.d_ks, g.d_innerNeumannBCWeir,
        g.d_dh, g.d_dqx, g.d_dqy, g.d_qxf, g.d_qyf,
        g.d_fe_1, g.d_fe_2, g.d_fe_3,
        g.d_fn_1, g.d_fn_2, g.d_fn_3,
        g.MROWS, g.NROWS, g.NCOLS,
        g.hdry, g.gacc, g.nuem, g.cvdef, g.dtfl, g.dxy, g.arbase);

    // Phase 3: Flux derivatives
    kernel_hydro_phase3<<<grid_dim, block>>>(
        g.d_fe_1, g.d_fe_2, g.d_fe_3,
        g.d_fn_1, g.d_fn_2, g.d_fn_3,
        g.d_innerNeumannBCWeir,
        g.d_dh, g.d_dqx, g.d_dqy, g.d_qxf, g.d_qyf,
        g.MROWS, g.NROWS, g.NCOLS, g.dtfl, inv_dxy);

    // Phase 4: State update + smoothing
    kernel_hydro_phase4<<<grid_dim, block>>>(
        g.d_z, g.d_zb, g.d_h, g.d_qx, g.d_qy, g.d_us, g.d_ldry,
        g.d_dh, g.d_dqx, g.d_dqy, g.d_qxf, g.d_qyf,
        g.d_innerNeumannBCWeir,
        g.MROWS, g.NROWS, g.NCOLS, g.hdry);

    // Ensure all kernels complete before returning
    CUDA_CHECK(cudaDeviceSynchronize());
}

// ============================================================================
// KERNEL: ADE concentration adjustment (parallelizable part)
// ============================================================================
__global__ void kernel_ade_adjust(
    double* __restrict__ d_conc_SW,
    const double* __restrict__ d_h,
    const double* __restrict__ d_h0,
    const float* __restrict__ d_ldry,
    const float* __restrict__ d_ldry_prev,
    unsigned int MROWS, unsigned int NROWS, unsigned int NCOLS)
{
    unsigned int irow = blockIdx.x * blockDim.x + threadIdx.x + 1;
    unsigned int icol = blockIdx.y * blockDim.y + threadIdx.y + 1;

    if (irow > NROWS || icol > NCOLS) return;

    unsigned int idx = IDX(irow, icol, MROWS);

    // Compute bounds from neighbors
    double c_w = d_conc_SW[IDX(irow - 1, icol, MROWS)];
    double c_e = d_conc_SW[IDX(irow + 1, icol, MROWS)];
    double c_s = d_conc_SW[IDX(irow, icol - 1, MROWS)];
    double c_n = d_conc_SW[IDX(irow, icol + 1, MROWS)];

    double cmax = fmax(c_w, fmax(c_e, fmax(c_s, c_n)));
    double cmin = fmin(c_w, fmin(c_e, fmin(c_s, c_n)));

    double hnew = d_h[idx];

    if (d_ldry[idx] == 0 && d_ldry_prev[idx] == 0) {
        double adjusted = d_conc_SW[idx] * d_h0[idx] / hnew;
        // Clamp to neighbor bounds
        adjusted = fmin(cmax, fmax(cmin, adjusted));
        d_conc_SW[idx] = adjusted;
    } else if (d_ldry[idx] == 1) {
        d_conc_SW[idx] = 0.0;
    }
}

void cuda_ade_adjust(CudaMemoryManager& cmem, int it) {
    if (it <= 1) return;

    CudaGridData& g = cmem.grid;
    dim3 block(BLOCK_X, BLOCK_Y);
    dim3 grid_dim(ceil_div(g.NROWS, BLOCK_X), ceil_div(g.NCOLS, BLOCK_Y));

    kernel_ade_adjust<<<grid_dim, block>>>(
        g.d_conc_SW, g.d_h, g.d_h0,
        g.d_ldry, g.d_ldry_prev,
        g.MROWS, g.NROWS, g.NCOLS);

    CUDA_CHECK(cudaDeviceSynchronize());
}

#undef IDX

#endif // USE_CUDA
