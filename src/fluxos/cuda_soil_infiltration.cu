// Copyright 2019, Diogo Costa
// CUDA kernel for soil infiltration loss (regular mesh)
// Horton decay: f(t) = fc + (f0 - fc) * exp(-k * t_wet)

#ifdef USE_CUDA

#include <cuda_runtime.h>
#include <cstdio>
#include "cuda_soil_infiltration.cuh"

#define SOIL_CUDA_CHECK(call)                                                  \
    do {                                                                        \
        cudaError_t err = call;                                                 \
        if (err != cudaSuccess) {                                               \
            fprintf(stderr, "CUDA error at %s:%d: %s\n", __FILE__, __LINE__,   \
                    cudaGetErrorString(err));                                    \
        }                                                                       \
    } while (0)

// ============================================================================
// Horton soil infiltration loss kernel (regular mesh, 2D grid)
// f(t) = fc + (f0 - fc) * exp(-k * t_wet)
// ============================================================================
__global__ void kernel_soil_infiltration(
    double* __restrict__ d_z,
    double* __restrict__ d_zb,
    double* __restrict__ d_h,
    float*  __restrict__ d_ldry,
    const float* __restrict__ d_innerNeumannBCWeir,
    double* __restrict__ d_soil_infil_rate,
    const double* __restrict__ d_soil_Ks,
    const double* __restrict__ d_soil_f0,
    const double* __restrict__ d_soil_k,
    double* __restrict__ d_soil_wetting_time,
    unsigned int NROWS, unsigned int NCOLS, size_t MCOLS,
    double hdry, double dt)
{
    unsigned int irow = blockIdx.y * blockDim.y + threadIdx.y + 1;
    unsigned int icol = blockIdx.x * blockDim.x + threadIdx.x + 1;

    if (irow > NROWS || icol > NCOLS) return;

    size_t idx = irow * MCOLS + icol;

    if (d_innerNeumannBCWeir[idx] == 1.0f) return;

    double Ks = d_soil_Ks[idx];
    if (Ks <= 0.0) return;

    double h_surface = fmax(d_z[idx] - d_zb[idx], 0.0);
    if (h_surface <= 0.0) return;

    // Horton decay: f(t) = fc + (f0 - fc) * exp(-k * t_wet)
    double f0 = d_soil_f0[idx];
    double k_decay = d_soil_k[idx];
    double t_wet = d_soil_wetting_time[idx];

    double f_rate = Ks + (f0 - Ks) * exp(-k_decay * t_wet);

    // Infiltration = min(f(t) * dt, available water)
    double infil = fmin(f_rate * dt, h_surface);

    // Update surface water
    d_z[idx] -= infil;
    d_h[idx] = fmax(d_z[idx] - d_zb[idx], 0.0);
    d_ldry[idx] = (d_h[idx] <= hdry) ? 1.0f : 0.0f;

    // Store rate for output
    d_soil_infil_rate[idx] = infil / dt;

    // Advance wetting time
    d_soil_wetting_time[idx] += dt;
}

// ============================================================================
// Host function to launch infiltration kernel
// ============================================================================
void cuda_apply_soil_infiltration(CudaGridData& g)
{
    dim3 block(16, 16);
    dim3 grid((g.NCOLS + block.x - 1) / block.x,
              (g.NROWS + block.y - 1) / block.y);

    kernel_soil_infiltration<<<grid, block>>>(
        g.d_z, g.d_zb, g.d_h, g.d_ldry,
        g.d_innerNeumannBCWeir,
        g.d_soil_infil_rate, g.d_soil_Ks,
        g.d_soil_f0, g.d_soil_k, g.d_soil_wetting_time,
        g.NROWS, g.NCOLS, g.MCOLS,
        g.hdry, g.dtfl);

    SOIL_CUDA_CHECK(cudaGetLastError());
}

#endif // USE_CUDA
