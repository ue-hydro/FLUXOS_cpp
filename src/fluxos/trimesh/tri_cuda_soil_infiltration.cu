// Copyright 2019, Diogo Costa
// CUDA kernel for soil infiltration loss (triangular mesh)
// Horton decay: f(t) = fc + (f0 - fc) * exp(-k * t_wet)

#ifdef USE_CUDA

#include <cuda_runtime.h>
#include <cstdio>
#include "tri_cuda_soil_infiltration.cuh"

#define TRI_SOIL_CUDA_CHECK(call)                                              \
    do {                                                                        \
        cudaError_t err = call;                                                 \
        if (err != cudaSuccess) {                                               \
            fprintf(stderr, "CUDA error at %s:%d: %s\n", __FILE__, __LINE__,   \
                    cudaGetErrorString(err));                                    \
        }                                                                       \
    } while (0)

// ============================================================================
// Horton soil infiltration loss kernel (triangular mesh, 1D cell array)
// f(t) = fc + (f0 - fc) * exp(-k * t_wet)
// ============================================================================
__global__ void kernel_tri_soil_infiltration(
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
    int num_cells, double hdry, double dt)
{
    int ci = blockIdx.x * blockDim.x + threadIdx.x;
    if (ci >= num_cells) return;

    if (d_innerNeumannBCWeir[ci] == 1.0f) return;

    double Ks = d_soil_Ks[ci];
    if (Ks <= 0.0) return;

    double h_surface = fmax(d_z[ci] - d_zb[ci], 0.0);
    if (h_surface <= 0.0) return;

    // Horton decay: f(t) = fc + (f0 - fc) * exp(-k * t_wet)
    double f0 = d_soil_f0[ci];
    double k_decay = d_soil_k[ci];
    double t_wet = d_soil_wetting_time[ci];

    double f_rate = Ks + (f0 - Ks) * exp(-k_decay * t_wet);

    double infil = fmin(f_rate * dt, h_surface);

    d_z[ci] -= infil;
    d_h[ci] = fmax(d_z[ci] - d_zb[ci], 0.0);
    d_ldry[ci] = (d_h[ci] <= hdry) ? 1.0f : 0.0f;

    d_soil_infil_rate[ci] = infil / dt;

    // Advance wetting time
    d_soil_wetting_time[ci] += dt;
}

// ============================================================================
// Host function to launch triangular mesh infiltration kernel
// ============================================================================
void tri_cuda_apply_soil_infiltration(TriCudaData& d)
{
    const int BLOCK_SIZE = 256;
    int num_blocks = (d.num_cells + BLOCK_SIZE - 1) / BLOCK_SIZE;

    kernel_tri_soil_infiltration<<<num_blocks, BLOCK_SIZE>>>(
        d.d_z, d.d_zb, d.d_h, d.d_ldry,
        d.d_innerNeumannBCWeir,
        d.d_soil_infil_rate, d.d_soil_Ks,
        d.d_soil_f0, d.d_soil_k, d.d_soil_wetting_time,
        d.num_cells, d.hdry, d.dtfl);

    TRI_SOIL_CUDA_CHECK(cudaGetLastError());
}

#endif // USE_CUDA
