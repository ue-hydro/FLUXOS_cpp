// Copyright 2019, Diogo Costa
// CUDA kernel declarations for triangular mesh
// Plain C++ declarations callable from .cpp files

#ifndef TRI_CUDA_KERNELS_H_INCLUDED
#define TRI_CUDA_KERNELS_H_INCLUDED

#ifdef USE_CUDA

#include "tri_cuda_memory.h"

// Courant condition: compute dtfl_min and hpall_max on GPU
void tri_cuda_courant_condition(TriCudaMemoryManager& cmem,
                                 double& dtfl_out,
                                 double& hpall_out);

// Full hydrodynamics timestep on GPU
// Launches 6 kernels: wetdry -> gradient -> limiter -> edge_flux -> accumulate -> update
void tri_cuda_hydrodynamics_calc(TriCudaMemoryManager& cmem);

// ADE concentration adjustment on GPU
void tri_cuda_ade_adjust(TriCudaMemoryManager& cmem);

#endif // USE_CUDA
#endif // TRI_CUDA_KERNELS_H_INCLUDED
