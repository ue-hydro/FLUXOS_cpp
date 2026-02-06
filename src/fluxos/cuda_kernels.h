// Copyright 2019, Diogo Costa
// CUDA GPU acceleration for FLUXOS - Kernel launch declarations
// These are plain C++ declarations callable from .cpp files

#ifndef CUDA_KERNELS_H_INCLUDED
#define CUDA_KERNELS_H_INCLUDED

#ifdef USE_CUDA

#include "cuda_memory.h"

// Courant condition: computes dtfl_min and hpall_max on GPU
void cuda_courant_condition(CudaMemoryManager& cmem,
                            double& dtfl_out,
                            double& hpall_out);

// Hydrodynamics: launches 4 phase kernels sequentially on GPU
void cuda_hydrodynamics_calc(CudaMemoryManager& cmem);

// ADE concentration adjustment (parallelizable part only)
void cuda_ade_adjust(CudaMemoryManager& cmem, int it);

#endif // USE_CUDA
#endif // CUDA_KERNELS_H_INCLUDED
