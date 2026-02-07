// Copyright 2019, Diogo Costa
// CUDA kernel for soil infiltration loss (triangular mesh)

#ifndef TRI_CUDA_SOIL_INFILTRATION_CUH
#define TRI_CUDA_SOIL_INFILTRATION_CUH

#ifdef USE_CUDA

#include "tri_cuda_memory.h"

// Launch soil infiltration loss kernel on GPU (triangular mesh)
void tri_cuda_apply_soil_infiltration(TriCudaData& d);

#endif // USE_CUDA
#endif // TRI_CUDA_SOIL_INFILTRATION_CUH
