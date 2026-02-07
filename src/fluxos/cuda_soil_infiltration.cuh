// Copyright 2019, Diogo Costa
// CUDA kernel for soil infiltration loss (regular mesh)

#ifndef CUDA_SOIL_INFILTRATION_CUH
#define CUDA_SOIL_INFILTRATION_CUH

#ifdef USE_CUDA

#include "cuda_memory.h"

// Launch soil infiltration loss kernel on GPU (regular mesh)
void cuda_apply_soil_infiltration(CudaGridData& g);

#endif // USE_CUDA
#endif // CUDA_SOIL_INFILTRATION_CUH
