// Copyright 2026, Diogo Costa, diogo.costa@uevora.pt
// This file is part of the FLUXOS model.

// This program, FLUXOS, is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
