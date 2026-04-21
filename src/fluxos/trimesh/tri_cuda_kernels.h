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
