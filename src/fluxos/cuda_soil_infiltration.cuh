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
// CUDA kernel for soil infiltration loss (regular mesh)

#ifndef CUDA_SOIL_INFILTRATION_CUH
#define CUDA_SOIL_INFILTRATION_CUH

#ifdef USE_CUDA

#include "cuda_memory.h"

// Launch soil infiltration loss kernel on GPU (regular mesh)
void cuda_apply_soil_infiltration(CudaGridData& g);

#endif // USE_CUDA
#endif // CUDA_SOIL_INFILTRATION_CUH
