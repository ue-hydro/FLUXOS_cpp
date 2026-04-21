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
// Multi-phase hydrodynamics driver for triangular mesh
// 5-phase approach: wet/dry -> gradients -> limiters -> edge fluxes -> accumulate+update

#ifndef TRI_HYDRODYNAMICS_H_INCLUDED
#define TRI_HYDRODYNAMICS_H_INCLUDED

#include "../GlobVar.h"
#include "TriMesh.h"
#include "TriSolution.h"

// Main hydrodynamics timestep driver
// Performs: wet/dry check -> gradient -> limiter -> edge fluxes -> accumulate -> state update
void tri_hydrodynamics_calc(
    GlobVar& ds,
    const TriMesh& mesh,
    TriSolution& sol);

// Compute Courant condition for triangular mesh
// Returns minimum dt from CFL = cfl * inradius / (|u| + c)
// Also computes maximum water depth hpall
void tri_courant_condition(
    GlobVar& ds,
    const TriMesh& mesh,
    TriSolution& sol,
    double& dtfl_out,
    double& hpall_out);

#endif // TRI_HYDRODYNAMICS_H_INCLUDED
