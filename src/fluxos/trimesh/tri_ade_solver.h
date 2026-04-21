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
// Edge-based advection-dispersion equation solver for triangular mesh
// Upwind MUSCL advection + Fick's law dispersion

#ifndef TRI_ADE_SOLVER_H_INCLUDED
#define TRI_ADE_SOLVER_H_INCLUDED

#include "../GlobVar.h"
#include "TriMesh.h"
#include "TriSolution.h"

// ADE solver for triangular mesh
// Edge-based formulation: loop over edges, compute advective + dispersive flux
void tri_adesolver_calc(
    GlobVar& ds,
    const TriMesh& mesh,
    TriSolution& sol,
    int it,
    int ichem);

#endif // TRI_ADE_SOLVER_H_INCLUDED
