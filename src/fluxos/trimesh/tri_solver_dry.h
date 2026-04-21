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
// Dry-front Ritter solution in edge-normal frame for triangular mesh

#ifndef TRI_SOLVER_DRY_H_INCLUDED
#define TRI_SOLVER_DRY_H_INCLUDED

#include "TriMesh.h"
#include "TriSolution.h"

// Compute Ritter dry-front flux for an edge where one side is dry
// This is called instead of the full Roe solver when one cell is wet
// and the neighbor is dry (wetting front)
void tri_edge_flux_dry(
    const TriMesh& mesh,
    TriSolution& sol,
    int edge_id,
    double gacc,
    double hdry,
    double dtfl);

#endif // TRI_SOLVER_DRY_H_INCLUDED
