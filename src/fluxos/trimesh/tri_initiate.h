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
// Initialization routines for triangular mesh

#ifndef TRI_INITIATE_H_INCLUDED
#define TRI_INITIATE_H_INCLUDED

#include "../GlobVar.h"
#include "TriMesh.h"
#include "TriSolution.h"
#include <fstream>

// Initialize triangular mesh solution:
// - Interpolate bed elevation from DEM to cell centroids
// - Set initial conditions
// - Apply boundary conditions from mesh tags
// - Find inflow cell
void tri_initiation(
    GlobVar& ds,
    TriMesh& mesh,
    TriSolution& sol,
    int nchem,
    std::ofstream& logFLUXOSfile);

// Apply boundary conditions based on mesh boundary tags
void tri_apply_boundary_conditions(
    GlobVar& ds,
    TriMesh& mesh,
    TriSolution& sol);

#endif // TRI_INITIATE_H_INCLUDED
