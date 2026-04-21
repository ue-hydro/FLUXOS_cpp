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
// Forcing routines for triangular mesh: meteo and inflow

#ifndef TRI_FORCING_H_INCLUDED
#define TRI_FORCING_H_INCLUDED

#include "../GlobVar.h"
#include "TriMesh.h"
#include "TriSolution.h"

// Add meteo forcing (rainfall/snowmelt) to all cells
bool tri_add_meteo(
    GlobVar& ds,
    const TriMesh& mesh,
    TriSolution& sol,
    int nchem);

// Add inflow at discharge point
bool tri_add_inflow(
    GlobVar& ds,
    const TriMesh& mesh,
    TriSolution& sol,
    int nchem);

#endif // TRI_FORCING_H_INCLUDED
