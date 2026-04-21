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
// Rotated Roe flux solver with MUSCL reconstruction + hydrostatic reconstruction
// for wet cells on unstructured triangular mesh

#ifndef TRI_SOLVER_WET_H_INCLUDED
#define TRI_SOLVER_WET_H_INCLUDED

#include "TriMesh.h"
#include "TriSolution.h"

// Compute Roe flux for a single edge between two wet cells
// Stores result in sol.flux_mass[edge_id], flux_momx, flux_momy
// Uses MUSCL reconstruction with Barth-Jespersen limiter
// and hydrostatic reconstruction for well-balancedness
void tri_edge_flux_wet(
    const TriMesh& mesh,
    TriSolution& sol,
    int edge_id,
    double gacc,
    double hdry,
    double cvdef,
    double nuem,
    double dtfl);

// Compute all edge fluxes (wet domain)
void tri_compute_edge_fluxes(
    const TriMesh& mesh,
    TriSolution& sol,
    double gacc,
    double hdry,
    double cvdef,
    double nuem,
    double dtfl);

#endif // TRI_SOLVER_WET_H_INCLUDED
