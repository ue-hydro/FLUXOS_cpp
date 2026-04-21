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
// Least-squares gradient reconstruction + Barth-Jespersen limiter
// for second-order MUSCL on unstructured triangular mesh

#ifndef TRI_GRADIENT_H_INCLUDED
#define TRI_GRADIENT_H_INCLUDED

#include "TriMesh.h"
#include "TriSolution.h"

// Compute least-squares gradients for z, qx, qy
// Uses precomputed LSQ inverse matrices from TriMesh
void tri_compute_gradients(
    const TriMesh& mesh,
    TriSolution& sol);

// Apply Barth-Jespersen limiter to gradients
// Ensures monotonicity: reconstructed values at edge midpoints
// do not exceed neighbor min/max values
void tri_apply_limiter(
    const TriMesh& mesh,
    TriSolution& sol);

#endif // TRI_GRADIENT_H_INCLUDED
