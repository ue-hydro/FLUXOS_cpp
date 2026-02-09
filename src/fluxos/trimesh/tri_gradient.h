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
