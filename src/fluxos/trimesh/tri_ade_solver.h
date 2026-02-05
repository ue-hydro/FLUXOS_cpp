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
