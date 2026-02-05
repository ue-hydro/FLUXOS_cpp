// Copyright 2019, Diogo Costa
// WINTRA soil release solver for triangular mesh
// Cell-based soil mass release using per-cell area

#ifndef TRI_WINTRA_SOLVER_H_INCLUDED
#define TRI_WINTRA_SOLVER_H_INCLUDED

#include "../GlobVar.h"
#include "TriMesh.h"
#include "TriSolution.h"

// WINTRA solver for triangular mesh
void tri_wintrasolver_calc(
    GlobVar& ds,
    const TriMesh& mesh,
    TriSolution& sol,
    int ichem);

#endif // TRI_WINTRA_SOLVER_H_INCLUDED
