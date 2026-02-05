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
