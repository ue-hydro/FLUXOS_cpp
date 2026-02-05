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
