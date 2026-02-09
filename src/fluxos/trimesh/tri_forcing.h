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
