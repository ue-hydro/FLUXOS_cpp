// Copyright 2019, Diogo Costa
// Soil infiltration loss for triangular mesh

#ifndef TRI_SOIL_INFILTRATION_H_INCLUDED
#define TRI_SOIL_INFILTRATION_H_INCLUDED

#include "../GlobVar.h"
#include "../soil_infiltration.h"
#include "TriMesh.h"
#include "TriSolution.h"

// ============================================================================
// Initialize soil parameters for triangular mesh cells
// Samples the soil map raster at cell centroids (nearest-neighbor)
// ============================================================================
bool tri_init_soil(
    GlobVar& ds,
    const TriMesh& mesh,
    TriSolution& sol,
    const SoilConfig& soil_config,
    std::ofstream& logFLUXOSfile);

// ============================================================================
// Apply soil infiltration loss for triangular mesh (one timestep)
// Subtracts min(Ks * dt, h_surface) from surface water
// Called after forcing, before hydrodynamic solver
// ============================================================================
void tri_apply_soil_infiltration(
    GlobVar& ds,
    const TriMesh& mesh,
    TriSolution& sol,
    const SoilConfig& soil_config);

#endif // TRI_SOIL_INFILTRATION_H_INCLUDED
