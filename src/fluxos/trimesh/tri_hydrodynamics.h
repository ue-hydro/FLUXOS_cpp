// Copyright 2019, Diogo Costa
// Multi-phase hydrodynamics driver for triangular mesh
// 5-phase approach: wet/dry -> gradients -> limiters -> edge fluxes -> accumulate+update

#ifndef TRI_HYDRODYNAMICS_H_INCLUDED
#define TRI_HYDRODYNAMICS_H_INCLUDED

#include "../GlobVar.h"
#include "TriMesh.h"
#include "TriSolution.h"

// Main hydrodynamics timestep driver
// Performs: wet/dry check -> gradient -> limiter -> edge fluxes -> accumulate -> state update
void tri_hydrodynamics_calc(
    GlobVar& ds,
    const TriMesh& mesh,
    TriSolution& sol);

// Compute Courant condition for triangular mesh
// Returns minimum dt from CFL = cfl * inradius / (|u| + c)
// Also computes maximum water depth hpall
void tri_courant_condition(
    GlobVar& ds,
    const TriMesh& mesh,
    TriSolution& sol,
    double& dtfl_out,
    double& hpall_out);

#endif // TRI_HYDRODYNAMICS_H_INCLUDED
