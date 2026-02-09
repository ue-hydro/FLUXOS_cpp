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
