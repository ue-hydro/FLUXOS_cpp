// Copyright 2019, Diogo Costa
// Flat-array solution storage for triangular mesh

#include "TriSolution.h"
#include "TriMesh.h"

// ============================================================================
// Allocate all arrays
// ============================================================================
void TriSolution::allocate(int ncells, int nedges, int num_chem)
{
    num_cells = ncells;
    num_edges = nedges;
    nchem = num_chem;

    // State variables (per cell)
    z.assign(ncells, 0.0);
    zb.assign(ncells, 0.0);
    h.assign(ncells, 0.0);
    qx.assign(ncells, 0.0);
    qy.assign(ncells, 0.0);
    ux.assign(ncells, 0.0);
    uy.assign(ncells, 0.0);
    us.assign(ncells, 0.0);
    ks.assign(ncells, 0.0);

    // Wet/dry flags
    ldry.assign(ncells, 1.0f);
    ldry_prev.assign(ncells, 1.0f);
    innerNeumannBCWeir.assign(ncells, 0.0f);

    // Gradient fields
    grad_z_x.assign(ncells, 0.0);
    grad_z_y.assign(ncells, 0.0);
    grad_qx_x.assign(ncells, 0.0);
    grad_qx_y.assign(ncells, 0.0);
    grad_qy_x.assign(ncells, 0.0);
    grad_qy_y.assign(ncells, 0.0);

    // Limiter values
    phi_z.assign(ncells, 1.0);
    phi_qx.assign(ncells, 1.0);
    phi_qy.assign(ncells, 1.0);

    // Edge fluxes
    flux_mass.assign(nedges, 0.0);
    flux_momx.assign(nedges, 0.0);
    flux_momy.assign(nedges, 0.0);

    // Accumulators
    dh.assign(ncells, 0.0);
    dqx.assign(ncells, 0.0);
    dqy.assign(ncells, 0.0);

    // ADE / WINTRA
    h0.assign(ncells, 0.0);
    soil_mass.assign(ncells, 0.0);
    twetimetracer.assign(ncells, 0.0);

    conc_SW.resize(num_chem);
    for (int ic = 0; ic < num_chem; ic++) {
        conc_SW[ic].assign(ncells, 0.0);
    }
}

// ============================================================================
// Zero state fields
// ============================================================================
void TriSolution::zero_state()
{
    std::fill(z.begin(), z.end(), 0.0);
    std::fill(h.begin(), h.end(), 0.0);
    std::fill(qx.begin(), qx.end(), 0.0);
    std::fill(qy.begin(), qy.end(), 0.0);
    std::fill(ux.begin(), ux.end(), 0.0);
    std::fill(uy.begin(), uy.end(), 0.0);
    std::fill(us.begin(), us.end(), 0.0);
    std::fill(ldry.begin(), ldry.end(), 1.0f);
    std::fill(ldry_prev.begin(), ldry_prev.end(), 1.0f);
}

// ============================================================================
// Zero accumulators and fluxes only
// ============================================================================
void TriSolution::zero_accumulators()
{
    std::fill(dh.begin(), dh.end(), 0.0);
    std::fill(dqx.begin(), dqx.end(), 0.0);
    std::fill(dqy.begin(), dqy.end(), 0.0);
    std::fill(flux_mass.begin(), flux_mass.end(), 0.0);
    std::fill(flux_momx.begin(), flux_momx.end(), 0.0);
    std::fill(flux_momy.begin(), flux_momy.end(), 0.0);
    std::fill(grad_z_x.begin(), grad_z_x.end(), 0.0);
    std::fill(grad_z_y.begin(), grad_z_y.end(), 0.0);
    std::fill(grad_qx_x.begin(), grad_qx_x.end(), 0.0);
    std::fill(grad_qx_y.begin(), grad_qx_y.end(), 0.0);
    std::fill(grad_qy_x.begin(), grad_qy_x.end(), 0.0);
    std::fill(grad_qy_y.begin(), grad_qy_y.end(), 0.0);
    std::fill(phi_z.begin(), phi_z.end(), 1.0);
    std::fill(phi_qx.begin(), phi_qx.end(), 1.0);
    std::fill(phi_qy.begin(), phi_qy.end(), 1.0);
}

// ============================================================================
// Clear all memory
// ============================================================================
void TriSolution::clear()
{
    z.clear(); zb.clear(); h.clear();
    qx.clear(); qy.clear();
    ux.clear(); uy.clear(); us.clear(); ks.clear();
    ldry.clear(); ldry_prev.clear(); innerNeumannBCWeir.clear();
    grad_z_x.clear(); grad_z_y.clear();
    grad_qx_x.clear(); grad_qx_y.clear();
    grad_qy_x.clear(); grad_qy_y.clear();
    phi_z.clear(); phi_qx.clear(); phi_qy.clear();
    flux_mass.clear(); flux_momx.clear(); flux_momy.clear();
    dh.clear(); dqx.clear(); dqy.clear();
    h0.clear(); soil_mass.clear(); twetimetracer.clear();
    conc_SW.clear();
    num_cells = 0;
    num_edges = 0;
    nchem = 0;
}
