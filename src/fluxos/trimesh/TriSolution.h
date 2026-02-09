// Copyright 2019, Diogo Costa
// Flat-array solution storage for triangular mesh
// All arrays indexed by cell ID (not arma::Mat)

#ifndef TRISOLUTION_H_INCLUDED
#define TRISOLUTION_H_INCLUDED

#include <vector>
#include <string>

// Forward declaration
class TriMesh;

// ============================================================================
// TriSolution: solution data on triangular mesh
// All fields stored as flat std::vector<double> indexed by cell ID
// ============================================================================
class TriSolution {
public:
    // ---- State variables (per cell) ----
    std::vector<double> z;          // water surface elevation [m]
    std::vector<double> zb;         // bed elevation [m]
    std::vector<double> h;          // water depth [m]
    std::vector<double> qx;         // x-momentum (discharge per unit width) [m^2/s]
    std::vector<double> qy;         // y-momentum (discharge per unit width) [m^2/s]
    std::vector<double> ux;         // x-velocity [m/s]
    std::vector<double> uy;         // y-velocity [m/s]
    std::vector<double> us;         // shear stress velocity [m/s]
    std::vector<double> ks;         // roughness height [m]

    // ---- Wet/dry flag (per cell) ----
    std::vector<float> ldry;        // 1.0 = dry, 0.0 = wet
    std::vector<float> ldry_prev;   // previous timestep dry flag
    std::vector<float> innerNeumannBCWeir; // inner boundary condition flag

    // ---- Gradient fields (per cell) ----
    std::vector<double> grad_z_x, grad_z_y;
    std::vector<double> grad_qx_x, grad_qx_y;
    std::vector<double> grad_qy_x, grad_qy_y;

    // ---- Limiter values (per cell) ----
    std::vector<double> phi_z;      // Barth-Jespersen limiter for z
    std::vector<double> phi_qx;     // Barth-Jespersen limiter for qx
    std::vector<double> phi_qy;     // Barth-Jespersen limiter for qy

    // ---- Edge fluxes (per edge) ----
    std::vector<double> flux_mass;  // mass flux through edge [m^3/s]
    std::vector<double> flux_momx;  // x-momentum flux [m^3/s * m/s]
    std::vector<double> flux_momy;  // y-momentum flux [m^3/s * m/s]

    // ---- Hydrostatic pressure source terms (per edge) ----
    // Bed-step correction from Audusse et al. hydrostatic reconstruction
    std::vector<double> src_pcorr_L;  // 0.5*g*(h_orig^2 - h_recon^2) left side
    std::vector<double> src_pcorr_R;  // 0.5*g*(h_orig^2 - h_recon^2) right side

    // ---- Accumulators (per cell) ----
    std::vector<double> dh;         // change in water depth
    std::vector<double> dqx;        // change in x-momentum
    std::vector<double> dqy;        // change in y-momentum

    // ---- ADE fields (per cell) ----
    std::vector<double> h0;             // depth at start of timestep
    std::vector<double> twetimetracer;  // wetting time tracer [hours]
    std::vector<std::vector<double>> conc_SW;  // concentration [nchem][ncells]

    // ---- Soil infiltration (per cell, Horton decay model) ----
    std::vector<double> soil_infil_rate;   // instantaneous infiltration rate [m/s]
    std::vector<double> soil_Ks;           // saturated hydraulic conductivity [m/s] = final rate (fc)
    std::vector<double> soil_f0;           // initial infiltration rate [m/s]
    std::vector<double> soil_k;            // Horton decay constant [1/s]
    std::vector<double> soil_wetting_time; // cumulative wetting time [s]
    std::vector<int>    soil_type_id;      // soil type ID (for output/diagnostics)

    // ---- Metadata ----
    int num_cells = 0;
    int num_edges = 0;
    int nchem = 0;

    // ---- Methods ----

    // Allocate all arrays for given mesh dimensions
    void allocate(int ncells, int nedges, int num_chem = 1);

    // Zero all solution fields (but not bed elevation or roughness)
    void zero_state();

    // Zero accumulators and fluxes only
    void zero_accumulators();

    // Clear all memory
    void clear();
};

#endif // TRISOLUTION_H_INCLUDED
