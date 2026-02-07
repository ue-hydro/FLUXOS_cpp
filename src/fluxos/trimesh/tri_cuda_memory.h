// Copyright 2019, Diogo Costa
// CUDA GPU memory management for triangular mesh
// Flat array device allocation and transfer

#ifndef TRI_CUDA_MEMORY_H_INCLUDED
#define TRI_CUDA_MEMORY_H_INCLUDED

#ifdef USE_CUDA

#include "TriMesh.h"
#include "TriSolution.h"

// ============================================================================
// Device data for triangular mesh: flat arrays
// ============================================================================
struct TriCudaData {
    int num_cells;
    int num_edges;

    // --- Mesh topology (read-only on device) ---
    // Edge connectivity
    int* d_edge_left;       // left cell for each edge
    int* d_edge_right;      // right cell for each edge (-1 = boundary)
    double* d_edge_nx;      // edge normal x
    double* d_edge_ny;      // edge normal y
    double* d_edge_length;  // edge length
    double* d_edge_mx;      // edge midpoint x
    double* d_edge_my;      // edge midpoint y
    double* d_edge_dist_LR; // distance between cell centroids
    int* d_edge_boundary_tag; // boundary tag

    // Cell geometry
    double* d_cell_cx;      // centroid x
    double* d_cell_cy;      // centroid y
    double* d_cell_area;    // cell area
    double* d_cell_inradius; // inscribed circle radius
    double* d_cell_lsq_inv; // LSQ inverse (4 doubles per cell: [a00,a01,a10,a11])

    // Cell connectivity: 3 edges and 3 neighbors per cell
    int* d_cell_edges;      // [3 * num_cells]: cell_edges[ci*3+e]
    int* d_cell_neighbors;  // [3 * num_cells]: cell_neighbors[ci*3+e]

    // --- Solution data (read-write on device) ---
    double* d_z;
    double* d_zb;
    double* d_h;
    double* d_qx;
    double* d_qy;
    double* d_ux;
    double* d_uy;
    double* d_us;
    double* d_ks;

    float* d_ldry;
    float* d_ldry_prev;
    float* d_innerNeumannBCWeir;

    // Gradients
    double* d_grad_z_x;
    double* d_grad_z_y;
    double* d_grad_qx_x;
    double* d_grad_qx_y;
    double* d_grad_qy_x;
    double* d_grad_qy_y;

    // Limiters
    double* d_phi_z;
    double* d_phi_qx;
    double* d_phi_qy;

    // Edge fluxes
    double* d_flux_mass;
    double* d_flux_momx;
    double* d_flux_momy;

    // Accumulators
    double* d_dh;
    double* d_dqx;
    double* d_dqy;

    // ADE
    double* d_h0;
    double* d_twetimetracer;
    double* d_conc_SW;      // single species buffer

    // Reduction buffer
    double* d_block_reduce;

    // Soil infiltration device arrays (Horton decay model)
    double* d_soil_infil_rate;
    double* d_soil_Ks;
    double* d_soil_f0;
    double* d_soil_k;
    double* d_soil_wetting_time;
    bool soil_allocated;

    // Scalar constants
    double hdry, gacc, cfl, cvdef, nuem, dtfl;
};

// ============================================================================
// TriCudaMemoryManager: manages device memory for triangular mesh
// ============================================================================
class TriCudaMemoryManager {
public:
    TriCudaData data;

    // Allocate device memory and copy mesh topology
    void allocate(const TriMesh& mesh, const TriSolution& sol);

    // Copy solution fields host -> device
    void copy_solution_to_device(const TriSolution& sol);

    // Copy solution fields device -> host
    void copy_solution_to_host(TriSolution& sol);

    // Copy output-only fields device -> host (z, h, qx, qy, ux, uy, us)
    void copy_output_to_host(TriSolution& sol);

    // Update scalar constants
    void update_scalars(double hdry, double gacc, double cfl,
                        double cvdef, double nuem, double dtfl);

    // Copy single concentration species to/from device
    void copy_conc_to_device(const TriSolution& sol, int ichem);
    void copy_conc_to_host(TriSolution& sol, int ichem);

    // Soil infiltration GPU memory
    void allocate_soil(int num_cells);
    void copy_soil_to_device(const TriSolution& sol);
    void copy_soil_to_host(TriSolution& sol);

    // Deallocate all device memory
    void deallocate();
};

#endif // USE_CUDA
#endif // TRI_CUDA_MEMORY_H_INCLUDED
