// Copyright 2019, Diogo Costa
// CUDA GPU memory management for triangular mesh

#ifdef USE_CUDA

#include "tri_cuda_memory.h"
#include <cuda_runtime.h>
#include <iostream>
#include <cstring>

#define TRI_CUDA_CHECK(call) \
    do { \
        cudaError_t err = call; \
        if (err != cudaSuccess) { \
            std::cerr << "CUDA error in " << __FILE__ << ":" << __LINE__ \
                      << ": " << cudaGetErrorString(err) << std::endl; \
            exit(EXIT_FAILURE); \
        } \
    } while(0)

// Helper macros
#define ALLOC_D(ptr, count, type) TRI_CUDA_CHECK(cudaMalloc(&(ptr), (count) * sizeof(type)))
#define FREE_D(ptr) if (ptr) { cudaFree(ptr); ptr = nullptr; }
#define H2D(dst, src, count, type) TRI_CUDA_CHECK(cudaMemcpy(dst, src, (count)*sizeof(type), cudaMemcpyHostToDevice))
#define D2H(dst, src, count, type) TRI_CUDA_CHECK(cudaMemcpy(dst, src, (count)*sizeof(type), cudaMemcpyDeviceToHost))

// ============================================================================
// Allocate device memory and copy mesh topology
// ============================================================================
void TriCudaMemoryManager::allocate(const TriMesh& mesh, const TriSolution& sol)
{
    int nc = mesh.num_cells;
    int ne = mesh.num_edges;
    data.num_cells = nc;
    data.num_edges = ne;
    data.soil_allocated = false;

    // Print GPU info
    int device;
    cudaGetDevice(&device);
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop, device);
    std::cout << "\n=== Triangular Mesh CUDA Initialization ===" << std::endl;
    std::cout << "  GPU: " << prop.name << std::endl;
    std::cout << "  Cells: " << nc << "  Edges: " << ne << std::endl;

    // --- Allocate mesh topology ---
    ALLOC_D(data.d_edge_left, ne, int);
    ALLOC_D(data.d_edge_right, ne, int);
    ALLOC_D(data.d_edge_nx, ne, double);
    ALLOC_D(data.d_edge_ny, ne, double);
    ALLOC_D(data.d_edge_length, ne, double);
    ALLOC_D(data.d_edge_mx, ne, double);
    ALLOC_D(data.d_edge_my, ne, double);
    ALLOC_D(data.d_edge_dist_LR, ne, double);
    ALLOC_D(data.d_edge_boundary_tag, ne, int);

    ALLOC_D(data.d_cell_cx, nc, double);
    ALLOC_D(data.d_cell_cy, nc, double);
    ALLOC_D(data.d_cell_area, nc, double);
    ALLOC_D(data.d_cell_inradius, nc, double);
    ALLOC_D(data.d_cell_lsq_inv, nc * 4, double);
    ALLOC_D(data.d_cell_edges, nc * 3, int);
    ALLOC_D(data.d_cell_neighbors, nc * 3, int);

    // --- Allocate solution arrays ---
    ALLOC_D(data.d_z, nc, double);
    ALLOC_D(data.d_zb, nc, double);
    ALLOC_D(data.d_h, nc, double);
    ALLOC_D(data.d_qx, nc, double);
    ALLOC_D(data.d_qy, nc, double);
    ALLOC_D(data.d_ux, nc, double);
    ALLOC_D(data.d_uy, nc, double);
    ALLOC_D(data.d_us, nc, double);
    ALLOC_D(data.d_ks, nc, double);

    ALLOC_D(data.d_ldry, nc, float);
    ALLOC_D(data.d_ldry_prev, nc, float);
    ALLOC_D(data.d_innerNeumannBCWeir, nc, float);

    ALLOC_D(data.d_grad_z_x, nc, double);
    ALLOC_D(data.d_grad_z_y, nc, double);
    ALLOC_D(data.d_grad_qx_x, nc, double);
    ALLOC_D(data.d_grad_qx_y, nc, double);
    ALLOC_D(data.d_grad_qy_x, nc, double);
    ALLOC_D(data.d_grad_qy_y, nc, double);

    ALLOC_D(data.d_phi_z, nc, double);
    ALLOC_D(data.d_phi_qx, nc, double);
    ALLOC_D(data.d_phi_qy, nc, double);

    ALLOC_D(data.d_flux_mass, ne, double);
    ALLOC_D(data.d_flux_momx, ne, double);
    ALLOC_D(data.d_flux_momy, ne, double);

    ALLOC_D(data.d_dh, nc, double);
    ALLOC_D(data.d_dqx, nc, double);
    ALLOC_D(data.d_dqy, nc, double);

    ALLOC_D(data.d_h0, nc, double);
    ALLOC_D(data.d_twetimetracer, nc, double);
    ALLOC_D(data.d_conc_SW, nc, double);

    // Reduction buffer (enough for 1 block per 256 threads)
    int num_blocks = (nc + 255) / 256;
    ALLOC_D(data.d_block_reduce, num_blocks * 2, double);

    // --- Copy mesh topology to device ---
    // Pack edge data
    std::vector<int> h_edge_left(ne), h_edge_right(ne), h_edge_btag(ne);
    std::vector<double> h_edge_nx(ne), h_edge_ny(ne), h_edge_len(ne);
    std::vector<double> h_edge_mx(ne), h_edge_my(ne), h_edge_dist(ne);

    for (int ei = 0; ei < ne; ei++) {
        h_edge_left[ei] = mesh.edges[ei].left_cell;
        h_edge_right[ei] = mesh.edges[ei].right_cell;
        h_edge_nx[ei] = mesh.edges[ei].nx;
        h_edge_ny[ei] = mesh.edges[ei].ny;
        h_edge_len[ei] = mesh.edges[ei].length;
        h_edge_mx[ei] = mesh.edges[ei].mx;
        h_edge_my[ei] = mesh.edges[ei].my;
        h_edge_dist[ei] = mesh.edges[ei].dist_LR;
        h_edge_btag[ei] = mesh.edges[ei].boundary_tag;
    }

    H2D(data.d_edge_left, h_edge_left.data(), ne, int);
    H2D(data.d_edge_right, h_edge_right.data(), ne, int);
    H2D(data.d_edge_nx, h_edge_nx.data(), ne, double);
    H2D(data.d_edge_ny, h_edge_ny.data(), ne, double);
    H2D(data.d_edge_length, h_edge_len.data(), ne, double);
    H2D(data.d_edge_mx, h_edge_mx.data(), ne, double);
    H2D(data.d_edge_my, h_edge_my.data(), ne, double);
    H2D(data.d_edge_dist_LR, h_edge_dist.data(), ne, double);
    H2D(data.d_edge_boundary_tag, h_edge_btag.data(), ne, int);

    // Pack cell data
    std::vector<double> h_cx(nc), h_cy(nc), h_area(nc), h_inrad(nc);
    std::vector<double> h_lsq(nc * 4);
    std::vector<int> h_cell_edges(nc * 3), h_cell_nbrs(nc * 3);

    for (int ci = 0; ci < nc; ci++) {
        h_cx[ci] = mesh.cells[ci].cx;
        h_cy[ci] = mesh.cells[ci].cy;
        h_area[ci] = mesh.cells[ci].area;
        h_inrad[ci] = mesh.cells[ci].inradius;
        h_lsq[ci * 4 + 0] = mesh.cells[ci].lsq_inv[0][0];
        h_lsq[ci * 4 + 1] = mesh.cells[ci].lsq_inv[0][1];
        h_lsq[ci * 4 + 2] = mesh.cells[ci].lsq_inv[1][0];
        h_lsq[ci * 4 + 3] = mesh.cells[ci].lsq_inv[1][1];
        for (int e = 0; e < 3; e++) {
            h_cell_edges[ci * 3 + e] = mesh.cells[ci].edges[e];
            h_cell_nbrs[ci * 3 + e] = mesh.cells[ci].neighbors[e];
        }
    }

    H2D(data.d_cell_cx, h_cx.data(), nc, double);
    H2D(data.d_cell_cy, h_cy.data(), nc, double);
    H2D(data.d_cell_area, h_area.data(), nc, double);
    H2D(data.d_cell_inradius, h_inrad.data(), nc, double);
    H2D(data.d_cell_lsq_inv, h_lsq.data(), nc * 4, double);
    H2D(data.d_cell_edges, h_cell_edges.data(), nc * 3, int);
    H2D(data.d_cell_neighbors, h_cell_nbrs.data(), nc * 3, int);

    size_t total_mem = (ne * (sizeof(int) * 3 + sizeof(double) * 6) +
                        nc * (sizeof(double) * 26 + sizeof(float) * 3 + sizeof(int) * 6) +
                        ne * sizeof(double) * 3);
    std::cout << "  Device memory: " << total_mem / (1024.0 * 1024.0) << " MB" << std::endl;
    std::cout << "========================================\n" << std::endl;
}

// ============================================================================
// Copy solution to device
// ============================================================================
void TriCudaMemoryManager::copy_solution_to_device(const TriSolution& sol)
{
    int nc = data.num_cells;
    H2D(data.d_z, sol.z.data(), nc, double);
    H2D(data.d_zb, sol.zb.data(), nc, double);
    H2D(data.d_h, sol.h.data(), nc, double);
    H2D(data.d_qx, sol.qx.data(), nc, double);
    H2D(data.d_qy, sol.qy.data(), nc, double);
    H2D(data.d_ux, sol.ux.data(), nc, double);
    H2D(data.d_uy, sol.uy.data(), nc, double);
    H2D(data.d_us, sol.us.data(), nc, double);
    H2D(data.d_ks, sol.ks.data(), nc, double);
    H2D(data.d_ldry, sol.ldry.data(), nc, float);
    H2D(data.d_ldry_prev, sol.ldry_prev.data(), nc, float);
    H2D(data.d_innerNeumannBCWeir, sol.innerNeumannBCWeir.data(), nc, float);
    H2D(data.d_h0, sol.h0.data(), nc, double);
    H2D(data.d_twetimetracer, sol.twetimetracer.data(), nc, double);
}

// ============================================================================
// Copy solution to host
// ============================================================================
void TriCudaMemoryManager::copy_solution_to_host(TriSolution& sol)
{
    int nc = data.num_cells;
    D2H(sol.z.data(), data.d_z, nc, double);
    D2H(sol.zb.data(), data.d_zb, nc, double);
    D2H(sol.h.data(), data.d_h, nc, double);
    D2H(sol.qx.data(), data.d_qx, nc, double);
    D2H(sol.qy.data(), data.d_qy, nc, double);
    D2H(sol.ux.data(), data.d_ux, nc, double);
    D2H(sol.uy.data(), data.d_uy, nc, double);
    D2H(sol.us.data(), data.d_us, nc, double);
    D2H(sol.ldry.data(), data.d_ldry, nc, float);
    D2H(sol.ldry_prev.data(), data.d_ldry_prev, nc, float);
    D2H(sol.h0.data(), data.d_h0, nc, double);
    D2H(sol.twetimetracer.data(), data.d_twetimetracer, nc, double);
}

// ============================================================================
// Copy output-only fields
// ============================================================================
void TriCudaMemoryManager::copy_output_to_host(TriSolution& sol)
{
    int nc = data.num_cells;
    D2H(sol.z.data(), data.d_z, nc, double);
    D2H(sol.h.data(), data.d_h, nc, double);
    D2H(sol.qx.data(), data.d_qx, nc, double);
    D2H(sol.qy.data(), data.d_qy, nc, double);
    D2H(sol.ux.data(), data.d_ux, nc, double);
    D2H(sol.uy.data(), data.d_uy, nc, double);
    D2H(sol.us.data(), data.d_us, nc, double);
    D2H(sol.ldry.data(), data.d_ldry, nc, float);
    D2H(sol.twetimetracer.data(), data.d_twetimetracer, nc, double);
}

// ============================================================================
// Update scalars
// ============================================================================
void TriCudaMemoryManager::update_scalars(double hdry, double gacc, double cfl,
                                           double cvdef, double nuem, double dtfl)
{
    data.hdry = hdry;
    data.gacc = gacc;
    data.cfl = cfl;
    data.cvdef = cvdef;
    data.nuem = nuem;
    data.dtfl = dtfl;
}

// ============================================================================
// Concentration transfer
// ============================================================================
void TriCudaMemoryManager::copy_conc_to_device(const TriSolution& sol, int ichem)
{
    H2D(data.d_conc_SW, sol.conc_SW[ichem].data(), data.num_cells, double);
}

void TriCudaMemoryManager::copy_conc_to_host(TriSolution& sol, int ichem)
{
    D2H(sol.conc_SW[ichem].data(), data.d_conc_SW, data.num_cells, double);
}

// ============================================================================
// Allocate soil infiltration device arrays
// ============================================================================
void TriCudaMemoryManager::allocate_soil(int num_cells)
{
    int nc = num_cells;
    ALLOC_D(data.d_soil_infil_rate, nc, double);
    ALLOC_D(data.d_soil_Ks, nc, double);
    ALLOC_D(data.d_soil_f0, nc, double);
    ALLOC_D(data.d_soil_k, nc, double);
    ALLOC_D(data.d_soil_wetting_time, nc, double);
    data.soil_allocated = true;
    std::cout << "  Triangular CUDA: soil arrays allocated" << std::endl;
}

void TriCudaMemoryManager::copy_soil_to_device(const TriSolution& sol)
{
    if (!data.soil_allocated) return;
    int nc = data.num_cells;
    H2D(data.d_soil_infil_rate, sol.soil_infil_rate.data(), nc, double);
    H2D(data.d_soil_Ks, sol.soil_Ks.data(), nc, double);
    H2D(data.d_soil_f0, sol.soil_f0.data(), nc, double);
    H2D(data.d_soil_k, sol.soil_k.data(), nc, double);
    H2D(data.d_soil_wetting_time, sol.soil_wetting_time.data(), nc, double);
}

void TriCudaMemoryManager::copy_soil_to_host(TriSolution& sol)
{
    if (!data.soil_allocated) return;
    int nc = data.num_cells;
    D2H(sol.soil_infil_rate.data(), data.d_soil_infil_rate, nc, double);
    D2H(sol.soil_wetting_time.data(), data.d_soil_wetting_time, nc, double);
}

// ============================================================================
// Deallocate
// ============================================================================
void TriCudaMemoryManager::deallocate()
{
    FREE_D(data.d_edge_left); FREE_D(data.d_edge_right);
    FREE_D(data.d_edge_nx); FREE_D(data.d_edge_ny);
    FREE_D(data.d_edge_length); FREE_D(data.d_edge_mx); FREE_D(data.d_edge_my);
    FREE_D(data.d_edge_dist_LR); FREE_D(data.d_edge_boundary_tag);
    FREE_D(data.d_cell_cx); FREE_D(data.d_cell_cy);
    FREE_D(data.d_cell_area); FREE_D(data.d_cell_inradius);
    FREE_D(data.d_cell_lsq_inv); FREE_D(data.d_cell_edges); FREE_D(data.d_cell_neighbors);
    FREE_D(data.d_z); FREE_D(data.d_zb); FREE_D(data.d_h);
    FREE_D(data.d_qx); FREE_D(data.d_qy);
    FREE_D(data.d_ux); FREE_D(data.d_uy); FREE_D(data.d_us); FREE_D(data.d_ks);
    FREE_D(data.d_ldry); FREE_D(data.d_ldry_prev); FREE_D(data.d_innerNeumannBCWeir);
    FREE_D(data.d_grad_z_x); FREE_D(data.d_grad_z_y);
    FREE_D(data.d_grad_qx_x); FREE_D(data.d_grad_qx_y);
    FREE_D(data.d_grad_qy_x); FREE_D(data.d_grad_qy_y);
    FREE_D(data.d_phi_z); FREE_D(data.d_phi_qx); FREE_D(data.d_phi_qy);
    FREE_D(data.d_flux_mass); FREE_D(data.d_flux_momx); FREE_D(data.d_flux_momy);
    FREE_D(data.d_dh); FREE_D(data.d_dqx); FREE_D(data.d_dqy);
    FREE_D(data.d_h0); FREE_D(data.d_twetimetracer);
    FREE_D(data.d_conc_SW); FREE_D(data.d_block_reduce);

    // Free soil arrays
    if (data.soil_allocated) {
        FREE_D(data.d_soil_infil_rate);
        FREE_D(data.d_soil_Ks);
        FREE_D(data.d_soil_f0);
        FREE_D(data.d_soil_k);
        FREE_D(data.d_soil_wetting_time);
        data.soil_allocated = false;
    }
}

#undef TRI_CUDA_CHECK
#undef ALLOC_D
#undef FREE_D
#undef H2D
#undef D2H

#endif // USE_CUDA
