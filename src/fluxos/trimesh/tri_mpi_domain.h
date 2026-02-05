// Copyright 2019, Diogo Costa
// MPI domain decomposition for unstructured triangular mesh
// METIS graph partitioning + halo exchange

#ifndef TRI_MPI_DOMAIN_H_INCLUDED
#define TRI_MPI_DOMAIN_H_INCLUDED

#include "TriMesh.h"
#include "TriSolution.h"
#include <vector>

#ifdef USE_MPI
#include <mpi.h>
#endif

// ============================================================================
// TriMPIDomain: manages MPI decomposition for unstructured mesh
// ============================================================================
class TriMPIDomain {
public:
    int rank = 0;
    int size = 1;
    bool is_root = true;

    // Partition info
    std::vector<int> cell_partition;    // partition ID for each cell
    int local_num_cells = 0;
    int local_num_edges = 0;

    // Local-to-global and global-to-local maps
    std::vector<int> local_to_global_cell;
    std::vector<int> global_to_local_cell;

    // Halo cells: cells owned by other ranks that share edges with local cells
    std::vector<int> halo_cells;        // global IDs of halo cells

    // Communication: neighbor ranks and their halo data
    struct NeighborComm {
        int rank;
        std::vector<int> send_cells;    // local cell IDs to send to this neighbor
        std::vector<int> recv_cells;    // local halo cell IDs to receive from this neighbor
    };
    std::vector<NeighborComm> neighbors;

    // ---- Methods ----

    // Initialize MPI (if USE_MPI is defined)
    void init();

    // Partition the mesh using METIS (or naive block partitioning if METIS unavailable)
    void partition_mesh(const TriMesh& mesh);

    // Exchange halo cell data for a given field
    void exchange_halo(std::vector<double>& field);
    void exchange_halo_float(std::vector<float>& field);

    // Global reductions
    double global_min(double local_val);
    double global_max(double local_val);
    double global_sum(double local_val);

    // Finalize MPI
    void finalize();
};

#endif // TRI_MPI_DOMAIN_H_INCLUDED
