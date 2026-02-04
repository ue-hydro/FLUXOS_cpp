// Copyright 2019, Diogo Costa
// MPI Domain Decomposition for FLUXOS
// Hybrid MPI+OpenMP implementation for HPC clusters

#ifndef MPI_DOMAIN_H_INCLUDED
#define MPI_DOMAIN_H_INCLUDED

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <armadillo>
#include <vector>
#include <memory>

// Domain decomposition class for 2D Cartesian grid
class MPIDomain {
public:
    // MPI info
    int world_rank;         // Global rank
    int world_size;         // Total number of processes
    int cart_rank;          // Rank in Cartesian communicator
    int coords[2];          // Process coordinates in 2D grid [row, col]
    int dims[2];            // Process grid dimensions
    int neighbors[4];       // Neighbor ranks: [NORTH, SOUTH, EAST, WEST]

#ifdef USE_MPI
    MPI_Comm cart_comm;     // Cartesian communicator
    MPI_Datatype col_type_double;  // MPI datatype for column exchange (double)
    MPI_Datatype col_type_float;   // MPI datatype for column exchange (float)
#endif

    // Local domain info
    unsigned int local_nrows;    // Local rows (excluding ghost cells)
    unsigned int local_ncols;    // Local cols (excluding ghost cells)
    unsigned int global_nrows;   // Global rows
    unsigned int global_ncols;   // Global cols
    unsigned int ghost_width;    // Ghost cell width (typically 1-2)

    // Global offsets for this process
    unsigned int row_offset;     // Starting row in global domain
    unsigned int col_offset;     // Starting col in global domain

    // Direction constants
    static constexpr int NORTH = 0;
    static constexpr int SOUTH = 1;
    static constexpr int EAST = 2;
    static constexpr int WEST = 3;

    MPIDomain();
    ~MPIDomain();

    // Initialize MPI and create domain decomposition
    void initialize(int argc, char* argv[], unsigned int global_rows, unsigned int global_cols, unsigned int ghost = 1);

    // Finalize MPI
    void finalize();

    // Check if this is the root process
    bool is_root() const { return world_rank == 0; }

    // Exchange ghost cells for a double matrix
    void exchange_ghost_cells(arma::Mat<double>& matrix);

    // Exchange ghost cells for a float matrix
    void exchange_ghost_cells(arma::Mat<float>& matrix);

    // Global reduction operations
    double global_min(double local_value);
    double global_max(double local_value);
    double global_sum(double local_value);

    // Gather full domain to root (for output)
    void gather_to_root(const arma::Mat<double>& local_mat, arma::Mat<double>& global_mat);

    // Scatter from root to all processes (for input)
    void scatter_from_root(const arma::Mat<double>& global_mat, arma::Mat<double>& local_mat);

    // Synchronization barrier
    void barrier();

    // Print domain info (debug)
    void print_info();

private:
    bool initialized;

    // Helper to compute local domain size
    void compute_local_sizes();

    // Create MPI datatypes for non-contiguous data
    void create_mpi_types();
};

// Global MPI domain instance
extern MPIDomain mpi_domain;

#endif // MPI_DOMAIN_H_INCLUDED
