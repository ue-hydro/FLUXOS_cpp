// Copyright 2019, Diogo Costa
// MPI Parallel I/O for FLUXOS
// Efficient parallel file output for HPC clusters

#ifndef MPI_IO_H_INCLUDED
#define MPI_IO_H_INCLUDED

#ifdef USE_MPI
#include <mpi.h>
#endif

#include <armadillo>
#include <string>
#include "GlobVar.h"
#include "mpi_domain.h"

// Parallel I/O helper class
class MPIParallelIO {
public:
    // Write results in parallel (each process writes its own portion)
    // Uses MPI-IO for efficient parallel file access
    static bool write_parallel_results(
        GlobVar& ds,
        unsigned int print_time,
        const std::string& output_folder);

    // Gather and write (traditional approach - root collects and writes)
    // Simpler but less scalable
    static bool write_gathered_results(
        GlobVar& ds,
        unsigned int print_time,
        const std::string& output_folder);

    // Write HDF5 output in parallel (if HDF5 parallel is available)
    static bool write_parallel_hdf5(
        GlobVar& ds,
        unsigned int print_time,
        const std::string& output_folder);

private:
    // Helper to write a single matrix using MPI-IO
    static bool write_matrix_mpi(
        const arma::Mat<double>& local_mat,
        const std::string& filename,
        unsigned int global_rows,
        unsigned int global_cols);
};

#endif // MPI_IO_H_INCLUDED
