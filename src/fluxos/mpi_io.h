// Copyright 2026, Diogo Costa, diogo.costa@uevora.pt
// This file is part of the FLUXOS model.

// This program, FLUXOS, is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
