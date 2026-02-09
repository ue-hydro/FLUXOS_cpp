// Copyright 2019, Diogo Costa
// MPI Parallel I/O for FLUXOS
// Efficient parallel file output for HPC clusters

#include "mpi_io.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>

bool MPIParallelIO::write_gathered_results(
    GlobVar& ds,
    unsigned int print_time,
    const std::string& output_folder)
{
#ifdef USE_MPI
    // Gather matrices to root process
    arma::Mat<double> global_h, global_qx, global_qy, global_z;

    // Gather water depth
    mpi_domain.gather_to_root(*ds.h, global_h);
    mpi_domain.gather_to_root(*ds.qx, global_qx);
    mpi_domain.gather_to_root(*ds.qy, global_qy);
    mpi_domain.gather_to_root(*ds.z, global_z);

    // Only root writes the file
    if (mpi_domain.is_root()) {
        std::stringstream filename;
        filename << output_folder << "/" << print_time << ".txt";

        std::ofstream outfile(filename.str());
        if (!outfile.is_open()) {
            std::cerr << "Error: Could not open output file: " << filename.str() << std::endl;
            return false;
        }

        outfile << std::fixed << std::setprecision(6);

        // Write header
        outfile << "# FLUXOS Output - Time: " << print_time << " seconds\n";
        outfile << "# Columns: row, col, h, qx, qy, z\n";

        // Write data (excluding ghost cells)
        for (unsigned int j = 1; j <= mpi_domain.global_ncols; j++) {
            for (unsigned int i = 1; i <= mpi_domain.global_nrows; i++) {
                outfile << i << " " << j << " "
                        << global_h(i, j) << " "
                        << global_qx(i, j) << " "
                        << global_qy(i, j) << " "
                        << global_z(i, j) << "\n";
            }
        }

        outfile.close();
    }

    // Synchronize all processes
    mpi_domain.barrier();
    return true;

#else
    // Non-MPI fallback - direct write
    std::stringstream filename;
    filename << output_folder << "/" << print_time << ".txt";

    std::ofstream outfile(filename.str());
    if (!outfile.is_open()) {
        return false;
    }

    outfile << std::fixed << std::setprecision(6);
    outfile << "# FLUXOS Output - Time: " << print_time << " seconds\n";
    outfile << "# Columns: row, col, h, qx, qy, z\n";

    for (unsigned int j = 1; j <= ds.NCOLS; j++) {
        for (unsigned int i = 1; i <= ds.NROWS; i++) {
            outfile << i << " " << j << " "
                    << (*ds.h)(i, j) << " "
                    << (*ds.qx)(i, j) << " "
                    << (*ds.qy)(i, j) << " "
                    << (*ds.z)(i, j) << "\n";
        }
    }

    outfile.close();
    return true;
#endif
}

bool MPIParallelIO::write_parallel_results(
    GlobVar& ds,
    unsigned int print_time,
    const std::string& output_folder)
{
#ifdef USE_MPI
    // Each process writes its own file portion
    // This is more scalable for very large domains

    std::stringstream filename;
    filename << output_folder << "/" << print_time << "_rank" << mpi_domain.cart_rank << ".txt";

    std::ofstream outfile(filename.str());
    if (!outfile.is_open()) {
        std::cerr << "Rank " << mpi_domain.cart_rank
                  << " Error: Could not open output file: " << filename.str() << std::endl;
        return false;
    }

    outfile << std::fixed << std::setprecision(6);

    // Write header with domain info
    outfile << "# FLUXOS Output - Time: " << print_time << " seconds\n";
    outfile << "# Rank: " << mpi_domain.cart_rank
            << " of " << mpi_domain.world_size << "\n";
    outfile << "# Local domain: rows [" << mpi_domain.row_offset << "-"
            << (mpi_domain.row_offset + mpi_domain.local_nrows - 1) << "], cols ["
            << mpi_domain.col_offset << "-"
            << (mpi_domain.col_offset + mpi_domain.local_ncols - 1) << "]\n";
    outfile << "# Columns: global_row, global_col, h, qx, qy, z\n";

    const unsigned int g = mpi_domain.ghost_width;

    // Write local data with global coordinates
    for (unsigned int j = 1; j <= mpi_domain.local_ncols; j++) {
        for (unsigned int i = 1; i <= mpi_domain.local_nrows; i++) {
            unsigned int global_row = mpi_domain.row_offset + i;
            unsigned int global_col = mpi_domain.col_offset + j;

            outfile << global_row << " " << global_col << " "
                    << (*ds.h)(i, j) << " "
                    << (*ds.qx)(i, j) << " "
                    << (*ds.qy)(i, j) << " "
                    << (*ds.z)(i, j) << "\n";
        }
    }

    outfile.close();

    // Synchronize
    mpi_domain.barrier();

    // Root process creates a manifest file listing all output files
    if (mpi_domain.is_root()) {
        std::stringstream manifest_filename;
        manifest_filename << output_folder << "/" << print_time << "_manifest.txt";

        std::ofstream manifest(manifest_filename.str());
        manifest << "# FLUXOS Parallel Output Manifest\n";
        manifest << "# Time: " << print_time << " seconds\n";
        manifest << "# Total processes: " << mpi_domain.world_size << "\n";
        manifest << "# Global domain: " << mpi_domain.global_nrows << " x "
                 << mpi_domain.global_ncols << "\n";
        manifest << "# Files:\n";

        for (int p = 0; p < mpi_domain.world_size; p++) {
            manifest << print_time << "_rank" << p << ".txt\n";
        }
        manifest.close();
    }

    return true;

#else
    // Non-MPI fallback
    return write_gathered_results(ds, print_time, output_folder);
#endif
}

bool MPIParallelIO::write_matrix_mpi(
    const arma::Mat<double>& local_mat,
    const std::string& filename,
    unsigned int global_rows,
    unsigned int global_cols)
{
#ifdef USE_MPI
    MPI_File fh;
    MPI_Status status;

    // Open file collectively
    int rc = MPI_File_open(mpi_domain.cart_comm, filename.c_str(),
                           MPI_MODE_CREATE | MPI_MODE_WRONLY,
                           MPI_INFO_NULL, &fh);

    if (rc != MPI_SUCCESS) {
        if (mpi_domain.is_root()) {
            std::cerr << "Error opening MPI file: " << filename << std::endl;
        }
        return false;
    }

    // Calculate file offset for this process
    // Each process writes its local data at the correct global position
    const unsigned int g = mpi_domain.ghost_width;
    const size_t elem_size = sizeof(double);

    // Write local data row by row to correct global positions
    for (unsigned int j = 0; j < mpi_domain.local_ncols; j++) {
        unsigned int global_col = mpi_domain.col_offset + j;
        MPI_Offset offset = (global_col * global_rows + mpi_domain.row_offset) * elem_size;

        // Extract column from local matrix (excluding ghost cells)
        std::vector<double> col_data(mpi_domain.local_nrows);
        for (unsigned int i = 0; i < mpi_domain.local_nrows; i++) {
            col_data[i] = local_mat(i + g, j + g);
        }

        MPI_File_write_at(fh, offset, col_data.data(),
                          mpi_domain.local_nrows, MPI_DOUBLE, &status);
    }

    MPI_File_close(&fh);
    return true;

#else
    // Non-MPI fallback
    std::ofstream outfile(filename, std::ios::binary);
    if (!outfile.is_open()) return false;

    outfile.write(reinterpret_cast<const char*>(local_mat.memptr()),
                  local_mat.n_elem * sizeof(double));
    outfile.close();
    return true;
#endif
}

bool MPIParallelIO::write_parallel_hdf5(
    GlobVar& ds,
    unsigned int print_time,
    const std::string& output_folder)
{
    // HDF5 parallel output would go here
    // Requires HDF5 compiled with parallel support
    // For now, fall back to gathered output
    return write_gathered_results(ds, print_time, output_folder);
}
