// Copyright 2019, Diogo Costa
// MPI Domain Decomposition for FLUXOS
// Hybrid MPI+OpenMP implementation for HPC clusters

#include "mpi_domain.h"
#include <iostream>
#include <cmath>
#include <algorithm>

// Global instance
MPIDomain mpi_domain;

MPIDomain::MPIDomain()
    : world_rank(0), world_size(1), cart_rank(0),
      local_nrows(0), local_ncols(0),
      global_nrows(0), global_ncols(0),
      ghost_width(1), row_offset(0), col_offset(0),
      initialized(false)
{
    coords[0] = coords[1] = 0;
    dims[0] = dims[1] = 1;
    neighbors[NORTH] = neighbors[SOUTH] = neighbors[EAST] = neighbors[WEST] = -1;
#ifdef USE_MPI
    cart_comm = MPI_COMM_NULL;
#endif
}

MPIDomain::~MPIDomain() {
    if (initialized) {
        finalize();
    }
}

void MPIDomain::initialize(int argc, char* argv[],
                           unsigned int global_rows,
                           unsigned int global_cols,
                           unsigned int ghost) {
#ifdef USE_MPI
    // Initialize MPI
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);

    if (provided < MPI_THREAD_FUNNELED) {
        std::cerr << "Warning: MPI does not support MPI_THREAD_FUNNELED" << std::endl;
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    global_nrows = global_rows;
    global_ncols = global_cols;
    ghost_width = ghost;

    // Let MPI choose optimal 2D decomposition
    dims[0] = dims[1] = 0;
    MPI_Dims_create(world_size, 2, dims);

    // Create 2D Cartesian communicator
    int periods[2] = {0, 0};  // Non-periodic boundaries
    int reorder = 1;          // Allow reordering for optimization
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, reorder, &cart_comm);

    // Get rank and coordinates in Cartesian grid
    MPI_Comm_rank(cart_comm, &cart_rank);
    MPI_Cart_coords(cart_comm, cart_rank, 2, coords);

    // Get neighbor ranks
    MPI_Cart_shift(cart_comm, 0, 1, &neighbors[WEST], &neighbors[EAST]);   // Row direction
    MPI_Cart_shift(cart_comm, 1, 1, &neighbors[SOUTH], &neighbors[NORTH]); // Col direction

    // Compute local domain sizes
    compute_local_sizes();

    // Create MPI datatypes for ghost exchange
    create_mpi_types();

    initialized = true;

    if (is_root()) {
        std::cout << "MPI initialized with " << world_size << " processes ("
                  << dims[0] << " x " << dims[1] << " grid)" << std::endl;
    }
#else
    // Non-MPI fallback
    world_rank = 0;
    world_size = 1;
    cart_rank = 0;
    coords[0] = coords[1] = 0;
    dims[0] = dims[1] = 1;
    global_nrows = global_rows;
    global_ncols = global_cols;
    ghost_width = ghost;
    local_nrows = global_rows;
    local_ncols = global_cols;
    row_offset = 0;
    col_offset = 0;
    initialized = true;
#endif
}

void MPIDomain::compute_local_sizes() {
#ifdef USE_MPI
    // Divide domain among processes
    // Row distribution
    unsigned int base_rows = global_nrows / dims[0];
    unsigned int extra_rows = global_nrows % dims[0];

    if (coords[0] < (int)extra_rows) {
        local_nrows = base_rows + 1;
        row_offset = coords[0] * (base_rows + 1);
    } else {
        local_nrows = base_rows;
        row_offset = extra_rows * (base_rows + 1) + (coords[0] - extra_rows) * base_rows;
    }

    // Column distribution
    unsigned int base_cols = global_ncols / dims[1];
    unsigned int extra_cols = global_ncols % dims[1];

    if (coords[1] < (int)extra_cols) {
        local_ncols = base_cols + 1;
        col_offset = coords[1] * (base_cols + 1);
    } else {
        local_ncols = base_cols;
        col_offset = extra_cols * (base_cols + 1) + (coords[1] - extra_cols) * base_cols;
    }
#endif
}

void MPIDomain::create_mpi_types() {
#ifdef USE_MPI
    // Create datatype for column exchange (non-contiguous in row-major)
    // Armadillo uses column-major storage, so columns are contiguous
    // but rows are not

    // For double matrices - a column vector of size local_nrows
    MPI_Type_vector(local_nrows + 2 * ghost_width, ghost_width, 1,
                    MPI_DOUBLE, &col_type_double);
    MPI_Type_commit(&col_type_double);

    // For float matrices
    MPI_Type_vector(local_nrows + 2 * ghost_width, ghost_width, 1,
                    MPI_FLOAT, &col_type_float);
    MPI_Type_commit(&col_type_float);
#endif
}

void MPIDomain::exchange_ghost_cells(arma::Mat<double>& matrix) {
#ifdef USE_MPI
    if (world_size == 1) return;

    MPI_Request requests[8];
    MPI_Status statuses[8];
    int req_count = 0;

    const unsigned int nrows = matrix.n_rows;
    const unsigned int ncols = matrix.n_cols;
    const unsigned int g = ghost_width;

    // Allocate send/receive buffers for row exchanges
    std::vector<double> send_north(nrows * g), recv_north(nrows * g);
    std::vector<double> send_south(nrows * g), recv_south(nrows * g);

    // Pack data for North/South (column direction) exchange
    // Send to North: last real column
    // Send to South: first real column
    if (neighbors[NORTH] >= 0) {
        for (unsigned int i = 0; i < nrows; i++) {
            for (unsigned int j = 0; j < g; j++) {
                send_north[i * g + j] = matrix(i, ncols - 2 * g + j);
            }
        }
        MPI_Isend(send_north.data(), nrows * g, MPI_DOUBLE,
                  neighbors[NORTH], 0, cart_comm, &requests[req_count++]);
        MPI_Irecv(recv_north.data(), nrows * g, MPI_DOUBLE,
                  neighbors[NORTH], 1, cart_comm, &requests[req_count++]);
    }

    if (neighbors[SOUTH] >= 0) {
        for (unsigned int i = 0; i < nrows; i++) {
            for (unsigned int j = 0; j < g; j++) {
                send_south[i * g + j] = matrix(i, g + j);
            }
        }
        MPI_Isend(send_south.data(), nrows * g, MPI_DOUBLE,
                  neighbors[SOUTH], 1, cart_comm, &requests[req_count++]);
        MPI_Irecv(recv_south.data(), nrows * g, MPI_DOUBLE,
                  neighbors[SOUTH], 0, cart_comm, &requests[req_count++]);
    }

    // Allocate send/receive buffers for East/West (row direction) exchange
    std::vector<double> send_east(ncols * g), recv_east(ncols * g);
    std::vector<double> send_west(ncols * g), recv_west(ncols * g);

    if (neighbors[EAST] >= 0) {
        for (unsigned int j = 0; j < ncols; j++) {
            for (unsigned int i = 0; i < g; i++) {
                send_east[j * g + i] = matrix(nrows - 2 * g + i, j);
            }
        }
        MPI_Isend(send_east.data(), ncols * g, MPI_DOUBLE,
                  neighbors[EAST], 2, cart_comm, &requests[req_count++]);
        MPI_Irecv(recv_east.data(), ncols * g, MPI_DOUBLE,
                  neighbors[EAST], 3, cart_comm, &requests[req_count++]);
    }

    if (neighbors[WEST] >= 0) {
        for (unsigned int j = 0; j < ncols; j++) {
            for (unsigned int i = 0; i < g; i++) {
                send_west[j * g + i] = matrix(g + i, j);
            }
        }
        MPI_Isend(send_west.data(), ncols * g, MPI_DOUBLE,
                  neighbors[WEST], 3, cart_comm, &requests[req_count++]);
        MPI_Irecv(recv_west.data(), ncols * g, MPI_DOUBLE,
                  neighbors[WEST], 2, cart_comm, &requests[req_count++]);
    }

    // Wait for all communications to complete
    MPI_Waitall(req_count, requests, statuses);

    // Unpack received data into ghost cells
    if (neighbors[NORTH] >= 0) {
        for (unsigned int i = 0; i < nrows; i++) {
            for (unsigned int j = 0; j < g; j++) {
                matrix(i, ncols - g + j) = recv_north[i * g + j];
            }
        }
    }

    if (neighbors[SOUTH] >= 0) {
        for (unsigned int i = 0; i < nrows; i++) {
            for (unsigned int j = 0; j < g; j++) {
                matrix(i, j) = recv_south[i * g + j];
            }
        }
    }

    if (neighbors[EAST] >= 0) {
        for (unsigned int j = 0; j < ncols; j++) {
            for (unsigned int i = 0; i < g; i++) {
                matrix(nrows - g + i, j) = recv_east[j * g + i];
            }
        }
    }

    if (neighbors[WEST] >= 0) {
        for (unsigned int j = 0; j < ncols; j++) {
            for (unsigned int i = 0; i < g; i++) {
                matrix(i, j) = recv_west[j * g + i];
            }
        }
    }
#endif
}

void MPIDomain::exchange_ghost_cells(arma::Mat<float>& matrix) {
#ifdef USE_MPI
    if (world_size == 1) return;

    MPI_Request requests[8];
    MPI_Status statuses[8];
    int req_count = 0;

    const unsigned int nrows = matrix.n_rows;
    const unsigned int ncols = matrix.n_cols;
    const unsigned int g = ghost_width;

    // Same logic as double version but with float
    std::vector<float> send_north(nrows * g), recv_north(nrows * g);
    std::vector<float> send_south(nrows * g), recv_south(nrows * g);
    std::vector<float> send_east(ncols * g), recv_east(ncols * g);
    std::vector<float> send_west(ncols * g), recv_west(ncols * g);

    // North/South exchanges
    if (neighbors[NORTH] >= 0) {
        for (unsigned int i = 0; i < nrows; i++) {
            for (unsigned int j = 0; j < g; j++) {
                send_north[i * g + j] = matrix(i, ncols - 2 * g + j);
            }
        }
        MPI_Isend(send_north.data(), nrows * g, MPI_FLOAT,
                  neighbors[NORTH], 0, cart_comm, &requests[req_count++]);
        MPI_Irecv(recv_north.data(), nrows * g, MPI_FLOAT,
                  neighbors[NORTH], 1, cart_comm, &requests[req_count++]);
    }

    if (neighbors[SOUTH] >= 0) {
        for (unsigned int i = 0; i < nrows; i++) {
            for (unsigned int j = 0; j < g; j++) {
                send_south[i * g + j] = matrix(i, g + j);
            }
        }
        MPI_Isend(send_south.data(), nrows * g, MPI_FLOAT,
                  neighbors[SOUTH], 1, cart_comm, &requests[req_count++]);
        MPI_Irecv(recv_south.data(), nrows * g, MPI_FLOAT,
                  neighbors[SOUTH], 0, cart_comm, &requests[req_count++]);
    }

    // East/West exchanges
    if (neighbors[EAST] >= 0) {
        for (unsigned int j = 0; j < ncols; j++) {
            for (unsigned int i = 0; i < g; i++) {
                send_east[j * g + i] = matrix(nrows - 2 * g + i, j);
            }
        }
        MPI_Isend(send_east.data(), ncols * g, MPI_FLOAT,
                  neighbors[EAST], 2, cart_comm, &requests[req_count++]);
        MPI_Irecv(recv_east.data(), ncols * g, MPI_FLOAT,
                  neighbors[EAST], 3, cart_comm, &requests[req_count++]);
    }

    if (neighbors[WEST] >= 0) {
        for (unsigned int j = 0; j < ncols; j++) {
            for (unsigned int i = 0; i < g; i++) {
                send_west[j * g + i] = matrix(g + i, j);
            }
        }
        MPI_Isend(send_west.data(), ncols * g, MPI_FLOAT,
                  neighbors[WEST], 3, cart_comm, &requests[req_count++]);
        MPI_Irecv(recv_west.data(), ncols * g, MPI_FLOAT,
                  neighbors[WEST], 2, cart_comm, &requests[req_count++]);
    }

    MPI_Waitall(req_count, requests, statuses);

    // Unpack
    if (neighbors[NORTH] >= 0) {
        for (unsigned int i = 0; i < nrows; i++) {
            for (unsigned int j = 0; j < g; j++) {
                matrix(i, ncols - g + j) = recv_north[i * g + j];
            }
        }
    }

    if (neighbors[SOUTH] >= 0) {
        for (unsigned int i = 0; i < nrows; i++) {
            for (unsigned int j = 0; j < g; j++) {
                matrix(i, j) = recv_south[i * g + j];
            }
        }
    }

    if (neighbors[EAST] >= 0) {
        for (unsigned int j = 0; j < ncols; j++) {
            for (unsigned int i = 0; i < g; i++) {
                matrix(nrows - g + i, j) = recv_east[j * g + i];
            }
        }
    }

    if (neighbors[WEST] >= 0) {
        for (unsigned int j = 0; j < ncols; j++) {
            for (unsigned int i = 0; i < g; i++) {
                matrix(i, j) = recv_west[j * g + i];
            }
        }
    }
#endif
}

double MPIDomain::global_min(double local_value) {
#ifdef USE_MPI
    double global_value;
    MPI_Allreduce(&local_value, &global_value, 1, MPI_DOUBLE, MPI_MIN, cart_comm);
    return global_value;
#else
    return local_value;
#endif
}

double MPIDomain::global_max(double local_value) {
#ifdef USE_MPI
    double global_value;
    MPI_Allreduce(&local_value, &global_value, 1, MPI_DOUBLE, MPI_MAX, cart_comm);
    return global_value;
#else
    return local_value;
#endif
}

double MPIDomain::global_sum(double local_value) {
#ifdef USE_MPI
    double global_value;
    MPI_Allreduce(&local_value, &global_value, 1, MPI_DOUBLE, MPI_SUM, cart_comm);
    return global_value;
#else
    return local_value;
#endif
}

void MPIDomain::gather_to_root(const arma::Mat<double>& local_mat,
                                arma::Mat<double>& global_mat) {
#ifdef USE_MPI
    if (world_size == 1) {
        global_mat = local_mat;
        return;
    }

    // Gather sizes from all processes
    std::vector<int> all_nrows(world_size), all_ncols(world_size);
    std::vector<int> all_row_offsets(world_size), all_col_offsets(world_size);

    int my_nrows = local_nrows;
    int my_ncols = local_ncols;
    int my_row_offset = row_offset;
    int my_col_offset = col_offset;

    MPI_Gather(&my_nrows, 1, MPI_INT, all_nrows.data(), 1, MPI_INT, 0, cart_comm);
    MPI_Gather(&my_ncols, 1, MPI_INT, all_ncols.data(), 1, MPI_INT, 0, cart_comm);
    MPI_Gather(&my_row_offset, 1, MPI_INT, all_row_offsets.data(), 1, MPI_INT, 0, cart_comm);
    MPI_Gather(&my_col_offset, 1, MPI_INT, all_col_offsets.data(), 1, MPI_INT, 0, cart_comm);

    // Extract local data (without ghost cells)
    arma::Mat<double> local_data = local_mat.submat(
        ghost_width, ghost_width,
        ghost_width + local_nrows - 1,
        ghost_width + local_ncols - 1
    );

    if (is_root()) {
        global_mat.zeros(global_nrows, global_ncols);

        // Copy root's data
        global_mat.submat(row_offset, col_offset,
                         row_offset + local_nrows - 1,
                         col_offset + local_ncols - 1) = local_data;

        // Receive from all other processes
        for (int p = 1; p < world_size; p++) {
            arma::Mat<double> recv_data(all_nrows[p], all_ncols[p]);
            MPI_Recv(recv_data.memptr(), all_nrows[p] * all_ncols[p], MPI_DOUBLE,
                     p, 0, cart_comm, MPI_STATUS_IGNORE);

            global_mat.submat(all_row_offsets[p], all_col_offsets[p],
                             all_row_offsets[p] + all_nrows[p] - 1,
                             all_col_offsets[p] + all_ncols[p] - 1) = recv_data;
        }
    } else {
        MPI_Send(local_data.memptr(), local_nrows * local_ncols, MPI_DOUBLE,
                 0, 0, cart_comm);
    }
#else
    global_mat = local_mat;
#endif
}

void MPIDomain::scatter_from_root(const arma::Mat<double>& global_mat,
                                   arma::Mat<double>& local_mat) {
#ifdef USE_MPI
    if (world_size == 1) {
        local_mat = global_mat;
        return;
    }

    // Gather sizes from all processes (same as gather)
    std::vector<int> all_nrows(world_size), all_ncols(world_size);
    std::vector<int> all_row_offsets(world_size), all_col_offsets(world_size);

    int my_nrows = local_nrows;
    int my_ncols = local_ncols;
    int my_row_offset = row_offset;
    int my_col_offset = col_offset;

    MPI_Gather(&my_nrows, 1, MPI_INT, all_nrows.data(), 1, MPI_INT, 0, cart_comm);
    MPI_Gather(&my_ncols, 1, MPI_INT, all_ncols.data(), 1, MPI_INT, 0, cart_comm);
    MPI_Gather(&my_row_offset, 1, MPI_INT, all_row_offsets.data(), 1, MPI_INT, 0, cart_comm);
    MPI_Gather(&my_col_offset, 1, MPI_INT, all_col_offsets.data(), 1, MPI_INT, 0, cart_comm);

    // Allocate local matrix with ghost cells
    local_mat.zeros(local_nrows + 2 * ghost_width, local_ncols + 2 * ghost_width);

    arma::Mat<double> local_data(local_nrows, local_ncols);

    if (is_root()) {
        // Extract and send to all processes
        for (int p = 1; p < world_size; p++) {
            arma::Mat<double> send_data = global_mat.submat(
                all_row_offsets[p], all_col_offsets[p],
                all_row_offsets[p] + all_nrows[p] - 1,
                all_col_offsets[p] + all_ncols[p] - 1
            );
            MPI_Send(send_data.memptr(), all_nrows[p] * all_ncols[p], MPI_DOUBLE,
                     p, 0, cart_comm);
        }

        // Copy root's portion
        local_data = global_mat.submat(row_offset, col_offset,
                                        row_offset + local_nrows - 1,
                                        col_offset + local_ncols - 1);
    } else {
        MPI_Recv(local_data.memptr(), local_nrows * local_ncols, MPI_DOUBLE,
                 0, 0, cart_comm, MPI_STATUS_IGNORE);
    }

    // Copy into local matrix (with ghost cell padding)
    local_mat.submat(ghost_width, ghost_width,
                     ghost_width + local_nrows - 1,
                     ghost_width + local_ncols - 1) = local_data;
#else
    local_mat = global_mat;
#endif
}

void MPIDomain::barrier() {
#ifdef USE_MPI
    MPI_Barrier(cart_comm);
#endif
}

void MPIDomain::finalize() {
#ifdef USE_MPI
    if (initialized) {
        MPI_Type_free(&col_type_double);
        MPI_Type_free(&col_type_float);
        if (cart_comm != MPI_COMM_NULL) {
            MPI_Comm_free(&cart_comm);
        }
        MPI_Finalize();
        initialized = false;
    }
#endif
}

void MPIDomain::print_info() {
#ifdef USE_MPI
    barrier();
    for (int p = 0; p < world_size; p++) {
        if (cart_rank == p) {
            std::cout << "Rank " << cart_rank << " (" << coords[0] << "," << coords[1] << "): "
                      << "local size = " << local_nrows << "x" << local_ncols
                      << ", offset = (" << row_offset << "," << col_offset << ")"
                      << ", neighbors = [N:" << neighbors[NORTH]
                      << " S:" << neighbors[SOUTH]
                      << " E:" << neighbors[EAST]
                      << " W:" << neighbors[WEST] << "]"
                      << std::endl;
        }
        barrier();
    }
#else
    std::cout << "MPI not enabled. Running in serial mode." << std::endl;
#endif
}
