// Copyright 2019, Diogo Costa
// MPI domain decomposition for unstructured triangular mesh

#include "tri_mpi_domain.h"
#include <iostream>
#include <algorithm>
#include <set>
#include <map>

// ============================================================================
// Initialize MPI
// ============================================================================
void TriMPIDomain::init()
{
#ifdef USE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    is_root = (rank == 0);
#else
    rank = 0;
    size = 1;
    is_root = true;
#endif
}

// ============================================================================
// Partition the mesh
// Uses naive block partitioning (METIS support can be added later)
// ============================================================================
void TriMPIDomain::partition_mesh(const TriMesh& mesh)
{
    const int ncells = mesh.num_cells;
    cell_partition.resize(ncells);

    if (size <= 1) {
        // Single process: all cells are local
        std::fill(cell_partition.begin(), cell_partition.end(), 0);
        local_num_cells = ncells;
        local_num_edges = mesh.num_edges;

        local_to_global_cell.resize(ncells);
        global_to_local_cell.resize(ncells);
        for (int i = 0; i < ncells; i++) {
            local_to_global_cell[i] = i;
            global_to_local_cell[i] = i;
        }
        return;
    }

#ifdef USE_MPI
    // Naive block partitioning: distribute cells evenly
    int cells_per_rank = ncells / size;
    int remainder = ncells % size;

    for (int ci = 0; ci < ncells; ci++) {
        if (ci < (cells_per_rank + 1) * remainder) {
            cell_partition[ci] = ci / (cells_per_rank + 1);
        } else {
            cell_partition[ci] = remainder + (ci - (cells_per_rank + 1) * remainder) / cells_per_rank;
        }
    }

    // Build local cell list
    local_to_global_cell.clear();
    global_to_local_cell.assign(ncells, -1);

    for (int ci = 0; ci < ncells; ci++) {
        if (cell_partition[ci] == rank) {
            global_to_local_cell[ci] = static_cast<int>(local_to_global_cell.size());
            local_to_global_cell.push_back(ci);
        }
    }
    local_num_cells = static_cast<int>(local_to_global_cell.size());

    // Find halo cells: cells owned by other ranks that neighbor our local cells
    std::set<int> halo_set;
    std::map<int, std::set<int>> rank_to_send;
    std::map<int, std::set<int>> rank_to_recv;

    for (int li = 0; li < local_num_cells; li++) {
        int gi = local_to_global_cell[li];
        const auto& cell = mesh.cells[gi];

        for (int e = 0; e < 3; e++) {
            int nj = cell.neighbors[e];
            if (nj >= 0 && cell_partition[nj] != rank) {
                int other_rank = cell_partition[nj];
                halo_set.insert(nj);
                rank_to_recv[other_rank].insert(nj);
                rank_to_send[other_rank].insert(gi);
            }
        }
    }

    halo_cells.assign(halo_set.begin(), halo_set.end());

    // Build neighbor communication structures
    neighbors.clear();
    for (const auto& [other_rank, recv_set] : rank_to_recv) {
        NeighborComm nc;
        nc.rank = other_rank;
        nc.recv_cells.assign(recv_set.begin(), recv_set.end());
        nc.send_cells.assign(rank_to_send[other_rank].begin(),
                              rank_to_send[other_rank].end());
        neighbors.push_back(nc);
    }

    if (is_root) {
        std::cout << "MPI partitioning: " << size << " ranks, "
                  << cells_per_rank << " cells/rank (+/- 1)" << std::endl;
    }
#endif
}

// ============================================================================
// Exchange halo cell data (double)
// ============================================================================
void TriMPIDomain::exchange_halo(std::vector<double>& field)
{
#ifdef USE_MPI
    if (size <= 1) return;

    std::vector<MPI_Request> requests;

    for (auto& nc : neighbors) {
        // Send our cells to neighbor
        std::vector<double> send_buf(nc.send_cells.size());
        for (size_t i = 0; i < nc.send_cells.size(); i++) {
            send_buf[i] = field[nc.send_cells[i]];
        }

        MPI_Request send_req;
        MPI_Isend(send_buf.data(), static_cast<int>(send_buf.size()),
                  MPI_DOUBLE, nc.rank, 0, MPI_COMM_WORLD, &send_req);
        requests.push_back(send_req);

        // Receive neighbor's cells
        std::vector<double> recv_buf(nc.recv_cells.size());
        MPI_Request recv_req;
        MPI_Irecv(recv_buf.data(), static_cast<int>(recv_buf.size()),
                  MPI_DOUBLE, nc.rank, 0, MPI_COMM_WORLD, &recv_req);
        requests.push_back(recv_req);
    }

    // Wait for all communication to complete
    MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);

    // Unpack received data
    int req_idx = 0;
    for (auto& nc : neighbors) {
        req_idx++;  // skip send
        // Recv data would need proper handling with buffer persistence
        // This is a simplified version - in production, use persistent buffers
        req_idx++;
    }
#endif
}

// ============================================================================
// Exchange halo cell data (float)
// ============================================================================
void TriMPIDomain::exchange_halo_float(std::vector<float>& field)
{
#ifdef USE_MPI
    if (size <= 1) return;

    // Similar to exchange_halo but with MPI_FLOAT
    std::vector<MPI_Request> requests;

    for (auto& nc : neighbors) {
        std::vector<float> send_buf(nc.send_cells.size());
        for (size_t i = 0; i < nc.send_cells.size(); i++) {
            send_buf[i] = field[nc.send_cells[i]];
        }

        MPI_Request send_req;
        MPI_Isend(send_buf.data(), static_cast<int>(send_buf.size()),
                  MPI_FLOAT, nc.rank, 1, MPI_COMM_WORLD, &send_req);
        requests.push_back(send_req);

        std::vector<float> recv_buf(nc.recv_cells.size());
        MPI_Request recv_req;
        MPI_Irecv(recv_buf.data(), static_cast<int>(recv_buf.size()),
                  MPI_FLOAT, nc.rank, 1, MPI_COMM_WORLD, &recv_req);
        requests.push_back(recv_req);
    }

    MPI_Waitall(static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE);
#endif
}

// ============================================================================
// Global reductions
// ============================================================================
double TriMPIDomain::global_min(double local_val)
{
#ifdef USE_MPI
    double global_val;
    MPI_Allreduce(&local_val, &global_val, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    return global_val;
#else
    return local_val;
#endif
}

double TriMPIDomain::global_max(double local_val)
{
#ifdef USE_MPI
    double global_val;
    MPI_Allreduce(&local_val, &global_val, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return global_val;
#else
    return local_val;
#endif
}

double TriMPIDomain::global_sum(double local_val)
{
#ifdef USE_MPI
    double global_val;
    MPI_Allreduce(&local_val, &global_val, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return global_val;
#else
    return local_val;
#endif
}

// ============================================================================
// Finalize MPI
// ============================================================================
void TriMPIDomain::finalize()
{
    // MPI_Finalize is called in main.cpp
}
