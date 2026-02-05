// Copyright 2019, Diogo Costa
// Initialization routines for triangular mesh

#include "tri_initiate.h"
#include "tri_mesh_io.h"
#include <iostream>
#include <cmath>
#include <algorithm>

// ============================================================================
// Initialize triangular mesh solution
// ============================================================================
void tri_initiation(
    GlobVar& ds,
    TriMesh& mesh,
    TriSolution& sol,
    int nchem,
    std::ofstream& logFLUXOSfile)
{
    const int ncells = mesh.num_cells;
    const int nedges = mesh.num_edges;

    // Allocate solution arrays
    sol.allocate(ncells, nedges, nchem);

    logFLUXOSfile << "\nTriangular mesh initialization:" << std::endl;
    logFLUXOSfile << "  Cells: " << ncells << std::endl;
    logFLUXOSfile << "  Edges: " << nedges << std::endl;

    // ---- Interpolate bed elevation from DEM to cell centroids ----
    // The DEM is in the regular grid stored in ds.zb
    // Bilinear interpolation from DEM grid to triangle centroids
    for (int ci = 0; ci < ncells; ci++) {
        double cx = mesh.cells[ci].cx;
        double cy = mesh.cells[ci].cy;

        // Convert to grid coordinates
        double gx = (cx - ds.XLLCORNER) / ds.dxy;
        double gy = (cy - ds.YLLCORNER) / ds.dxy;

        // Grid indices
        int ix = static_cast<int>(std::floor(gx));
        int iy = static_cast<int>(std::floor(gy));

        // Clamp to valid range
        ix = std::max(1, std::min(ix, static_cast<int>(ds.NROWS)));
        iy = std::max(1, std::min(iy, static_cast<int>(ds.NCOLS)));

        // Use nearest-neighbor for simplicity (can upgrade to bilinear)
        double zb_val = (*ds.zb)(ix, iy);

        // Check for NODATA
        if (zb_val == ds.NODATA_VALUE || zb_val == 0.0) {
            sol.innerNeumannBCWeir[ci] = 1.0f;
            // Average from valid neighbors
            double sum_zb = 0.0;
            int count = 0;
            for (int e = 0; e < 3; e++) {
                int nj = mesh.cells[ci].neighbors[e];
                if (nj >= 0) {
                    double ncx = mesh.cells[nj].cx;
                    double ncy = mesh.cells[nj].cy;
                    double ngx = (ncx - ds.XLLCORNER) / ds.dxy;
                    double ngy = (ncy - ds.YLLCORNER) / ds.dxy;
                    int nix = std::max(1, std::min(static_cast<int>(std::floor(ngx)),
                                                    static_cast<int>(ds.NROWS)));
                    int niy = std::max(1, std::min(static_cast<int>(std::floor(ngy)),
                                                    static_cast<int>(ds.NCOLS)));
                    double nzb = (*ds.zb)(nix, niy);
                    if (nzb != ds.NODATA_VALUE && nzb != 0.0) {
                        sum_zb += std::fabs(nzb);
                        count++;
                    }
                }
            }
            if (count > 0) {
                zb_val = sum_zb / count;
            }
        } else {
            zb_val = std::fabs(zb_val);
        }

        sol.zb[ci] = zb_val;
        sol.z[ci] = zb_val;  // Initially dry (water surface = bed)
        sol.h[ci] = 0.0;
        sol.ks[ci] = ds.hdry;  // Same as regular mesh: hdry from ks
        sol.ldry[ci] = 1.0f;   // Start dry

        // Soil mass initialization for WINTRA
        if (ds.wintra) {
            sol.soil_mass[ci] = ds.soil_conc_bckgrd;
        }
    }

    // ---- Find inflow cell ----
    if (ds.inflow_xcoord != 0 || ds.inflow_ycoord != 0) {
        double inflow_x = static_cast<double>(ds.inflow_xcoord);
        double inflow_y = static_cast<double>(ds.inflow_ycoord);

        int inflow_cell = mesh.find_cell(inflow_x, inflow_y);
        if (inflow_cell >= 0) {
            logFLUXOSfile << "  Inflow cell found: " << inflow_cell
                          << " at (" << mesh.cells[inflow_cell].cx << ", "
                          << mesh.cells[inflow_cell].cy << ")" << std::endl;
            // Store inflow cell index in ds for later use
            // We use inflow_nrow to store the cell index for triangular mesh
            ds.inflow_nrow = inflow_cell;
            ds.inflow_ncol = 0;  // Not used for triangular mesh
        } else {
            logFLUXOSfile << "  WARNING: Inflow location not found in mesh!" << std::endl;
        }
    }

    // Apply boundary conditions
    tri_apply_boundary_conditions(ds, mesh, sol);

    logFLUXOSfile << "  Initialization complete." << std::endl;
}

// ============================================================================
// Apply boundary conditions based on mesh boundary tags
// ============================================================================
void tri_apply_boundary_conditions(
    GlobVar& ds,
    TriMesh& mesh,
    TriSolution& sol)
{
    // Mark boundary cells based on boundary_conditions config
    // Default: wall boundaries (zero-flux, handled naturally by edge solver)

    for (int ei = 0; ei < mesh.num_edges; ei++) {
        auto& edge = mesh.edges[ei];
        if (edge.is_boundary && edge.left_cell >= 0) {
            mesh.cells[edge.left_cell].is_boundary = true;

            // Boundary type is determined at runtime during flux computation
            // based on edge.boundary_tag and the JSON config
        }
    }
}
