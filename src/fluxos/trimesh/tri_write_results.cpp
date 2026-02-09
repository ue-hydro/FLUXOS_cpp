// Copyright 2019, Diogo Costa
// VTK XML Unstructured Grid (.vtu) output for triangular mesh

#include "tri_write_results.h"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>

// ============================================================================
// Write a single timestep to VTK .vtu file
// ============================================================================
bool tri_write_results(
    GlobVar& ds,
    const TriMesh& mesh,
    const TriSolution& sol,
    int print_tag,
    std::chrono::duration<double> elapsed_seconds)
{
    std::string filename = ds.output_folder + "/" + std::to_string(print_tag) + ".vtu";

    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "ERROR: Cannot write VTK file: " << filename << std::endl;
        return false;
    }

    const int ncells = mesh.num_cells;
    const int nverts = mesh.num_vertices;

    file << std::scientific << std::setprecision(8);

    // VTK XML header
    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
    file << "  <UnstructuredGrid>\n";
    file << "    <Piece NumberOfPoints=\"" << nverts
         << "\" NumberOfCells=\"" << ncells << "\">\n";

    // ---- Point data (none for now, cell-centered data below) ----

    // ---- Cell data ----
    file << "      <CellData>\n";

    // Water depth
    file << "        <DataArray type=\"Float64\" Name=\"h\" format=\"ascii\">\n";
    for (int ci = 0; ci < ncells; ci++) {
        if (sol.h[ci] > ds.h_min_print) {
            file << "          " << sol.h[ci] << "\n";
        } else {
            file << "          " << 0.0 << "\n";
        }
    }
    file << "        </DataArray>\n";

    // Water surface elevation
    file << "        <DataArray type=\"Float64\" Name=\"z\" format=\"ascii\">\n";
    for (int ci = 0; ci < ncells; ci++) {
        file << "          " << sol.z[ci] << "\n";
    }
    file << "        </DataArray>\n";

    // Bed elevation
    file << "        <DataArray type=\"Float64\" Name=\"zb\" format=\"ascii\">\n";
    for (int ci = 0; ci < ncells; ci++) {
        file << "          " << sol.zb[ci] << "\n";
    }
    file << "        </DataArray>\n";

    // Velocity (2-component vector)
    file << "        <DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int ci = 0; ci < ncells; ci++) {
        file << "          " << sol.ux[ci] << " " << sol.uy[ci] << " 0.0\n";
    }
    file << "        </DataArray>\n";

    // Discharge
    file << "        <DataArray type=\"Float64\" Name=\"discharge\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int ci = 0; ci < ncells; ci++) {
        file << "          " << sol.qx[ci] << " " << sol.qy[ci] << " 0.0\n";
    }
    file << "        </DataArray>\n";

    // Shear velocity
    file << "        <DataArray type=\"Float64\" Name=\"us\" format=\"ascii\">\n";
    for (int ci = 0; ci < ncells; ci++) {
        file << "          " << sol.us[ci] << "\n";
    }
    file << "        </DataArray>\n";

    // Concentration (first chemical species)
    if (sol.nchem > 0 && !sol.conc_SW.empty()) {
        file << "        <DataArray type=\"Float64\" Name=\"conc_SW\" format=\"ascii\">\n";
        for (int ci = 0; ci < ncells; ci++) {
            file << "          " << sol.conc_SW[0][ci] << "\n";
        }
        file << "        </DataArray>\n";
    }

    // Wetting time tracer
    file << "        <DataArray type=\"Float64\" Name=\"twetimetracer\" format=\"ascii\">\n";
    for (int ci = 0; ci < ncells; ci++) {
        file << "          " << sol.twetimetracer[ci] << "\n";
    }
    file << "        </DataArray>\n";

    // Cell area (useful for debugging)
    file << "        <DataArray type=\"Float64\" Name=\"cell_area\" format=\"ascii\">\n";
    for (int ci = 0; ci < ncells; ci++) {
        file << "          " << mesh.cells[ci].area << "\n";
    }
    file << "        </DataArray>\n";

    // Soil infiltration fields (only output if soil module is active)
    if (!sol.soil_Ks.empty() && sol.soil_Ks.size() == static_cast<size_t>(ncells)) {
        bool has_soil_data = false;
        for (int ci = 0; ci < ncells; ci++) {
            if (sol.soil_Ks[ci] > 0.0) { has_soil_data = true; break; }
        }
        if (has_soil_data) {
            // Instantaneous infiltration rate
            file << "        <DataArray type=\"Float64\" Name=\"soil_infil_rate\" format=\"ascii\">\n";
            for (int ci = 0; ci < ncells; ci++) {
                file << "          " << sol.soil_infil_rate[ci] << "\n";
            }
            file << "        </DataArray>\n";

            // Soil type ID
            file << "        <DataArray type=\"Int32\" Name=\"soil_type\" format=\"ascii\">\n";
            for (int ci = 0; ci < ncells; ci++) {
                file << "          " << sol.soil_type_id[ci] << "\n";
            }
            file << "        </DataArray>\n";
        }
    }

    file << "      </CellData>\n";

    // ---- Points (vertices) ----
    file << "      <Points>\n";
    file << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (int vi = 0; vi < nverts; vi++) {
        file << "          " << mesh.vertices[vi].x << " "
             << mesh.vertices[vi].y << " 0.0\n";
    }
    file << "        </DataArray>\n";
    file << "      </Points>\n";

    // ---- Cells (triangles) ----
    file << "      <Cells>\n";

    // Connectivity
    file << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (int ci = 0; ci < ncells; ci++) {
        file << "          " << mesh.cells[ci].verts[0] << " "
             << mesh.cells[ci].verts[1] << " "
             << mesh.cells[ci].verts[2] << "\n";
    }
    file << "        </DataArray>\n";

    // Offsets
    file << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (int ci = 0; ci < ncells; ci++) {
        file << "          " << (ci + 1) * 3 << "\n";
    }
    file << "        </DataArray>\n";

    // Types (5 = VTK_TRIANGLE)
    file << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n";
    for (int ci = 0; ci < ncells; ci++) {
        file << "          5\n";
    }
    file << "        </DataArray>\n";

    file << "      </Cells>\n";

    file << "    </Piece>\n";
    file << "  </UnstructuredGrid>\n";
    file << "</VTKFile>\n";

    file.close();
    return true;
}

// ============================================================================
// Write .pvd time series collection file
// ============================================================================
void tri_write_pvd(
    const std::string& output_folder,
    const std::vector<std::pair<double, std::string>>& time_files)
{
    std::string pvd_file = output_folder + "/timeseries.pvd";
    std::ofstream file(pvd_file);
    if (!file.is_open()) {
        std::cerr << "WARNING: Cannot write PVD file: " << pvd_file << std::endl;
        return;
    }

    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
    file << "  <Collection>\n";

    for (const auto& tf : time_files) {
        file << "    <DataSet timestep=\"" << tf.first
             << "\" file=\"" << tf.second << "\"/>\n";
    }

    file << "  </Collection>\n";
    file << "</VTKFile>\n";

    file.close();
}
