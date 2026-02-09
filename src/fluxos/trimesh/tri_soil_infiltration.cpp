// Copyright 2019, Diogo Costa
// Soil infiltration loss for triangular mesh
// Horton decay: f(t) = fc + (f0 - fc) * exp(-k * t_wet)

#include "tri_soil_infiltration.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>

// ============================================================================
// Initialize soil parameters for triangular mesh cells
// ============================================================================
bool tri_init_soil(
    GlobVar& ds,
    const TriMesh& mesh,
    TriSolution& sol,
    const SoilConfig& soil_config,
    std::ofstream& logFLUXOSfile)
{
    if (!soil_config.enabled) return false;

    const int ncells = mesh.num_cells;

    // Allocate soil fields
    sol.soil_infil_rate.assign(ncells, 0.0);
    sol.soil_Ks.assign(ncells, 0.0);
    sol.soil_f0.assign(ncells, 0.0);
    sol.soil_k.assign(ncells, 0.0);
    sol.soil_wetting_time.assign(ncells, 0.0);
    sol.soil_type_id.assign(ncells, 0);

    // If no soil map, use uniform soil type 1
    if (soil_config.soil_map_file.empty()) {
        logFLUXOSfile << "  TriMesh soil: no SOIL_MAP; using uniform soil type 1\n";

        auto it = soil_config.soil_types.find(1);
        if (it == soil_config.soil_types.end()) {
            std::cout << "ERROR: No soil type 1 defined for uniform fallback" << std::endl;
            return true;
        }
        const SoilParams& sp = it->second;

        for (int ci = 0; ci < ncells; ci++) {
            if (sol.innerNeumannBCWeir[ci] == 1.0f) continue;
            sol.soil_type_id[ci] = 1;
            sol.soil_Ks[ci] = sp.Ks;
            sol.soil_f0[ci] = sp.f0;
            sol.soil_k[ci] = sp.k;
        }
        return false;
    }

    // ---- Read soil map raster and sample at cell centroids ----
    std::ifstream mapfile(soil_config.soil_map_file);
    if (!mapfile.is_open()) {
        std::cout << "ERROR: Cannot open soil map file: "
                  << soil_config.soil_map_file << std::endl;
        return true;
    }

    std::string line;
    int map_ncols = 0, map_nrows = 0;
    double map_xll = 0.0, map_yll = 0.0, map_cellsize = 0.0;
    int map_nodata = -9999;

    for (int i = 0; i < 6; i++) {
        getline(mapfile, line);
        std::stringstream ss(line);
        std::string key;
        ss >> key;
        std::transform(key.begin(), key.end(), key.begin(), ::toupper);

        if (key == "NCOLS") ss >> map_ncols;
        else if (key == "NROWS") ss >> map_nrows;
        else if (key == "XLLCORNER") ss >> map_xll;
        else if (key == "YLLCORNER") ss >> map_yll;
        else if (key == "CELLSIZE") ss >> map_cellsize;
        else if (key == "NODATA_VALUE") ss >> map_nodata;
    }

    std::vector<std::vector<int>> soil_grid(map_nrows, std::vector<int>(map_ncols, map_nodata));
    for (int r = map_nrows - 1; r >= 0; r--) {
        getline(mapfile, line);
        std::stringstream ss(line);
        for (int c = 0; c < map_ncols; c++) {
            ss >> soil_grid[r][c];
        }
    }
    mapfile.close();

    int assigned = 0;
    for (int ci = 0; ci < ncells; ci++) {
        if (sol.innerNeumannBCWeir[ci] == 1.0f) continue;

        double cx = mesh.cells[ci].cx;
        double cy = mesh.cells[ci].cy;

        int gc = static_cast<int>(std::floor((cx - map_xll) / map_cellsize));
        int gr = static_cast<int>(std::floor((cy - map_yll) / map_cellsize));

        gc = std::max(0, std::min(gc, map_ncols - 1));
        gr = std::max(0, std::min(gr, map_nrows - 1));

        int soil_id = soil_grid[gr][gc];

        if (soil_id == map_nodata) {
            sol.soil_type_id[ci] = 0;
            continue;
        }

        sol.soil_type_id[ci] = soil_id;

        auto it = soil_config.soil_types.find(soil_id);
        if (it != soil_config.soil_types.end()) {
            sol.soil_Ks[ci] = it->second.Ks;
            sol.soil_f0[ci] = it->second.f0;
            sol.soil_k[ci] = it->second.k;
            assigned++;
        }
    }

    logFLUXOSfile << "  TriMesh soil: " << assigned << " of " << ncells
                  << " cells assigned soil parameters from map\n";
    return false;
}

// ============================================================================
// Apply Horton soil infiltration loss for triangular mesh (one timestep)
// f(t) = fc + (f0 - fc) * exp(-k * t_wet)
// ============================================================================
void tri_apply_soil_infiltration(
    GlobVar& ds,
    const TriMesh& mesh,
    TriSolution& sol,
    const SoilConfig& soil_config)
{
    if (!soil_config.enabled) return;

    const double dt = ds.dtfl;
    const double hdry = ds.hdry;
    const int ncells = mesh.num_cells;

    #pragma omp parallel for schedule(static)
    for (int ci = 0; ci < ncells; ci++) {
        if (sol.innerNeumannBCWeir[ci] == 1.0f) continue;

        double Ks = sol.soil_Ks[ci];
        if (Ks <= 0.0) continue;

        double h_surface = std::fmax(sol.z[ci] - sol.zb[ci], 0.0);
        if (h_surface <= 0.0) continue;

        // Horton decay: f(t) = fc + (f0 - fc) * exp(-k * t_wet)
        double f0 = sol.soil_f0[ci];
        double k_decay = sol.soil_k[ci];
        double t_wet = sol.soil_wetting_time[ci];

        double f_rate = Ks + (f0 - Ks) * std::exp(-k_decay * t_wet);

        double infil = std::fmin(f_rate * dt, h_surface);

        sol.z[ci] -= infil;
        sol.h[ci] = std::fmax(sol.z[ci] - sol.zb[ci], 0.0);
        sol.ldry[ci] = (sol.h[ci] <= hdry) ? 1.0f : 0.0f;
        sol.soil_infil_rate[ci] = infil / dt;

        // Advance wetting time
        sol.soil_wetting_time[ci] += dt;
    }
}
