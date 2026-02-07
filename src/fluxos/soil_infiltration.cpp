// Copyright 2019, Diogo Costa
// Soil infiltration loss model for FLUXOS - Regular mesh implementation
// Horton decay: f(t) = fc + (f0 - fc) * exp(-k * t_wet)

#include "soil_infiltration.h"
#include "GlobVar.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>

// ============================================================================
// Parse soil infiltration configuration from JSON
// ============================================================================
bool read_soil_config(
    const nlohmann::json& master_json,
    SoilConfig& soil_config,
    std::ofstream& logFLUXOSfile)
{
    auto it = master_json.find("SOIL_INFILTRATION");
    if (it == master_json.end()) {
        soil_config.enabled = false;
        logFLUXOSfile << "Soil infiltration: DISABLED (no SOIL_INFILTRATION section)\n";
        return false; // Not an error - just not configured
    }

    const auto& soil_json = *it;

    // Check STATUS
    auto status_it = soil_json.find("STATUS");
    if (status_it == soil_json.end() || !status_it->get<bool>()) {
        soil_config.enabled = false;
        logFLUXOSfile << "Soil infiltration: DISABLED\n";
        return false;
    }

    soil_config.enabled = true;

    // Soil map file
    auto map_it = soil_json.find("SOIL_MAP");
    if (map_it != soil_json.end()) {
        soil_config.soil_map_file = map_it->get<std::string>();
    }

    // Get USDA lookup table for reference
    auto usda_table = get_usda_soil_table();

    // Parse soil types
    auto types_it = soil_json.find("SOIL_TYPES");
    if (types_it != soil_json.end()) {
        for (auto& [key, val] : types_it->items()) {
            int type_id = std::stoi(key);
            SoilParams params;

            // Check if it references a USDA class by name
            auto usda_it = val.find("USDA_CLASS");
            if (usda_it != val.end()) {
                std::string usda_class = usda_it->get<std::string>();
                std::transform(usda_class.begin(), usda_class.end(),
                               usda_class.begin(), ::tolower);
                auto lookup = usda_table.find(usda_class);
                if (lookup != usda_table.end()) {
                    params = lookup->second;
                } else {
                    std::cout << "WARNING: Unknown USDA class '" << usda_class
                              << "' for soil type " << type_id << std::endl;
                    continue;
                }
            } else {
                // Manual parameter specification
                params.name = val.value("NAME", "Unknown");
                params.Ks = val.value("Ks", 0.0);
                params.f0 = val.value("f0", params.Ks * 5.0); // default: 5x Ks
                params.k = val.value("k", 2.78e-4);            // default: ~1/hr
            }

            // Allow override of individual params even with USDA class
            if (val.contains("NAME")) params.name = val["NAME"].get<std::string>();
            if (val.contains("Ks")) params.Ks = val["Ks"].get<double>();
            if (val.contains("f0")) params.f0 = val["f0"].get<double>();
            if (val.contains("k")) params.k = val["k"].get<double>();

            soil_config.soil_types[type_id] = params;

            logFLUXOSfile << "  Soil type " << type_id << ": " << params.name
                          << " (Ks=" << params.Ks << " m/s"
                          << ", f0=" << params.f0 << " m/s"
                          << ", k=" << params.k << " 1/s)\n";
        }
    }

    logFLUXOSfile << "Soil infiltration: ENABLED (Horton decay model)\n";
    logFLUXOSfile << "  Soil map: " << soil_config.soil_map_file << "\n";
    logFLUXOSfile << "  Number of soil types: " << soil_config.soil_types.size() << "\n";

    return false; // no error
}

// ============================================================================
// Read soil type map (ESRI ASCII grid) and initialize soil fields
// for regular mesh. Assigns per-cell Ks, f0, k from the lookup table.
// ============================================================================
bool read_soil_map_regular(
    GlobVar& ds,
    const SoilConfig& soil_config,
    std::ofstream& logFLUXOSfile)
{
    if (!soil_config.enabled) return false;

    // If no soil map file, use a default soil type (ID=1) for all cells
    if (soil_config.soil_map_file.empty()) {
        logFLUXOSfile << "  No SOIL_MAP file; using uniform soil type 1\n";

        auto it = soil_config.soil_types.find(1);
        if (it == soil_config.soil_types.end()) {
            std::cout << "ERROR: No soil type 1 defined for uniform fallback" << std::endl;
            return true; // error
        }
        const SoilParams& sp = it->second;

        for (unsigned int icol = 1; icol <= ds.NCOLS; icol++) {
            for (unsigned int irow = 1; irow <= ds.NROWS; irow++) {
                if ((*ds.innerNeumannBCWeir)(irow, icol) == 1.0f) continue;

                (*ds.soil_type)(irow, icol) = 1;
                (*ds.soil_Ks)(irow, icol) = sp.Ks;
                (*ds.soil_f0)(irow, icol) = sp.f0;
                (*ds.soil_k)(irow, icol) = sp.k;
            }
        }
        return false;
    }

    // Read soil map from ESRI ASCII file
    std::ifstream mapfile(soil_config.soil_map_file);
    if (!mapfile.is_open()) {
        std::cout << "ERROR: Cannot open soil map file: "
                  << soil_config.soil_map_file << std::endl;
        return true;
    }

    // Skip header lines (same format as DEM)
    std::string line;
    for (int i = 0; i < 6; i++) {
        getline(mapfile, line);
    }

    // Read grid data (bottom-to-top row order, same as DEM)
    std::string temp_str;
    int temp_int;
    int irow = ds.NROWS;
    while (!mapfile.eof() && irow >= 1) {
        std::stringstream str_strm;
        getline(mapfile, line);
        str_strm << line;

        unsigned int icol = 1;
        while (!str_strm.eof() && icol <= ds.NCOLS) {
            str_strm >> temp_str;
            if (std::stringstream(temp_str) >> temp_int) {
                int soil_id = temp_int;

                if ((*ds.innerNeumannBCWeir)(irow, icol) == 1.0f) {
                    icol++;
                    temp_str = "";
                    continue;
                }

                (*ds.soil_type)(irow, icol) = soil_id;

                auto it = soil_config.soil_types.find(soil_id);
                if (it != soil_config.soil_types.end()) {
                    (*ds.soil_Ks)(irow, icol) = it->second.Ks;
                    (*ds.soil_f0)(irow, icol) = it->second.f0;
                    (*ds.soil_k)(irow, icol) = it->second.k;
                } else {
                    (*ds.soil_Ks)(irow, icol) = 0.0;
                    (*ds.soil_f0)(irow, icol) = 0.0;
                    (*ds.soil_k)(irow, icol) = 0.0;
                }
            }
            icol++;
            temp_str = "";
        }
        irow--;
    }
    mapfile.close();

    logFLUXOSfile << "  Soil map loaded successfully: "
                  << soil_config.soil_map_file << "\n";
    return false;
}

// ============================================================================
// Apply Horton soil infiltration loss to regular mesh (one timestep)
// f(t) = fc + (f0 - fc) * exp(-k * t_wet)
// where fc = Ks, f0 = initial rate, k = decay constant, t_wet = wetting time
// ============================================================================
void apply_soil_infiltration(
    GlobVar& ds,
    const SoilConfig& soil_config)
{
    if (!soil_config.enabled) return;

    const double dt = ds.dtfl;
    const double hdry = ds.hdry;

    #pragma omp parallel for collapse(2) schedule(static)
    for (unsigned int icol = 1; icol <= ds.NCOLS; icol++) {
        for (unsigned int irow = 1; irow <= ds.NROWS; irow++) {
            if ((*ds.innerNeumannBCWeir)(irow, icol) == 1.0f) continue;

            double Ks = (*ds.soil_Ks)(irow, icol);
            if (Ks <= 0.0) continue;

            double h_surface = std::fmax((*ds.z)(irow, icol) - (*ds.zb)(irow, icol), 0.0);
            if (h_surface <= 0.0) continue;

            // Horton decay: f(t) = fc + (f0 - fc) * exp(-k * t_wet)
            double f0 = (*ds.soil_f0)(irow, icol);
            double k_decay = (*ds.soil_k)(irow, icol);
            double t_wet = (*ds.soil_wetting_time)(irow, icol);

            double f_rate = Ks + (f0 - Ks) * std::exp(-k_decay * t_wet);

            // Infiltration = min(f(t) * dt, available water)
            double infil = std::fmin(f_rate * dt, h_surface);

            // Update surface water
            (*ds.z)(irow, icol) -= infil;
            (*ds.h)(irow, icol) = std::fmax((*ds.z)(irow, icol) - (*ds.zb)(irow, icol), 0.0);

            // Update dry flag
            (*ds.ldry)(irow, icol) = ((*ds.h)(irow, icol) <= hdry) ? 1.0f : 0.0f;

            // Store instantaneous infiltration rate for output
            (*ds.soil_infil_rate)(irow, icol) = infil / dt;

            // Advance wetting time (cell is wet since h_surface > 0)
            (*ds.soil_wetting_time)(irow, icol) += dt;
        }
    }
}
