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
// Soil infiltration loss model for FLUXOS
// Horton decay: f(t) = fc + (f0 - fc) * exp(-k * t_wet)
// f0 = initial infiltration rate, fc = final rate (≈ Ks), k = decay constant
// t_wet = cumulative wetting time per cell

#ifndef SOIL_INFILTRATION_H_INCLUDED
#define SOIL_INFILTRATION_H_INCLUDED

#include <string>
#include <vector>
#include <map>
#include "jnlohmann/json.hpp"

class GlobVar;

// ============================================================================
// SoilParams: Horton infiltration parameters for a single soil type
// ============================================================================
struct SoilParams {
    std::string name;       // Soil type name (e.g., "Sand", "Clay")
    double Ks;              // Saturated hydraulic conductivity [m/s] = final rate (fc)
    double f0;              // Initial (maximum) infiltration rate [m/s]
    double k;               // Horton decay constant [1/s]
};

// ============================================================================
// SoilConfig: runtime configuration for the soil infiltration module
// ============================================================================
struct SoilConfig {
    bool enabled = false;
    std::string soil_map_file;          // Path to soil type raster (.asc)
    std::map<int, SoilParams> soil_types; // Lookup: soil type ID -> parameters
};

// ============================================================================
// USDA soil texture classes with Horton infiltration parameters
// Ks (fc): Rawls, Brakensiek & Miller (1983), J. Hydraulic Eng.
// f0: estimated as ~5x Ks (typical for initially dry soil)
// k: decay constant calibrated per texture (coarser soils decay faster)
// ============================================================================
inline std::map<std::string, SoilParams> get_usda_soil_table()
{
    std::map<std::string, SoilParams> table;

    //                                 name              Ks=fc [m/s]   f0 [m/s]     k [1/s]
    table["sand"]            = {"Sand",              1.176e-4, 5.880e-4, 5.56e-4};
    table["loamy_sand"]      = {"Loamy Sand",        2.988e-5, 1.494e-4, 4.17e-4};
    table["sandy_loam"]      = {"Sandy Loam",        1.089e-5, 5.445e-5, 2.78e-4};
    table["loam"]            = {"Loam",              3.396e-6, 1.698e-5, 1.67e-4};
    table["silt_loam"]       = {"Silt Loam",         1.889e-6, 9.445e-6, 1.11e-4};
    table["sandy_clay_loam"] = {"Sandy Clay Loam",   1.186e-6, 5.930e-6, 8.33e-5};
    table["clay_loam"]       = {"Clay Loam",         6.389e-7, 3.195e-6, 5.56e-5};
    table["silty_clay_loam"] = {"Silty Clay Loam",   4.167e-7, 2.084e-6, 4.17e-5};
    table["sandy_clay"]      = {"Sandy Clay",        3.333e-7, 1.667e-6, 3.33e-5};
    table["silty_clay"]      = {"Silty Clay",        2.500e-7, 1.250e-6, 2.78e-5};
    table["clay"]            = {"Clay",              1.667e-7, 8.335e-7, 1.67e-5};
    table["impervious"]      = {"Impervious",        0.0,      0.0,      0.0};

    return table;
}

// ============================================================================
// Configuration parsing
// ============================================================================
bool read_soil_config(
    const nlohmann::json& master_json,
    SoilConfig& soil_config,
    std::ofstream& logFLUXOSfile);

// ============================================================================
// Soil map reading (ESRI ASCII grid format, same as DEM)
// Returns per-cell soil type IDs for the regular mesh
// ============================================================================
bool read_soil_map_regular(
    GlobVar& ds,
    const SoilConfig& soil_config,
    std::ofstream& logFLUXOSfile);

// ============================================================================
// Regular mesh: apply Horton soil infiltration loss for one timestep
// f(t) = fc + (f0 - fc) * exp(-k * t_wet); infil = min(f(t) * dt, h_surface)
// Called after forcing (rainfall/inflow), before hydrodynamic solver
// ============================================================================
void apply_soil_infiltration(
    GlobVar& ds,
    const SoilConfig& soil_config);

#endif // SOIL_INFILTRATION_H_INCLUDED
