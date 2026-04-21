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
// VTK XML Unstructured Grid (.vtu) output for triangular mesh
// + .pvd time series file for ParaView

#ifndef TRI_WRITE_RESULTS_H_INCLUDED
#define TRI_WRITE_RESULTS_H_INCLUDED

#include "../GlobVar.h"
#include "TriMesh.h"
#include "TriSolution.h"
#include <chrono>
#include <string>
#include <vector>

// Write a single timestep to VTK .vtu file
bool tri_write_results(
    GlobVar& ds,
    const TriMesh& mesh,
    const TriSolution& sol,
    int print_tag,
    std::chrono::duration<double> elapsed_seconds);

// Write .pvd time series collection file
void tri_write_pvd(
    const std::string& output_folder,
    const std::vector<std::pair<double, std::string>>& time_files);

#endif // TRI_WRITE_RESULTS_H_INCLUDED
