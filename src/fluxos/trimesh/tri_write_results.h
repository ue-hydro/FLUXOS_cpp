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
