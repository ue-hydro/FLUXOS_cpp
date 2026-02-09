// Copyright 2019, Diogo Costa
// Mesh I/O for triangular meshes: Gmsh .msh v2.2 and Triangle .node/.ele/.edge

#ifndef TRI_MESH_IO_H_INCLUDED
#define TRI_MESH_IO_H_INCLUDED

#include <string>
#include "TriMesh.h"

// Read Gmsh .msh v2.2 format
// Returns true on success
bool read_gmsh_msh(const std::string& filename, TriMesh& mesh);

// Read Triangle format (.node + .ele files)
// base_filename should be without extension (e.g., "domain" for domain.node, domain.ele)
bool read_triangle_format(const std::string& base_filename, TriMesh& mesh);

// Auto-detect format and read
bool read_mesh_file(const std::string& filename, const std::string& format, TriMesh& mesh);

#endif // TRI_MESH_IO_H_INCLUDED
