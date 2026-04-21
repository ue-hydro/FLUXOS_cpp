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
