// Copyright 2019, Diogo Costa
// Mesh I/O: Gmsh .msh v2.2 and Triangle .node/.ele readers

#include "tri_mesh_io.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <map>
#include <cmath>

// ============================================================================
// Read Gmsh .msh v2.2 format
// ============================================================================
bool read_gmsh_msh(const std::string& filename, TriMesh& mesh)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "ERROR: Cannot open mesh file: " << filename << std::endl;
        return false;
    }

    mesh.clear();

    std::string line;
    // Map from gmsh node ID (1-based, possibly non-contiguous) to internal index
    std::map<int, int> node_id_map;

    // Temporary storage for boundary line elements
    struct BoundaryLine {
        int v0, v1;
        int physical_tag;
    };
    std::vector<BoundaryLine> boundary_lines;

    while (std::getline(file, line)) {
        // Trim whitespace
        line.erase(0, line.find_first_not_of(" \t\r\n"));

        // ---- Parse $Nodes section ----
        if (line == "$Nodes") {
            int num_nodes;
            file >> num_nodes;
            mesh.vertices.resize(num_nodes);
            mesh.num_vertices = num_nodes;

            bool any_nonzero_z = false;
            for (int i = 0; i < num_nodes; i++) {
                int node_id;
                double x, y, z_val;
                file >> node_id >> x >> y >> z_val;

                node_id_map[node_id] = i;
                mesh.vertices[i].x = x;
                mesh.vertices[i].y = y;
                mesh.vertices[i].z = z_val;

                if (std::fabs(z_val) > 1e-10) {
                    any_nonzero_z = true;
                }
            }
            // If the mesh file contained non-zero z values, mark that vertex
            // elevations are available (Python preprocessing wrote DEM z into vertices)
            mesh.has_vertex_elevations = any_nonzero_z;
            std::getline(file, line); // consume rest of line
            std::getline(file, line); // $EndNodes
            continue;
        }

        // ---- Parse $Elements section ----
        if (line == "$Elements") {
            int num_elements;
            file >> num_elements;

            for (int i = 0; i < num_elements; i++) {
                int elem_id, elem_type, num_tags;
                file >> elem_id >> elem_type >> num_tags;

                std::vector<int> tags(num_tags);
                for (int t = 0; t < num_tags; t++) {
                    file >> tags[t];
                }

                if (elem_type == 1) {
                    // Line element (2 nodes) - boundary edge
                    int n0, n1;
                    file >> n0 >> n1;
                    BoundaryLine bl;
                    bl.v0 = node_id_map[n0];
                    bl.v1 = node_id_map[n1];
                    bl.physical_tag = (num_tags > 0) ? tags[0] : 0;
                    boundary_lines.push_back(bl);
                }
                else if (elem_type == 2) {
                    // Triangle element (3 nodes)
                    int n0, n1, n2;
                    file >> n0 >> n1 >> n2;
                    TriCell cell;
                    cell.verts[0] = node_id_map[n0];
                    cell.verts[1] = node_id_map[n1];
                    cell.verts[2] = node_id_map[n2];
                    cell.boundary_tag = (num_tags > 0) ? tags[0] : 0;
                    mesh.cells.push_back(cell);
                }
                else {
                    // Skip other element types (points, quads, etc.)
                    // Read remaining nodes for this element
                    // Element node counts: 15=1node, 1=2nodes, 2=3nodes, 3=4nodes, etc.
                    int nnodes = 0;
                    switch (elem_type) {
                        case 15: nnodes = 1; break; // point
                        case 1:  nnodes = 2; break; // line
                        case 2:  nnodes = 3; break; // triangle
                        case 3:  nnodes = 4; break; // quad
                        case 4:  nnodes = 4; break; // tet
                        case 5:  nnodes = 8; break; // hex
                        default: nnodes = 0; break;
                    }
                    // Already read for types 1 and 2, skip for others
                    if (elem_type != 1 && elem_type != 2) {
                        for (int n = 0; n < nnodes; n++) {
                            int dummy;
                            file >> dummy;
                        }
                    }
                }
            }
            std::getline(file, line); // consume rest
            std::getline(file, line); // $EndElements
            continue;
        }
    }

    file.close();

    mesh.num_cells = static_cast<int>(mesh.cells.size());

    if (mesh.num_cells == 0) {
        std::cerr << "ERROR: No triangular elements found in " << filename << std::endl;
        return false;
    }

    // Build edges
    mesh.build_edges();

    // Apply boundary tags from line elements to edges
    // Create map from sorted vertex pair to boundary tag
    std::map<std::pair<int, int>, int> boundary_tag_map;
    for (const auto& bl : boundary_lines) {
        auto key = std::make_pair(std::min(bl.v0, bl.v1), std::max(bl.v0, bl.v1));
        boundary_tag_map[key] = bl.physical_tag;
    }

    for (int ei = 0; ei < mesh.num_edges; ei++) {
        auto& edge = mesh.edges[ei];
        auto key = std::make_pair(std::min(edge.v0, edge.v1), std::max(edge.v0, edge.v1));
        auto it = boundary_tag_map.find(key);
        if (it != boundary_tag_map.end()) {
            edge.boundary_tag = it->second;
        }
    }

    // Compute geometry
    mesh.compute_geometry();

    std::cout << "Successfully loaded Gmsh mesh: " << filename << std::endl;
    return true;
}

// ============================================================================
// Read Triangle format (.node + .ele files)
// ============================================================================
bool read_triangle_format(const std::string& base_filename, TriMesh& mesh)
{
    mesh.clear();

    // ---- Read .node file ----
    std::string node_file = base_filename + ".node";
    std::ifstream nf(node_file);
    if (!nf.is_open()) {
        std::cerr << "ERROR: Cannot open node file: " << node_file << std::endl;
        return false;
    }

    int num_nodes, dim, num_attr, num_boundary_markers;
    nf >> num_nodes >> dim >> num_attr >> num_boundary_markers;

    mesh.vertices.resize(num_nodes);
    mesh.num_vertices = num_nodes;

    std::map<int, int> node_id_map;  // external ID -> internal index

    for (int i = 0; i < num_nodes; i++) {
        int node_id;
        double x, y;
        nf >> node_id >> x >> y;

        // Skip attributes and boundary markers
        for (int a = 0; a < num_attr; a++) {
            double dummy;
            nf >> dummy;
        }
        for (int b = 0; b < num_boundary_markers; b++) {
            int dummy;
            nf >> dummy;
        }

        node_id_map[node_id] = i;
        mesh.vertices[i].x = x;
        mesh.vertices[i].y = y;
    }
    nf.close();

    // ---- Read .ele file ----
    std::string ele_file = base_filename + ".ele";
    std::ifstream ef(ele_file);
    if (!ef.is_open()) {
        std::cerr << "ERROR: Cannot open element file: " << ele_file << std::endl;
        return false;
    }

    int num_triangles, nodes_per_tri, num_tri_attr;
    ef >> num_triangles >> nodes_per_tri >> num_tri_attr;

    if (nodes_per_tri != 3) {
        std::cerr << "ERROR: Only 3-node triangles supported, got " << nodes_per_tri << std::endl;
        return false;
    }

    mesh.cells.resize(num_triangles);
    mesh.num_cells = num_triangles;

    for (int i = 0; i < num_triangles; i++) {
        int tri_id, n0, n1, n2;
        ef >> tri_id >> n0 >> n1 >> n2;

        mesh.cells[i].verts[0] = node_id_map[n0];
        mesh.cells[i].verts[1] = node_id_map[n1];
        mesh.cells[i].verts[2] = node_id_map[n2];

        // Read optional attributes
        if (num_tri_attr > 0) {
            double attr;
            ef >> attr;
            mesh.cells[i].boundary_tag = static_cast<int>(attr);
        }
    }
    ef.close();

    // Build edges and compute geometry
    mesh.build_edges();

    // ---- Read optional .edge file for boundary tags ----
    std::string edge_file = base_filename + ".edge";
    std::ifstream edf(edge_file);
    if (edf.is_open()) {
        int num_file_edges, has_boundary_markers;
        edf >> num_file_edges >> has_boundary_markers;

        std::map<std::pair<int, int>, int> boundary_tag_map;
        for (int i = 0; i < num_file_edges; i++) {
            int edge_id, ev0, ev1;
            int bmarker = 0;
            edf >> edge_id >> ev0 >> ev1;
            if (has_boundary_markers) {
                edf >> bmarker;
            }
            if (bmarker != 0) {
                int iv0 = node_id_map[ev0];
                int iv1 = node_id_map[ev1];
                auto key = std::make_pair(std::min(iv0, iv1), std::max(iv0, iv1));
                boundary_tag_map[key] = bmarker;
            }
        }
        edf.close();

        for (int ei = 0; ei < mesh.num_edges; ei++) {
            auto& edge = mesh.edges[ei];
            auto key = std::make_pair(std::min(edge.v0, edge.v1), std::max(edge.v0, edge.v1));
            auto it = boundary_tag_map.find(key);
            if (it != boundary_tag_map.end()) {
                edge.boundary_tag = it->second;
            }
        }
    }

    mesh.compute_geometry();

    std::cout << "Successfully loaded Triangle mesh: " << base_filename << std::endl;
    return true;
}

// ============================================================================
// Auto-detect format and read
// ============================================================================
bool read_mesh_file(const std::string& filename, const std::string& format, TriMesh& mesh)
{
    if (format == "gmsh" || format == "msh") {
        return read_gmsh_msh(filename, mesh);
    } else if (format == "triangle" || format == "tri") {
        return read_triangle_format(filename, mesh);
    } else {
        // Try to detect from extension
        if (filename.size() > 4 && filename.substr(filename.size() - 4) == ".msh") {
            return read_gmsh_msh(filename, mesh);
        } else if (filename.size() > 5 && filename.substr(filename.size() - 5) == ".node") {
            std::string base = filename.substr(0, filename.size() - 5);
            return read_triangle_format(base, mesh);
        } else {
            // Default to gmsh
            return read_gmsh_msh(filename, mesh);
        }
    }
}
