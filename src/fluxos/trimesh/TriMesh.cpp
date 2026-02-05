// Copyright 2019, Diogo Costa
// Triangular mesh geometry computation and connectivity

#include "TriMesh.h"
#include <map>
#include <algorithm>
#include <cassert>
#include <cmath>

// ============================================================================
// Compute all geometry: areas, centroids, normals, inradii
// ============================================================================
void TriMesh::compute_geometry()
{
    min_area = std::numeric_limits<double>::max();
    max_area = 0.0;
    min_inradius = std::numeric_limits<double>::max();

    // Bounding box
    xmin = std::numeric_limits<double>::max();
    ymin = std::numeric_limits<double>::max();
    xmax = -std::numeric_limits<double>::max();
    ymax = -std::numeric_limits<double>::max();

    for (int i = 0; i < num_vertices; i++) {
        xmin = std::min(xmin, vertices[i].x);
        xmax = std::max(xmax, vertices[i].x);
        ymin = std::min(ymin, vertices[i].y);
        ymax = std::max(ymax, vertices[i].y);
    }

    // Cell geometry
    for (int i = 0; i < num_cells; i++) {
        auto& cell = cells[i];
        const Vertex& v0 = vertices[cell.verts[0]];
        const Vertex& v1 = vertices[cell.verts[1]];
        const Vertex& v2 = vertices[cell.verts[2]];

        // Centroid
        cell.cx = (v0.x + v1.x + v2.x) / 3.0;
        cell.cy = (v0.y + v1.y + v2.y) / 3.0;

        // Area (signed, ensures CCW orientation)
        double signed_area = 0.5 * ((v1.x - v0.x) * (v2.y - v0.y) -
                                     (v2.x - v0.x) * (v1.y - v0.y));
        cell.area = std::fabs(signed_area);

        // If area is negative, vertices are CW -> swap to CCW
        if (signed_area < 0.0) {
            std::swap(cell.verts[1], cell.verts[2]);
        }

        // Perimeter
        double e0 = std::sqrt((v1.x - v0.x) * (v1.x - v0.x) +
                               (v1.y - v0.y) * (v1.y - v0.y));
        double e1 = std::sqrt((v2.x - v1.x) * (v2.x - v1.x) +
                               (v2.y - v1.y) * (v2.y - v1.y));
        double e2 = std::sqrt((v0.x - v2.x) * (v0.x - v2.x) +
                               (v0.y - v2.y) * (v0.y - v2.y));
        double perimeter = e0 + e1 + e2;

        // Inradius = 2 * area / perimeter
        cell.inradius = (perimeter > 0.0) ? (2.0 * cell.area / perimeter) : 0.0;

        min_area = std::min(min_area, cell.area);
        max_area = std::max(max_area, cell.area);
        min_inradius = std::min(min_inradius, cell.inradius);
    }

    // Edge geometry
    double total_edge_length = 0.0;
    for (int i = 0; i < num_edges; i++) {
        auto& edge = edges[i];
        const Vertex& va = vertices[edge.v0];
        const Vertex& vb = vertices[edge.v1];

        // Length
        double dx = vb.x - va.x;
        double dy = vb.y - va.y;
        edge.length = std::sqrt(dx * dx + dy * dy);
        total_edge_length += edge.length;

        // Midpoint
        edge.mx = 0.5 * (va.x + vb.x);
        edge.my = 0.5 * (va.y + vb.y);

        // Normal: rotate edge vector 90 degrees CCW, then normalize
        // Normal points from left_cell to right_cell
        double len = edge.length;
        if (len > 0.0) {
            edge.nx = dy / len;    // outward normal x
            edge.ny = -dx / len;   // outward normal y
        }

        // Ensure normal points from left_cell outward
        if (edge.left_cell >= 0) {
            const TriCell& lc = cells[edge.left_cell];
            double dot = edge.nx * (edge.mx - lc.cx) + edge.ny * (edge.my - lc.cy);
            if (dot < 0.0) {
                edge.nx = -edge.nx;
                edge.ny = -edge.ny;
            }
        }

        // Distance between cell centroids
        if (edge.left_cell >= 0 && edge.right_cell >= 0) {
            const TriCell& lc = cells[edge.left_cell];
            const TriCell& rc = cells[edge.right_cell];
            double dcx = rc.cx - lc.cx;
            double dcy = rc.cy - lc.cy;
            edge.dist_LR = std::sqrt(dcx * dcx + dcy * dcy);
        } else if (edge.left_cell >= 0) {
            // Boundary edge: use distance from centroid to edge midpoint times 2
            const TriCell& lc = cells[edge.left_cell];
            double dcx = edge.mx - lc.cx;
            double dcy = edge.my - lc.cy;
            edge.dist_LR = 2.0 * std::sqrt(dcx * dcx + dcy * dcy);
        }
    }

    avg_edge_length = (num_edges > 0) ? (total_edge_length / num_edges) : 0.0;

    // Compute LSQ gradient matrices
    compute_lsq_matrices();
}

// ============================================================================
// Build edge-cell connectivity from cell-vertex data
// ============================================================================
void TriMesh::build_edges()
{
    // Map from sorted vertex pair -> edge index
    std::map<std::pair<int, int>, int> edge_map;

    edges.clear();
    num_edges = 0;

    for (int ci = 0; ci < num_cells; ci++) {
        auto& cell = cells[ci];

        for (int e = 0; e < 3; e++) {
            int va = cell.verts[e];
            int vb = cell.verts[(e + 1) % 3];

            // Sort vertex pair for consistent lookup
            auto key = std::make_pair(std::min(va, vb), std::max(va, vb));

            auto it = edge_map.find(key);
            if (it == edge_map.end()) {
                // New edge
                Edge edge;
                edge.v0 = va;
                edge.v1 = vb;
                edge.left_cell = ci;
                edge.right_cell = -1;  // boundary until proven otherwise

                int eidx = num_edges;
                edge_map[key] = eidx;
                edges.push_back(edge);
                cell.edges[e] = eidx;
                num_edges++;
            } else {
                // Existing edge - this is the right cell
                int eidx = it->second;
                edges[eidx].right_cell = ci;
                cell.edges[e] = eidx;
            }
        }
    }

    // Build neighbor connectivity from edges
    for (int ci = 0; ci < num_cells; ci++) {
        cells[ci].neighbors = {-1, -1, -1};
    }

    num_boundary_edges = 0;
    num_internal_edges = 0;

    for (int ei = 0; ei < num_edges; ei++) {
        auto& edge = edges[ei];

        if (edge.right_cell >= 0) {
            // Internal edge
            num_internal_edges++;
            edge.is_boundary = false;

            // Set neighbor pointers
            int lc = edge.left_cell;
            int rc = edge.right_cell;

            // Find which local edge slot to use
            for (int e = 0; e < 3; e++) {
                if (cells[lc].edges[e] == ei) {
                    cells[lc].neighbors[e] = rc;
                }
                if (cells[rc].edges[e] == ei) {
                    cells[rc].neighbors[e] = lc;
                }
            }
        } else {
            // Boundary edge
            num_boundary_edges++;
            edge.is_boundary = true;
            cells[edge.left_cell].is_boundary = true;
        }
    }
}

// ============================================================================
// Compute least-squares gradient matrices for all cells
// Uses inverse-distance weighted least-squares
// ============================================================================
void TriMesh::compute_lsq_matrices()
{
    for (int ci = 0; ci < num_cells; ci++) {
        auto& cell = cells[ci];

        // Build 2x2 matrix: A^T * W * A where A has rows [dx_j, dy_j]
        // and W is inverse-distance weighting
        double a00 = 0.0, a01 = 0.0, a11 = 0.0;

        for (int e = 0; e < 3; e++) {
            int nj = cell.neighbors[e];
            double dx, dy, w;

            if (nj >= 0) {
                // Internal neighbor
                dx = cells[nj].cx - cell.cx;
                dy = cells[nj].cy - cell.cy;
                double dist = std::sqrt(dx * dx + dy * dy);
                w = (dist > 1e-15) ? (1.0 / dist) : 0.0;
            } else {
                // Boundary edge: use edge midpoint as ghost
                int eidx = cell.edges[e];
                dx = edges[eidx].mx - cell.cx;
                dy = edges[eidx].my - cell.cy;
                double dist = std::sqrt(dx * dx + dy * dy);
                w = (dist > 1e-15) ? (1.0 / dist) : 0.0;
            }

            a00 += w * dx * dx;
            a01 += w * dx * dy;
            a11 += w * dy * dy;
        }

        // Invert 2x2 symmetric matrix: [a00, a01; a01, a11]
        double det = a00 * a11 - a01 * a01;
        if (std::fabs(det) > 1e-30) {
            double inv_det = 1.0 / det;
            cell.lsq_inv[0][0] = a11 * inv_det;
            cell.lsq_inv[0][1] = -a01 * inv_det;
            cell.lsq_inv[1][0] = -a01 * inv_det;
            cell.lsq_inv[1][1] = a00 * inv_det;
        } else {
            // Degenerate - zero gradient
            cell.lsq_inv[0][0] = 0.0;
            cell.lsq_inv[0][1] = 0.0;
            cell.lsq_inv[1][0] = 0.0;
            cell.lsq_inv[1][1] = 0.0;
        }
    }
}

// ============================================================================
// Validate mesh integrity
// ============================================================================
bool TriMesh::validate() const
{
    bool valid = true;

    // Check vertex indices in cells
    for (int ci = 0; ci < num_cells; ci++) {
        for (int e = 0; e < 3; e++) {
            if (cells[ci].verts[e] < 0 || cells[ci].verts[e] >= num_vertices) {
                std::cerr << "TriMesh::validate: Cell " << ci
                          << " has invalid vertex index " << cells[ci].verts[e] << std::endl;
                valid = false;
            }
        }
        if (cells[ci].area <= 0.0) {
            std::cerr << "TriMesh::validate: Cell " << ci
                      << " has non-positive area " << cells[ci].area << std::endl;
            valid = false;
        }
    }

    // Check edge-cell connectivity
    for (int ei = 0; ei < num_edges; ei++) {
        if (edges[ei].left_cell < 0 || edges[ei].left_cell >= num_cells) {
            std::cerr << "TriMesh::validate: Edge " << ei
                      << " has invalid left_cell " << edges[ei].left_cell << std::endl;
            valid = false;
        }
        if (edges[ei].right_cell >= num_cells) {
            std::cerr << "TriMesh::validate: Edge " << ei
                      << " has invalid right_cell " << edges[ei].right_cell << std::endl;
            valid = false;
        }
    }

    return valid;
}

// ============================================================================
// Print mesh statistics
// ============================================================================
void TriMesh::print_info(std::ostream& os) const
{
    os << "\n=== Triangular Mesh Info ===" << std::endl;
    os << "  Vertices       : " << num_vertices << std::endl;
    os << "  Cells          : " << num_cells << std::endl;
    os << "  Edges          : " << num_edges << std::endl;
    os << "    Internal     : " << num_internal_edges << std::endl;
    os << "    Boundary     : " << num_boundary_edges << std::endl;
    os << "  Domain         : [" << xmin << ", " << xmax << "] x ["
       << ymin << ", " << ymax << "]" << std::endl;
    os << "  Cell area      : min=" << min_area << "  max=" << max_area << std::endl;
    os << "  Min inradius   : " << min_inradius << std::endl;
    os << "  Avg edge length: " << avg_edge_length << std::endl;
    os << "==========================\n" << std::endl;
}

// ============================================================================
// Point-in-triangle: find which cell contains point (px, py)
// Brute force O(N) - suitable for initialization, not per-timestep
// ============================================================================
int TriMesh::find_cell(double px, double py) const
{
    for (int ci = 0; ci < num_cells; ci++) {
        if (point_in_triangle(px, py, ci)) {
            return ci;
        }
    }
    return -1;  // not found
}

// ============================================================================
// Get local edge index
// ============================================================================
int TriMesh::local_edge_index(int cell_id, int edge_id) const
{
    for (int e = 0; e < 3; e++) {
        if (cells[cell_id].edges[e] == edge_id) return e;
    }
    return -1;
}

// ============================================================================
// Clear
// ============================================================================
void TriMesh::clear()
{
    vertices.clear();
    cells.clear();
    edges.clear();
    num_vertices = 0;
    num_cells = 0;
    num_edges = 0;
    num_boundary_edges = 0;
    num_internal_edges = 0;
}

// ============================================================================
// Triangle area from 3 vertex indices
// ============================================================================
double TriMesh::triangle_area(int v0, int v1, int v2) const
{
    const Vertex& a = vertices[v0];
    const Vertex& b = vertices[v1];
    const Vertex& c = vertices[v2];
    return 0.5 * std::fabs((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y));
}

// ============================================================================
// Point-in-triangle test using barycentric coordinates
// ============================================================================
bool TriMesh::point_in_triangle(double px, double py, int cell_id) const
{
    const TriCell& cell = cells[cell_id];
    const Vertex& v0 = vertices[cell.verts[0]];
    const Vertex& v1 = vertices[cell.verts[1]];
    const Vertex& v2 = vertices[cell.verts[2]];

    double d1 = (px - v1.x) * (v0.y - v1.y) - (v0.x - v1.x) * (py - v1.y);
    double d2 = (px - v2.x) * (v1.y - v2.y) - (v1.x - v2.x) * (py - v2.y);
    double d3 = (px - v0.x) * (v2.y - v0.y) - (v2.x - v0.x) * (py - v0.y);

    bool has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
    bool has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

    return !(has_neg && has_pos);
}
