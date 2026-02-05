// Copyright 2019, Diogo Costa
// Triangular mesh data structures for FLUXOS unstructured mesh support
// Edge-based finite volume method on unstructured triangular grids

#ifndef TRIMESH_H_INCLUDED
#define TRIMESH_H_INCLUDED

#include <vector>
#include <array>
#include <string>
#include <cmath>
#include <limits>
#include <iostream>

// ============================================================================
// Vertex: a point in the mesh
// ============================================================================
struct Vertex {
    double x = 0.0;
    double y = 0.0;
};

// ============================================================================
// TriCell: a triangular cell
// ============================================================================
struct TriCell {
    std::array<int, 3> verts = {-1, -1, -1};       // vertex indices (CCW order)
    std::array<int, 3> edges = {-1, -1, -1};        // edge indices
    std::array<int, 3> neighbors = {-1, -1, -1};    // neighbor cell indices (-1 = boundary)

    // Geometry (precomputed)
    double cx = 0.0, cy = 0.0;     // centroid
    double area = 0.0;             // cell area
    double inradius = 0.0;         // inscribed circle radius = 2*area/perimeter

    // Least-squares gradient matrix (precomputed 2x2 inverse)
    // Used for gradient reconstruction: [grad_x, grad_y]^T = lsq_inv * rhs
    double lsq_inv[2][2] = {{0.0, 0.0}, {0.0, 0.0}};

    // Boundary flags
    bool is_boundary = false;       // true if cell has at least one boundary edge
    int boundary_tag = 0;           // boundary condition tag (from mesh file)
};

// ============================================================================
// Edge: connects two cells (or one cell + boundary)
// ============================================================================
struct Edge {
    int left_cell = -1;     // cell index on left side (always valid)
    int right_cell = -1;    // cell index on right side (-1 = boundary)

    int v0 = -1, v1 = -1;  // vertex indices defining the edge

    // Normal vector (unit, points from left to right cell)
    double nx = 0.0, ny = 0.0;

    // Edge properties
    double length = 0.0;
    double mx = 0.0, my = 0.0;     // midpoint

    // Distance between cell centroids (for diffusion)
    double dist_LR = 0.0;

    // Boundary
    int boundary_tag = 0;   // 0 = internal, >0 = boundary tag from mesh file
    bool is_boundary = false;
};

// ============================================================================
// TriMesh: complete unstructured triangular mesh
// ============================================================================
class TriMesh {
public:
    // Mesh data
    std::vector<Vertex> vertices;
    std::vector<TriCell> cells;
    std::vector<Edge> edges;

    // Counts
    int num_vertices = 0;
    int num_cells = 0;
    int num_edges = 0;
    int num_boundary_edges = 0;
    int num_internal_edges = 0;

    // Domain bounding box
    double xmin = 0.0, xmax = 0.0;
    double ymin = 0.0, ymax = 0.0;

    // Mesh quality metrics
    double min_area = 0.0;
    double max_area = 0.0;
    double min_inradius = 0.0;
    double avg_edge_length = 0.0;

    // ---- Methods ----

    // Compute all geometry: areas, centroids, normals, inradii, LSQ matrices
    void compute_geometry();

    // Build edge-cell connectivity from cell-vertex data
    void build_edges();

    // Compute least-squares gradient matrices for all cells
    void compute_lsq_matrices();

    // Validate mesh integrity
    bool validate() const;

    // Print mesh statistics
    void print_info(std::ostream& os = std::cout) const;

    // Point-in-triangle lookup: find which cell contains point (px, py)
    // Returns cell index or -1 if not found
    int find_cell(double px, double py) const;

    // Get the local edge index (0,1,2) of a given global edge in a cell
    int local_edge_index(int cell_id, int edge_id) const;

    // Clear all data
    void clear();

private:
    // Compute area of triangle from 3 vertices
    double triangle_area(int v0, int v1, int v2) const;

    // Check if point is inside triangle
    bool point_in_triangle(double px, double py, int cell_id) const;
};

#endif // TRIMESH_H_INCLUDED
