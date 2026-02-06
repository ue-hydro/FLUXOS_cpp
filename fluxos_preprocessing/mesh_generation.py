#!/usr/bin/env python3
"""
Slope-based adaptive triangular mesh generation for FLUXOS.
Uses Gmsh (via pygmsh) to create meshes where higher slopes produce finer resolution.
DEM elevations are embedded as vertex z-coordinates in the output .msh file.
"""

import os
import sys
import tempfile
import numpy as np


def compute_slope(elevation, cellsize, nodata_value):
    """
    Compute slope magnitude |grad(z)| = sqrt((dz/dx)^2 + (dz/dy)^2).

    Parameters
    ----------
    elevation : numpy.ndarray
        2D elevation array.
    cellsize : float
        Cell size in CRS units (meters).
    nodata_value : float
        NODATA value.

    Returns
    -------
    slope : numpy.ndarray
        Slope magnitude array (same shape as elevation).
    """
    elev = elevation.astype(np.float64).copy()
    mask = np.isclose(elev, nodata_value)
    elev[mask] = np.nan

    # numpy.gradient returns [dz/dy, dz/dx] for 2D array
    dz_dy, dz_dx = np.gradient(elev, cellsize)

    slope = np.sqrt(dz_dx**2 + dz_dy**2)
    slope[np.isnan(slope)] = 0.0

    return slope


def build_domain_boundary(elevation, transform, nodata_value):
    """
    Extract the outer boundary polygon of valid (non-NODATA) DEM cells.

    Parameters
    ----------
    elevation : numpy.ndarray
        2D elevation array.
    transform : affine.Affine
        Affine transform (pixel to CRS coordinates).
    nodata_value : float
        NODATA value.

    Returns
    -------
    boundary_coords : list of (float, float)
        Ordered list of (x, y) boundary coordinates in CRS units.
    """
    import rasterio.features
    from shapely.geometry import shape, MultiPolygon
    from shapely.ops import unary_union

    # Create binary mask: 1 = valid data, 0 = nodata
    valid_mask = (~np.isclose(elevation, nodata_value) & (elevation != 0.0)).astype(np.uint8)

    # Vectorize the mask to polygons
    shapes_gen = rasterio.features.shapes(valid_mask, transform=transform)
    polygons = []
    for geom, value in shapes_gen:
        if value == 1:
            poly = shape(geom)
            if poly.is_valid:
                polygons.append(poly)

    if not polygons:
        print("ERROR: No valid DEM cells found for boundary extraction")
        sys.exit(1)

    # Merge all valid polygons into one
    merged = unary_union(polygons)

    # Get the largest polygon if MultiPolygon
    if isinstance(merged, MultiPolygon):
        merged = max(merged.geoms, key=lambda p: p.area)
        print("  WARNING: Multiple disconnected valid regions found, using largest")

    # Simplify boundary to reduce vertex count (tolerance = 1 pixel)
    pixel_size = abs(transform.a)
    simplified = merged.simplify(pixel_size * 0.5, preserve_topology=True)

    # Extract exterior ring coordinates
    coords = list(simplified.exterior.coords)

    # Remove duplicate last point (closed ring)
    if len(coords) > 1 and coords[0] == coords[-1]:
        coords = coords[:-1]

    print(f"  Domain boundary: {len(coords)} vertices")
    print(f"  Domain area: {simplified.area:.0f} m^2")

    return coords


def create_size_field(slope, min_size, max_size, slope_factor):
    """
    Map normalized slope to mesh element size.

    Higher slope -> smaller elements (finer mesh).
    Lower slope -> larger elements (coarser mesh).

    size(x,y) = max_size - slope_norm(x,y) * (max_size - min_size) * slope_factor

    Parameters
    ----------
    slope : numpy.ndarray
        Slope magnitude array.
    min_size : float
        Minimum element size (finest mesh, at steepest slopes).
    max_size : float
        Maximum element size (coarsest mesh, at flat areas).
    slope_factor : float
        Aggressiveness of slope-based refinement (0 = uniform, 1 = full range).

    Returns
    -------
    size_field : numpy.ndarray
        Element size field (same shape as slope).
    """
    slope_max = np.nanmax(slope)
    if slope_max > 0:
        slope_norm = slope / slope_max
    else:
        slope_norm = np.zeros_like(slope)
        print("  WARNING: Slope is zero everywhere, using uniform mesh size")

    factor = min(max(slope_factor, 0.0), 1.0)
    size_field = max_size - slope_norm * (max_size - min_size) * factor
    size_field = np.clip(size_field, min_size, max_size)

    print(f"  Size field range: {size_field.min():.2f} to {size_field.max():.2f}")

    return size_field


def _write_pos_file(filepath, size_field, transform):
    """
    Write a Gmsh .pos file (text-based post-processing view) containing the size field.
    This can be loaded as a background mesh field in Gmsh.

    Parameters
    ----------
    filepath : str
        Output .pos file path.
    size_field : numpy.ndarray
        2D size field array.
    transform : affine.Affine
        Affine transform for the grid.
    """
    nrows, ncols = size_field.shape

    with open(filepath, 'w') as f:
        f.write('View "size_field" {\n')

        # Write as scalar triangles (ST) covering each pixel pair
        # Each pixel (i,j) creates a quad, split into 2 triangles
        for i in range(nrows - 1):
            for j in range(ncols - 1):
                # Four corners of the pixel quad
                x0 = transform.c + j * transform.a
                y0 = transform.f + i * transform.e
                x1 = x0 + transform.a
                y1 = y0 + transform.e

                s00 = size_field[i, j]
                s10 = size_field[i, j + 1]
                s01 = size_field[i + 1, j]
                s11 = size_field[i + 1, j + 1]

                # Skip if any value is invalid
                if np.isnan(s00) or np.isnan(s10) or np.isnan(s01) or np.isnan(s11):
                    continue

                # Triangle 1: (x0,y0), (x1,y0), (x0,y1)
                f.write(f'ST({x0},{y0},0,{x1},{y0},0,{x0},{y1},0){{{s00},{s10},{s01}}};\n')
                # Triangle 2: (x1,y0), (x1,y1), (x0,y1)
                f.write(f'ST({x1},{y0},0,{x1},{y1},0,{x0},{y1},0){{{s10},{s11},{s01}}};\n')

        f.write('};\n')


def generate_adaptive_mesh(boundary_coords, size_field, transform, elevation,
                           nodata_value, min_size, max_size, output_path):
    """
    Generate an adaptive triangular mesh using Gmsh.

    Steps:
    1. Create boundary polygon in Gmsh
    2. Set up size field from slope data
    3. Generate 2D triangular mesh
    4. Interpolate DEM elevations to mesh vertices
    5. Write .msh v2.2 with z-coordinates

    Parameters
    ----------
    boundary_coords : list of (float, float)
        Domain boundary polygon vertices.
    size_field : numpy.ndarray
        Element size field array.
    transform : affine.Affine
        Affine transform for the size field / DEM grid.
    elevation : numpy.ndarray
        2D elevation array.
    nodata_value : float
        NODATA value.
    min_size : float
        Minimum element size.
    max_size : float
        Maximum element size.
    output_path : str
        Output .msh file path.
    """
    import gmsh
    from scipy.interpolate import RegularGridInterpolator

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.model.add("fluxos_domain")

    print("\n  Creating Gmsh geometry...")

    # 1. Create boundary polygon
    lc = max_size  # default characteristic length (overridden by background field)
    points = []
    for x, y in boundary_coords:
        tag = gmsh.model.geo.addPoint(x, y, 0.0, lc)
        points.append(tag)

    lines = []
    for i in range(len(points)):
        j = (i + 1) % len(points)
        tag = gmsh.model.geo.addLine(points[i], points[j])
        lines.append(tag)

    loop = gmsh.model.geo.addCurveLoop(lines)
    surface = gmsh.model.geo.addPlaneSurface([loop])
    gmsh.model.geo.synchronize()

    # Add physical groups for boundary and surface
    # All boundary lines get tag 1 (wall by default)
    gmsh.model.addPhysicalGroup(1, lines, tag=1, name="boundary")
    gmsh.model.addPhysicalGroup(2, [surface], tag=1, name="domain")

    # 2. Set up background mesh size field from the slope-based size data
    print("  Writing background size field...")

    pos_tmpfile = tempfile.mktemp(suffix='.pos')
    try:
        _write_pos_file(pos_tmpfile, size_field, transform)

        # Load the .pos file as a view
        gmsh.merge(pos_tmpfile)

        # Create PostView background field
        bg_field = gmsh.model.mesh.field.add("PostView")
        gmsh.model.mesh.field.setNumber(bg_field, "ViewIndex", 0)
        gmsh.model.mesh.field.setAsBackgroundMesh(bg_field)

    except Exception as e:
        print(f"  WARNING: Could not create background size field ({e})")
        print(f"  Falling back to uniform mesh size: {max_size}")
        # Fallback: set uniform size
        gmsh.model.mesh.field.remove(bg_field) if 'bg_field' in dir() else None

    finally:
        if os.path.exists(pos_tmpfile):
            os.remove(pos_tmpfile)

    # Disable default size constraints so background field controls everything
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    # Set algorithm (Delaunay = 5, Frontal-Delaunay = 6)
    gmsh.option.setNumber("Mesh.Algorithm", 6)

    # 3. Generate 2D mesh
    print("  Generating mesh...")
    gmsh.model.mesh.generate(2)

    # 4. Extract mesh data and interpolate DEM elevations
    print("  Interpolating DEM elevations to mesh vertices...")

    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    # node_coords is flat: [x0, y0, z0, x1, y1, z1, ...]

    # Build interpolator for elevation
    nrows, ncols = elevation.shape
    x_coords = np.array([transform.c + j * transform.a + transform.a / 2.0
                         for j in range(ncols)])
    y_coords = np.array([transform.f + i * transform.e + transform.e / 2.0
                         for i in range(nrows)])

    # Rasterio convention: transform.e is negative (y decreases downward)
    # RegularGridInterpolator needs monotonically increasing axes
    # y_coords goes from top to bottom (decreasing), so reverse both
    elev_clean = elevation.astype(np.float64).copy()
    elev_clean[np.isclose(elev_clean, nodata_value)] = np.nan
    elev_clean[elev_clean == 0.0] = np.nan

    # Reverse y so it's increasing
    y_increasing = y_coords[::-1]
    elev_increasing = elev_clean[::-1, :]

    interp = RegularGridInterpolator(
        (y_increasing, x_coords),
        elev_increasing,
        method='linear',
        bounds_error=False,
        fill_value=np.nan
    )

    # Set z for each vertex
    num_nodes = len(node_tags)
    node_coords_array = np.array(node_coords).reshape(-1, 3)

    for i in range(num_nodes):
        x = node_coords_array[i, 0]
        y = node_coords_array[i, 1]
        z = float(interp((y, x)))

        if np.isnan(z) or z == 0.0:
            # Try nearest neighbor for boundary/edge vertices
            # Find nearest valid cell
            col = int((x - transform.c) / transform.a)
            row = int((y - transform.f) / transform.e)
            col = max(0, min(col, ncols - 1))
            row = max(0, min(row, nrows - 1))
            z_nn = elevation[row, col]
            if not np.isclose(z_nn, nodata_value) and z_nn != 0.0:
                z = abs(float(z_nn))
            else:
                z = 0.0  # will be handled by C++ NODATA logic

        node_coords_array[i, 2] = abs(z)

    # Update gmsh node coordinates with z values
    gmsh.model.mesh.setNodes(
        dim=2,
        tag=surface,
        nodeTags=node_tags,
        coord=node_coords_array.flatten()
    )

    # 5. Write .msh v2.2
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 0)
    gmsh.write(output_path)

    # Print mesh statistics
    elem_types, elem_tags, elem_node_tags = gmsh.model.mesh.getElements(dim=2)
    num_triangles = sum(len(tags) for tags in elem_tags)

    print(f"\n  Mesh written to: {output_path}")
    print(f"    Vertices:  {num_nodes}")
    print(f"    Triangles: {num_triangles}")
    print(f"    Format:    Gmsh .msh v2.2 (with DEM z in vertices)")

    gmsh.finalize()


def cmd_mesh(args):
    """
    CLI handler for the 'mesh' subcommand.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed CLI arguments.
    """
    from dem_operations import read_geotiff, write_esri_ascii

    # Read GeoTIFF
    print(f"\nReading GeoTIFF: {args.input}")
    elevation, transform, crs, nodata = read_geotiff(args.input)
    nrows, ncols = elevation.shape
    cellsize = abs(transform.a)
    print(f"  Loaded: {ncols} cols x {nrows} rows, cellsize={cellsize:.4f}")

    # Compute slope
    print("\nComputing slope...")
    slope = compute_slope(elevation, cellsize, nodata)
    slope_max = np.nanmax(slope)
    print(f"  Max slope: {slope_max:.4f}")

    # Create size field
    print("\nCreating adaptive size field...")
    size_field = create_size_field(slope, args.min_size, args.max_size, args.slope_factor)

    # Extract domain boundary
    print("\nExtracting domain boundary...")
    boundary_coords = build_domain_boundary(elevation, transform, nodata)

    # Generate mesh
    print(f"\nGenerating adaptive triangular mesh...")
    print(f"  Min element size: {args.min_size}")
    print(f"  Max element size: {args.max_size}")
    print(f"  Slope factor: {args.slope_factor}")

    generate_adaptive_mesh(
        boundary_coords, size_field, transform, elevation,
        nodata, args.min_size, args.max_size, args.output
    )

    # Generate companion .asc for DEM_FILE in modset.json
    # (main.cpp always reads DEM via get_domain_size/read_geo, even for triangular mesh)
    asc_path = os.path.splitext(args.output)[0] + '_dem.asc'
    print(f"\nGenerating companion ASC for modset.json DEM_FILE...")

    xllcorner = transform.c
    yllcorner = transform.f + transform.e * nrows

    write_esri_ascii(asc_path, elevation, xllcorner, yllcorner, cellsize, nodata)

    print(f"\n=== Summary ===")
    print(f"  Mesh file:   {args.output} (use as MESH_FILE in modset.json)")
    print(f"  DEM file:    {asc_path} (use as DEM_FILE in modset.json)")
    print(f"\nDone!")
