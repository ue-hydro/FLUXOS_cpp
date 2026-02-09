#!/usr/bin/env python3
"""
FLUXOS KML Exporter
====================

Export FLUXOS simulation results as lightweight KML + PNG files for animated
visualization in Google Earth.

Each timestep produces one small transparent PNG raster (GroundOverlay) instead
of thousands of polygon placemarks, keeping files compact (<1 MB total).

Supports both regular (Cartesian) and triangular mesh output formats.

Usage:
    # Regular mesh
    python fluxos_viewer.py --results-dir ./Results --dem ./terrain.asc --utm-zone 10

    # Triangular mesh
    python fluxos_viewer.py --results-dir ./Results_tri --dem ./terrain.asc \\
        --mesh-type triangular --utm-zone 10

    # Custom variable and colour range
    python fluxos_viewer.py --results-dir ./Results --dem ./terrain.asc \\
        --utm-zone 10 --variable velocity --clim 0 1.5

Requirements:
    pip install numpy pyproj
"""

import argparse
import datetime
import glob
import io
import os
import re
import struct
import sys
import xml.etree.ElementTree as ET
import zlib

import numpy as np

try:
    from pyproj import Proj
except ImportError:
    print("ERROR: pyproj is required.  Install with: pip install pyproj")
    sys.exit(1)


# ═══════════════════════════════════════════════════════════════════════════════
#  MINIMAL PNG WRITER  (no Pillow / PIL dependency)
# ═══════════════════════════════════════════════════════════════════════════════

def _write_png(rgba_array):
    """
    Encode an (H, W, 4) uint8 RGBA array as a PNG file in memory.

    Returns bytes of the PNG file.
    """
    h, w = rgba_array.shape[:2]

    def _chunk(chunk_type, data):
        c = chunk_type + data
        crc = struct.pack(">I", zlib.crc32(c) & 0xFFFFFFFF)
        return struct.pack(">I", len(data)) + c + crc

    # Signature
    sig = b"\x89PNG\r\n\x1a\n"

    # IHDR
    ihdr_data = struct.pack(">IIBBBBB", w, h, 8, 6, 0, 0, 0)  # 8-bit RGBA
    ihdr = _chunk(b"IHDR", ihdr_data)

    # IDAT — raw image data with filter byte 0 (None) per row
    raw_rows = []
    for y in range(h):
        raw_rows.append(b"\x00")  # filter: None
        raw_rows.append(rgba_array[y].tobytes())
    raw = b"".join(raw_rows)
    compressed = zlib.compress(raw, 9)
    idat = _chunk(b"IDAT", compressed)

    # IEND
    iend = _chunk(b"IEND", b"")

    return sig + ihdr + idat + iend


def _build_legend_png(clim, variable, opacity=180, width=300, height=40):
    """
    Create a colour-bar legend as a PNG image (returned as bytes).
    """
    var_labels = {
        "h": "Water Depth [m]",
        "velocity": "Velocity [m/s]",
        "conc_SW": "Concentration [mg/L]",
    }
    # Gradient bar matching the blue->cyan->yellow->red ramp (no text rendering)
    vals = np.linspace(clim[0] + 1e-9, clim[1], width).reshape(1, width)
    vals = np.broadcast_to(vals, (height, width)).copy()
    rgba = _value_to_rgba(vals, clim[0], clim[1], opacity)
    return _write_png(rgba)


# ═══════════════════════════════════════════════════════════════════════════════
#  DATA I/O
# ═══════════════════════════════════════════════════════════════════════════════

def read_dem_asc(filepath):
    """Read an ESRI ASCII Grid (.asc) DEM file."""
    meta = {}
    header_lines = 6
    with open(filepath, "r") as f:
        for _ in range(header_lines):
            line = f.readline().strip()
            parts = line.split()
            key = parts[0].lower()
            val = parts[1]
            if key == "ncols":
                meta["ncols"] = int(val)
            elif key == "nrows":
                meta["nrows"] = int(val)
            elif key == "xllcorner":
                meta["xllcorner"] = float(val)
            elif key == "yllcorner":
                meta["yllcorner"] = float(val)
            elif key == "cellsize":
                meta["cellsize"] = float(val)
            elif key == "nodata_value":
                meta["nodata"] = float(val)

    elevation = np.loadtxt(filepath, skiprows=header_lines)
    if elevation.shape != (meta["nrows"], meta["ncols"]):
        raise ValueError(
            f"DEM shape {elevation.shape} != header "
            f"({meta['nrows']}, {meta['ncols']})"
        )
    nodata = meta.get("nodata", -99999.0)
    elevation[elevation == nodata] = np.nan

    valid = elevation[~np.isnan(elevation)]
    print(f"  DEM loaded: {meta['ncols']}x{meta['nrows']} cells, "
          f"cellsize={meta['cellsize']}m")
    if len(valid) > 0:
        print(f"  Elevation range: {valid.min():.1f} - {valid.max():.1f} m")
    return elevation, meta


def read_regular_results(results_dir):
    """Read regular mesh .txt output files."""
    pattern = os.path.join(results_dir, "*.txt")
    files = glob.glob(pattern)
    timesteps = []
    for f in files:
        basename = os.path.splitext(os.path.basename(f))[0]
        try:
            t = float(basename)
            timesteps.append((t, f))
        except ValueError:
            continue
    timesteps.sort(key=lambda x: x[0])
    if not timesteps:
        print(f"  WARNING: No numeric .txt files in {results_dir}")
        return []
    print(f"  Found {len(timesteps)} regular-mesh timesteps")
    results = []
    for t, filepath in timesteps:
        try:
            data = np.genfromtxt(filepath, delimiter=",", skip_header=1)
            if data.ndim == 1:
                data = data.reshape(1, -1)
            results.append({"time": t, "data": data})
        except Exception as e:
            print(f"  WARNING: Could not read {filepath}: {e}")
    if results:
        print(f"  Loaded {len(results)} timesteps "
              f"({results[0]['time']:.0f}s - {results[-1]['time']:.0f}s)")
    return results


def read_triangular_results(results_dir):
    """Read triangular mesh .vtu output files (XML VTK format)."""
    pvd_files = glob.glob(os.path.join(results_dir, "*.pvd"))
    timesteps = []
    if pvd_files:
        pvd_path = pvd_files[0]
        print(f"  Reading time series from: {os.path.basename(pvd_path)}")
        tree = ET.parse(pvd_path)
        root = tree.getroot()
        for dataset in root.iter("DataSet"):
            t = float(dataset.get("timestep", 0))
            vtu_file = dataset.get("file", "")
            if not os.path.isabs(vtu_file):
                vtu_file = os.path.join(results_dir, vtu_file)
            timesteps.append((t, vtu_file))
    else:
        vtu_files = glob.glob(os.path.join(results_dir, "*.vtu"))
        for f in vtu_files:
            basename = os.path.splitext(os.path.basename(f))[0]
            try:
                t = float(basename)
                timesteps.append((t, f))
            except ValueError:
                continue
    timesteps.sort(key=lambda x: x[0])
    if not timesteps:
        print(f"  WARNING: No .vtu files in {results_dir}")
        return []
    print(f"  Found {len(timesteps)} triangular-mesh timesteps")
    results = []
    for t, filepath in timesteps:
        try:
            pts, cells, cdata = _parse_vtu(filepath)
            results.append({"time": t, "points": pts, "cells": cells,
                            "cell_data": cdata})
        except Exception as e:
            print(f"  WARNING: Could not read {filepath}: {e}")
    if results:
        print(f"  Loaded {len(results)} timesteps "
              f"({results[0]['time']:.0f}s - {results[-1]['time']:.0f}s)")
    return results


def _parse_vtu(filepath):
    """Minimal VTU (XML) parser -- no PyVista dependency."""
    tree = ET.parse(filepath)
    root = tree.getroot()
    piece = root.find(".//Piece")
    n_points = int(piece.get("NumberOfPoints"))
    n_cells = int(piece.get("NumberOfCells"))

    pts_da = piece.find(".//Points/DataArray")
    pts_flat = np.fromstring(pts_da.text.strip(), sep=" ")
    points = pts_flat.reshape(n_points, 3)

    cells_section = piece.find(".//Cells")
    conn_da = offsets_da = None
    for da in cells_section.findall("DataArray"):
        name = da.get("Name", "")
        if name == "connectivity":
            conn_da = da
        elif name == "offsets":
            offsets_da = da
    conn = np.fromstring(conn_da.text.strip(), dtype=int, sep=" ")
    offsets = np.fromstring(offsets_da.text.strip(), dtype=int, sep=" ")

    triangles = []
    prev = 0
    for off in offsets:
        cell_verts = conn[prev:off]
        if len(cell_verts) == 3:
            triangles.append(tuple(cell_verts))
        prev = off

    cell_data = {}
    cd_section = piece.find(".//CellData")
    if cd_section is not None:
        for da in cd_section.findall("DataArray"):
            name = da.get("Name", "unknown")
            n_comp = int(da.get("NumberOfComponents", "1"))
            vals = np.fromstring(da.text.strip(), sep=" ")
            if n_comp == 1:
                cell_data[name] = vals[:n_cells]
            else:
                cell_data[name] = vals.reshape(n_cells, n_comp)
    return points, triangles, cell_data


# ═══════════════════════════════════════════════════════════════════════════════
#  COORDINATE TRANSFORM
# ═══════════════════════════════════════════════════════════════════════════════

def make_utm_to_latlon(utm_zone, northern=True):
    """Return a function (easting, northing) -> (lon, lat)."""
    p = Proj(proj="utm", zone=utm_zone, ellps="WGS84", datum="WGS84",
             south=not northern)

    def transform(x, y):
        lon, lat = p(x, y, inverse=True)
        return lon, lat

    # Vectorised version for arrays
    def transform_arrays(x_arr, y_arr):
        return p(x_arr, y_arr, inverse=True)

    transform.arrays = transform_arrays
    return transform


def detect_utm_zone(xll, yll):
    """
    UTM zone cannot be uniquely determined from (easting, northing) alone.
    Defaults to zone 10 with a warning.
    """
    northern = yll < 10_000_000
    print("  WARNING: UTM zone cannot be auto-detected reliably from "
          "easting/northing alone.")
    print("           Defaulting to zone 10 (western Canada / US Pacific NW).")
    print("           Use --utm-zone to override if this is incorrect.")
    return 10, northern


# ═══════════════════════════════════════════════════════════════════════════════
#  RASTERIZE TO PNG
# ═══════════════════════════════════════════════════════════════════════════════

def _value_to_rgba(values, vmin, vmax, opacity=180):
    """
    Vectorised: map an array of scalar values to RGBA (H, W, 4) or (N, 4).
    Transparent where value <= 0 or NaN.
    """
    # Replace NaN with 0 before computing colours (masked out below)
    safe = np.where(np.isnan(values), 0.0, values)
    t = np.clip((safe - vmin) / max(vmax - vmin, 1e-12), 0.0, 1.0)
    r = (200 * (1 - t)).astype(np.uint8)
    g = (220 * (1 - t) + 50 * t).astype(np.uint8)
    b = (255 * (1 - t * 0.3)).astype(np.uint8)
    a = (opacity * (0.3 + 0.7 * t)).astype(np.uint8)

    # Transparent for dry / NaN cells
    mask = np.isnan(values) | (values <= 0)
    r[mask] = 0
    g[mask] = 0
    b[mask] = 0
    a[mask] = 0

    rgba = np.stack([r, g, b, a], axis=-1)
    return rgba


def rasterize_regular(dem, meta, data, variable, clim, h_min, opacity):
    """
    Rasterize one regular-mesh timestep into an RGBA image (nrows x ncols x 4).
    """
    nrows, ncols = dem.shape
    cs = meta["cellsize"]
    xll = meta["xllcorner"]
    yll = meta["yllcorner"]

    x_coords = data[:, 2]
    y_coords = data[:, 3]
    h_depth = data[:, 5]

    if variable == "h":
        var_values = h_depth.copy()
    elif variable == "velocity":
        ux, uy = data[:, 6], data[:, 7]
        var_values = np.sqrt(ux**2 + uy**2)
    elif variable == "conc_SW":
        var_values = data[:, 11].copy()
    else:
        var_values = h_depth.copy()

    # Map scattered points to grid
    col_idx = np.round((x_coords - xll) / cs).astype(int)
    row_idx = np.round((nrows - 1) - (y_coords - yll) / cs).astype(int)
    col_idx = np.clip(col_idx, 0, ncols - 1)
    row_idx = np.clip(row_idx, 0, nrows - 1)

    grid = np.full((nrows, ncols), np.nan)
    wet = h_depth > h_min
    grid[row_idx[wet], col_idx[wet]] = var_values[wet]

    rgba = _value_to_rgba(grid, clim[0], clim[1], opacity)
    wet_count = int(np.sum(wet))
    return rgba.astype(np.uint8), wet_count


def rasterize_triangular(points, cells, cell_data, variable, clim,
                          h_min, opacity, nrows, ncols, xll, yll, cs):
    """
    Rasterize one triangular-mesh timestep into an RGBA image (nrows x ncols x 4).

    Uses vertex-based interpolation for smooth rendering:
    1. Average cell values to vertices (area-weighted)
    2. Rasterize each triangle using barycentric interpolation of vertex values
    3. This produces smooth gradients instead of flat per-cell shading.
    """
    if "h" not in cell_data:
        return np.zeros((nrows, ncols, 4), dtype=np.uint8), 0

    h = cell_data["h"]
    if variable == "h":
        var_values = h
    elif variable == "velocity":
        if "velocity" in cell_data:
            vel = cell_data["velocity"]
            if vel.ndim == 2:
                var_values = np.sqrt(vel[:, 0]**2 + vel[:, 1]**2)
            else:
                var_values = np.abs(vel)
        else:
            var_values = h
    elif variable == "conc_SW":
        var_values = cell_data.get("conc_SW", h)
    else:
        var_values = h

    # ── Step 1: Interpolate cell values to vertices (area-weighted average) ──
    npts = len(points)
    vertex_val_sum = np.zeros(npts)
    vertex_weight = np.zeros(npts)
    for ci, (i0, i1, i2) in enumerate(cells):
        if h[ci] <= h_min:
            continue
        val = var_values[ci]
        # Use 1.0 weight (uniform); could use cell area for area-weighting
        vertex_val_sum[i0] += val
        vertex_val_sum[i1] += val
        vertex_val_sum[i2] += val
        vertex_weight[i0] += 1.0
        vertex_weight[i1] += 1.0
        vertex_weight[i2] += 1.0

    # Vertex values: average of surrounding wet cells
    with np.errstate(invalid='ignore', divide='ignore'):
        vertex_vals = np.where(vertex_weight > 0, vertex_val_sum / vertex_weight, 0.0)
    vertex_wet = vertex_weight > 0

    # ── Step 2: Rasterize with barycentric interpolation ──
    grid = np.full((nrows, ncols), np.nan)
    wet_count = 0

    for ci, (i0, i1, i2) in enumerate(cells):
        if h[ci] <= h_min:
            continue
        wet_count += 1

        # Vertex values for this triangle
        v0_val = vertex_vals[i0]
        v1_val = vertex_vals[i1]
        v2_val = vertex_vals[i2]

        # Triangle vertices in grid coordinates
        p0, p1, p2 = points[i0], points[i1], points[i2]
        c0g = (p0[0] - xll) / cs
        c1g = (p1[0] - xll) / cs
        c2g = (p2[0] - xll) / cs
        r0g = (nrows - 1) - (p0[1] - yll) / cs
        r1g = (nrows - 1) - (p1[1] - yll) / cs
        r2g = (nrows - 1) - (p2[1] - yll) / cs

        # Bounding box of triangle
        r_min = max(0, int(np.floor(min(r0g, r1g, r2g))))
        r_max = min(nrows - 1, int(np.ceil(max(r0g, r1g, r2g))))
        c_min = max(0, int(np.floor(min(c0g, c1g, c2g))))
        c_max = min(ncols - 1, int(np.ceil(max(c0g, c1g, c2g))))

        # Precompute barycentric denominator
        denom = (r1g - r2g) * (c0g - c2g) + (c2g - c1g) * (r0g - r2g)
        if abs(denom) < 1e-12:
            continue  # Degenerate triangle

        inv_denom = 1.0 / denom

        # Scan-fill with barycentric interpolation
        for r in range(r_min, r_max + 1):
            for c in range(c_min, c_max + 1):
                px, py = float(c), float(r)
                # Barycentric coordinates
                w0 = ((r1g - r2g) * (px - c2g) + (c2g - c1g) * (py - r2g)) * inv_denom
                w1 = ((r2g - r0g) * (px - c2g) + (c0g - c2g) * (py - r2g)) * inv_denom
                w2 = 1.0 - w0 - w1

                if w0 >= -1e-6 and w1 >= -1e-6 and w2 >= -1e-6:
                    # Interpolated value
                    grid[r, c] = w0 * v0_val + w1 * v1_val + w2 * v2_val

    rgba = _value_to_rgba(grid, clim[0], clim[1], opacity)
    return rgba.astype(np.uint8), wet_count


# ═══════════════════════════════════════════════════════════════════════════════
#  GRID RASTERIZER  (mesh wireframe for the first frame)
# ═══════════════════════════════════════════════════════════════════════════════

def _draw_line(img, r0, c0, r1, c1, rgba):
    """Bresenham's line into an (H, W, 4) image."""
    h, w = img.shape[:2]
    dr = abs(r1 - r0)
    dc = abs(c1 - c0)
    sr = 1 if r0 < r1 else -1
    sc = 1 if c0 < c1 else -1
    err = dc - dr
    while True:
        if 0 <= r0 < h and 0 <= c0 < w:
            img[r0, c0] = rgba
        if r0 == r1 and c0 == c1:
            break
        e2 = 2 * err
        if e2 > -dr:
            err -= dr
            c0 += sc
        if e2 < dc:
            err += dc
            r0 += sr


def rasterize_regular_grid(dem, meta, grid_color=(255, 255, 0, 160), spacing=1):
    """
    Draw the regular mesh grid lines onto a transparent RGBA image.
    Only draws lines where the DEM has valid (non-NaN) cells.
    """
    nrows, ncols = dem.shape
    img = np.zeros((nrows, ncols, 4), dtype=np.uint8)
    rgba = np.array(grid_color, dtype=np.uint8)

    # Mark domain outline and cell boundaries
    valid = ~np.isnan(dem)

    # Draw cell boundaries: highlight every Nth row and column
    step = max(1, spacing)
    # Horizontal lines
    for r in range(0, nrows, step):
        in_run = False
        for c in range(ncols):
            if valid[r, c]:
                if not in_run:
                    in_run = True
                img[r, c] = rgba
            else:
                in_run = False

    # Vertical lines
    for c in range(0, ncols, step):
        in_run = False
        for r in range(nrows):
            if valid[r, c]:
                if not in_run:
                    in_run = True
                img[r, c] = rgba
            else:
                in_run = False

    # Domain boundary: draw edges between valid and invalid
    boundary_rgba = np.array((255, 200, 0, 220), dtype=np.uint8)
    for r in range(nrows):
        for c in range(ncols):
            if not valid[r, c]:
                continue
            # Check 4-neighbours — if any is invalid or out-of-bounds, this is boundary
            is_edge = False
            for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                nr, nc = r + dr, c + dc
                if nr < 0 or nr >= nrows or nc < 0 or nc >= ncols:
                    is_edge = True
                elif np.isnan(dem[nr, nc]):
                    is_edge = True
            if is_edge:
                img[r, c] = boundary_rgba

    return img


def rasterize_triangular_grid(points, cells, nrows, ncols, xll, yll, cs,
                                grid_color=(255, 255, 0, 160)):
    """
    Draw the triangular mesh wireframe onto a transparent RGBA image.
    """
    img = np.zeros((nrows, ncols, 4), dtype=np.uint8)
    rgba = np.array(grid_color, dtype=np.uint8)

    for i0, i1, i2 in cells:
        p0, p1, p2 = points[i0], points[i1], points[i2]
        # Convert to pixel coordinates
        verts = []
        for p in [p0, p1, p2]:
            c = int(round((p[0] - xll) / cs))
            r = int(round((nrows - 1) - (p[1] - yll) / cs))
            verts.append((r, c))

        # Draw three edges
        _draw_line(img, verts[0][0], verts[0][1], verts[1][0], verts[1][1], rgba)
        _draw_line(img, verts[1][0], verts[1][1], verts[2][0], verts[2][1], rgba)
        _draw_line(img, verts[2][0], verts[2][1], verts[0][0], verts[0][1], rgba)

    return img


# ═══════════════════════════════════════════════════════════════════════════════
#  KML GENERATION  (one KML per timestep + master KML)
# ═══════════════════════════════════════════════════════════════════════════════

def format_time(seconds):
    """Format simulation time to readable string."""
    if seconds < 60:
        return f"{seconds:.0f}s"
    elif seconds < 3600:
        return f"{seconds/60:.1f} min"
    elif seconds < 86400:
        return f"{seconds/3600:.1f} hr"
    else:
        return f"{seconds/86400:.1f} days"


def _ground_overlay_kml(png_filename, name, north, south, east, west,
                         ts_begin, ts_end):
    """Create a KML GroundOverlay referencing a PNG image file."""
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
<Document>
  <name>{name}</name>
  <GroundOverlay>
    <name>{name}</name>
    <TimeSpan>
      <begin>{ts_begin}</begin>
      <end>{ts_end}</end>
    </TimeSpan>
    <Icon>
      <href>{png_filename}</href>
    </Icon>
    <LatLonBox>
      <north>{north:.8f}</north>
      <south>{south:.8f}</south>
      <east>{east:.8f}</east>
      <west>{west:.8f}</west>
    </LatLonBox>
    <color>ffffffff</color>
  </GroundOverlay>
</Document>
</kml>
"""


def _master_kml(timestep_entries, title="FLUXOS Simulation",
                grid_kml_file=None):
    """
    Create a master KML that uses NetworkLinks to load per-timestep KML files.
    timestep_entries: list of (kml_filename, name, ts_begin, ts_end)
    grid_kml_file: optional always-visible grid overlay
    """
    grid_link = ""
    if grid_kml_file:
        grid_link = f"""  <NetworkLink>
    <name>Mesh Grid</name>
    <visibility>1</visibility>
    <Link>
      <href>{grid_kml_file}</href>
    </Link>
  </NetworkLink>
"""

    links = ""
    for kml_file, name, ts_begin, ts_end in timestep_entries:
        links += f"""  <NetworkLink>
    <name>{name}</name>
    <TimeSpan>
      <begin>{ts_begin}</begin>
      <end>{ts_end}</end>
    </TimeSpan>
    <Link>
      <href>{kml_file}</href>
    </Link>
  </NetworkLink>
"""

    return f"""<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
<Document>
  <name>{title}</name>
  <description>FLUXOS flood simulation results.
Use the time slider in Google Earth to animate.</description>
  <open>1</open>
{grid_link}{links}</Document>
</kml>
"""


# ═══════════════════════════════════════════════════════════════════════════════
#  EXPORT
# ═══════════════════════════════════════════════════════════════════════════════

def export_kml(dem, meta, results, mesh_type, variable, clim, utm_to_ll,
               sim_start, output_dir, opacity=180, h_min=0.001):
    """
    Export results as individual KML + PNG files (one per timestep).
    """
    nrows, ncols = dem.shape
    cs = meta["cellsize"]
    xll = meta["xllcorner"]
    yll = meta["yllcorner"]

    # Compute bounding box in lat/lon
    sw_lon, sw_lat = utm_to_ll(xll, yll)
    ne_lon, ne_lat = utm_to_ll(xll + ncols * cs, yll + nrows * cs)
    # Also check NW and SE corners for better bounds
    nw_lon, nw_lat = utm_to_ll(xll, yll + nrows * cs)
    se_lon, se_lat = utm_to_ll(xll + ncols * cs, yll)

    north = max(sw_lat, ne_lat, nw_lat, se_lat)
    south = min(sw_lat, ne_lat, nw_lat, se_lat)
    east = max(sw_lon, ne_lon, nw_lon, se_lon)
    west = min(sw_lon, ne_lon, nw_lon, se_lon)

    print(f"  Bounding box: lat [{south:.6f}, {north:.6f}], "
          f"lon [{west:.6f}, {east:.6f}]")

    os.makedirs(output_dir, exist_ok=True)

    timestep_entries = []
    total_size = 0

    # ── Frame 0: mesh grid ───────────────────────────────────
    if mesh_type == "regular":
        # Grid lines every 10 cells for visibility
        grid_step = max(1, min(nrows, ncols) // 60)
        grid_img = rasterize_regular_grid(dem, meta, spacing=grid_step)
    else:
        grid_img = rasterize_triangular_grid(
            results[0]["points"], results[0]["cells"],
            nrows, ncols, xll, yll, cs)

    grid_png_name = "grid.png"
    grid_png_path = os.path.join(output_dir, grid_png_name)
    grid_png_bytes = _write_png(grid_img)
    with open(grid_png_path, "wb") as f:
        f.write(grid_png_bytes)

    # Grid frame is always visible (no TimeSpan)
    grid_kml_name = "grid.kml"
    grid_kml_path = os.path.join(output_dir, grid_kml_name)
    grid_kml_content = f"""<?xml version="1.0" encoding="UTF-8"?>
<kml xmlns="http://www.opengis.net/kml/2.2">
<Document>
  <name>Mesh Grid</name>
  <GroundOverlay>
    <name>Mesh Grid ({mesh_type})</name>
    <Icon>
      <href>{grid_png_name}</href>
    </Icon>
    <LatLonBox>
      <north>{north:.8f}</north>
      <south>{south:.8f}</south>
      <east>{east:.8f}</east>
      <west>{west:.8f}</west>
    </LatLonBox>
    <color>ffffffff</color>
  </GroundOverlay>
</Document>
</kml>
"""
    with open(grid_kml_path, "w") as f:
        f.write(grid_kml_content)

    total_size += len(grid_png_bytes)
    print(f"  Grid frame: {mesh_type} mesh, PNG {len(grid_png_bytes)/1024:.0f} KB")

    for ri, r in enumerate(results):
        t_sec = r["time"]

        # Time span
        t_start = sim_start + datetime.timedelta(seconds=t_sec)
        if ri + 1 < len(results):
            t_end = sim_start + datetime.timedelta(
                seconds=results[ri + 1]["time"])
        else:
            t_end = t_start + datetime.timedelta(seconds=t_sec * 0.5)

        ts_begin = t_start.strftime("%Y-%m-%dT%H:%M:%SZ")
        ts_end = t_end.strftime("%Y-%m-%dT%H:%M:%SZ")

        # Rasterize
        if mesh_type == "regular":
            rgba, wet = rasterize_regular(
                dem, meta, r["data"], variable, clim, h_min, opacity)
        else:
            rgba, wet = rasterize_triangular(
                r["points"], r["cells"], r["cell_data"],
                variable, clim, h_min, opacity, nrows, ncols, xll, yll, cs)

        if wet == 0:
            print(f"  Frame {ri+1}/{len(results)}: t={format_time(t_sec)} "
                  f"- no wet cells, skipped")
            continue

        # Write PNG
        png_name = f"frame_{ri:04d}.png"
        png_path = os.path.join(output_dir, png_name)
        png_bytes = _write_png(rgba)
        with open(png_path, "wb") as f:
            f.write(png_bytes)

        # Write per-timestep KML
        kml_name = f"frame_{ri:04d}.kml"
        kml_path = os.path.join(output_dir, kml_name)
        frame_label = f"t = {format_time(t_sec)}"
        kml_content = _ground_overlay_kml(
            png_name, frame_label, north, south, east, west,
            ts_begin, ts_end)
        with open(kml_path, "w") as f:
            f.write(kml_content)

        timestep_entries.append((kml_name, frame_label, ts_begin, ts_end))
        png_size = len(png_bytes)
        total_size += png_size

        print(f"  Frame {ri+1}/{len(results)}: t={format_time(t_sec)}, "
              f"{wet} wet cells, PNG {png_size/1024:.0f} KB")

    # Write legend PNG
    legend_bytes = _build_legend_png(clim, variable, opacity)
    legend_path = os.path.join(output_dir, "legend.png")
    with open(legend_path, "wb") as f:
        f.write(legend_bytes)

    # Write master KML (with grid as always-visible layer)
    master_path = os.path.join(output_dir, "fluxos_animation.kml")
    master_kml = _master_kml(timestep_entries, f"FLUXOS {mesh_type.title()} Mesh",
                              grid_kml_file=grid_kml_name)
    with open(master_path, "w") as f:
        f.write(master_kml)

    n_frames = len(timestep_entries)
    print(f"\n  Exported {n_frames} frames to: {output_dir}/")
    print(f"  Total PNG size: {total_size/1024:.0f} KB")
    print(f"  Master KML: {master_path}")
    print(f"\n  Open in Google Earth:")
    print(f"    open {master_path}")
    print(f"  Then use the time slider to animate through timesteps.")


# ═══════════════════════════════════════════════════════════════════════════════
#  MAIN / CLI
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="FLUXOS KML Exporter -- export simulation results as "
                    "KML + PNG for Google Earth with time animation.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
  # Regular mesh
  python fluxos_viewer.py --results-dir ./Results --dem ./terrain.asc --utm-zone 10

  # Triangular mesh
  python fluxos_viewer.py --results-dir ./Results_tri --dem ./terrain.asc \\
      --mesh-type triangular --utm-zone 10

  # Custom variable and colour range
  python fluxos_viewer.py --results-dir ./Results --dem ./terrain.asc \\
      --utm-zone 10 --variable velocity --clim 0 1.5

Output:
  A directory of .kml + .png files.  Open fluxos_animation.kml in Google Earth.
  Use the time slider to animate through timesteps.
""",
    )
    parser.add_argument(
        "--results-dir", required=True,
        help="Directory containing simulation output files (.txt or .vtu)",
    )
    parser.add_argument(
        "--dem", required=True,
        help="Path to ESRI ASCII Grid DEM file (.asc)",
    )
    parser.add_argument(
        "--mesh-type", choices=["regular", "triangular"], default="regular",
        help="Mesh type: 'regular' (.txt) or 'triangular' (.vtu). Default: regular",
    )
    parser.add_argument(
        "--variable", choices=["h", "velocity", "conc_SW"], default="h",
        help="Variable to visualize. Default: h (water depth)",
    )
    parser.add_argument(
        "--clim", type=float, nargs=2, metavar=("MIN", "MAX"), default=None,
        help="Colour limits for the variable (e.g. --clim 0 2.0)",
    )
    parser.add_argument(
        "--utm-zone", type=int, default=None,
        help="UTM zone number (1-60). Auto-detected from DEM if omitted.",
    )
    parser.add_argument(
        "--southern-hemisphere", action="store_true",
        help="Set if the site is in the southern hemisphere.",
    )
    parser.add_argument(
        "--sim-start", type=str, default=None,
        help="Simulation start date-time (ISO format, e.g. '2009-01-01T00:00:00'). "
             "Defaults to value in modset.json if found, else 2000-01-01.",
    )
    parser.add_argument(
        "--output-dir", "-o", type=str, default=None,
        help="Output directory for KML + PNG files. "
             "Default: fluxos_kml_<mesh_type>",
    )
    parser.add_argument(
        "--opacity", type=int, default=140,
        help="Water overlay opacity (0-255). Default: 140",
    )
    parser.add_argument(
        "--h-min", type=float, default=0.001,
        help="Minimum water depth to include (m). Default: 0.001",
    )

    args = parser.parse_args()

    # ── Validate ─────────────────────────────────────────────
    if not os.path.isfile(args.dem):
        print(f"ERROR: DEM file not found: {args.dem}")
        sys.exit(1)
    if not os.path.isdir(args.results_dir):
        print(f"ERROR: Results directory not found: {args.results_dir}")
        sys.exit(1)

    if args.output_dir is None:
        args.output_dir = f"fluxos_kml_{args.mesh_type}"

    # ── Step 1: Read DEM ─────────────────────────────────────
    print("\n[1/4] Reading DEM...")
    dem, meta = read_dem_asc(args.dem)

    # ── Step 2: Coordinate transform ─────────────────────────
    print("[2/4] Setting up coordinate transform...")
    if args.utm_zone is not None:
        zone = args.utm_zone
        northern = not args.southern_hemisphere
    else:
        zone, northern = detect_utm_zone(meta["xllcorner"], meta["yllcorner"])

    utm_to_ll = make_utm_to_latlon(zone, northern)
    test_lon, test_lat = utm_to_ll(meta["xllcorner"], meta["yllcorner"])
    print(f"  UTM Zone {zone}{'N' if northern else 'S'}")
    print(f"  DEM SW corner -> lat={test_lat:.6f}, lon={test_lon:.6f}")

    # ── Simulation start time ────────────────────────────────
    if args.sim_start:
        sim_start = datetime.datetime.fromisoformat(args.sim_start)
    else:
        sim_start = _read_sim_start()
    print(f"  Simulation start: {sim_start.isoformat()}")

    # ── Step 3: Read results ─────────────────────────────────
    print(f"[3/4] Reading {args.mesh_type} mesh results...")
    if args.mesh_type == "regular":
        results = read_regular_results(args.results_dir)
    else:
        results = read_triangular_results(args.results_dir)

    if not results:
        print("ERROR: No results found.")
        sys.exit(1)

    # ── Compute colour limits ────────────────────────────────
    if args.clim is not None:
        clim = tuple(args.clim)
    else:
        clim = _auto_clim(results, args.mesh_type, args.variable, args.h_min)
    print(f"  Colour range: {clim[0]:.4f} - {clim[1]:.4f}")

    # ── Step 4: Export ───────────────────────────────────────
    print(f"[4/4] Exporting KML + PNG ({args.mesh_type} mesh)...")
    export_kml(
        dem, meta, results, args.mesh_type, args.variable, clim,
        utm_to_ll, sim_start, args.output_dir, args.opacity, args.h_min,
    )


def _read_sim_start():
    """Try to read SIM_DATETIME_START from modset.json."""
    import json
    for path in ["modset.json", "bin/modset.json"]:
        if os.path.isfile(path):
            try:
                with open(path) as f:
                    cfg = json.load(f)
                dt_str = cfg.get("SIM_DATETIME_START", "")
                if dt_str:
                    return datetime.datetime.fromisoformat(dt_str)
            except Exception:
                pass
    return datetime.datetime(2000, 1, 1)


def _auto_clim(results, mesh_type, variable, h_min):
    """Auto-detect colour limits from the data."""
    all_vals = []
    for r in results:
        if mesh_type == "regular":
            data = r["data"]
            h = data[:, 5]
            wet = h > h_min
            if not np.any(wet):
                continue
            if variable == "h":
                vals = h[wet]
            elif variable == "velocity":
                vals = np.sqrt(data[:, 6][wet]**2 + data[:, 7][wet]**2)
            elif variable == "conc_SW":
                vals = data[:, 11][wet]
            else:
                vals = h[wet]
        else:
            cdata = r["cell_data"]
            if "h" not in cdata:
                continue
            h = cdata["h"]
            wet = h > h_min
            if not np.any(wet):
                continue
            if variable == "h":
                vals = h[wet]
            elif variable == "velocity":
                if "velocity" in cdata:
                    vel = cdata["velocity"]
                    if vel.ndim == 2:
                        vals = np.sqrt(vel[wet, 0]**2 + vel[wet, 1]**2)
                    else:
                        vals = np.abs(vel[wet])
                else:
                    vals = h[wet]
            elif variable == "conc_SW":
                vals = cdata.get("conc_SW", h)[wet]
            else:
                vals = h[wet]
        valid = vals[~np.isnan(vals)]
        if len(valid) > 0:
            all_vals.append(valid)
    if all_vals:
        all_vals = np.concatenate(all_vals)
        return (0.0, float(np.percentile(all_vals, 98)))
    return (0.0, 1.0)


if __name__ == "__main__":
    main()
