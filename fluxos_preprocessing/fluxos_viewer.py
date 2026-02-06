#!/usr/bin/env python3
"""
FLUXOS 3D Flood Viewer
======================

Interactive 3D visualization of FLUXOS simulation results on top of terrain.
Supports both regular (Cartesian) and triangular mesh output formats.

Usage:
    # Regular mesh
    python fluxos_viewer.py --results-dir ./Results --dem ./terrain.asc --mesh-type regular

    # Triangular mesh
    python fluxos_viewer.py --results-dir ./fluxos_out --dem ./terrain.asc --mesh-type triangular

    # With options
    python fluxos_viewer.py --results-dir ./Results --dem ./terrain.asc \\
        --variable h --exaggeration 5.0 --water-opacity 0.7 --clim 0.0 2.0

Requirements:
    pip install pyvista numpy
"""

import argparse
import glob
import os
import re
import sys
import xml.etree.ElementTree as ET

import numpy as np

try:
    import pyvista as pv
except ImportError:
    print("ERROR: pyvista is required. Install with: pip install pyvista")
    sys.exit(1)


# ═══════════════════════════════════════════════════════════════════════════════
#  DATA I/O
# ═══════════════════════════════════════════════════════════════════════════════

def read_dem_asc(filepath):
    """
    Read an ESRI ASCII Grid (.asc) DEM file.

    Returns:
        elevation (np.ndarray): 2D array of elevations (nrows x ncols), NaN for NODATA.
        meta (dict): Header metadata {ncols, nrows, xllcorner, yllcorner, cellsize, nodata}.
    """
    meta = {}
    header_lines = 6
    header_keys = ["ncols", "nrows", "xllcorner", "yllcorner", "cellsize", "nodata_value"]

    with open(filepath, "r") as f:
        for i in range(header_lines):
            line = f.readline().strip()
            parts = line.split()
            key = parts[0].lower()
            val = parts[1]
            # Map to standard names
            if key in ("ncols",):
                meta["ncols"] = int(val)
            elif key in ("nrows",):
                meta["nrows"] = int(val)
            elif key in ("xllcorner",):
                meta["xllcorner"] = float(val)
            elif key in ("yllcorner",):
                meta["yllcorner"] = float(val)
            elif key in ("cellsize",):
                meta["cellsize"] = float(val)
            elif key in ("nodata_value",):
                meta["nodata"] = float(val)

    # Read elevation data (north-to-south, row 0 = northernmost)
    elevation = np.loadtxt(filepath, skiprows=header_lines)

    if elevation.shape != (meta["nrows"], meta["ncols"]):
        raise ValueError(
            f"DEM data shape {elevation.shape} doesn't match header "
            f"({meta['nrows']}, {meta['ncols']})"
        )

    # Replace NODATA with NaN
    nodata = meta.get("nodata", -99999.0)
    elevation[elevation == nodata] = np.nan

    print(f"  DEM loaded: {meta['ncols']}x{meta['nrows']} cells, "
          f"cellsize={meta['cellsize']}m")
    valid = elevation[~np.isnan(elevation)]
    if len(valid) > 0:
        print(f"  Elevation range: {valid.min():.1f} - {valid.max():.1f} m")

    return elevation, meta


def read_regular_results(results_dir):
    """
    Read regular mesh .txt output files.

    Returns:
        list of dict: [{"time": float, "data": np.ndarray}, ...] sorted by time.
            data columns: Xcoord, Ycoord, z, h, ux, uy, conc_SW
    """
    # Find all .txt files with numeric names (timestep in seconds)
    pattern = os.path.join(results_dir, "*.txt")
    files = glob.glob(pattern)

    timesteps = []
    for f in files:
        basename = os.path.splitext(os.path.basename(f))[0]
        # Only include files with numeric names
        try:
            t = float(basename)
            timesteps.append((t, f))
        except ValueError:
            continue

    timesteps.sort(key=lambda x: x[0])

    if not timesteps:
        print(f"  WARNING: No numeric .txt files found in {results_dir}")
        return []

    print(f"  Found {len(timesteps)} regular mesh timesteps")

    results = []
    for t, filepath in timesteps:
        try:
            # Read CSV with header
            data = np.genfromtxt(filepath, delimiter=",", skip_header=1)
            if data.ndim == 1:
                data = data.reshape(1, -1)
            # Extract relevant columns:
            # 0=irow, 1=icol, 2=Xcoord, 3=Ycoord, 4=z, 5=h, 6=ux, 7=uy,
            # 8=qx*dxy, 9=qy*dxy, 10=us, 11=conc_SW, 12=fe_1, 13=fn_1, 14=twetimetracer
            results.append({"time": t, "data": data})
        except Exception as e:
            print(f"  WARNING: Could not read {filepath}: {e}")
            continue

    print(f"  Loaded {len(results)} timesteps "
          f"({results[0]['time']:.0f}s - {results[-1]['time']:.0f}s)")
    return results


def read_triangular_results(results_dir):
    """
    Read triangular mesh .vtu output files (with optional .pvd time series).

    Returns:
        list of dict: [{"time": float, "mesh": pyvista.UnstructuredGrid}, ...] sorted by time.
    """
    # Try to find .pvd file first
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
        # Fallback: glob for .vtu files with numeric names
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
        print(f"  WARNING: No .vtu files found in {results_dir}")
        return []

    print(f"  Found {len(timesteps)} triangular mesh timesteps")

    results = []
    for t, filepath in timesteps:
        try:
            mesh = pv.read(filepath)
            results.append({"time": t, "mesh": mesh})
        except Exception as e:
            print(f"  WARNING: Could not read {filepath}: {e}")
            continue

    print(f"  Loaded {len(results)} timesteps "
          f"({results[0]['time']:.0f}s - {results[-1]['time']:.0f}s)")
    return results


# ═══════════════════════════════════════════════════════════════════════════════
#  MESH CONSTRUCTION
# ═══════════════════════════════════════════════════════════════════════════════

def build_terrain_mesh(dem, meta, exaggeration=5.0):
    """
    Build a PyVista StructuredGrid for the terrain surface.

    Args:
        dem (np.ndarray): 2D elevation array (nrows x ncols), NaN for NODATA.
        meta (dict): DEM metadata.
        exaggeration (float): Vertical exaggeration factor.

    Returns:
        pyvista.StructuredGrid: Terrain surface mesh.
    """
    nrows, ncols = dem.shape
    cellsize = meta["cellsize"]
    xll = meta["xllcorner"]
    yll = meta["yllcorner"]

    # Build coordinate arrays
    # X increases with column index (west to east)
    x = xll + np.arange(ncols) * cellsize
    # Y: ASC row 0 is northernmost, so flip for Y coordinates
    y = yll + (nrows - 1 - np.arange(nrows)) * cellsize

    # Create 2D meshgrids
    X, Y = np.meshgrid(x, y)

    # Handle NaN for Z: replace with minimum valid elevation for rendering
    Z = dem.copy()
    valid_mask = ~np.isnan(Z)
    if np.any(valid_mask):
        z_min = np.nanmin(Z)
        Z[~valid_mask] = z_min
    else:
        Z[~valid_mask] = 0.0

    Z_exag = Z * exaggeration

    # Build StructuredGrid (needs 3D arrays)
    grid = pv.StructuredGrid(X, Y, Z_exag)

    # Add elevation as point data for coloring
    grid.point_data["elevation"] = dem.flatten(order="C")

    # Store NODATA mask
    nodata_mask = np.isnan(dem).flatten(order="C").astype(float)
    grid.point_data["nodata_mask"] = nodata_mask

    return grid


def build_water_regular(dem, meta, timestep_data, variable="h", exaggeration=5.0):
    """
    Build a water surface mesh for one regular-mesh timestep.

    Args:
        dem (np.ndarray): 2D elevation array.
        meta (dict): DEM metadata.
        timestep_data (np.ndarray): Output data array (rows of CSV columns).
        variable (str): Variable to visualize ("h", "velocity", "conc_SW").
        exaggeration (float): Vertical exaggeration.

    Returns:
        pyvista.PolyData or None: Water surface mesh, or None if no wet cells.
    """
    if timestep_data is None or len(timestep_data) == 0:
        return None

    nrows, ncols = dem.shape
    cellsize = meta["cellsize"]
    xll = meta["xllcorner"]
    yll = meta["yllcorner"]

    # Extract columns from output data
    # 0=irow, 1=icol, 2=Xcoord, 3=Ycoord, 4=z, 5=h, 6=ux, 7=uy,
    # 8=qx*dxy, 9=qy*dxy, 10=us, 11=conc_SW, ...
    x_coords = timestep_data[:, 2]
    y_coords = timestep_data[:, 3]
    z_surface = timestep_data[:, 4]  # water surface elevation
    h_depth = timestep_data[:, 5]    # water depth

    # Select the visualization variable
    if variable == "h":
        var_values = h_depth
        var_name = "Water Depth [m]"
    elif variable == "velocity":
        ux = timestep_data[:, 6]
        uy = timestep_data[:, 7]
        var_values = np.sqrt(ux**2 + uy**2)
        var_name = "Velocity [m/s]"
    elif variable == "conc_SW":
        var_values = timestep_data[:, 11]
        var_name = "Concentration [mg/L]"
    else:
        var_values = h_depth
        var_name = "Water Depth [m]"

    # Filter for wet cells only (h > 0)
    wet_mask = h_depth > 1e-6
    if not np.any(wet_mask):
        return None

    x_wet = x_coords[wet_mask]
    y_wet = y_coords[wet_mask]
    z_wet = z_surface[wet_mask] * exaggeration
    var_wet = var_values[wet_mask]

    # Build a surface from scattered wet cell points
    # Map to grid indices and create a structured water surface
    col_idx = np.round((x_wet - xll) / cellsize).astype(int)
    row_idx = np.round((nrows - 1) - (y_wet - yll) / cellsize).astype(int)

    # Clamp indices
    col_idx = np.clip(col_idx, 0, ncols - 1)
    row_idx = np.clip(row_idx, 0, nrows - 1)

    # Create full grids initialized with NaN
    water_z = np.full((nrows, ncols), np.nan)
    water_var = np.full((nrows, ncols), np.nan)

    for i in range(len(x_wet)):
        r, c = row_idx[i], col_idx[i]
        water_z[r, c] = z_wet[i]
        water_var[r, c] = var_wet[i]

    # Build coordinate arrays (same as terrain)
    x = xll + np.arange(ncols) * cellsize
    y = yll + (nrows - 1 - np.arange(nrows)) * cellsize
    X, Y = np.meshgrid(x, y)

    # Create points and faces for wet cells only
    # Collect quads for cells that are wet and have wet neighbors
    points = []
    faces = []
    scalars = []
    pt_idx = 0

    for r in range(nrows - 1):
        for c in range(ncols - 1):
            # Check if this cell and its right/bottom neighbors have water
            corners_z = [water_z[r, c], water_z[r, c+1],
                         water_z[r+1, c+1], water_z[r+1, c]]
            if any(np.isnan(z) for z in corners_z):
                continue

            # All 4 corners are wet: create a quad
            p0 = [X[r, c],     Y[r, c],     corners_z[0]]
            p1 = [X[r, c+1],   Y[r, c+1],   corners_z[1]]
            p2 = [X[r+1, c+1], Y[r+1, c+1], corners_z[2]]
            p3 = [X[r+1, c],   Y[r+1, c],   corners_z[3]]

            points.extend([p0, p1, p2, p3])
            faces.append([4, pt_idx, pt_idx+1, pt_idx+2, pt_idx+3])
            # Cell scalar = average of 4 corner values
            cell_val = np.nanmean([water_var[r, c], water_var[r, c+1],
                                   water_var[r+1, c+1], water_var[r+1, c]])
            scalars.append(cell_val)
            pt_idx += 4

    if not points:
        return None

    points = np.array(points)
    faces = np.array([item for sublist in faces for item in sublist])
    scalars = np.array(scalars)

    water_mesh = pv.PolyData(points, faces)
    water_mesh.cell_data[var_name] = scalars
    water_mesh._var_name = var_name  # Store for later reference

    return water_mesh


def build_water_triangular(timestep_mesh, variable="h", exaggeration=5.0):
    """
    Build a water surface mesh for one triangular-mesh timestep.

    Args:
        timestep_mesh (pyvista.UnstructuredGrid): VTU mesh with cell data.
        variable (str): Variable to visualize ("h", "velocity", "conc_SW").
        exaggeration (float): Vertical exaggeration.

    Returns:
        pyvista.UnstructuredGrid or None: Water surface with selected variable.
    """
    mesh = timestep_mesh.copy()

    # Get available arrays
    cell_arrays = list(mesh.cell_data.keys())

    # Get water depth
    if "h" not in cell_arrays:
        print("  WARNING: 'h' array not found in VTU cell data")
        return None

    h = mesh.cell_data["h"]

    # Get bed elevation
    if "zb" in cell_arrays:
        zb = mesh.cell_data["zb"]
    else:
        zb = np.zeros(mesh.n_cells)

    # Select visualization variable
    if variable == "h":
        var_values = h.copy()
        var_name = "Water Depth [m]"
    elif variable == "velocity":
        if "velocity" in cell_arrays:
            vel = mesh.cell_data["velocity"]
            if vel.ndim == 2 and vel.shape[1] >= 2:
                var_values = np.sqrt(vel[:, 0]**2 + vel[:, 1]**2)
            else:
                var_values = np.abs(vel)
        else:
            var_values = h.copy()
            var_name = "Water Depth [m]"
            print("  WARNING: 'velocity' not found, falling back to 'h'")
        var_name = "Velocity [m/s]"
    elif variable == "conc_SW":
        if "conc_SW" in cell_arrays:
            var_values = mesh.cell_data["conc_SW"]
        else:
            var_values = h.copy()
            print("  WARNING: 'conc_SW' not found, falling back to 'h'")
        var_name = "Concentration [mg/L]"
    else:
        var_values = h.copy()
        var_name = "Water Depth [m]"

    # Convert cell data to point data for smooth elevation interpolation
    mesh_with_points = mesh.cell_data_to_point_data()

    # Set Z coordinates: water surface = zb + h, exaggerated
    if "zb" in mesh_with_points.point_data and "h" in mesh_with_points.point_data:
        zb_pts = mesh_with_points.point_data["zb"]
        h_pts = mesh_with_points.point_data["h"]
        new_z = (zb_pts + h_pts) * exaggeration
    elif "z" in mesh_with_points.point_data:
        new_z = mesh_with_points.point_data["z"] * exaggeration
    else:
        # Fallback: use existing points z
        new_z = mesh.points[:, 2] * exaggeration

    # Update Z coordinates
    pts = mesh.points.copy()
    pts[:, 2] = new_z if len(new_z) == len(pts) else pts[:, 2] * exaggeration
    mesh.points = pts

    # Add the selected variable
    mesh.cell_data[var_name] = var_values

    # Threshold to keep only wet cells
    eps = 1e-4
    wet_mask = h > eps
    if not np.any(wet_mask):
        return None

    # Extract wet cells
    cell_ids = np.where(wet_mask)[0]
    if len(cell_ids) == 0:
        return None

    water_mesh = mesh.extract_cells(cell_ids)
    water_mesh._var_name = var_name

    return water_mesh


# ═══════════════════════════════════════════════════════════════════════════════
#  TERRAIN MESH FOR TRIANGULAR (from DEM)
# ═══════════════════════════════════════════════════════════════════════════════

def build_terrain_for_triangular(dem, meta, exaggeration=5.0):
    """
    Build terrain surface for use with triangular mesh results.
    Same as build_terrain_mesh but returns a surface suitable for overlay.
    """
    return build_terrain_mesh(dem, meta, exaggeration)


# ═══════════════════════════════════════════════════════════════════════════════
#  VIEWER
# ═══════════════════════════════════════════════════════════════════════════════

def format_time(seconds):
    """Format simulation time in seconds to a readable string."""
    if seconds < 60:
        return f"{seconds:.0f}s"
    elif seconds < 3600:
        return f"{seconds/60:.1f} min"
    elif seconds < 86400:
        return f"{seconds/3600:.1f} hr"
    else:
        return f"{seconds/86400:.1f} days"


def launch_viewer(terrain, water_frames, times, variable, water_opacity=0.7,
                  clim=None):
    """
    Create an interactive PyVista 3D viewer with terrain and animated water.

    Args:
        terrain (pyvista.StructuredGrid): Terrain surface mesh.
        water_frames (list): List of water surface meshes (one per timestep).
        times (list[float]): Timestep values in seconds.
        variable (str): Name of the variable being visualized.
        water_opacity (float): Water surface opacity (0-1).
        clim (tuple or None): Color limits (min, max) for the water variable.
    """
    if not water_frames:
        print("No water frames to display. Showing terrain only.")
        pl = pv.Plotter()
        pl.add_mesh(terrain, scalars="elevation", cmap="gist_earth",
                    nan_opacity=0.0, show_scalar_bar=True,
                    scalar_bar_args={"title": "Elevation [m]"})
        pl.add_text("FLUXOS Terrain (no simulation data)", font_size=12)
        pl.show()
        return

    # Determine variable name from first water frame
    var_name = getattr(water_frames[0], "_var_name", "Water Depth [m]")

    # Compute global color limits if not specified
    if clim is None:
        all_vals = []
        for wf in water_frames:
            if wf is not None and var_name in wf.cell_data:
                vals = wf.cell_data[var_name]
                valid = vals[~np.isnan(vals)]
                if len(valid) > 0:
                    all_vals.append(valid)
        if all_vals:
            all_vals = np.concatenate(all_vals)
            clim = (0.0, np.percentile(all_vals, 98))
        else:
            clim = (0.0, 1.0)

    # ── Create plotter ──────────────────────────────────────────
    pl = pv.Plotter()
    pl.set_background("white")

    # Add terrain (static)
    pl.add_mesh(
        terrain,
        scalars="elevation",
        cmap="gist_earth",
        nan_opacity=0.0,
        show_scalar_bar=True,
        scalar_bar_args={
            "title": "Elevation [m]",
            "position_x": 0.05,
            "position_y": 0.05,
            "width": 0.3,
            "height": 0.05,
        },
    )

    # Find first non-None water frame
    first_idx = 0
    for i, wf in enumerate(water_frames):
        if wf is not None:
            first_idx = i
            break

    # Add first water frame
    water_actor = None
    if water_frames[first_idx] is not None:
        water_actor = pl.add_mesh(
            water_frames[first_idx],
            scalars=var_name,
            cmap="Blues",
            clim=clim,
            opacity=water_opacity,
            show_scalar_bar=True,
            scalar_bar_args={
                "title": var_name,
                "position_x": 0.05,
                "position_y": 0.15,
                "width": 0.3,
                "height": 0.05,
            },
            name="water_surface",
        )

    # Time label
    time_label = pl.add_text(
        f"Time: {format_time(times[first_idx])}",
        position="upper_right",
        font_size=14,
        color="black",
        name="time_label",
    )

    # Title
    pl.add_text(
        "FLUXOS 3D Flood Viewer",
        position="upper_left",
        font_size=10,
        color="gray",
    )

    # ── State for animation ────────────────────────────────────
    state = {"current_idx": first_idx, "playing": False}

    def update_timestep(idx):
        """Update the displayed water surface to timestep idx."""
        idx = int(round(idx))
        idx = max(0, min(idx, len(water_frames) - 1))
        state["current_idx"] = idx

        # Update time label
        pl.add_text(
            f"Time: {format_time(times[idx])}",
            position="upper_right",
            font_size=14,
            color="black",
            name="time_label",
        )

        # Update water surface
        wf = water_frames[idx]
        if wf is not None and var_name in wf.cell_data:
            pl.add_mesh(
                wf,
                scalars=var_name,
                cmap="Blues",
                clim=clim,
                opacity=water_opacity,
                show_scalar_bar=False,
                name="water_surface",
            )
        else:
            # Remove water surface if this timestep has no data
            pl.remove_actor("water_surface")

    # ── Slider widget ──────────────────────────────────────────
    if len(water_frames) > 1:
        pl.add_slider_widget(
            update_timestep,
            rng=[0, len(water_frames) - 1],
            value=first_idx,
            title="Timestep",
            pointa=(0.25, 0.92),
            pointb=(0.75, 0.92),
            style="modern",
            fmt="%.0f",
        )

    # ── Key bindings for animation ─────────────────────────────
    def next_frame():
        idx = min(state["current_idx"] + 1, len(water_frames) - 1)
        update_timestep(idx)
        pl.render()

    def prev_frame():
        idx = max(state["current_idx"] - 1, 0)
        update_timestep(idx)
        pl.render()

    def toggle_play():
        state["playing"] = not state["playing"]
        if state["playing"]:
            print("  Animation: PLAY (press 'p' to pause)")
        else:
            print("  Animation: PAUSE")

    pl.add_key_event("Right", next_frame)
    pl.add_key_event("Left", prev_frame)
    pl.add_key_event("p", toggle_play)

    # ── Camera setup ───────────────────────────────────────────
    pl.reset_camera()
    # Set an elevated oblique view
    bounds = terrain.bounds
    cx = (bounds[0] + bounds[1]) / 2
    cy = (bounds[2] + bounds[3]) / 2
    cz = (bounds[4] + bounds[5]) / 2
    dx = bounds[1] - bounds[0]
    dy = bounds[3] - bounds[2]
    dist = max(dx, dy) * 0.8

    pl.camera_position = [
        (cx - dist * 0.5, cy - dist * 0.5, cz + dist * 0.6),  # camera
        (cx, cy, cz),                                            # focal point
        (0, 0, 1),                                                # view up
    ]

    print("\n  Controls:")
    print("    Mouse:  Orbit (left), Pan (middle/shift+left), Zoom (scroll)")
    print("    Keys:   Left/Right arrows = prev/next frame")
    print("            p = play/pause animation")
    print("            q = quit")
    print()

    # ── Auto-play timer ────────────────────────────────────────
    def timer_callback(step):
        if state["playing"]:
            idx = state["current_idx"] + 1
            if idx >= len(water_frames):
                idx = 0  # loop
            update_timestep(idx)

    if len(water_frames) > 1:
        pl.add_callback(timer_callback, interval=500)

    # ── Show ───────────────────────────────────────────────────
    pl.show()


# ═══════════════════════════════════════════════════════════════════════════════
#  MAIN / CLI
# ═══════════════════════════════════════════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="FLUXOS 3D Flood Viewer - Interactive visualization of "
                    "simulation results on 3D terrain.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Examples:
  # View regular mesh results
  python fluxos_viewer.py --results-dir ./Results --dem ./terrain.asc

  # View triangular mesh results
  python fluxos_viewer.py --results-dir ./fluxos_out --dem ./terrain.asc --mesh-type triangular

  # Custom settings
  python fluxos_viewer.py --results-dir ./Results --dem ./terrain.asc \\
      --variable velocity --exaggeration 3.0 --water-opacity 0.8 --clim 0 1.5
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
        help="Mesh type: 'regular' for Cartesian (.txt), 'triangular' for unstructured (.vtu). Default: regular",
    )
    parser.add_argument(
        "--variable", choices=["h", "velocity", "conc_SW"], default="h",
        help="Variable to visualize. Default: h (water depth)",
    )
    parser.add_argument(
        "--exaggeration", type=float, default=5.0,
        help="Vertical exaggeration factor. Default: 5.0",
    )
    parser.add_argument(
        "--water-opacity", type=float, default=0.7,
        help="Water surface opacity (0=transparent, 1=opaque). Default: 0.7",
    )
    parser.add_argument(
        "--clim", type=float, nargs=2, metavar=("MIN", "MAX"), default=None,
        help="Color limits for the visualization variable (e.g., --clim 0 2.0)",
    )
    parser.add_argument(
        "--terrain-only", action="store_true",
        help="Show terrain DEM only (no simulation results)",
    )

    args = parser.parse_args()

    # ── Validate inputs ──────────────────────────────────────
    if not os.path.isfile(args.dem):
        print(f"ERROR: DEM file not found: {args.dem}")
        sys.exit(1)

    if not args.terrain_only and not os.path.isdir(args.results_dir):
        print(f"ERROR: Results directory not found: {args.results_dir}")
        sys.exit(1)

    # ── Step 1: Read DEM ─────────────────────────────────────
    print("\n[1/4] Reading DEM...")
    dem, meta = read_dem_asc(args.dem)

    # ── Step 2: Build terrain mesh ───────────────────────────
    print("[2/4] Building terrain mesh...")
    terrain = build_terrain_mesh(dem, meta, args.exaggeration)
    print(f"  Terrain mesh: {terrain.n_points} points")

    if args.terrain_only:
        print("\n  Displaying terrain only...")
        pl = pv.Plotter()
        pl.add_mesh(terrain, scalars="elevation", cmap="gist_earth",
                    nan_opacity=0.0, show_scalar_bar=True,
                    scalar_bar_args={"title": "Elevation [m]"})
        pl.add_text("FLUXOS Terrain", font_size=12, position="upper_left")
        pl.show()
        return

    # ── Step 3: Read simulation results ──────────────────────
    print(f"[3/4] Reading {args.mesh_type} mesh results...")

    if args.mesh_type == "regular":
        results = read_regular_results(args.results_dir)
    else:
        results = read_triangular_results(args.results_dir)

    if not results:
        print("  No results found. Showing terrain only.")
        launch_viewer(terrain, [], [], args.variable)
        return

    # ── Step 4: Build water frames ───────────────────────────
    print("[4/4] Building water surface frames...")
    water_frames = []
    times = []

    if args.mesh_type == "regular":
        for i, r in enumerate(results):
            wf = build_water_regular(
                dem, meta, r["data"], args.variable, args.exaggeration
            )
            water_frames.append(wf)
            times.append(r["time"])
            wet_cells = 0 if wf is None else wf.n_cells
            if (i + 1) % max(1, len(results) // 5) == 0 or i == 0:
                print(f"  Frame {i+1}/{len(results)}: "
                      f"t={format_time(r['time'])}, {wet_cells} wet cells")
    else:
        for i, r in enumerate(results):
            wf = build_water_triangular(
                r["mesh"], args.variable, args.exaggeration
            )
            water_frames.append(wf)
            times.append(r["time"])
            wet_cells = 0 if wf is None else wf.n_cells
            if (i + 1) % max(1, len(results) // 5) == 0 or i == 0:
                print(f"  Frame {i+1}/{len(results)}: "
                      f"t={format_time(r['time'])}, {wet_cells} wet cells")

    # Count valid frames
    valid_frames = sum(1 for wf in water_frames if wf is not None)
    print(f"\n  {valid_frames}/{len(water_frames)} frames have water data")

    # ── Launch viewer ────────────────────────────────────────
    print("\nLaunching 3D viewer...")
    launch_viewer(
        terrain,
        water_frames,
        times,
        args.variable,
        water_opacity=args.water_opacity,
        clim=tuple(args.clim) if args.clim else None,
    )


if __name__ == "__main__":
    main()
