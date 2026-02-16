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
import shutil
import struct
import subprocess
import sys
import xml.etree.ElementTree as ET
import math
import zlib

import numpy as np

try:
    from pyproj import Proj
except ImportError:
    print("ERROR: pyproj is required.  Install with: pip install pyproj")
    sys.exit(1)

# Video export (optional -- needed only for --export-video)
try:
    from PIL import Image, ImageDraw, ImageFont
    _HAS_PIL = True
except ImportError:
    _HAS_PIL = False


def _check_ffmpeg():
    """Return True if ffmpeg is on PATH."""
    return shutil.which("ffmpeg") is not None


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

def _extract_velocity_regular(data, meta):
    """Extract (vx, vy) from regular mesh data, with discharge fallback."""
    ux = data[:, 6]
    uy = data[:, 7]
    if np.any(ux != 0) or np.any(uy != 0):
        return ux, uy
    h = data[:, 5]
    cs = meta["cellsize"]
    qx_dxy = data[:, 8]
    qy_dxy = data[:, 9]
    vx = np.where(h > 1e-6, qx_dxy / (h * cs), 0.0)
    vy = np.where(h > 1e-6, qy_dxy / (h * cs), 0.0)
    return vx, vy


def _extract_velocity_triangular(cell_data):
    """Extract (vx, vy) from triangular mesh cell data, with discharge fallback."""
    if "velocity" in cell_data:
        vel = cell_data["velocity"]
        if vel.ndim == 2:
            vx, vy = vel[:, 0], vel[:, 1]
        else:
            vx, vy = vel, np.zeros_like(vel)
        if np.any(vx != 0) or np.any(vy != 0):
            return vx, vy
    if "discharge" in cell_data and "h" in cell_data:
        disc = cell_data["discharge"]
        h = cell_data["h"]
        if disc.ndim == 2:
            vx = np.where(h > 1e-6, disc[:, 0] / h, 0.0)
            vy = np.where(h > 1e-6, disc[:, 1] / h, 0.0)
            return vx, vy
    n = len(cell_data.get("h", []))
    return np.zeros(n), np.zeros(n)


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


# ═══════════════════════════════════════════════════════════════════════════════
#  HILLSHADE (for video background)
# ═══════════════════════════════════════════════════════════════════════════════

def _compute_hillshade(dem, cellsize, azimuth=315, altitude=45):
    """Compute a hillshade from a DEM.

    Returns an (H, W) uint8 grayscale image (0-255).
    NaN cells in the DEM are rendered as neutral gray (128).
    """
    filled = dem.copy()
    nan_mask = np.isnan(filled)
    if np.any(nan_mask):
        # Fill NaN with nearest valid value to avoid gradient edge artifacts
        from scipy.ndimage import distance_transform_edt
        try:
            idx = distance_transform_edt(nan_mask, return_distances=False,
                                          return_indices=True)
            filled = filled[tuple(idx)]
        except Exception:
            filled[nan_mask] = np.nanmean(dem)

    dy, dx = np.gradient(filled, cellsize)
    slope = np.arctan(np.sqrt(dx**2 + dy**2))
    aspect = np.arctan2(-dy, dx)

    az_rad = np.radians(azimuth)
    alt_rad = np.radians(altitude)

    hs = (np.sin(alt_rad) * np.cos(slope) +
          np.cos(alt_rad) * np.sin(slope) * np.cos(az_rad - aspect))
    hs = np.clip(hs, 0, 1)
    hs_uint8 = (hs * 255).astype(np.uint8)
    hs_uint8[nan_mask] = 128
    return hs_uint8


# ═══════════════════════════════════════════════════════════════════════════════
#  FLOW VISUALIZATION UTILITIES
# ═══════════════════════════════════════════════════════════════════════════════

def _velocity_colormap(t):
    """Map normalised velocity t∈[0,1] to (R,G,B).  Warm ramp: white→yellow→orange→red."""
    t = max(0.0, min(1.0, t))
    # Piecewise linear through 5 stops
    stops = [
        (0.00, 255, 255, 255),
        (0.25, 255, 255, 100),
        (0.50, 255, 180,  50),
        (0.75, 255,  80,  20),
        (1.00, 200,  20,  20),
    ]
    for i in range(len(stops) - 1):
        t0, r0, g0, b0 = stops[i]
        t1, r1, g1, b1 = stops[i + 1]
        if t <= t1:
            f = (t - t0) / (t1 - t0) if t1 > t0 else 0.0
            return (int(r0 + (r1 - r0) * f),
                    int(g0 + (g1 - g0) * f),
                    int(b0 + (b1 - b0) * f))
    return (200, 20, 20)


def _bilinear_interpolate(grid, r, c):
    """Bilinear interpolation of a 2D grid at fractional (row, col)."""
    h, w = grid.shape
    if r < 0 or r >= h - 1 or c < 0 or c >= w - 1:
        # Clamp to nearest edge value for edge pixels
        ri = max(0, min(h - 1, int(round(r))))
        ci = max(0, min(w - 1, int(round(c))))
        return float(grid[ri, ci])
    r0 = int(math.floor(r))
    c0 = int(math.floor(c))
    r1 = min(r0 + 1, h - 1)
    c1 = min(c0 + 1, w - 1)
    fr = r - r0
    fc = c - c0
    val = (grid[r0, c0] * (1 - fr) * (1 - fc) +
           grid[r0, c1] * (1 - fr) * fc +
           grid[r1, c0] * fr * (1 - fc) +
           grid[r1, c1] * fr * fc)
    return float(val)


# ═══════════════════════════════════════════════════════════════════════════════
#  ADVECTED PARTICLE TRAIL ENGINE  (earth.nullschool-style)
# ═══════════════════════════════════════════════════════════════════════════════

def _bilinear_interpolate_vec(grid, rows, cols):
    """Vectorized bilinear interpolation for arrays of (row, col) positions.

    grid : float64[H, W]
    rows, cols : float64[N]
    Returns : float64[N]
    """
    H, W = grid.shape
    r0 = np.clip(np.floor(rows).astype(np.intp), 0, H - 2)
    c0 = np.clip(np.floor(cols).astype(np.intp), 0, W - 2)
    r1 = r0 + 1
    c1 = c0 + 1
    fr = np.clip(rows - r0, 0.0, 1.0)
    fc = np.clip(cols - c0, 0.0, 1.0)
    return (grid[r0, c0] * (1.0 - fr) * (1.0 - fc) +
            grid[r0, c1] * (1.0 - fr) * fc +
            grid[r1, c0] * fr * (1.0 - fc) +
            grid[r1, c1] * fr * fc)


def _rk4_advect(grid_vx, grid_vy, rows, cols, dt):
    """Vectorized RK4 integration step for particle advection.

    Returns (new_rows, new_cols, speed) where speed = |v| at the new position.
    Velocity convention: dr = -vy (north = row-decreasing), dc = vx.
    """
    def vel_at(r, c):
        vx = _bilinear_interpolate_vec(grid_vx, r, c)
        vy = _bilinear_interpolate_vec(grid_vy, r, c)
        return -vy, vx  # (dr, dc) in pixel space

    k1r, k1c = vel_at(rows, cols)
    k2r, k2c = vel_at(rows + 0.5 * dt * k1r, cols + 0.5 * dt * k1c)
    k3r, k3c = vel_at(rows + 0.5 * dt * k2r, cols + 0.5 * dt * k2c)
    k4r, k4c = vel_at(rows + dt * k3r, cols + dt * k3c)

    new_rows = rows + (dt / 6.0) * (k1r + 2.0 * k2r + 2.0 * k3r + k4r)
    new_cols = cols + (dt / 6.0) * (k1c + 2.0 * k2c + 2.0 * k3c + k4c)

    vx_f = _bilinear_interpolate_vec(grid_vx, new_rows, new_cols)
    vy_f = _bilinear_interpolate_vec(grid_vy, new_rows, new_cols)
    speed = np.sqrt(vx_f ** 2 + vy_f ** 2)
    return new_rows, new_cols, speed


# ── Advected-mode colormap ────────────────────────────────────────────────────

_ADVECTED_STOPS = np.array([
    # t,     R,    G,    B
    [0.00,  15,   20,   80],   # deep navy
    [0.15,  30,   70,  170],   # blue
    [0.40,  60,  190,  255],   # cyan
    [0.70, 200,  255,  255],   # light cyan / near-white
    [0.90, 255,  255,  180],   # warm white
    [1.00, 255,  230,  80],    # golden yellow
], dtype=np.float64)


def _advected_speed_colormap(speeds, vel_scale):
    """Vectorized: map speed array -> float32[N, 3] RGB in [0, 1]."""
    t = np.clip(speeds / max(vel_scale, 1e-12), 0.0, 1.0)
    stops = _ADVECTED_STOPS
    ts = stops[:, 0]
    rgb = np.zeros((len(t), 3), dtype=np.float32)
    for i in range(len(ts) - 1):
        mask = (t >= ts[i]) & (t <= ts[i + 1])
        if not np.any(mask):
            continue
        frac = (t[mask] - ts[i]) / max(ts[i + 1] - ts[i], 1e-12)
        for ch in range(3):
            v0 = stops[i, 1 + ch] / 255.0
            v1 = stops[i + 1, 1 + ch] / 255.0
            rgb[mask, ch] = v0 + (v1 - v0) * frac
    return rgb


def _advected_speed_colormap_scalar(speed, vel_scale):
    """Scalar version for colorbar drawing. Returns (R, G, B) ints 0-255."""
    t = max(0.0, min(1.0, speed / max(vel_scale, 1e-12)))
    stops = _ADVECTED_STOPS
    ts = stops[:, 0]
    for i in range(len(ts) - 1):
        if t <= ts[i + 1]:
            frac = (t - ts[i]) / max(ts[i + 1] - ts[i], 1e-12)
            r = int(stops[i, 1] + (stops[i + 1, 1] - stops[i, 1]) * frac)
            g = int(stops[i, 2] + (stops[i + 1, 2] - stops[i, 2]) * frac)
            b = int(stops[i, 3] + (stops[i + 1, 3] - stops[i, 3]) * frac)
            return (r, g, b)
    return (int(stops[-1, 1]), int(stops[-1, 2]), int(stops[-1, 3]))


# ── Velocity field temporal interpolation ─────────────────────────────────────

def _lerp_velocity_fields(vx0, vy0, vx1, vy1, frac):
    """Linear interpolation between two velocity fields."""
    return vx0 * (1.0 - frac) + vx1 * frac, vy0 * (1.0 - frac) + vy1 * frac


# ── Anti-aliased particle splatting ───────────────────────────────────────────

def _draw_particles_antialiased(canvas, rows, cols, colors, brightness):
    """Splat particles onto float32[H, W, 3] canvas with sub-pixel AA.

    canvas    : float32[H, W, 3]  -- modified in-place (values in [0, ~2])
    rows, cols: float64[N]
    colors    : float32[N, 3]     -- RGB in [0, 1]
    brightness: float64[N]        -- [0, 1]
    """
    H, W = canvas.shape[:2]
    r_int = np.floor(rows).astype(np.intp)
    c_int = np.floor(cols).astype(np.intp)
    fr = (rows - r_int).astype(np.float32)
    fc = (cols - c_int).astype(np.float32)

    # Bilinear weight distribution across 2×2 neighbourhood
    for dr in range(2):
        for dc in range(2):
            rr = r_int + dr
            cc = c_int + dc
            wr = fr if dr == 1 else (1.0 - fr)
            wc = fc if dc == 1 else (1.0 - fc)
            w = wr * wc * brightness.astype(np.float32)
            valid = (rr >= 0) & (rr < H) & (cc >= 0) & (cc < W)
            idx_r = rr[valid]
            idx_c = cc[valid]
            wv = w[valid]
            for ch in range(3):
                np.add.at(canvas[:, :, ch], (idx_r, idx_c),
                          colors[valid, ch] * wv)


# ── Particle System ───────────────────────────────────────────────────────────

class ParticleSystem:
    """Stateful particle system for earth.nullschool-style advection trails."""

    def __init__(self, n_particles, wet_mask, rng_seed=42):
        self.n = n_particles
        self.row = np.zeros(n_particles, dtype=np.float64)
        self.col = np.zeros(n_particles, dtype=np.float64)
        self.age = np.zeros(n_particles, dtype=np.int32)
        self.max_age = np.zeros(n_particles, dtype=np.int32)
        self.speed = np.zeros(n_particles, dtype=np.float64)
        self.stagnant_count = np.zeros(n_particles, dtype=np.int32)
        self.alive = np.zeros(n_particles, dtype=bool)
        self.rng = np.random.default_rng(rng_seed)
        self.wet_mask = wet_mask
        self._wet_indices = np.argwhere(wet_mask)
        self._spawn_at(np.arange(n_particles))

    def _spawn_at(self, indices):
        """Respawn specific particles at random wet cell locations."""
        n = len(indices)
        if len(self._wet_indices) == 0:
            self.alive[indices] = False
            return
        chosen = self.rng.choice(len(self._wet_indices), size=n,
                                  replace=(n > len(self._wet_indices)))
        base = self._wet_indices[chosen].astype(np.float64)
        base += self.rng.uniform(-0.4, 0.4, base.shape)
        self.row[indices] = base[:, 0]
        self.col[indices] = base[:, 1]
        self.age[indices] = 0
        self.max_age[indices] = self.rng.integers(60, 200, size=n)
        self.speed[indices] = 0.0
        self.stagnant_count[indices] = 0
        self.alive[indices] = True

    def advect(self, grid_vx, grid_vy, dt_pixels):
        """Advance all alive particles using RK4."""
        mask = self.alive
        if not np.any(mask):
            return
        idx = np.where(mask)[0]
        r = self.row[idx]
        c = self.col[idx]
        new_r, new_c, spd = _rk4_advect(grid_vx, grid_vy, r, c, dt_pixels)
        self.row[idx] = new_r
        self.col[idx] = new_c
        self.speed[idx] = spd
        self.age[idx] += 1

        H, W = grid_vx.shape
        # Kill: out of bounds
        oob = (new_r < 0) | (new_r >= H - 1) | (new_c < 0) | (new_c >= W - 1)
        # Kill: in dry cell
        r_int = np.clip(np.round(new_r).astype(np.intp), 0, H - 1)
        c_int = np.clip(np.round(new_c).astype(np.intp), 0, W - 1)
        dry = ~self.wet_mask[r_int, c_int]
        # Kill: exceeded max age
        old = self.age[idx] > self.max_age[idx]
        # Kill: stagnant
        slow = spd < 0.005
        self.stagnant_count[idx] = np.where(slow,
                                             self.stagnant_count[idx] + 1, 0)
        stagnant = self.stagnant_count[idx] > 15
        dead = oob | dry | old | stagnant
        self.alive[idx[dead]] = False

    def respawn_dead(self):
        """Respawn all dead particles at new random wet locations."""
        dead_idx = np.where(~self.alive)[0]
        if len(dead_idx) > 0:
            self._spawn_at(dead_idx)

    def update_wet_mask(self, new_mask):
        """Update wet mask (e.g., when flooding extent changes)."""
        self.wet_mask = new_mask
        self._wet_indices = np.argwhere(new_mask)

    def get_draw_data(self):
        """Return (rows, cols, speeds, ages, max_ages) for alive particles."""
        mask = self.alive
        return (self.row[mask], self.col[mask], self.speed[mask],
                self.age[mask], self.max_age[mask])


# ── Velocity Field Cache (sliding window) ────────────────────────────────────

class VelocityFieldCache:
    """Lazy-loading cache with a sliding window of velocity fields."""

    def __init__(self, results, mesh_type, dem, meta, h_min, scale):
        self.results = results
        self.mesh_type = mesh_type
        self.dem = dem
        self.meta = meta
        self.h_min = h_min
        self.scale = scale
        self._cache = {}        # idx -> (vx_hr, vy_hr, wet_hr)
        self._times = np.array([r["time"] for r in results], dtype=np.float64)
        self._nrows, self._ncols = dem.shape
        self._cs = meta["cellsize"]
        self._xll = meta["xllcorner"]
        self._yll = meta["yllcorner"]

    def get(self, idx):
        """Get upscaled velocity field for timestep index."""
        if idx not in self._cache:
            self._evict_if_full()
            self._load(idx)
        return self._cache[idx]

    def _load(self, idx):
        r = self.results[idx]
        if self.mesh_type == "regular":
            vx, vy, wet = _build_velocity_grids_regular(
                self.dem, self.meta, r["data"], self.h_min)
        else:
            vx, vy, wet = _build_velocity_grids_triangular(
                r["points"], r["cells"], r["cell_data"],
                self.h_min, self._nrows, self._ncols,
                self._xll, self._yll, self._cs)
        s = self.scale
        vx_hr = _upscale_velocity_grid(vx, s)
        vy_hr = _upscale_velocity_grid(vy, s)
        wet_hr = _upscale_velocity_grid(
            wet.astype(np.uint8), s).astype(bool)
        self._cache[idx] = (vx_hr, vy_hr, wet_hr)

    def _evict_if_full(self):
        while len(self._cache) > 4:
            oldest = min(self._cache.keys())
            del self._cache[oldest]

    def find_bracket(self, sim_time):
        """Return (idx_lo, idx_hi, frac) for temporal interpolation."""
        times = self._times
        if sim_time <= times[0]:
            return 0, 0, 0.0
        if sim_time >= times[-1]:
            n = len(times) - 1
            return n, n, 0.0
        idx = int(np.searchsorted(times, sim_time, side='right')) - 1
        idx = max(0, min(idx, len(times) - 2))
        t0, t1 = times[idx], times[idx + 1]
        frac = (sim_time - t0) / max(t1 - t0, 1e-12)
        return idx, idx + 1, float(frac)


def _build_velocity_grids_regular(dem, meta, data, h_min):
    """Build (grid_vx, grid_vy, wet_mask) from regular mesh data."""
    nrows, ncols = dem.shape
    cs = meta["cellsize"]
    xll = meta["xllcorner"]
    yll = meta["yllcorner"]

    vx, vy = _extract_velocity_regular(data, meta)
    h_depth = data[:, 5]
    x_coords = data[:, 2]
    y_coords = data[:, 3]

    col_idx = np.round((x_coords - xll) / cs).astype(int)
    row_idx = np.round((nrows - 1) - (y_coords - yll) / cs).astype(int)
    col_idx = np.clip(col_idx, 0, ncols - 1)
    row_idx = np.clip(row_idx, 0, nrows - 1)

    grid_vx = np.zeros((nrows, ncols))
    grid_vy = np.zeros((nrows, ncols))
    wet_mask = np.zeros((nrows, ncols), dtype=bool)
    wet = h_depth > h_min
    grid_vx[row_idx[wet], col_idx[wet]] = vx[wet]
    grid_vy[row_idx[wet], col_idx[wet]] = vy[wet]
    wet_mask[row_idx[wet], col_idx[wet]] = True
    return grid_vx, grid_vy, wet_mask


def _build_velocity_grids_triangular(points, cells, cell_data, h_min,
                                      nrows, ncols, xll, yll, cs):
    """Build (grid_vx, grid_vy, wet_mask) from triangular mesh data."""
    if "h" not in cell_data:
        return (np.zeros((nrows, ncols)), np.zeros((nrows, ncols)),
                np.zeros((nrows, ncols), dtype=bool))

    h = cell_data["h"]
    vx, vy = _extract_velocity_triangular(cell_data)
    n_cells = len(cells)

    # Bin cells into pixel-sized blocks and average
    grid_vx = np.zeros((nrows, ncols))
    grid_vy = np.zeros((nrows, ncols))
    grid_cnt = np.zeros((nrows, ncols))

    for ci in range(n_cells):
        if h[ci] <= h_min:
            continue
        i0, i1, i2 = cells[ci]
        p0, p1, p2 = points[i0], points[i1], points[i2]
        cx = (p0[0] + p1[0] + p2[0]) / 3.0
        cy = (p0[1] + p1[1] + p2[1]) / 3.0
        col = int(round((cx - xll) / cs))
        row = int(round((nrows - 1) - (cy - yll) / cs))
        if 0 <= row < nrows and 0 <= col < ncols:
            grid_vx[row, col] += vx[ci]
            grid_vy[row, col] += vy[ci]
            grid_cnt[row, col] += 1

    valid = grid_cnt > 0
    grid_vx[valid] /= grid_cnt[valid]
    grid_vy[valid] /= grid_cnt[valid]
    return grid_vx, grid_vy, valid


def _draw_thick_line(img, r0, c0, r1, c1, rgba, width=1):
    """Draw a line with variable thickness using parallel Bresenham lines."""
    if width <= 1:
        _draw_line(img, r0, c0, r1, c1, rgba)
        return
    # Perpendicular direction
    dr = r1 - r0
    dc = c1 - c0
    length = math.sqrt(dr * dr + dc * dc)
    if length < 0.5:
        return
    # Normalised perpendicular
    pr = -dc / length
    pc = dr / length
    half = width / 2.0
    for i in range(width):
        offset = -half + 0.5 + i
        ro = int(round(offset * pr))
        co = int(round(offset * pc))
        _draw_line(img, r0 + ro, c0 + co, r1 + ro, c1 + co, rgba)


def _draw_segment(img, r0, c0, r1, c1, r_col, g_col, b_col, alpha, width=1):
    """Draw a coloured line segment with specified RGBA and width."""
    a = max(0, min(255, int(alpha)))
    if a < 5:
        return
    rgba = np.array([r_col, g_col, b_col, a], dtype=np.uint8)
    _draw_thick_line(img, int(round(r0)), int(round(c0)),
                     int(round(r1)), int(round(c1)), rgba, width)


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
#  VELOCITY ARROW OVERLAY
# ═══════════════════════════════════════════════════════════════════════════════

_ARROW_FG = np.array([255, 255, 255, 220], dtype=np.uint8)
_ARROW_BG = np.array([0, 0, 0, 160], dtype=np.uint8)


def _draw_arrow_outlined(img, r_from, c_from, r_to, c_to, head_size=3):
    """Draw an arrow with a dark outline for contrast."""
    for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
        _draw_arrow(img, r_from + dr, c_from + dc,
                    r_to + dr, c_to + dc, _ARROW_BG, head_size)
    _draw_arrow(img, r_from, c_from, r_to, c_to, _ARROW_FG, head_size)


def draw_velocity_arrows_regular(img, dem, meta, data, h_min=0.001,
                                  arrow_spacing=15, max_arrow_px=12,
                                  vel_scale=None):
    """Overlay velocity arrows onto an RGBA image for regular mesh (in-place)."""
    nrows, ncols = dem.shape
    cs = meta["cellsize"]
    xll = meta["xllcorner"]
    yll = meta["yllcorner"]

    vx, vy = _extract_velocity_regular(data, meta)
    h_depth = data[:, 5]
    x_coords = data[:, 2]
    y_coords = data[:, 3]

    col_idx = np.round((x_coords - xll) / cs).astype(int)
    row_idx = np.round((nrows - 1) - (y_coords - yll) / cs).astype(int)
    col_idx = np.clip(col_idx, 0, ncols - 1)
    row_idx = np.clip(row_idx, 0, nrows - 1)

    grid_vx = np.zeros((nrows, ncols))
    grid_vy = np.zeros((nrows, ncols))
    wet = h_depth > h_min
    grid_vx[row_idx[wet], col_idx[wet]] = vx[wet]
    grid_vy[row_idx[wet], col_idx[wet]] = vy[wet]

    if vel_scale is None or vel_scale < 1e-9:
        vmag_all = np.sqrt(vx[wet] ** 2 + vy[wet] ** 2)
        valid = vmag_all[vmag_all > 1e-6]
        vel_scale = float(np.percentile(valid, 95)) if len(valid) > 0 else 1.0

    half = arrow_spacing // 2
    for r in range(half, nrows, arrow_spacing):
        for c in range(half, ncols, arrow_spacing):
            gvx = grid_vx[r, c]
            gvy = grid_vy[r, c]
            vmag = math.sqrt(gvx * gvx + gvy * gvy)
            if vmag < 1e-6:
                continue
            length_px = min(max_arrow_px, max_arrow_px * vmag / vel_scale)
            if length_px < 2:
                continue
            # vy positive = northward = pixel row decreasing
            dr = -gvy / vmag * length_px
            dc = gvx / vmag * length_px
            r_to = int(round(r + dr))
            c_to = int(round(c + dc))
            hs = max(2, int(length_px * 0.35))
            _draw_arrow_outlined(img, r, c, r_to, c_to, hs)


def draw_velocity_arrows_triangular(img, points, cells, cell_data, h_min,
                                     nrows, ncols, xll, yll, cs,
                                     arrow_spacing=15, max_arrow_px=12,
                                     vel_scale=None):
    """Overlay velocity arrows onto an RGBA image for triangular mesh (in-place)."""
    if "h" not in cell_data:
        return

    h = cell_data["h"]
    vx, vy = _extract_velocity_triangular(cell_data)

    # Compute cell centroids in pixel coordinates
    n_cells = len(cells)
    cx_px = np.zeros(n_cells)
    cy_px = np.zeros(n_cells)
    for ci, (i0, i1, i2) in enumerate(cells):
        p0, p1, p2 = points[i0], points[i1], points[i2]
        cx_utm = (p0[0] + p1[0] + p2[0]) / 3.0
        cy_utm = (p0[1] + p1[1] + p2[1]) / 3.0
        cx_px[ci] = (cx_utm - xll) / cs
        cy_px[ci] = (nrows - 1) - (cy_utm - yll) / cs

    wet = h > h_min

    if vel_scale is None or vel_scale < 1e-9:
        vmag_all = np.sqrt(vx[wet] ** 2 + vy[wet] ** 2)
        valid = vmag_all[vmag_all > 1e-6]
        vel_scale = float(np.percentile(valid, 95)) if len(valid) > 0 else 1.0

    # Bin cells into arrow_spacing blocks and average velocities
    n_bins_r = max(1, nrows // arrow_spacing + 1)
    n_bins_c = max(1, ncols // arrow_spacing + 1)
    bin_vx = np.zeros((n_bins_r, n_bins_c))
    bin_vy = np.zeros((n_bins_r, n_bins_c))
    bin_cnt = np.zeros((n_bins_r, n_bins_c))

    for ci in range(n_cells):
        if not wet[ci]:
            continue
        br = int(cy_px[ci] / arrow_spacing)
        bc = int(cx_px[ci] / arrow_spacing)
        if 0 <= br < n_bins_r and 0 <= bc < n_bins_c:
            bin_vx[br, bc] += vx[ci]
            bin_vy[br, bc] += vy[ci]
            bin_cnt[br, bc] += 1

    half = arrow_spacing // 2
    for br in range(n_bins_r):
        for bc in range(n_bins_c):
            if bin_cnt[br, bc] < 1:
                continue
            avg_vx = bin_vx[br, bc] / bin_cnt[br, bc]
            avg_vy = bin_vy[br, bc] / bin_cnt[br, bc]
            vmag = math.sqrt(avg_vx * avg_vx + avg_vy * avg_vy)
            if vmag < 1e-6:
                continue
            length_px = min(max_arrow_px, max_arrow_px * vmag / vel_scale)
            if length_px < 2:
                continue
            r_orig = br * arrow_spacing + half
            c_orig = bc * arrow_spacing + half
            dr = -avg_vy / vmag * length_px
            dc = avg_vx / vmag * length_px
            r_to = int(round(r_orig + dr))
            c_to = int(round(c_orig + dc))
            hs = max(2, int(length_px * 0.35))
            _draw_arrow_outlined(img, r_orig, c_orig, r_to, c_to, hs)


# ═══════════════════════════════════════════════════════════════════════════════
#  STREAMLINE VISUALIZATION
# ═══════════════════════════════════════════════════════════════════════════════

def _integrate_streamline(grid_vx, grid_vy, r_start, c_start,
                           max_steps=20, step_size=1.5):
    """Euler integration of a streamline.  Returns [(r, c, vmag), ...]."""
    nrows, ncols = grid_vx.shape
    path = []
    r, c = float(r_start), float(c_start)
    for _ in range(max_steps):
        vx = _bilinear_interpolate(grid_vx, r, c)
        vy = _bilinear_interpolate(grid_vy, r, c)
        vmag = math.sqrt(vx * vx + vy * vy)
        if vmag < 1e-6:
            break
        path.append((r, c, vmag))
        # vy positive = northward = pixel row decreasing
        dr = -vy / vmag * step_size
        dc = vx / vmag * step_size
        r_new = r + dr
        c_new = c + dc
        if r_new < 0 or r_new >= nrows or c_new < 0 or c_new >= ncols:
            break
        r, c = r_new, c_new
    # Append final point
    if path:
        vx = _bilinear_interpolate(grid_vx, r, c)
        vy = _bilinear_interpolate(grid_vy, r, c)
        vmag = math.sqrt(vx * vx + vy * vy)
        path.append((r, c, vmag))
    return path


def _generate_seeds(wet_mask, density=15, jitter=True, seed=42):
    """Generate seed points on a regular grid across wet cells with optional jitter."""
    nrows, ncols = wet_mask.shape
    rng = np.random.default_rng(seed)
    seeds = []
    half = density // 2
    for row in range(half, nrows, density):
        for col in range(half, ncols, density):
            if wet_mask[row, col]:
                if jitter:
                    jr = rng.integers(-density // 4, density // 4 + 1)
                    jc = rng.integers(-density // 4, density // 4 + 1)
                    r = max(0, min(nrows - 1, row + jr))
                    c = max(0, min(ncols - 1, col + jc))
                    if wet_mask[r, c]:
                        seeds.append((float(r), float(c)))
                    else:
                        seeds.append((float(row), float(col)))
                else:
                    seeds.append((float(row), float(col)))
    return seeds


def _draw_streamline_path(img, path, vel_scale):
    """Draw a single streamline path with tapered opacity, velocity color, and outline."""
    n = len(path)
    if n < 2:
        return
    for i in range(n - 1):
        r0, c0, vmag0 = path[i]
        r1, c1, vmag1 = path[i + 1]
        vmag = (vmag0 + vmag1) * 0.5
        t = min(vmag / vel_scale, 1.0) if vel_scale > 1e-9 else 0.5
        # Taper: bright at head (i=0), fading toward tail
        progress = i / max(n - 1, 1)
        alpha = 220 * (1.0 - progress * 0.75)
        # Width: 1px slow → 3px fast
        w = 1 + int(t * 2)
        # Color
        rc, gc, bc = _velocity_colormap(t)
        # Dark outline
        outline_a = max(0, int(alpha * 0.6))
        for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            _draw_segment(img, r0 + dr, c0 + dc, r1 + dr, c1 + dc,
                          0, 0, 0, outline_a, w)
        # Foreground
        _draw_segment(img, r0, c0, r1, c1, rc, gc, bc, alpha, w)


def draw_streamlines(img, grid_vx, grid_vy, wet_mask, vel_scale,
                      density=15, max_steps=20):
    """Draw streamlines onto an RGBA image (in-place)."""
    seeds = _generate_seeds(wet_mask, density)
    for r_seed, c_seed in seeds:
        path = _integrate_streamline(grid_vx, grid_vy, r_seed, c_seed,
                                      max_steps=max_steps, step_size=1.5)
        _draw_streamline_path(img, path, vel_scale)


# ═══════════════════════════════════════════════════════════════════════════════
#  PARTICLE TRACE VISUALIZATION
# ═══════════════════════════════════════════════════════════════════════════════

def _generate_particle_positions(wet_mask, count=2000, seed=42):
    """Generate random particle positions across wet cells."""
    wet_idx = np.argwhere(wet_mask)
    if len(wet_idx) == 0:
        return np.zeros((0, 2))
    rng = np.random.default_rng(seed)
    n = min(count, len(wet_idx))
    chosen = rng.choice(len(wet_idx), size=n, replace=(count > len(wet_idx)))
    positions = wet_idx[chosen].astype(float)
    # Sub-pixel jitter
    positions += rng.uniform(-0.4, 0.4, positions.shape)
    return positions


def draw_particles(img, grid_vx, grid_vy, wet_mask, vel_scale,
                    count=2000, tail_length=8, frame_index=0):
    """Draw particle streaks onto an RGBA image (in-place).

    frame_index shifts the random seed so particles move between frames,
    creating an animation effect in the KML time slider.
    """
    nrows, ncols = grid_vx.shape
    seed = 1000 + frame_index * 7919
    positions = _generate_particle_positions(wet_mask, count, seed)
    if len(positions) == 0:
        return

    step_size = 0.8
    for pi in range(len(positions)):
        r, c = float(positions[pi, 0]), float(positions[pi, 1])
        trail = [(r, c)]
        for _ in range(tail_length):
            vx = _bilinear_interpolate(grid_vx, r, c)
            vy = _bilinear_interpolate(grid_vy, r, c)
            vmag = math.sqrt(vx * vx + vy * vy)
            if vmag < 1e-6:
                break
            dr = -vy / vmag * step_size
            dc = vx / vmag * step_size
            r_new = r + dr
            c_new = c + dc
            if r_new < 0 or r_new >= nrows or c_new < 0 or c_new >= ncols:
                break
            r, c = r_new, c_new
            trail.append((r, c))

        if len(trail) < 2:
            continue

        # Draw the streak: head is bright, tail fades
        for i in range(len(trail) - 1):
            r0, c0 = trail[i]
            r1, c1 = trail[i + 1]
            vx = _bilinear_interpolate(grid_vx, r0, c0)
            vy = _bilinear_interpolate(grid_vy, r0, c0)
            vmag = math.sqrt(vx * vx + vy * vy)
            t = min(vmag / vel_scale, 1.0) if vel_scale > 1e-9 else 0.5
            # Head (i=0) is bright, tail fades
            progress = i / max(len(trail) - 1, 1)
            alpha = 200 * (1.0 - progress * 0.8)
            rc, gc, bc = _velocity_colormap(t)
            _draw_segment(img, r0, c0, r1, c1, rc, gc, bc, alpha, 1)


# ═══════════════════════════════════════════════════════════════════════════════
#  FLOW OVERLAY DISPATCHER
# ═══════════════════════════════════════════════════════════════════════════════

def _draw_flow_overlay(img, grid_vx, grid_vy, wet_mask, vel_scale,
                        flow_config, frame_index=0,
                        dem=None, meta=None, data=None,
                        points=None, cells=None, cell_data=None,
                        h_min=0.001, mesh_type="regular",
                        nrows=0, ncols=0, xll=0, yll=0, cs=1):
    """Dispatch to the appropriate flow visualization mode."""
    mode = flow_config.get("MODE", "arrows")

    if mode == "arrows":
        spacing = flow_config.get("ARROW_SPACING", 15)
        max_px = flow_config.get("ARROW_MAX_PX", 12)
        if mesh_type == "regular":
            draw_velocity_arrows_regular(img, dem, meta, data, h_min,
                                          spacing, max_px, vel_scale)
        else:
            draw_velocity_arrows_triangular(img, points, cells, cell_data,
                                             h_min, nrows, ncols, xll, yll, cs,
                                             spacing, max_px, vel_scale)

    elif mode == "streamlines":
        density = flow_config.get("STREAMLINE_DENSITY", 15)
        length = flow_config.get("STREAMLINE_LENGTH", 20)
        draw_streamlines(img, grid_vx, grid_vy, wet_mask, vel_scale,
                          density=density, max_steps=length)

    elif mode == "particles":
        count = flow_config.get("PARTICLE_COUNT", 2000)
        tail = flow_config.get("PARTICLE_TAIL", 8)
        draw_particles(img, grid_vx, grid_vy, wet_mask, vel_scale,
                        count=count, tail_length=tail, frame_index=frame_index)


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


def _draw_arrow(img, r_from, c_from, r_to, c_to, rgba, head_size=3):
    """Draw an arrow (shaft + two barbs) onto an RGBA image."""
    dr = r_to - r_from
    dc = c_to - c_from
    length = math.sqrt(dr * dr + dc * dc)
    if length < 2:
        return
    _draw_line(img, r_from, c_from, r_to, c_to, rgba)
    angle = math.atan2(dr, dc)
    hs = max(2, int(head_size))
    for barb_offset in (5 * math.pi / 6, -5 * math.pi / 6):
        ba = angle + barb_offset
        rb = int(round(r_to + hs * math.sin(ba)))
        cb = int(round(c_to + hs * math.cos(ba)))
        _draw_line(img, r_to, c_to, rb, cb, rgba)


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

def _upscale_nearest(img, scale):
    """Upscale an (H, W, 4) RGBA image by integer factor using nearest neighbour."""
    if scale <= 1:
        return img
    return np.repeat(np.repeat(img, scale, axis=0), scale, axis=1)


def _upscale_velocity_grid(grid, scale):
    """Upscale a 2D velocity grid by repeating each cell.

    Values are preserved (not interpolated) so that bilinear interpolation
    during streamline/particle integration operates on a finer pixel grid.
    """
    if scale <= 1:
        return grid
    return np.repeat(np.repeat(grid, scale, axis=0), scale, axis=1)


def _generate_frames(dem, meta, results, mesh_type, variable, clim,
                     opacity, h_min, draw_arrows, flow_config, vel_scale,
                     scale=1):
    """Yield (frame_index, time_seconds, rgba, wet_count) per timestep.

    The yielded rgba is an (H*scale, W*scale, 4) uint8 RGBA image
    with flow overlay already composited.  Dry frames are skipped.
    """
    nrows, ncols = dem.shape
    cs = meta["cellsize"]
    xll = meta["xllcorner"]
    yll = meta["yllcorner"]

    if flow_config is None:
        flow_config = {"MODE": "arrows", "ARROW_SPACING": 15, "ARROW_MAX_PX": 12,
                       "STREAMLINE_DENSITY": 15, "STREAMLINE_LENGTH": 20,
                       "PARTICLE_COUNT": 2000, "PARTICLE_TAIL": 8}

    for ri, r in enumerate(results):
        t_sec = r["time"]

        # Rasterize at native DEM resolution
        if mesh_type == "regular":
            rgba, wet = rasterize_regular(
                dem, meta, r["data"], variable, clim, h_min, opacity)
        else:
            rgba, wet = rasterize_triangular(
                r["points"], r["cells"], r["cell_data"],
                variable, clim, h_min, opacity, nrows, ncols, xll, yll, cs)

        if wet == 0:
            continue

        # Upscale depth colours to high-res canvas
        rgba = _upscale_nearest(rgba, scale)

        # Overlay flow visualisation (drawn at high-res)
        if draw_arrows and vel_scale is not None:
            if mesh_type == "regular":
                gvx, gvy, wmask = _build_velocity_grids_regular(
                    dem, meta, r["data"], h_min)
            else:
                gvx, gvy, wmask = _build_velocity_grids_triangular(
                    r["points"], r["cells"], r["cell_data"],
                    h_min, nrows, ncols, xll, yll, cs)
            gvx_hr = _upscale_velocity_grid(gvx, scale)
            gvy_hr = _upscale_velocity_grid(gvy, scale)
            wmask_hr = _upscale_velocity_grid(
                wmask.astype(np.uint8), scale).astype(bool)

            hr_flow_config = dict(flow_config)
            hr_flow_config["ARROW_SPACING"] = flow_config.get("ARROW_SPACING", 15) * scale
            hr_flow_config["ARROW_MAX_PX"] = flow_config.get("ARROW_MAX_PX", 12) * scale
            hr_flow_config["STREAMLINE_DENSITY"] = flow_config.get("STREAMLINE_DENSITY", 15) * scale
            hr_flow_config["STREAMLINE_LENGTH"] = flow_config.get("STREAMLINE_LENGTH", 20) * scale
            hr_flow_config["PARTICLE_TAIL"] = flow_config.get("PARTICLE_TAIL", 8) * scale
            hr_flow_config["PARTICLE_COUNT"] = flow_config.get("PARTICLE_COUNT", 2000) * scale

            _draw_flow_overlay(
                rgba, gvx_hr, gvy_hr, wmask_hr, vel_scale, hr_flow_config,
                frame_index=ri, dem=dem, meta=meta,
                data=r.get("data"), points=r.get("points"),
                cells=r.get("cells"), cell_data=r.get("cell_data"),
                h_min=h_min, mesh_type=mesh_type,
                nrows=nrows * scale, ncols=ncols * scale,
                xll=xll, yll=yll, cs=cs / scale)

        yield ri, t_sec, rgba, wet


def export_kml(dem, meta, results, mesh_type, variable, clim, utm_to_ll,
               sim_start, output_dir, opacity=180, h_min=0.001,
               draw_arrows=False, flow_config=None, scale=1):
    """
    Export results as individual KML + PNG files (one per timestep).

    scale: integer supersampling factor (1 = native DEM resolution,
           3 = 3× resolution in each dimension → 9× pixels).
    """
    nrows, ncols = dem.shape
    cs = meta["cellsize"]
    xll = meta["xllcorner"]
    yll = meta["yllcorner"]

    # Compute bounding box in lat/lon
    sw_lon, sw_lat = utm_to_ll(xll, yll)
    ne_lon, ne_lat = utm_to_ll(xll + ncols * cs, yll + nrows * cs)
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
        grid_step = max(1, min(nrows, ncols) // 60)
        grid_img = rasterize_regular_grid(dem, meta, spacing=grid_step)
    else:
        grid_img = rasterize_triangular_grid(
            results[0]["points"], results[0]["cells"],
            nrows, ncols, xll, yll, cs)

    grid_img = _upscale_nearest(grid_img, scale)

    grid_png_name = "grid.png"
    grid_png_path = os.path.join(output_dir, grid_png_name)
    grid_png_bytes = _write_png(grid_img)
    with open(grid_png_path, "wb") as f:
        f.write(grid_png_bytes)

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
    out_h, out_w = grid_img.shape[:2]
    if scale > 1:
        print(f"  Grid frame: {mesh_type} mesh, {out_w}x{out_h}px ({scale}x upscale), "
              f"PNG {len(grid_png_bytes)/1024:.0f} KB")
    else:
        print(f"  Grid frame: {mesh_type} mesh, PNG {len(grid_png_bytes)/1024:.0f} KB")

    # Flow config defaults
    if flow_config is None:
        flow_config = {"MODE": "arrows", "ARROW_SPACING": 15, "ARROW_MAX_PX": 12,
                       "STREAMLINE_DENSITY": 15, "STREAMLINE_LENGTH": 20,
                       "PARTICLE_COUNT": 2000, "PARTICLE_TAIL": 8}

    # Pre-compute velocity scale
    vel_scale = None
    if draw_arrows:
        vel_scale = _auto_vel_scale(results, mesh_type, meta, h_min)
        mode = flow_config.get("MODE", "arrows")
        print(f"  Flow mode: {mode} (vel_scale={vel_scale:.4f} m/s)")

    # ── Generate frames using shared generator ────────────────
    n_total = len(results)
    for ri, t_sec, rgba, wet in _generate_frames(
            dem, meta, results, mesh_type, variable, clim,
            opacity, h_min, draw_arrows, flow_config, vel_scale, scale):

        # Time span
        t_start = sim_start + datetime.timedelta(seconds=t_sec)
        if ri + 1 < n_total:
            t_end = sim_start + datetime.timedelta(
                seconds=results[ri + 1]["time"])
        else:
            t_end = t_start + datetime.timedelta(seconds=t_sec * 0.5)

        ts_begin = t_start.strftime("%Y-%m-%dT%H:%M:%SZ")
        ts_end = t_end.strftime("%Y-%m-%dT%H:%M:%SZ")

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

        print(f"  Frame {ri+1}/{n_total}: t={format_time(t_sec)}, "
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
#  VIDEO EXPORT (GIF / MP4)
# ═══════════════════════════════════════════════════════════════════════════════

def _subsample_indices(n_total, max_frames):
    """Return evenly-spaced indices when n_total exceeds max_frames."""
    if n_total <= max_frames:
        return set(range(n_total))
    step = n_total / max_frames
    return set(int(round(i * step)) for i in range(max_frames))


def _composite_frame(hillshade_rgb, water_rgba):
    """Alpha-blend water RGBA over hillshade RGB background.

    Returns an (H, W, 3) uint8 RGB image.
    """
    alpha = water_rgba[:, :, 3:4].astype(np.float32) / 255.0
    water_rgb = water_rgba[:, :, :3].astype(np.float32)
    bg = hillshade_rgb.astype(np.float32)
    out = water_rgb * alpha + bg * (1.0 - alpha)
    return np.clip(out, 0, 255).astype(np.uint8)


def _get_font(size=20):
    """Get a TrueType font, falling back to default bitmap font."""
    if not _HAS_PIL:
        return None
    paths = [
        "/System/Library/Fonts/Helvetica.ttc",
        "/System/Library/Fonts/SFNSMono.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
    ]
    for p in paths:
        if os.path.isfile(p):
            try:
                return ImageFont.truetype(p, size)
            except Exception:
                continue
    try:
        return ImageFont.load_default()
    except Exception:
        return None


def _draw_text_on_image(pil_img, text, x, y, font, bg_color=(0, 0, 0, 160),
                         fg_color=(255, 255, 255, 255)):
    """Draw text with a semi-transparent background rectangle."""
    draw = ImageDraw.Draw(pil_img, "RGBA")
    bbox = draw.textbbox((x, y), text, font=font)
    pad = 6
    draw.rectangle(
        [bbox[0] - pad, bbox[1] - pad, bbox[2] + pad, bbox[3] + pad],
        fill=bg_color,
    )
    draw.text((x, y), text, fill=fg_color, font=font)


def _draw_colorbar(pil_img, clim, variable, width=200, height=18, margin=12):
    """Draw a horizontal colorbar in the bottom-right corner."""
    draw = ImageDraw.Draw(pil_img, "RGBA")
    img_w, img_h = pil_img.size
    font = _get_font(14)

    # Bar position
    x0 = img_w - width - margin
    y0 = img_h - height - margin - 20  # leave room for labels
    x1 = x0 + width
    y1 = y0 + height

    # Background panel
    draw.rectangle([x0 - 6, y0 - 22, x1 + 6, img_h - margin + 6],
                    fill=(0, 0, 0, 140))

    # Title
    var_label = {"h": "Water depth (m)", "velocity": "Velocity (m/s)",
                 "conc_SW": "Concentration"}.get(variable, variable)
    draw.text((x0, y0 - 18), var_label, fill=(255, 255, 255, 220), font=font)

    # Gradient bar (pixel by pixel along width)
    vmin, vmax = clim
    for px in range(width):
        t = px / max(width - 1, 1)
        val = vmin + t * (vmax - vmin)
        # Reuse same colour logic as _value_to_rgba
        tc = max(0.0, min(1.0, (val - vmin) / max(vmax - vmin, 1e-12)))
        r = int(200 * (1 - tc))
        g = int(220 * (1 - tc) + 50 * tc)
        b = int(255 * (1 - tc * 0.3))
        draw.line([(x0 + px, y0), (x0 + px, y1)], fill=(r, g, b, 220))

    # Min/max labels
    draw.text((x0, y1 + 2), f"{vmin:.2f}", fill=(255, 255, 255, 200),
              font=_get_font(12))
    draw.text((x1 - 30, y1 + 2), f"{vmax:.2f}", fill=(255, 255, 255, 200),
              font=_get_font(12))


def export_video(dem, meta, results, mesh_type, variable, clim,
                 output_path, fmt="gif", opacity=180, h_min=0.001,
                 draw_arrows=False, flow_config=None, scale=1,
                 fps=15, max_frames=200):
    """Export simulation results as an animated video (GIF or MP4).

    Parameters
    ----------
    fmt : str
        "gif" (Pillow, no external deps) or "mp4" (requires ffmpeg).
    scale : int
        Resolution multiplier (1 = native DEM resolution).
    fps : int
        Frames per second in the output video.
    max_frames : int
        Subsample if total frames exceed this.
    """
    if fmt == "gif" and not _HAS_PIL:
        print("ERROR: GIF export requires Pillow.  Install with: pip install Pillow")
        sys.exit(1)
    if fmt == "mp4" and not _check_ffmpeg():
        print("ERROR: MP4 export requires ffmpeg.  Install with:")
        print("  brew install ffmpeg    # macOS")
        print("  apt install ffmpeg     # Ubuntu/Debian")
        print("  Or use --export-video gif instead.")
        sys.exit(1)

    nrows, ncols = dem.shape
    cs = meta["cellsize"]

    # ── Hillshade background ──────────────────────────────────
    print("  Computing hillshade...")
    hs = _compute_hillshade(dem, cs)
    hs_up = _upscale_nearest(
        np.stack([hs, hs, hs, np.full_like(hs, 255)], axis=-1), scale)
    hillshade_rgb = hs_up[:, :, :3]  # (H*scale, W*scale, 3)
    out_h, out_w = hillshade_rgb.shape[:2]
    print(f"  Video resolution: {out_w}x{out_h} ({scale}x)")

    # ── Flow config defaults ──────────────────────────────────
    if flow_config is None:
        flow_config = {"MODE": "arrows", "ARROW_SPACING": 15, "ARROW_MAX_PX": 12,
                       "STREAMLINE_DENSITY": 15, "STREAMLINE_LENGTH": 20,
                       "PARTICLE_COUNT": 2000, "PARTICLE_TAIL": 8}

    vel_scale = None
    if draw_arrows:
        vel_scale = _auto_vel_scale(results, mesh_type, meta, h_min)
        mode = flow_config.get("MODE", "arrows")
        print(f"  Flow mode: {mode} (vel_scale={vel_scale:.4f} m/s)")

    # ── Subsample if too many frames ──────────────────────────
    n_total = len(results)
    keep = _subsample_indices(n_total, max_frames)
    print(f"  Processing {min(n_total, max_frames)} of {n_total} timesteps "
          f"(fps={fps})...")

    font = _get_font(max(16, int(out_h / 40)))

    # ── Generate and composite frames ─────────────────────────
    if fmt == "gif":
        pil_frames = []
        frame_count = 0
        for ri, t_sec, rgba, wet in _generate_frames(
                dem, meta, results, mesh_type, variable, clim,
                opacity, h_min, draw_arrows, flow_config, vel_scale, scale):

            if ri not in keep:
                continue

            # Composite over hillshade
            rgb = _composite_frame(hillshade_rgb, rgba)
            pil_img = Image.fromarray(rgb, "RGB")

            # Text overlay: timestamp
            _draw_text_on_image(pil_img, f"t = {format_time(t_sec)}",
                                12, 12, font)

            # Colorbar
            _draw_colorbar(pil_img, clim, variable,
                            width=min(200, out_w // 4), margin=12)

            # Quantize to 256 colours for GIF
            pil_frames.append(pil_img.quantize(colors=256,
                                                method=Image.Quantize.FASTOCTREE))

            frame_count += 1
            if frame_count % 20 == 0 or frame_count == 1:
                print(f"    Frame {frame_count}: t={format_time(t_sec)}, "
                      f"{wet} wet cells")

        if not pil_frames:
            print("  WARNING: No wet frames — video not created.")
            return

        print(f"  Writing GIF ({frame_count} frames, {fps} fps)...")
        duration_ms = int(1000 / fps)
        pil_frames[0].save(
            output_path,
            format="GIF",
            save_all=True,
            append_images=pil_frames[1:],
            duration=duration_ms,
            loop=0,
            optimize=False,
            disposal=2,
        )

        fsize = os.path.getsize(output_path)
        print(f"\n  Video saved: {output_path}")
        print(f"  Size: {fsize / 1024 / 1024:.1f} MB  |  "
              f"{frame_count} frames  |  {fps} fps  |  "
              f"{frame_count / fps:.1f}s duration")
        print(f"\n  Open with:  open {output_path}")

    elif fmt == "mp4":
        # Pipe raw RGB frames to ffmpeg
        cmd = [
            "ffmpeg", "-y",
            "-f", "rawvideo",
            "-pix_fmt", "rgb24",
            "-s", f"{out_w}x{out_h}",
            "-r", str(fps),
            "-i", "-",
            "-c:v", "libx264",
            "-pix_fmt", "yuv420p",
            "-crf", "20",
            "-preset", "medium",
            output_path,
        ]
        proc = subprocess.Popen(cmd, stdin=subprocess.PIPE,
                                 stdout=subprocess.DEVNULL,
                                 stderr=subprocess.DEVNULL)

        frame_count = 0
        for ri, t_sec, rgba, wet in _generate_frames(
                dem, meta, results, mesh_type, variable, clim,
                opacity, h_min, draw_arrows, flow_config, vel_scale, scale):

            if ri not in keep:
                continue

            rgb = _composite_frame(hillshade_rgb, rgba)
            pil_img = Image.fromarray(rgb, "RGB")
            _draw_text_on_image(pil_img, f"t = {format_time(t_sec)}",
                                12, 12, font)
            _draw_colorbar(pil_img, clim, variable,
                            width=min(200, out_w // 4), margin=12)

            proc.stdin.write(np.array(pil_img).tobytes())
            frame_count += 1
            if frame_count % 20 == 0 or frame_count == 1:
                print(f"    Frame {frame_count}: t={format_time(t_sec)}, "
                      f"{wet} wet cells")

        proc.stdin.close()
        proc.wait()

        if proc.returncode != 0:
            print("  ERROR: ffmpeg exited with error.")
            sys.exit(1)

        fsize = os.path.getsize(output_path)
        print(f"\n  Video saved: {output_path}")
        print(f"  Size: {fsize / 1024 / 1024:.1f} MB  |  "
              f"{frame_count} frames  |  {fps} fps  |  "
              f"{frame_count / fps:.1f}s duration")
        print(f"\n  Play with:  open {output_path}")


# ═══════════════════════════════════════════════════════════════════════════════
#  ADVECTED PARTICLE TRAIL VIDEO EXPORT  (earth.nullschool style)
# ═══════════════════════════════════════════════════════════════════════════════

def _draw_colorbar_advected(pil_img, vel_scale, width=200, height=18, margin=12):
    """Draw a speed colorbar (blue→cyan→white→yellow) for advected mode."""
    draw = ImageDraw.Draw(pil_img, "RGBA")
    img_w, img_h = pil_img.size
    font = _get_font(14)

    x0 = img_w - width - margin
    y0 = img_h - height - margin - 20
    x1 = x0 + width

    # Background panel
    draw.rectangle([x0 - 6, y0 - 22, x1 + 6, img_h - margin + 6],
                    fill=(0, 0, 0, 180))

    # Title
    draw.text((x0, y0 - 18), "Flow speed (m/s)",
              fill=(180, 210, 255, 220), font=font)

    # Gradient bar
    for px in range(width):
        t = px / max(width - 1, 1)
        speed = t * vel_scale
        r, g, b = _advected_speed_colormap_scalar(speed, vel_scale)
        draw.line([(x0 + px, y0), (x0 + px, y0 + height)],
                  fill=(r, g, b, 220))

    # Labels
    small_font = _get_font(12)
    draw.text((x0, y0 + height + 2), "0",
              fill=(180, 210, 255, 200), font=small_font)
    draw.text((x1 - 40, y0 + height + 2), f"{vel_scale:.2f}",
              fill=(180, 210, 255, 200), font=small_font)


def export_video_advected(dem, meta, results, mesh_type, variable, clim,
                          output_path, opacity=180, h_min=0.001,
                          scale=2, fps=30, max_frames=900,
                          n_particles=10000, fade_factor=0.96,
                          dark_factor=0.25):
    """Export earth.nullschool-style advected particle trail animation as MP4.

    Particles persist across frames, are advected through the velocity field
    via RK4 integration, and drawn onto a fade-texture canvas that produces
    smooth comet-like trails on a dark background.
    """
    if not _HAS_PIL:
        print("ERROR: Advected video export requires Pillow. "
              "Install with: pip install Pillow")
        sys.exit(1)
    if not _check_ffmpeg():
        print("ERROR: Advected video export requires ffmpeg. Install with:")
        print("  brew install ffmpeg    # macOS")
        print("  apt install ffmpeg     # Ubuntu/Debian")
        sys.exit(1)

    nrows, ncols = dem.shape
    cs = meta["cellsize"]

    # ── Hillshade background (dark) ──────────────────────────────
    print("  Computing dark hillshade background...")
    hs = _compute_hillshade(dem, cs)
    hs_rgba = np.stack([hs, hs, hs, np.full_like(hs, 255)], axis=-1)
    hs_up = _upscale_nearest(hs_rgba, scale)
    hillshade_rgb = hs_up[:, :, :3].astype(np.float32)
    # Darken and add a slight blue tint for cinematic look
    dark_hs = hillshade_rgb * dark_factor
    dark_hs[:, :, 2] = np.clip(dark_hs[:, :, 2] * 1.3, 0, 255)  # blue tint
    dark_hs = dark_hs.astype(np.uint8)
    out_h, out_w = dark_hs.shape[:2]
    print(f"  Video resolution: {out_w}x{out_h} ({scale}x)")

    # ── Velocity field cache ──────────────────────────────────────
    vel_cache = VelocityFieldCache(results, mesh_type, dem, meta, h_min, scale)
    vel_scale_val = _auto_vel_scale(results, mesh_type, meta, h_min)
    print(f"  Velocity scale: {vel_scale_val:.4f} m/s (95th percentile)")
    print(f"  Particles: {n_particles}, fade: {fade_factor}, dark: {dark_factor}")

    # ── Compute timing ────────────────────────────────────────────
    times = [r["time"] for r in results]
    sim_start_t = times[0]
    sim_end_t = times[-1]
    sim_duration = sim_end_t - sim_start_t
    if sim_duration <= 0:
        print("  ERROR: No time range in results.")
        return

    # Total video frames
    total_frames = min(max_frames, int(sim_duration / 10.0 * 3))  # ~3 frames per 10s step
    total_frames = max(total_frames, 60)  # at least 60 frames
    sim_dt_per_frame = sim_duration / total_frames  # seconds per video frame

    # Pixel size in metres
    pixel_size_m = cs / scale
    # dt_pixels: how many pixel-units a 1 m/s flow moves per video frame
    dt_pixels = sim_dt_per_frame / pixel_size_m
    # Sub-step if too large
    max_step = 1.5
    n_substeps = max(1, int(math.ceil(dt_pixels / max_step)))
    substep_dt = dt_pixels / n_substeps

    print(f"  Simulation: {format_time(sim_start_t)} → {format_time(sim_end_t)} "
          f"({len(results)} timesteps)")
    print(f"  Video: {total_frames} frames @ {fps} fps = "
          f"{total_frames / fps:.1f}s | {n_substeps} RK4 substeps/frame")

    # ── Initialize particle system ────────────────────────────────
    # Use wet mask from a mid-simulation timestep for initial spawn
    mid_idx = min(len(results) - 1, len(results) // 4)
    vx_init, vy_init, wet_init = vel_cache.get(mid_idx)
    particles = ParticleSystem(n_particles, wet_init, rng_seed=42)

    # ── Trail canvas ──────────────────────────────────────────────
    trail_canvas = np.zeros((out_h, out_w, 3), dtype=np.float32)

    # ── Open ffmpeg ───────────────────────────────────────────────
    cmd = [
        "ffmpeg", "-y",
        "-f", "rawvideo",
        "-pix_fmt", "rgb24",
        "-s", f"{out_w}x{out_h}",
        "-r", str(fps),
        "-i", "-",
        "-c:v", "libx264",
        "-pix_fmt", "yuv420p",
        "-crf", "18",
        "-preset", "medium",
        output_path,
    ]
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE,
                             stdout=subprocess.DEVNULL,
                             stderr=subprocess.DEVNULL)

    font = _get_font(max(16, int(out_h / 40)))

    print(f"  Rendering frames...")
    for frame_i in range(total_frames):
        sim_time = sim_start_t + frame_i * sim_dt_per_frame

        # (a) Get bracketing velocity fields
        idx_lo, idx_hi, frac = vel_cache.find_bracket(sim_time)
        vx_lo, vy_lo, wet_lo = vel_cache.get(idx_lo)
        vx_hi, vy_hi, wet_hi = vel_cache.get(idx_hi)

        # (b) Interpolate velocity field
        if idx_lo == idx_hi:
            vx_now, vy_now = vx_lo, vy_lo
            wet_now = wet_lo
        else:
            vx_now, vy_now = _lerp_velocity_fields(
                vx_lo, vy_lo, vx_hi, vy_hi, frac)
            wet_now = wet_lo | wet_hi

        # (c) Update particle wet mask
        particles.update_wet_mask(wet_now)

        # (d) Advect particles (sub-stepped RK4)
        for _ in range(n_substeps):
            particles.advect(vx_now, vy_now, substep_dt)

        # (e) Respawn dead particles
        particles.respawn_dead()

        # (f) Fade trail canvas
        trail_canvas *= fade_factor

        # (g) Draw particles onto trail canvas
        rows, cols, speeds, ages, max_ages = particles.get_draw_data()
        if len(rows) > 0:
            colors = _advected_speed_colormap(speeds, vel_scale_val)
            # Brightness: bright when young, fade near max_age
            age_frac = ages.astype(np.float64) / np.maximum(
                max_ages.astype(np.float64), 1.0)
            brightness = np.clip(1.0 - age_frac ** 1.5, 0.1, 1.0)
            _draw_particles_antialiased(
                trail_canvas, rows, cols, colors, brightness)

        # (h) Composite: dark hillshade + trail (additive)
        trail_rgb = np.clip(trail_canvas * 255.0, 0.0, 255.0)
        frame_rgb = np.clip(
            dark_hs.astype(np.float32) + trail_rgb, 0.0, 255.0
        ).astype(np.uint8)

        # (i) HUD overlays
        pil_img = Image.fromarray(frame_rgb, "RGB")
        _draw_text_on_image(pil_img, f"t = {format_time(sim_time)}",
                            12, 12, font)
        _draw_colorbar_advected(pil_img, vel_scale_val,
                                 width=min(200, out_w // 4), margin=12)

        # (j) Pipe to ffmpeg
        proc.stdin.write(np.array(pil_img).tobytes())

        if frame_i % 50 == 0 or frame_i == total_frames - 1:
            alive_count = int(np.sum(particles.alive))
            print(f"    Frame {frame_i + 1}/{total_frames}: "
                  f"t={format_time(sim_time)}, "
                  f"{alive_count} alive particles")

    proc.stdin.close()
    proc.wait()

    if proc.returncode != 0:
        print("  ERROR: ffmpeg exited with error.")
        sys.exit(1)

    fsize = os.path.getsize(output_path)
    print(f"\n  Video saved: {output_path}")
    print(f"  Size: {fsize / 1024 / 1024:.1f} MB  |  "
          f"{total_frames} frames  |  {fps} fps  |  "
          f"{total_frames / fps:.1f}s duration")
    print(f"\n  Play with:  open {output_path}")


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
    parser.add_argument(
        "--velocity-arrows", action="store_true", default=None,
        help="Overlay velocity direction arrows on the flood map.",
    )
    parser.add_argument(
        "--no-velocity-arrows", action="store_true", default=False,
        help="Disable velocity arrows (overrides JSON config).",
    )
    parser.add_argument(
        "--arrow-spacing", type=int, default=15,
        help="Spacing between velocity arrows in pixels. Default: 15",
    )
    parser.add_argument(
        "--arrow-max-px", type=int, default=12,
        help="Maximum arrow length in pixels. Default: 12",
    )
    parser.add_argument(
        "--modset", type=str, default=None,
        help="Path to modset JSON config file (reads VELOCITY_ARROWS section).",
    )
    parser.add_argument(
        "--flow-style", choices=["arrows", "streamlines", "particles", "advected"],
        default=None,
        help="Flow visualization style: 'arrows' (classic), 'streamlines' "
             "(curved flow lines), 'particles' (animated streaks), or "
             "'advected' (earth.nullschool-style persistent particle trails "
             "with fading tails — requires --export-video mp4). "
             "Default: read from JSON config, or 'streamlines'.",
    )
    parser.add_argument(
        "--streamline-density", type=int, default=None,
        help="Spacing between streamline seed points (pixels). Default: 15",
    )
    parser.add_argument(
        "--streamline-length", type=int, default=None,
        help="Max integration steps per streamline. Default: 20",
    )
    parser.add_argument(
        "--particle-count", type=int, default=None,
        help="Number of particles per frame. Default: 2000",
    )
    parser.add_argument(
        "--particle-tail", type=int, default=None,
        help="Particle streak length (integration steps). Default: 8",
    )
    parser.add_argument(
        "--scale", type=int, default=3, choices=[1, 2, 3, 4, 5],
        help="Supersampling factor for PNG output. Higher = sharper but larger "
             "files. scale=3 renders at 3x DEM resolution. Default: 3",
    )

    # ── Video export ──────────────────────────────────────────
    parser.add_argument(
        "--export-video", choices=["gif", "mp4"], default=None,
        help="Export as video instead of KML.  'gif' uses Pillow (no extra "
             "deps); 'mp4' requires ffmpeg.",
    )
    parser.add_argument(
        "--video-fps", type=int, default=15,
        help="Video frame rate. Default: 15",
    )
    parser.add_argument(
        "--video-scale", type=int, default=None, choices=[1, 2, 3, 4, 5],
        help="Resolution scale for video (default: 2). Separate from "
             "--scale which controls KML export.",
    )
    parser.add_argument(
        "--video-max-frames", type=int, default=200,
        help="Max frames in video; evenly subsamples if exceeded. Default: 200",
    )
    parser.add_argument(
        "--video-output", type=str, default=None,
        help="Output path for video file.  Default: fluxos_<mesh>.<ext>",
    )

    # ── Advected particle trail options ──────────────────────────
    parser.add_argument(
        "--advected-particles", type=int, default=10000,
        help="Number of persistent particles for advected mode. Default: 10000",
    )
    parser.add_argument(
        "--advected-fade", type=float, default=0.96,
        help="Trail fade factor per frame (0.90=short, 0.99=long). Default: 0.96",
    )
    parser.add_argument(
        "--advected-dark", type=float, default=0.25,
        help="Hillshade darkening factor (0=black, 1=full bright). Default: 0.25",
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

    # ── Flow visualisation config (JSON + CLI overrides) ────
    draw_arrows = False
    va_cfg = _read_velocity_arrows_config(args.modset)
    if va_cfg.get("STATUS", False):
        draw_arrows = True

    # Build flow_config from JSON defaults
    flow_config = {
        "MODE": va_cfg.get("MODE", "streamlines"),
        "ARROW_SPACING": va_cfg.get("ARROW_SPACING", args.arrow_spacing),
        "ARROW_MAX_PX": va_cfg.get("ARROW_MAX_PX", args.arrow_max_px),
        "STREAMLINE_DENSITY": va_cfg.get("STREAMLINE_DENSITY", 15),
        "STREAMLINE_LENGTH": va_cfg.get("STREAMLINE_LENGTH", 20),
        "PARTICLE_COUNT": va_cfg.get("PARTICLE_COUNT", 2000),
        "PARTICLE_TAIL": va_cfg.get("PARTICLE_TAIL", 8),
    }

    # CLI flags override JSON
    if args.velocity_arrows:
        draw_arrows = True
    if args.flow_style:
        draw_arrows = True
        flow_config["MODE"] = args.flow_style
    if args.no_velocity_arrows:
        draw_arrows = False
    if args.streamline_density is not None:
        flow_config["STREAMLINE_DENSITY"] = args.streamline_density
    if args.streamline_length is not None:
        flow_config["STREAMLINE_LENGTH"] = args.streamline_length
    if args.particle_count is not None:
        flow_config["PARTICLE_COUNT"] = args.particle_count
    if args.particle_tail is not None:
        flow_config["PARTICLE_TAIL"] = args.particle_tail

    # ── Step 4: Export ───────────────────────────────────────
    # Advected mode requires --export-video mp4
    if (flow_config and flow_config.get("MODE") == "advected"
            and not args.export_video):
        print("NOTE: Advected mode requires --export-video mp4. Adding it automatically.")
        args.export_video = "mp4"

    if args.export_video:
        vscale = args.video_scale if args.video_scale else 2
        ext = args.export_video

        # Advected mode: earth.nullschool-style particle trails
        is_advected = (flow_config and
                       flow_config.get("MODE") == "advected")
        if is_advected:
            voutput = (args.video_output or
                       f"fluxos_{args.mesh_type}_advected.mp4")
            suffix = f" @ {vscale}x"
            print(f"[4/4] Exporting advected particle trail MP4 "
                  f"({args.mesh_type} mesh{suffix})...")
            export_video_advected(
                dem, meta, results, args.mesh_type, args.variable, clim,
                voutput, opacity=args.opacity, h_min=args.h_min,
                scale=vscale, fps=args.video_fps,
                max_frames=args.video_max_frames,
                n_particles=args.advected_particles,
                fade_factor=args.advected_fade,
                dark_factor=args.advected_dark,
            )
        else:
            voutput = args.video_output or f"fluxos_{args.mesh_type}.{ext}"
            suffix = f" @ {vscale}x" if vscale > 1 else ""
            print(f"[4/4] Exporting {ext.upper()} video "
                  f"({args.mesh_type} mesh{suffix})...")
            export_video(
                dem, meta, results, args.mesh_type, args.variable, clim,
                voutput, fmt=ext, opacity=args.opacity, h_min=args.h_min,
                draw_arrows=draw_arrows, flow_config=flow_config,
                scale=vscale, fps=args.video_fps,
                max_frames=args.video_max_frames,
            )
    else:
        scale = args.scale
        suffix = f" @ {scale}x" if scale > 1 else ""
        print(f"[4/4] Exporting KML + PNG ({args.mesh_type} mesh{suffix})...")
        export_kml(
            dem, meta, results, args.mesh_type, args.variable, clim,
            utm_to_ll, sim_start, args.output_dir, args.opacity, args.h_min,
            draw_arrows=draw_arrows, flow_config=flow_config, scale=scale,
        )


def _read_velocity_arrows_config(modset_path=None):
    """Read VELOCITY_ARROWS section from modset JSON config."""
    import json
    paths = [modset_path] if modset_path else ["modset.json", "bin/modset.json"]
    for path in paths:
        if path and os.path.isfile(path):
            try:
                with open(path) as f:
                    cfg = json.load(f)
                return cfg.get("VELOCITY_ARROWS", {})
            except Exception:
                pass
    return {}


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


def _auto_vel_scale(results, mesh_type, meta, h_min):
    """Compute a global velocity scale (95th percentile) for consistent arrows."""
    all_vmag = []
    for r in results:
        if mesh_type == "regular":
            vx, vy = _extract_velocity_regular(r["data"], meta)
            h = r["data"][:, 5]
        else:
            vx, vy = _extract_velocity_triangular(r["cell_data"])
            h = r["cell_data"].get("h", np.array([]))
        wet = h > h_min
        if np.any(wet):
            vmag = np.sqrt(vx[wet] ** 2 + vy[wet] ** 2)
            valid = vmag[vmag > 1e-6]
            if len(valid) > 0:
                all_vmag.append(valid)
    if all_vmag:
        return float(np.percentile(np.concatenate(all_vmag), 95))
    return 1.0


if __name__ == "__main__":
    main()
