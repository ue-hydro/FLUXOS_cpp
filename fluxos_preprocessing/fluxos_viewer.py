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
            vx, vy, wet, _h, _conc = _build_velocity_grids_regular(
                self.dem, self.meta, r["data"], self.h_min)
        else:
            vx, vy, wet, _h, _conc = _build_velocity_grids_triangular(
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
    """Build (grid_vx, grid_vy, wet_mask, grid_h, grid_conc) from regular mesh data."""
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
    grid_h = np.zeros((nrows, ncols))
    grid_conc = np.zeros((nrows, ncols))
    wet_mask = np.zeros((nrows, ncols), dtype=bool)
    wet = h_depth > h_min
    grid_vx[row_idx[wet], col_idx[wet]] = vx[wet]
    grid_vy[row_idx[wet], col_idx[wet]] = vy[wet]
    grid_h[row_idx[wet], col_idx[wet]] = h_depth[wet]
    wet_mask[row_idx[wet], col_idx[wet]] = True

    # Extract concentration (column 11) if available
    if data.shape[1] > 11:
        conc = data[:, 11]
        conc_valid = np.where(conc >= 0, conc, 0.0)  # -1 = inactive
        grid_conc[row_idx[wet], col_idx[wet]] = conc_valid[wet]
    return grid_vx, grid_vy, wet_mask, grid_h, grid_conc


def _build_velocity_grids_triangular(points, cells, cell_data, h_min,
                                      nrows, ncols, xll, yll, cs):
    """Build (grid_vx, grid_vy, wet_mask, grid_h, grid_conc) from triangular mesh data."""
    if "h" not in cell_data:
        return (np.zeros((nrows, ncols)), np.zeros((nrows, ncols)),
                np.zeros((nrows, ncols), dtype=bool),
                np.zeros((nrows, ncols)),
                np.zeros((nrows, ncols)))

    h = cell_data["h"]
    vx, vy = _extract_velocity_triangular(cell_data)
    conc_data = cell_data.get("conc_SW", np.zeros_like(h))
    n_cells = len(cells)

    # Bin cells into pixel-sized blocks and average
    grid_vx = np.zeros((nrows, ncols))
    grid_vy = np.zeros((nrows, ncols))
    grid_h = np.zeros((nrows, ncols))
    grid_conc = np.zeros((nrows, ncols))
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
            grid_h[row, col] += h[ci]
            c = conc_data[ci] if ci < len(conc_data) else 0.0
            grid_conc[row, col] += max(0.0, c)
            grid_cnt[row, col] += 1

    valid = grid_cnt > 0
    grid_vx[valid] /= grid_cnt[valid]
    grid_vy[valid] /= grid_cnt[valid]
    grid_h[valid] /= grid_cnt[valid]
    grid_conc[valid] /= grid_cnt[valid]
    return grid_vx, grid_vy, valid, grid_h, grid_conc


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
                gvx, gvy, wmask, _h, _conc = _build_velocity_grids_regular(
                    dem, meta, r["data"], h_min)
            else:
                gvx, gvy, wmask, _h, _conc = _build_velocity_grids_triangular(
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
#  SATELLITE IMAGERY DOWNLOAD
# ═══════════════════════════════════════════════════════════════════════════════

def _download_satellite(sw_lat, sw_lon, ne_lat, ne_lon,
                        img_w, img_h, output_path,
                        basin_mask=None):
    """Download satellite imagery from Esri World Imagery and reproject to
    match a DEM grid.

    The DEM lives in UTM (locally equirectangular), so its pixel grid maps
    linearly to lat/lon over small areas.  Satellite tiles are in Web
    Mercator, where the Y-axis is nonlinear in latitude.  To align them we
    resample every output pixel individually: compute the lat/lon for that
    DEM cell, convert lat to Mercator-Y, and bilinearly sample the mosaic.

    If *basin_mask* is provided (a bool array, shape == (img_h, img_w),
    True = valid DEM cell), pixels outside the basin are darkened.
    """
    import urllib.request
    import io

    # Choose zoom level: aim for ~1 pixel per DEM cell
    lat_mid = (sw_lat + ne_lat) / 2.0
    lon_span = ne_lon - sw_lon
    meters_per_pixel_z0 = 156543.03392 * math.cos(math.radians(lat_mid))
    target_mpp = (lon_span * 111320 * math.cos(math.radians(lat_mid))) / img_w
    zoom = max(1, min(18, int(round(math.log2(meters_per_pixel_z0 / target_mpp)))))
    print(f"  Satellite: zoom level {zoom}")

    def lat_lon_to_tile(lat, lon, z):
        n = 2 ** z
        x = int((lon + 180.0) / 360.0 * n)
        lat_r = math.radians(lat)
        y = int((1.0 - math.log(math.tan(lat_r) + 1.0 / math.cos(lat_r))
                 / math.pi) / 2.0 * n)
        return max(0, min(n - 1, x)), max(0, min(n - 1, y))

    def tile_bounds(tx, ty, z):
        """Return (west, south, east, north) in degrees for tile (tx, ty)."""
        n = 2 ** z
        west = tx / n * 360.0 - 180.0
        east = (tx + 1) / n * 360.0 - 180.0
        north = math.degrees(math.atan(math.sinh(
            math.pi * (1 - 2 * ty / n))))
        south = math.degrees(math.atan(math.sinh(
            math.pi * (1 - 2 * (ty + 1) / n))))
        return west, south, east, north

    tx_min, ty_min = lat_lon_to_tile(ne_lat, sw_lon, zoom)  # NW corner
    tx_max, ty_max = lat_lon_to_tile(sw_lat, ne_lon, zoom)  # SE corner

    n_tiles_x = tx_max - tx_min + 1
    n_tiles_y = ty_max - ty_min + 1
    print(f"  Satellite: {n_tiles_x}x{n_tiles_y} = "
          f"{n_tiles_x * n_tiles_y} tiles")

    # Download all tiles into a big mosaic (256 px per tile)
    tile_px = 256
    mosaic_w = n_tiles_x * tile_px
    mosaic_h = n_tiles_y * tile_px
    mosaic = np.full((mosaic_h, mosaic_w, 3), 128, dtype=np.uint8)

    url_template = ("https://server.arcgisonline.com/ArcGIS/rest/services/"
                    "World_Imagery/MapServer/tile/{z}/{y}/{x}")

    for ty in range(ty_min, ty_max + 1):
        for tx in range(tx_min, tx_max + 1):
            url = url_template.format(z=zoom, y=ty, x=tx)
            try:
                req = urllib.request.Request(url, headers={
                    "User-Agent": "FLUXOS-viewer/1.0"})
                with urllib.request.urlopen(req, timeout=15) as resp:
                    tile_data = resp.read()
                try:
                    from PIL import Image as PILImage
                    tile_img = PILImage.open(io.BytesIO(tile_data))
                    tile_arr = np.array(tile_img.convert("RGB"))
                except ImportError:
                    continue
                oy = (ty - ty_min) * tile_px
                ox = (tx - tx_min) * tile_px
                th, tw = tile_arr.shape[:2]
                mosaic[oy:oy+th, ox:ox+tw] = tile_arr[:, :, :3]
            except Exception as e:
                print(f"    Tile {tx},{ty}: {e}")

    # ── Reproject: Mercator mosaic → equirectangular DEM grid ──
    # Compute the geographic bounds of the full tile mosaic
    w0, _, _, n0 = tile_bounds(tx_min, ty_min, zoom)
    _, s1, e1, _ = tile_bounds(tx_max, ty_max, zoom)
    mosaic_west, mosaic_east = w0, e1
    mosaic_north, mosaic_south = n0, s1

    def lat_to_merc_y(lat):
        return math.log(math.tan(math.pi / 4 + math.radians(lat) / 2))

    merc_n = lat_to_merc_y(mosaic_north)
    merc_s = lat_to_merc_y(mosaic_south)
    merc_range = merc_n - merc_s
    lon_range = mosaic_east - mosaic_west

    # Build output image by sampling the Mercator mosaic per-pixel.
    # For each output row r (0=north=top), compute the latitude, convert
    # to Mercator-Y, and find the mosaic row.  Longitude is linear.
    out = np.full((img_h, img_w, 3), 40, dtype=np.uint8)

    # Pre-compute mosaic column for each output column (lon is linear)
    col_lons = sw_lon + (np.arange(img_w) + 0.5) / img_w * (ne_lon - sw_lon)
    mosaic_cols = ((col_lons - mosaic_west) / lon_range * mosaic_w).astype(
        np.float64)

    for r in range(img_h):
        # DEM row 0 = north (top), row img_h-1 = south (bottom)
        frac = (r + 0.5) / img_h          # 0..1 from north to south
        lat = ne_lat - frac * (ne_lat - sw_lat)  # linear in lat (equirect)

        # Convert this latitude to a mosaic row (Mercator Y)
        merc_y = lat_to_merc_y(lat)
        mosaic_row_f = (merc_n - merc_y) / merc_range * mosaic_h

        # Bilinear Y interpolation
        ry0 = max(0, min(mosaic_h - 2, int(mosaic_row_f)))
        ry1 = ry0 + 1
        fy = mosaic_row_f - ry0

        for c in range(img_w):
            mx = mosaic_cols[c]
            cx0 = max(0, min(mosaic_w - 2, int(mx)))
            cx1 = cx0 + 1
            fx = mx - cx0

            # Bilinear sample
            p00 = mosaic[ry0, cx0].astype(np.float64)
            p10 = mosaic[ry0, cx1].astype(np.float64)
            p01 = mosaic[ry1, cx0].astype(np.float64)
            p11 = mosaic[ry1, cx1].astype(np.float64)
            val = (p00 * (1 - fx) * (1 - fy) + p10 * fx * (1 - fy) +
                   p01 * (1 - fx) * fy + p11 * fx * fy)
            out[r, c] = np.clip(val, 0, 255).astype(np.uint8)

    # ── Clip to basin (darken outside valid DEM cells) ─────────
    if basin_mask is not None:
        # Darken pixels outside the basin
        outside = ~basin_mask
        out[outside] = (out[outside].astype(np.float32) * 0.25).astype(
            np.uint8)

    # ── Save ───────────────────────────────────────────────────
    try:
        from PIL import Image as PILImage
        pil_img = PILImage.fromarray(out)
        pil_img.save(output_path, "JPEG", quality=85)
        fsize = os.path.getsize(output_path)
        print(f"  Satellite: {fsize/1024:.0f} KB "
              f"({img_w}x{img_h} px)")
    except ImportError:
        sat_rgba = np.stack([out[:, :, 0], out[:, :, 1], out[:, :, 2],
                             np.full(out.shape[:2], 255, np.uint8)], axis=-1)
        sat_bytes = _write_png(sat_rgba)
        out_png = output_path.replace(".jpg", ".png")
        with open(out_png, "wb") as f:
            f.write(sat_bytes)
        print(f"  Satellite (PNG fallback): {len(sat_bytes)/1024:.0f} KB")


# ═══════════════════════════════════════════════════════════════════════════════
#  INTERACTIVE WEBGL VIEWER EXPORT
# ═══════════════════════════════════════════════════════════════════════════════

def export_webgl(dem, meta, results, mesh_type, variable, clim,
                 output_dir="fluxos_web", h_min=0.001, step=5,
                 n_particles=65536, utm_to_ll=None):
    """Export FLUXOS results as an interactive WebGL particle viewer.

    Creates a directory with index.html + data files that can be served
    via any HTTP server for an earth.nullschool-style visualization.
    """
    import json as _json

    nrows, ncols = dem.shape
    cs = meta["cellsize"]

    os.makedirs(output_dir, exist_ok=True)
    data_dir = os.path.join(output_dir, "data")
    vel_dir = os.path.join(data_dir, "velocity")
    os.makedirs(vel_dir, exist_ok=True)

    # ── Subsample timesteps ───────────────────────────────────
    indices = list(range(0, len(results), step))
    n_frames = len(indices)
    print(f"  Subsampling: {n_frames} of {len(results)} timesteps "
          f"(every {step}th)")

    # ── Pass 1: compute global velocity range + max depth ─────
    print("  Computing global velocity range + max depth...")
    all_vx = []
    all_vy = []
    all_depths = []
    all_concs = []
    has_concentration = False
    for idx in indices:
        r = results[idx]
        if mesh_type == "regular":
            gvx, gvy, wet, grid_h, grid_conc = _build_velocity_grids_regular(
                dem, meta, r["data"], h_min)
        else:
            gvx, gvy, wet, grid_h, grid_conc = _build_velocity_grids_triangular(
                r["points"], r["cells"], r["cell_data"],
                h_min, nrows, ncols,
                meta["xllcorner"], meta["yllcorner"], cs)
        if np.any(wet):
            all_vx.append(gvx[wet])
            all_vy.append(gvy[wet])
            all_depths.append(grid_h[wet])
            conc_wet = grid_conc[wet]
            if np.any(conc_wet > 0):
                has_concentration = True
            all_concs.append(conc_wet)

    # Use 99.5th percentile for velocity range to clip outlier spikes
    if all_vx:
        cat_vx = np.concatenate(all_vx)
        cat_vy = np.concatenate(all_vy)
        gvx_min = float(np.percentile(cat_vx, 0.5))
        gvx_max = float(np.percentile(cat_vx, 99.5))
        gvy_min = float(np.percentile(cat_vy, 0.5))
        gvy_max = float(np.percentile(cat_vy, 99.5))
    else:
        gvx_min, gvx_max = -0.1, 0.1
        gvy_min, gvy_max = -0.1, 0.1

    # Ensure non-zero range
    if gvx_max - gvx_min < 1e-9:
        gvx_min, gvx_max = -0.1, 0.1
    if gvy_max - gvy_min < 1e-9:
        gvy_min, gvy_max = -0.1, 0.1

    # Use 99.9th percentile for h_max to avoid outlier artefacts
    if all_depths:
        all_h = np.concatenate(all_depths)
        gh_max = float(np.percentile(all_h, 99.9))
    else:
        gh_max = 1.0
    if gh_max < 1e-6:
        gh_max = 1.0

    # Concentration max: use 95th percentile of significant values
    # The ADE solver produces very skewed distributions, so filter to
    # values above 1% of the overall max to focus on the actual plume.
    conc_max = 0.0
    if has_concentration and all_concs:
        cat_conc = np.concatenate(all_concs)
        cat_conc_pos = cat_conc[cat_conc > 0]
        if len(cat_conc_pos) > 0:
            raw_max = float(cat_conc_pos.max())
            # Keep only values above 1% of max (the actual plume)
            threshold = raw_max * 0.01
            plume = cat_conc_pos[cat_conc_pos >= threshold]
            if len(plume) > 0:
                conc_max = float(np.percentile(plume, 95))
            else:
                conc_max = raw_max
    if conc_max < 1e-6:
        conc_max = 1.0  # fallback

    print(f"  Velocity range: vx=[{gvx_min:.4f}, {gvx_max:.4f}], "
          f"vy=[{gvy_min:.4f}, {gvy_max:.4f}] m/s")
    print(f"  Max water depth: {gh_max:.4f} m")
    if has_concentration:
        print(f"  Max concentration: {conc_max:.4f} mg/L")

    # ── Pass 2: export velocity + depth + concentration PNGs ──
    print(f"  Exporting {n_frames} velocity+depth PNGs...")
    if has_concentration:
        conc_dir = os.path.join(data_dir, "conc")
        os.makedirs(conc_dir, exist_ok=True)
        print(f"  Also exporting {n_frames} concentration PNGs...")
    times = []
    total_bytes = 0
    for fi, idx in enumerate(indices):
        r = results[idx]
        times.append(float(r["time"]))

        if mesh_type == "regular":
            gvx, gvy, wet, grid_h, grid_conc = _build_velocity_grids_regular(
                dem, meta, r["data"], h_min)
        else:
            gvx, gvy, wet, grid_h, grid_conc = _build_velocity_grids_triangular(
                r["points"], r["cells"], r["cell_data"],
                h_min, nrows, ncols,
                meta["xllcorner"], meta["yllcorner"], cs)

        # Encode vx → R, vy → G, depth → B, wet → A
        r_ch = ((gvx - gvx_min) / (gvx_max - gvx_min) * 255.0
                ).clip(0, 255).astype(np.uint8)
        g_ch = ((gvy - gvy_min) / (gvy_max - gvy_min) * 255.0
                ).clip(0, 255).astype(np.uint8)
        b_ch = (grid_h / gh_max * 255.0).clip(0, 255).astype(np.uint8)
        a_ch = np.where(wet, np.uint8(255), np.uint8(0))

        rgba = np.stack([r_ch, g_ch, b_ch, a_ch], axis=-1)
        png_bytes = _write_png(rgba)
        fname = f"v_{fi:04d}.png"
        with open(os.path.join(vel_dir, fname), "wb") as f:
            f.write(png_bytes)
        total_bytes += len(png_bytes)

        # Export concentration PNG: R=concentration, G=0, B=0, A=wet
        if has_concentration:
            c_r = (grid_conc / conc_max * 255.0).clip(0, 255).astype(np.uint8)
            c_g = np.zeros_like(c_r)
            c_b = np.zeros_like(c_r)
            c_a = np.where(wet & (grid_conc > 0), np.uint8(255), np.uint8(0))
            c_rgba = np.stack([c_r, c_g, c_b, c_a], axis=-1)
            c_bytes = _write_png(c_rgba)
            with open(os.path.join(conc_dir, f"c_{fi:04d}.png"), "wb") as f:
                f.write(c_bytes)
            total_bytes += len(c_bytes)

        if fi % 50 == 0 or fi == n_frames - 1:
            print(f"    [{fi+1}/{n_frames}] t={format_time(r['time'])}, "
                  f"PNG {len(png_bytes)/1024:.0f} KB")

    print(f"  Total data: {total_bytes / 1024 / 1024:.1f} MB")

    # ── Export hillshade PNG ──────────────────────────────────
    print("  Exporting hillshade...")
    hs = _compute_hillshade(dem, cs)
    hs_rgba = np.stack([hs, hs, hs, np.full_like(hs, 255)], axis=-1)
    hs_bytes = _write_png(hs_rgba)
    with open(os.path.join(data_dir, "hillshade.png"), "wb") as f:
        f.write(hs_bytes)
    print(f"  Hillshade: {len(hs_bytes)/1024:.0f} KB")

    # ── Export heightmap PNG (16-bit encoded) ─────────────────
    print("  Exporting heightmap...")
    nan_mask = np.isnan(dem)
    valid_dem = dem.copy()
    valid_dem[nan_mask] = 0.0  # fill NaN for encoding

    z_min = float(np.nanmin(dem)) if not np.all(nan_mask) else 0.0
    z_max = float(np.nanmax(dem)) if not np.all(nan_mask) else 1.0
    if z_max - z_min < 0.01:
        z_max = z_min + 1.0  # avoid division by zero

    z_norm = np.clip((valid_dem - z_min) / (z_max - z_min), 0.0, 1.0)
    z_16bit = (z_norm * 65535.0).astype(np.uint16)
    r_ch = (z_16bit >> 8).astype(np.uint8)      # high byte
    g_ch = (z_16bit & 0xFF).astype(np.uint8)     # low byte
    b_ch = np.where(nan_mask, np.uint8(0), np.uint8(255))  # validity mask
    a_ch = np.full_like(r_ch, 255)               # alpha = 255 always

    hm_rgba = np.stack([r_ch, g_ch, b_ch, a_ch], axis=-1)
    hm_bytes = _write_png(hm_rgba)
    with open(os.path.join(data_dir, "heightmap.png"), "wb") as f:
        f.write(hm_bytes)
    print(f"  Heightmap: {len(hm_bytes)/1024:.0f} KB "
          f"(z_min={z_min:.2f}, z_max={z_max:.2f})")

    # ── Download satellite image ────────────────────────────
    bbox_ll = None
    if utm_to_ll is not None:
        xll = meta["xllcorner"]
        yll = meta["yllcorner"]
        sw_lon, sw_lat = utm_to_ll(xll, yll)
        ne_lon, ne_lat = utm_to_ll(xll + ncols * cs, yll + nrows * cs)
        bbox_ll = {
            "sw_lat": float(sw_lat), "sw_lon": float(sw_lon),
            "ne_lat": float(ne_lat), "ne_lon": float(ne_lon),
        }
        print(f"  Bounding box: ({sw_lat:.6f},{sw_lon:.6f}) - "
              f"({ne_lat:.6f},{ne_lon:.6f})")

        # Download satellite tiles, reproject, and clip to basin
        sat_path = os.path.join(data_dir, "satellite.jpg")
        basin_mask = ~np.isnan(dem)   # True where DEM is valid
        _download_satellite(sw_lat, sw_lon, ne_lat, ne_lon,
                            ncols, nrows, sat_path,
                            basin_mask=basin_mask)

    # ── Export metadata JSON ──────────────────────────────────
    metadata = {
        "width": int(ncols),
        "height": int(nrows),
        "n_timesteps": n_frames,
        "times": times,
        "vx_min": float(gvx_min),
        "vx_max": float(gvx_max),
        "vy_min": float(gvy_min),
        "vy_max": float(gvy_max),
        "cellsize": float(cs),
        "default_particles": n_particles,
        "z_min": z_min,
        "z_max": z_max,
        "h_max": float(gh_max),
        "height_exaggeration": 0.1,
    }
    if has_concentration:
        metadata["has_concentration"] = True
        metadata["conc_max"] = float(conc_max)
    if bbox_ll:
        metadata["bbox"] = bbox_ll
        metadata["has_satellite"] = True
    meta_path = os.path.join(data_dir, "metadata.json")
    with open(meta_path, "w") as f:
        _json.dump(metadata, f, indent=2)

    # ── Generate HTML viewer ──────────────────────────────────
    print("  Generating WebGL viewer...")
    html = _generate_webgl_html()
    html_path = os.path.join(output_dir, "index.html")
    with open(html_path, "w") as f:
        f.write(html)

    print(f"\n  WebGL viewer exported to: {output_dir}/")
    print(f"  {n_frames} velocity frames + hillshade + metadata + viewer")
    print(f"\n  To view, start a local server:")
    print(f"    python3 -m http.server 8080 -d {output_dir}")
    print(f"  Then open:  http://localhost:8080")


def _generate_webgl_html():
    """Generate the complete self-contained 3D WebGL particle viewer HTML.

    Reads the HTML from the canonical index.html template in the fluxos_web
    directory (sibling to fluxos_preprocessing).  Falls back to the embedded
    copy below if the file doesn't exist (e.g. standalone script usage).
    """
    import pathlib as _pathlib
    # Try canonical location first
    _html_path = _pathlib.Path(__file__).parent.parent / 'fluxos_web' / 'index.html'
    if _html_path.exists():
        return _html_path.read_text()
    # Fallback: embedded (may be out-of-date)
    return '''\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>FLUXOS Flood Simulation — Interactive Flow Viewer</title>
<style>
* { margin: 0; padding: 0; box-sizing: border-box; }
html, body { width: 100%; height: 100%; overflow: hidden; background: #0a0a12; }
canvas#gl { width: 100%; height: 100%; display: block; }
#controls {
    position: fixed; bottom: 16px; left: 16px;
    background: rgba(5,8,20,0.82); color: #b0c4de;
    padding: 16px 20px; border-radius: 10px;
    font-family: 'SF Mono', 'Fira Code', 'Consolas', monospace;
    font-size: 12px; min-width: 300px;
    border: 1px solid rgba(60,90,140,0.3);
    backdrop-filter: blur(8px); -webkit-backdrop-filter: blur(8px);
    user-select: none;
}
#controls h2 { font-size: 15px; color: #e0ecff; margin-bottom: 12px;
    font-weight: 600; letter-spacing: 0.5px; }
.ctrl-row { margin-bottom: 8px; display: flex; align-items: center;
    justify-content: space-between; }
.ctrl-row label { flex: 0 0 100px; color: #8899bb; }
.ctrl-row input[type=range] { flex: 1; margin: 0 8px; }
.ctrl-row .val { flex: 0 0 60px; text-align: right; color: #ccddef;
    font-variant-numeric: tabular-nums; }
.btn { background: rgba(40,60,100,0.5); color: #b0c4de; border: 1px solid
    rgba(80,120,180,0.4); border-radius: 5px; padding: 6px 14px;
    cursor: pointer; font-family: inherit; font-size: 12px;
    transition: background 0.15s; }
.btn:hover { background: rgba(60,90,140,0.6); }
.btn.active { background: rgba(80,140,220,0.4); border-color: rgba(100,160,240,0.6); }
#info { position: fixed; top: 16px; right: 16px;
    background: rgba(5,8,20,0.7); color: #8899bb; padding: 10px 14px;
    border-radius: 8px; font-family: 'SF Mono', monospace; font-size: 11px;
    border: 1px solid rgba(60,90,140,0.2); pointer-events: none; }
#info div { margin-bottom: 3px; }
#info .hl { color: #ccddef; }
#colorbar-wrap { position: fixed; bottom: 20px; right: 16px;
    background: rgba(5,8,20,0.7); padding: 8px 12px 6px;
    border-radius: 8px; border: 1px solid rgba(60,90,140,0.2); }
#colorbar-wrap .cb-title { font-family: monospace; font-size: 11px;
    color: #8899bb; margin-bottom: 4px; }
#colorbar-wrap .cb-labels { display: flex; justify-content: space-between;
    font-family: monospace; font-size: 10px; color: #8899bb; margin-top: 2px; }
canvas#colorbar { display: block; border-radius: 3px; }
#loading { position: fixed; top: 50%; left: 50%; transform: translate(-50%,-50%);
    color: #5577aa; font-family: monospace; font-size: 16px; }
input[type=range] { -webkit-appearance: none; height: 4px;
    background: rgba(60,90,140,0.4); border-radius: 2px; outline: none; }
input[type=range]::-webkit-slider-thumb { -webkit-appearance: none;
    width: 14px; height: 14px; border-radius: 50%;
    background: #5588cc; cursor: pointer; border: 2px solid #334466; }
</style>
</head>
<body>
<canvas id="gl"></canvas>

<div id="loading">Loading simulation data...</div>

<div id="controls" style="display:none;">
    <h2>FLUXOS Flood Simulation</h2>
    <div class="ctrl-row">
        <label>Time</label>
        <input type="range" id="sl-time" min="0" max="100" value="0">
        <span class="val" id="v-time">0s</span>
    </div>
    <div class="ctrl-row" style="justify-content:center; gap:8px; margin:10px 0;">
        <button class="btn active" id="btn-play">&#9654; Play</button>
        <button class="btn" id="btn-reset">&#8634; Reset</button>
        <button class="btn" id="btn-fs">&#x26F6; Fullscreen</button>
    </div>
    <div class="ctrl-row">
        <label>Speed</label>
        <input type="range" id="sl-speed" min="5" max="300" value="80">
        <span class="val" id="v-speed">0.8x</span>
    </div>
    <div class="ctrl-row">
        <label>Particles</label>
        <input type="range" id="sl-particles" min="10" max="18" value="14">
        <span class="val" id="v-particles">65K</span>
    </div>
    <div class="ctrl-row">
        <label>Trail length</label>
        <input type="range" id="sl-trail" min="880" max="998" value="970">
        <span class="val" id="v-trail">0.96</span>
    </div>
</div>

<div id="info" style="display:none;">
    <div>t = <span class="hl" id="i-time">0s</span></div>
    <div><span class="hl" id="i-fps">--</span> fps</div>
    <div><span class="hl" id="i-particles">--</span> particles</div>
</div>

<div id="colorbar-wrap" style="display:none;">
    <div class="cb-title">Flow speed (m/s)</div>
    <canvas id="colorbar" width="180" height="14"></canvas>
    <div class="cb-labels"><span id="cb-min">0</span><span id="cb-max">--</span></div>
</div>

<!-- ═══════ GLSL SHADERS ═══════ -->

<script id="quad-vs" type="x-shader/x-vertex">
precision mediump float;
attribute vec2 a_pos;
varying vec2 v_tex;
void main() {
    v_tex = a_pos;
    gl_Position = vec4(1.0 - 2.0 * a_pos, 0.0, 1.0);
}
</script>

<script id="update-fs" type="x-shader/x-fragment">
precision highp float;
uniform sampler2D u_particles;
uniform sampler2D u_velocity;
uniform sampler2D u_velocity_next;
uniform sampler2D u_spawn_tex;    // 1D texture (N×1): RG=x coord, BA=y coord of wet cells
uniform float u_spawn_count;      // number of entries in spawn texture
uniform float u_time_frac;
uniform float u_rand_seed;
uniform float u_speed_factor;
uniform float u_drop_rate;
uniform float u_drop_rate_bump;
uniform vec2 u_vel_range_x;
uniform vec2 u_vel_range_y;
uniform vec2 u_grid_res;
varying vec2 v_tex;

float rand(vec2 co) {
    return fract(sin(dot(co, vec2(12.9898, 78.233))) * 43758.5453);
}

vec2 lookup_vel(sampler2D tex, vec2 uv) {
    vec4 t = texture2D(tex, uv);
    float vx = mix(u_vel_range_x.x, u_vel_range_x.y, t.r);
    float vy = mix(u_vel_range_y.x, u_vel_range_y.y, t.g);
    float wet = step(0.5, t.a);
    return vec2(vx, vy) * wet;
}

// Decode a spawn position from the spawn lookup texture
vec2 getSpawnPos(float idx) {
    float u = (idx + 0.5) / u_spawn_count;
    vec4 s = texture2D(u_spawn_tex, vec2(u, 0.5));
    return vec2(s.r + s.b / 255.0, s.g + s.a / 255.0);
}

void main() {
    vec4 enc = texture2D(u_particles, v_tex);
    vec2 pos = vec2(enc.r / 255.0 + enc.b, enc.g / 255.0 + enc.a);

    // Check current wetness
    vec4 curVel = texture2D(u_velocity, pos);
    bool curWet = curVel.a > 0.5;

    vec2 v0 = lookup_vel(u_velocity, pos);
    vec2 v1 = lookup_vel(u_velocity_next, pos);
    vec2 vel = mix(v0, v1, u_time_frac);
    float spd = length(vel);

    // Convert velocity (m/s) to texture-coordinate offset
    vec2 off = vec2(
        vel.x * u_speed_factor / u_grid_res.x,
       -vel.y * u_speed_factor / u_grid_res.y
    );
    vec2 np2 = pos + off;

    bool oob = np2.x < 0.0 || np2.x > 1.0 || np2.y < 0.0 || np2.y > 1.0;
    bool destDry = texture2D(u_velocity, np2).a < 0.5;

    float s1 = rand(pos + u_rand_seed);
    float drop = u_drop_rate + spd * u_drop_rate_bump;

    bool needRespawn = oob || !curWet || destDry || s1 > (1.0 - drop);

    if (needRespawn) {
        // Sample a random wet cell from the spawn lookup texture
        float ridx = floor(rand(v_tex + u_rand_seed + 1.3) * u_spawn_count);
        np2 = getSpawnPos(ridx);
    }

    // encode position as 16-bit per axis
    vec2 f = fract(np2 * 255.0);
    vec2 i = floor(np2 * 255.0) / 255.0;
    gl_FragColor = vec4(f.x, f.y, i.x, i.y);
}
</script>

<script id="draw-vs" type="x-shader/x-vertex">
precision mediump float;
attribute float a_index;
uniform sampler2D u_particles;
uniform float u_particles_res;
uniform sampler2D u_velocity;
uniform sampler2D u_velocity_next;
uniform float u_time_frac;
uniform vec2 u_vel_range_x;
uniform vec2 u_vel_range_y;
uniform float u_point_size;
varying float v_speed_t;

void main() {
    float col = mod(a_index, u_particles_res);
    float row = floor(a_index / u_particles_res);
    vec2 uv = (vec2(col, row) + 0.5) / u_particles_res;
    vec4 enc = texture2D(u_particles, uv);
    vec2 pos = vec2(enc.r / 255.0 + enc.b, enc.g / 255.0 + enc.a);

    vec4 t0 = texture2D(u_velocity, pos);
    vec4 t1 = texture2D(u_velocity_next, pos);
    vec4 vt = mix(t0, t1, u_time_frac);

    float vx = mix(u_vel_range_x.x, u_vel_range_x.y, vt.r);
    float vy = mix(u_vel_range_y.x, u_vel_range_y.y, vt.g);
    float spd = length(vec2(vx, vy));
    float mx = max(abs(u_vel_range_x.y), abs(u_vel_range_y.y));
    v_speed_t = clamp(spd / mx, 0.0, 1.0);

    float wet = step(0.5, vt.a);
    gl_PointSize = u_point_size * wet;
    gl_Position = vec4(2.0 * pos.x - 1.0, 1.0 - 2.0 * pos.y, 0.0, 1.0);
}
</script>

<script id="draw-fs" type="x-shader/x-fragment">
precision mediump float;
uniform sampler2D u_color_ramp;
varying float v_speed_t;
void main() {
    gl_FragColor = texture2D(u_color_ramp, vec2(v_speed_t, 0.5));
}
</script>

<script id="screen-fs" type="x-shader/x-fragment">
precision mediump float;
uniform sampler2D u_screen;
uniform float u_opacity;
varying vec2 v_tex;
void main() {
    vec4 c = texture2D(u_screen, 1.0 - v_tex);
    gl_FragColor = vec4(floor(255.0 * c * u_opacity) / 255.0);
}
</script>

<script id="composite-fs" type="x-shader/x-fragment">
precision mediump float;
uniform sampler2D u_hillshade;
uniform sampler2D u_trails;
uniform float u_dark;
varying vec2 v_tex;
void main() {
    vec3 hs = texture2D(u_hillshade, v_tex).rgb * u_dark;
    hs.b = min(1.0, hs.b * 1.3);
    vec3 tr = texture2D(u_trails, v_tex).rgb;
    gl_FragColor = vec4(min(vec3(1.0), hs + tr), 1.0);
}
</script>

<!-- ═══════ JAVASCRIPT ═══════ -->
<script>
'use strict';

// ── Color ramp: deep navy → blue → cyan → white → yellow ──
function buildColorRamp() {
    const stops = [
        [0.00, 15, 20, 80], [0.15, 30, 70, 170], [0.40, 60, 190, 255],
        [0.70, 200, 255, 255], [0.90, 255, 255, 180], [1.00, 255, 230, 80]
    ];
    const data = new Uint8Array(256 * 4);
    for (let p = 0; p < 256; p++) {
        const t = p / 255;
        let si = 0;
        for (let i = 0; i < stops.length - 1; i++) {
            if (t >= stops[i][0] && t <= stops[i+1][0]) { si = i; break; }
            if (i === stops.length - 2) si = i;
        }
        const f = (t - stops[si][0]) / Math.max(stops[si+1][0] - stops[si][0], 1e-6);
        data[p*4+0] = Math.round(stops[si][1] + (stops[si+1][1] - stops[si][1]) * f);
        data[p*4+1] = Math.round(stops[si][2] + (stops[si+1][2] - stops[si][2]) * f);
        data[p*4+2] = Math.round(stops[si][3] + (stops[si+1][3] - stops[si][3]) * f);
        data[p*4+3] = 255;
    }
    return data;
}

function drawColorbar() {
    const c = document.getElementById('colorbar');
    const ctx = c.getContext('2d');
    const ramp = buildColorRamp();
    for (let x = 0; x < c.width; x++) {
        const i = Math.floor(x / c.width * 255) * 4;
        ctx.fillStyle = `rgb(${ramp[i]},${ramp[i+1]},${ramp[i+2]})`;
        ctx.fillRect(x, 0, 1, c.height);
    }
}

// ── WebGL helpers ──
function createShader(gl, type, src) {
    const s = gl.createShader(type);
    gl.shaderSource(s, src);
    gl.compileShader(s);
    if (!gl.getShaderParameter(s, gl.COMPILE_STATUS))
        throw new Error('Shader: ' + gl.getShaderInfoLog(s));
    return s;
}
function createProgram(gl, vs, fs) {
    const p = gl.createProgram();
    gl.attachShader(p, vs);
    gl.attachShader(p, fs);
    gl.linkProgram(p);
    if (!gl.getProgramParameter(p, gl.LINK_STATUS))
        throw new Error('Program: ' + gl.getProgramInfoLog(p));
    return p;
}
function createTex(gl, filter, data, w, h) {
    const t = gl.createTexture();
    gl.bindTexture(gl.TEXTURE_2D, t);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, filter);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, filter);
    if (data instanceof HTMLImageElement || data instanceof ImageData) {
        // Critical: prevent browser from premultiplying alpha, which would corrupt velocity data
        gl.pixelStorei(gl.UNPACK_PREMULTIPLY_ALPHA_WEBGL, false);
        gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, false);
        gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, data);
    } else {
        gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, w, h, 0, gl.RGBA, gl.UNSIGNED_BYTE, data);
    }
    return t;
}
function createFBO(gl, tex) {
    const fb = gl.createFramebuffer();
    gl.bindFramebuffer(gl.FRAMEBUFFER, fb);
    gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, tex, 0);
    return fb;
}
function bindAttr(gl, buf, attr, size) {
    gl.bindBuffer(gl.ARRAY_BUFFER, buf);
    gl.enableVertexAttribArray(attr);
    gl.vertexAttribPointer(attr, size, gl.FLOAT, false, 0, 0);
}

// ── Viewer class ──
class FloodViewer {
    constructor(canvas) {
        this.gl = canvas.getContext('webgl', { antialias: false, premultipliedAlpha: false, alpha: false });
        if (!this.gl) throw new Error('WebGL not supported');
        this.canvas = canvas;
        this.meta = null;
        this.currentFrame = 0;
        this.timeFrac = 0;
        this.playing = true;
        this.speed = 0.8;
        this.fadeOpacity = 0.97;
        this.dropRate = 0.003;
        this.dropRateBump = 0.01;
        this.darkFactor = 0.25;
        this.pointSize = 1.5;
        this.particleRes = 128;
        this.numParticles = 128 * 128;  // 16K particles — good for sparse flood domains
        this.frameCache = new Map();
        this.frameCacheMax = 12;
        this.speedFactor = 0.5;
        this._lastT = 0;
        this._fpsFrames = 0;
        this._fpsTime = 0;
        this._fps = 0;
    }

    async init() {
        const resp = await fetch('data/metadata.json');
        this.meta = await resp.json();
        const m = this.meta;

        document.getElementById('sl-time').max = m.n_timesteps - 1;
        document.getElementById('cb-max').textContent =
            Math.max(Math.abs(m.vx_max), Math.abs(m.vy_max)).toFixed(2);

        const simDt = m.times.length > 1 ? m.times[1] - m.times[0] : 10;
        // speedFactor: velocity (m/s) × speedFactor / grid_res → texture offset per frame
        // = dt_sec_per_render_frame / cellsize_m
        // (shader divides by grid_res for each axis separately)
        const dtPerRenderFrame = simDt / 50.0;  // ~50 render frames per data timestep
        this.speedFactor = dtPerRenderFrame / m.cellsize;  // cells per frame for 1 m/s

        const gl = this.gl;
        this.compileShaders(gl);
        this.createBuffers(gl);
        this.initParticleState(gl);
        this.initScreenTextures(gl);

        // Color ramp texture
        const rampData = buildColorRamp();
        this.colorRampTex = createTex(gl, gl.LINEAR, rampData, 256, 1);

        // Load hillshade
        this.hillshadeTex = await this.loadImage(gl, 'data/hillshade.png', gl.LINEAR);

        // Find first frame with water and load it
        let startFrame = 0;
        for (let i = 0; i < Math.min(m.n_timesteps, 20); i++) {
            const f = await this.loadVelFrame(i);
            console.log(`Frame ${i}: ${f.spawnCount} wet spawn cells`);
            if (f.spawnCount > 10) { startFrame = i; break; }
        }
        this.currentFrame = startFrame;
        const nextF = Math.min(startFrame + 1, m.n_timesteps - 1);
        if (!this.frameCache.has(nextF)) await this.loadVelFrame(nextF);
        console.log(`Starting at frame ${startFrame}, speedFactor: ${this.speedFactor.toFixed(4)}, simDt: ${simDt}s`);

        // Re-initialize particles at wet cell positions from the start frame
        const startEntry = this.frameCache.get(startFrame);
        if (startEntry && startEntry.spawnCount > 1) {
            this.initParticlesFromSpawn(gl, startEntry);
        }

        document.getElementById('loading').style.display = 'none';
        document.getElementById('controls').style.display = '';
        document.getElementById('info').style.display = '';
        document.getElementById('colorbar-wrap').style.display = '';
        drawColorbar();

        this._lastT = performance.now();
        this.render();
    }

    compileShaders(gl) {
        const quadVS = createShader(gl, gl.VERTEX_SHADER,
            document.getElementById('quad-vs').textContent);
        const drawVS = createShader(gl, gl.VERTEX_SHADER,
            document.getElementById('draw-vs').textContent);

        const mkFrag = id => createShader(gl, gl.FRAGMENT_SHADER,
            document.getElementById(id).textContent);

        this.updateProg = createProgram(gl, quadVS, mkFrag('update-fs'));
        this.drawProg = createProgram(gl, drawVS, mkFrag('draw-fs'));
        this.screenProg = createProgram(gl, quadVS, mkFrag('screen-fs'));
        this.compositeProg = createProgram(gl, quadVS, mkFrag('composite-fs'));
    }

    createBuffers(gl) {
        // Fullscreen quad
        this.quadBuf = gl.createBuffer();
        gl.bindBuffer(gl.ARRAY_BUFFER, this.quadBuf);
        gl.bufferData(gl.ARRAY_BUFFER,
            new Float32Array([0,0, 1,0, 0,1, 0,1, 1,0, 1,1]), gl.STATIC_DRAW);

        // Particle index buffer
        this.indexBuf = gl.createBuffer();
        const idx = new Float32Array(this.numParticles);
        for (let i = 0; i < this.numParticles; i++) idx[i] = i;
        gl.bindBuffer(gl.ARRAY_BUFFER, this.indexBuf);
        gl.bufferData(gl.ARRAY_BUFFER, idx, gl.STATIC_DRAW);
    }

    initParticleState(gl) {
        const n = this.particleRes * this.particleRes * 4;
        const data = new Uint8Array(n);
        for (let i = 0; i < this.particleRes * this.particleRes; i++) {
            const x = Math.random(), y = Math.random();
            data[i*4+0] = Math.floor(256 * (x * 255 - Math.floor(x * 255)));
            data[i*4+1] = Math.floor(256 * (y * 255 - Math.floor(y * 255)));
            data[i*4+2] = Math.floor(x * 255);
            data[i*4+3] = Math.floor(y * 255);
        }
        this.particleTex0 = createTex(gl, gl.NEAREST, data,
            this.particleRes, this.particleRes);
        this.particleTex1 = createTex(gl, gl.NEAREST,
            new Uint8Array(n), this.particleRes, this.particleRes);
        this.particleFB0 = createFBO(gl, this.particleTex0);
        this.particleFB1 = createFBO(gl, this.particleTex1);
    }

    // Initialize particles at wet cell positions from a frame entry
    initParticlesFromSpawn(gl, frameEntry) {
        const wp = frameEntry.wetPositions;
        const nWet = Math.floor(wp.length / 2);
        if (nWet < 1) return;

        const n = this.particleRes * this.particleRes;
        const data = new Uint8Array(n * 4);
        for (let i = 0; i < n; i++) {
            const ri = Math.floor(Math.random() * nWet);
            const x = wp[ri * 2];
            const y = wp[ri * 2 + 1];
            data[i*4+0] = Math.floor(256 * (x * 255 - Math.floor(x * 255)));
            data[i*4+1] = Math.floor(256 * (y * 255 - Math.floor(y * 255)));
            data[i*4+2] = Math.floor(x * 255);
            data[i*4+3] = Math.floor(y * 255);
        }
        gl.bindTexture(gl.TEXTURE_2D, this.particleTex0);
        gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, this.particleRes, this.particleRes,
            0, gl.RGBA, gl.UNSIGNED_BYTE, data);
        console.log(`Initialized ${n} particles at ${nWet} wet positions`);
    }

    initScreenTextures(gl) {
        const w = gl.canvas.width, h = gl.canvas.height;
        const empty = new Uint8Array(w * h * 4);
        if (this.screenTex0) gl.deleteTexture(this.screenTex0);
        if (this.screenTex1) gl.deleteTexture(this.screenTex1);
        if (this.screenFB0) gl.deleteFramebuffer(this.screenFB0);
        if (this.screenFB1) gl.deleteFramebuffer(this.screenFB1);
        this.screenTex0 = createTex(gl, gl.NEAREST, empty, w, h);
        this.screenTex1 = createTex(gl, gl.NEAREST, empty, w, h);
        this.screenFB0 = createFBO(gl, this.screenTex0);
        this.screenFB1 = createFBO(gl, this.screenTex1);
    }

    async loadImage(gl, url, filter) {
        return new Promise((res, rej) => {
            const img = new Image();
            img.onload = () => {
                const t = createTex(gl, filter, img, 0, 0);
                res(t);
            };
            img.onerror = rej;
            img.src = url;
        });
    }

    // Load a velocity frame PNG and extract wet cell positions for spawn texture
    async loadVelFrame(idx) {
        if (this.frameCache.has(idx)) return this.frameCache.get(idx);
        const gl = this.gl;
        const m = this.meta;

        // Load as Image
        const img = await new Promise((res, rej) => {
            const im = new Image();
            im.onload = () => res(im);
            im.onerror = rej;
            im.src = `data/velocity/v_${String(idx).padStart(4,'0')}.png`;
        });

        // Create velocity WebGL texture
        const tex = createTex(gl, gl.LINEAR, img, 0, 0);

        // Extract wet cell positions by drawing to an offscreen canvas
        const oc = document.createElement('canvas');
        oc.width = img.width; oc.height = img.height;
        const ctx = oc.getContext('2d');
        ctx.drawImage(img, 0, 0);
        const pixels = ctx.getImageData(0, 0, img.width, img.height).data;

        // Collect wet cell positions (subsample if too many)
        const wetPositions = [];
        const step = Math.max(1, Math.floor(img.width * img.height / 50000)); // check at most 50K pixels
        for (let i = 0; i < img.width * img.height; i += step) {
            if (pixels[i * 4 + 3] > 128) { // alpha > 0.5 → wet
                const col = i % img.width;
                const row = Math.floor(i / img.width);
                wetPositions.push((col + 0.5) / img.width, (row + 0.5) / img.height);
            }
        }

        // Build spawn texture (Nx1 RGBA, positions encoded in RG + BA for precision)
        const maxSpawn = 2048;  // keep within WebGL 1 min texture dimension limit
        let nSpawn = Math.floor(wetPositions.length / 2);
        if (nSpawn > maxSpawn) nSpawn = maxSpawn;
        if (nSpawn < 1) nSpawn = 1; // at least one entry (will be 0.5, 0.5)

        const spawnData = new Uint8Array(nSpawn * 4);
        for (let i = 0; i < nSpawn; i++) {
            let si = i;
            if (wetPositions.length / 2 > maxSpawn) {
                // subsample uniformly
                si = Math.floor(i * (wetPositions.length / 2) / maxSpawn);
            }
            const x = si < wetPositions.length / 2 ? wetPositions[si * 2] : 0.5;
            const y = si < wetPositions.length / 2 ? wetPositions[si * 2 + 1] : 0.5;
            // Encode: R=high byte of x, G=high byte of y, B=low byte of x, A=low byte of y
            spawnData[i * 4 + 0] = Math.floor(x * 255);         // x high
            spawnData[i * 4 + 1] = Math.floor(y * 255);         // y high
            spawnData[i * 4 + 2] = Math.floor((x * 255 - Math.floor(x * 255)) * 255); // x low
            spawnData[i * 4 + 3] = Math.floor((y * 255 - Math.floor(y * 255)) * 255); // y low
        }

        const spawnTex = createTex(gl, gl.NEAREST, spawnData, nSpawn, 1);

        this.frameCache.set(idx, { velTex: tex, spawnTex: spawnTex, spawnCount: nSpawn, wetPositions: wetPositions });
        if (this.frameCache.size > this.frameCacheMax) {
            const oldest = this.frameCache.keys().next().value;
            const entry = this.frameCache.get(oldest);
            gl.deleteTexture(entry.velTex);
            gl.deleteTexture(entry.spawnTex);
            this.frameCache.delete(oldest);
        }
        return this.frameCache.get(idx);
    }

    preload() {
        for (let i = 1; i <= 5; i++) {
            const fi = (this.currentFrame + i) % this.meta.n_timesteps;
            if (!this.frameCache.has(fi)) this.loadVelFrame(fi);
        }
    }

    setParticleCount(power) {
        const gl = this.gl;
        this.particleRes = Math.pow(2, Math.ceil(power / 2));
        this.numParticles = this.particleRes * this.particleRes;
        if (this.particleFB0) gl.deleteFramebuffer(this.particleFB0);
        if (this.particleFB1) gl.deleteFramebuffer(this.particleFB1);
        if (this.particleTex0) gl.deleteTexture(this.particleTex0);
        if (this.particleTex1) gl.deleteTexture(this.particleTex1);
        this.initParticleState(gl);
        // Re-seed at wet positions if available
        const entry = this.frameCache.get(this.currentFrame);
        if (entry && entry.spawnCount > 1) this.initParticlesFromSpawn(gl, entry);
        const idx = new Float32Array(this.numParticles);
        for (let i = 0; i < this.numParticles; i++) idx[i] = i;
        gl.bindBuffer(gl.ARRAY_BUFFER, this.indexBuf);
        gl.bufferData(gl.ARRAY_BUFFER, idx, gl.STATIC_DRAW);
    }

    // ── Render loop ──
    render() {
        const now = performance.now();
        const dt = (now - this._lastT) / 1000;
        this._lastT = now;

        // FPS counter
        this._fpsFrames++;
        this._fpsTime += dt;
        if (this._fpsTime >= 0.5) {
            this._fps = Math.round(this._fpsFrames / this._fpsTime);
            this._fpsFrames = 0;
            this._fpsTime = 0;
            document.getElementById('i-fps').textContent = this._fps;
        }

        if (this.playing) {
            this.timeFrac += dt * this.speed * 1.5;
            while (this.timeFrac >= 1.0) {
                this.timeFrac -= 1.0;
                this.currentFrame = (this.currentFrame + 1) % this.meta.n_timesteps;
                this.preload();
            }
        }

        // Update UI
        const fi = this.currentFrame;
        const tSec = this.meta.times[fi] || 0;
        document.getElementById('sl-time').value = fi;
        document.getElementById('v-time').textContent = fmtTime(tSec);
        document.getElementById('i-time').textContent = fmtTime(tSec);
        document.getElementById('i-particles').textContent =
            (this.numParticles >= 1000 ?
                Math.round(this.numParticles/1024) + 'K' : this.numParticles);

        const gl = this.gl;
        const m = this.meta;
        const curEntry = this.frameCache.get(fi);
        const nextFi = (fi + 1) % m.n_timesteps;
        const nextEntry = this.frameCache.get(nextFi) || curEntry;
        if (!curEntry) { requestAnimationFrame(() => this.render()); return; }

        const velTex = curEntry.velTex;
        const velTexNext = nextEntry.velTex;
        const spawnTex = curEntry.spawnTex;
        const spawnCount = curEntry.spawnCount;

        // (1) Update particles
        gl.bindFramebuffer(gl.FRAMEBUFFER, this.particleFB1);
        gl.viewport(0, 0, this.particleRes, this.particleRes);
        const up = this.updateProg;
        gl.useProgram(up);

        bindAttr(gl, this.quadBuf, gl.getAttribLocation(up, 'a_pos'), 2);

        gl.activeTexture(gl.TEXTURE0);
        gl.bindTexture(gl.TEXTURE_2D, this.particleTex0);
        gl.uniform1i(gl.getUniformLocation(up, 'u_particles'), 0);

        gl.activeTexture(gl.TEXTURE1);
        gl.bindTexture(gl.TEXTURE_2D, velTex);
        gl.uniform1i(gl.getUniformLocation(up, 'u_velocity'), 1);

        gl.activeTexture(gl.TEXTURE2);
        gl.bindTexture(gl.TEXTURE_2D, velTexNext);
        gl.uniform1i(gl.getUniformLocation(up, 'u_velocity_next'), 2);

        gl.activeTexture(gl.TEXTURE3);
        gl.bindTexture(gl.TEXTURE_2D, spawnTex);
        gl.uniform1i(gl.getUniformLocation(up, 'u_spawn_tex'), 3);
        gl.uniform1f(gl.getUniformLocation(up, 'u_spawn_count'), spawnCount);

        gl.uniform1f(gl.getUniformLocation(up, 'u_time_frac'), this.timeFrac);
        gl.uniform1f(gl.getUniformLocation(up, 'u_rand_seed'), Math.random());
        gl.uniform1f(gl.getUniformLocation(up, 'u_speed_factor'), this.speedFactor * this.speed);
        gl.uniform1f(gl.getUniformLocation(up, 'u_drop_rate'), this.dropRate);
        gl.uniform1f(gl.getUniformLocation(up, 'u_drop_rate_bump'), this.dropRateBump);
        gl.uniform2f(gl.getUniformLocation(up, 'u_vel_range_x'), m.vx_min, m.vx_max);
        gl.uniform2f(gl.getUniformLocation(up, 'u_vel_range_y'), m.vy_min, m.vy_max);
        gl.uniform2f(gl.getUniformLocation(up, 'u_grid_res'), m.width, m.height);

        gl.drawArrays(gl.TRIANGLES, 0, 6);

        // Swap particle textures
        [this.particleTex0, this.particleTex1] = [this.particleTex1, this.particleTex0];
        [this.particleFB0, this.particleFB1] = [this.particleFB1, this.particleFB0];

        // (2) Fade screen + draw particles
        const cw = gl.canvas.width, ch = gl.canvas.height;
        gl.bindFramebuffer(gl.FRAMEBUFFER, this.screenFB1);
        gl.viewport(0, 0, cw, ch);

        // Fade previous screen
        const sp = this.screenProg;
        gl.useProgram(sp);
        bindAttr(gl, this.quadBuf, gl.getAttribLocation(sp, 'a_pos'), 2);
        gl.activeTexture(gl.TEXTURE0);
        gl.bindTexture(gl.TEXTURE_2D, this.screenTex0);
        gl.uniform1i(gl.getUniformLocation(sp, 'u_screen'), 0);
        gl.uniform1f(gl.getUniformLocation(sp, 'u_opacity'), this.fadeOpacity);
        gl.drawArrays(gl.TRIANGLES, 0, 6);

        // Draw particles on top (additive)
        gl.enable(gl.BLEND);
        gl.blendFunc(gl.SRC_ALPHA, gl.ONE);

        const dp = this.drawProg;
        gl.useProgram(dp);
        bindAttr(gl, this.indexBuf, gl.getAttribLocation(dp, 'a_index'), 1);

        gl.activeTexture(gl.TEXTURE0);
        gl.bindTexture(gl.TEXTURE_2D, this.particleTex0);
        gl.uniform1i(gl.getUniformLocation(dp, 'u_particles'), 0);
        gl.uniform1f(gl.getUniformLocation(dp, 'u_particles_res'), this.particleRes);

        gl.activeTexture(gl.TEXTURE1);
        gl.bindTexture(gl.TEXTURE_2D, velTex);
        gl.uniform1i(gl.getUniformLocation(dp, 'u_velocity'), 1);

        gl.activeTexture(gl.TEXTURE2);
        gl.bindTexture(gl.TEXTURE_2D, velTexNext);
        gl.uniform1i(gl.getUniformLocation(dp, 'u_velocity_next'), 2);

        gl.uniform1f(gl.getUniformLocation(dp, 'u_time_frac'), this.timeFrac);
        gl.uniform2f(gl.getUniformLocation(dp, 'u_vel_range_x'), m.vx_min, m.vx_max);
        gl.uniform2f(gl.getUniformLocation(dp, 'u_vel_range_y'), m.vy_min, m.vy_max);
        gl.uniform1f(gl.getUniformLocation(dp, 'u_point_size'), this.pointSize);

        gl.activeTexture(gl.TEXTURE3);
        gl.bindTexture(gl.TEXTURE_2D, this.colorRampTex);
        gl.uniform1i(gl.getUniformLocation(dp, 'u_color_ramp'), 3);

        gl.drawArrays(gl.POINTS, 0, this.numParticles);
        gl.disable(gl.BLEND);

        // Swap screen textures
        [this.screenTex0, this.screenTex1] = [this.screenTex1, this.screenTex0];
        [this.screenFB0, this.screenFB1] = [this.screenFB1, this.screenFB0];

        // (3) Composite to canvas
        gl.bindFramebuffer(gl.FRAMEBUFFER, null);
        gl.viewport(0, 0, cw, ch);

        const cp = this.compositeProg;
        gl.useProgram(cp);
        bindAttr(gl, this.quadBuf, gl.getAttribLocation(cp, 'a_pos'), 2);

        gl.activeTexture(gl.TEXTURE0);
        gl.bindTexture(gl.TEXTURE_2D, this.hillshadeTex);
        gl.uniform1i(gl.getUniformLocation(cp, 'u_hillshade'), 0);

        gl.activeTexture(gl.TEXTURE1);
        gl.bindTexture(gl.TEXTURE_2D, this.screenTex0);
        gl.uniform1i(gl.getUniformLocation(cp, 'u_trails'), 1);

        gl.uniform1f(gl.getUniformLocation(cp, 'u_dark'), this.darkFactor);
        gl.drawArrays(gl.TRIANGLES, 0, 6);

        requestAnimationFrame(() => this.render());
    }
}

function fmtTime(s) {
    if (s < 60) return Math.round(s) + 's';
    if (s < 3600) return (s/60).toFixed(1) + ' min';
    return (s/3600).toFixed(1) + ' hr';
}

// ── Init ──
const canvas = document.getElementById('gl');
function resize() {
    canvas.width = window.innerWidth;
    canvas.height = window.innerHeight;
}
resize();
window.addEventListener('resize', () => {
    resize();
    if (viewer) viewer.initScreenTextures(viewer.gl);
});

let viewer;
(async () => {
    viewer = new FloodViewer(canvas);
    await viewer.init();

    // Wire controls
    document.getElementById('sl-time').addEventListener('input', e => {
        viewer.currentFrame = parseInt(e.target.value);
        viewer.timeFrac = 0;
        viewer.preload();
    });
    document.getElementById('btn-play').addEventListener('click', () => {
        viewer.playing = !viewer.playing;
        const b = document.getElementById('btn-play');
        b.innerHTML = viewer.playing ? '&#9646;&#9646; Pause' : '&#9654; Play';
        b.classList.toggle('active', viewer.playing);
    });
    document.getElementById('btn-reset').addEventListener('click', () => {
        viewer.currentFrame = 0;
        viewer.timeFrac = 0;
        viewer.initScreenTextures(viewer.gl);
        viewer.initParticleState(viewer.gl);
    });
    document.getElementById('btn-fs').addEventListener('click', () => {
        document.documentElement.requestFullscreen().catch(() => {});
    });
    document.getElementById('sl-speed').addEventListener('input', e => {
        viewer.speed = parseInt(e.target.value) / 100;
        document.getElementById('v-speed').textContent = viewer.speed.toFixed(1) + 'x';
    });
    document.getElementById('sl-particles').addEventListener('input', e => {
        const pw = parseInt(e.target.value);
        viewer.setParticleCount(pw);
        const n = viewer.numParticles;
        document.getElementById('v-particles').textContent =
            n >= 1024 ? Math.round(n/1024) + 'K' : n;
    });
    document.getElementById('sl-trail').addEventListener('input', e => {
        viewer.fadeOpacity = parseInt(e.target.value) / 1000;
        document.getElementById('v-trail').textContent = viewer.fadeOpacity.toFixed(2);
    });
})();
</script>
</body>
</html>
'''


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

    # ── WebGL interactive viewer export ──────────────────────
    parser.add_argument(
        "--export-webgl", action="store_true", default=False,
        help="Export an interactive WebGL particle viewer (HTML + data files).",
    )
    parser.add_argument(
        "--webgl-output", type=str, default="fluxos_web",
        help="Output directory for WebGL viewer. Default: fluxos_web/",
    )
    parser.add_argument(
        "--webgl-step", type=int, default=5,
        help="Subsample every Nth timestep for WebGL export. Default: 5",
    )
    parser.add_argument(
        "--webgl-particles", type=int, default=65536,
        help="Default particle count in viewer. Default: 65536",
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
    if args.export_webgl:
        print(f"[4/4] Exporting interactive WebGL viewer "
              f"({args.mesh_type} mesh)...")
        export_webgl(
            dem, meta, results, args.mesh_type, args.variable, clim,
            output_dir=args.webgl_output, h_min=args.h_min,
            step=args.webgl_step, n_particles=args.webgl_particles,
            utm_to_ll=utm_to_ll,
        )
        return

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
