# Copyright 2026, Diogo Costa, diogo.costa@uevora.pt
# This file is part of the FLUXOS model.

# This program, FLUXOS, is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Readers for FLUXOS simulation output.

Two formats are supported:
- Regular (Cartesian) mesh: per-timestep ``<seconds>.txt`` CSVs with columns
  ``irow,icol,Xcoord,Ycoord,z,h,ux,uy,qx*dxy,qy*dxy,us,conc_SW,fe_1,fn_1,
  twetimetracer,soil_infil_rate``. Only wet cells (h > H_MIN_TO_PRINT) are written.
- Triangular mesh: per-timestep ``<seconds>.vtu`` files (VTK XML Unstructured
  Grid) + a ``.pvd`` collection. Cell-data arrays usually include ``h``, ``ux``,
  ``uy``, optionally ``conc_SW``.
"""

from __future__ import annotations

import glob
import json
import os
import re
from typing import Iterator

import numpy as np


# --- Modset ---------------------------------------------------------------

def read_modset(path: str) -> dict:
    with open(path) as f:
        return json.load(f)


def detect_mesh_type(modset: dict) -> str:
    return modset.get("MESH_TYPE", "regular")


# --- ESRI ASCII DEM header ------------------------------------------------

def read_asc_header(path: str) -> dict:
    """Read the 6-line ESRI ASCII header. Returns keys lowercased."""
    hdr = {}
    with open(path) as f:
        for _ in range(6):
            line = f.readline().strip()
            if not line:
                break
            parts = line.split()
            if len(parts) < 2:
                continue
            k, v = parts[0].lower(), parts[1]
            try:
                hdr[k] = int(v) if re.fullmatch(r"-?\d+", v) else float(v)
            except ValueError:
                hdr[k] = v
    return hdr


# --- Regular mesh readers -------------------------------------------------

# Canonical column order written by FLUXOS (see src/write_functions.cpp).
_REG_COLS = [
    "irow", "icol", "x", "y", "z", "h", "ux", "uy",
    "qx_dxy", "qy_dxy", "us", "conc_SW",
    "fe_1", "fn_1", "twetime", "soil_infil_rate",
]


def list_regular_timesteps(results_dir: str) -> list[tuple[int, str]]:
    """Return ``[(time_s, filepath), ...]`` sorted by time."""
    out = []
    for fp in glob.glob(os.path.join(results_dir, "*.txt")):
        stem = os.path.basename(fp)[:-4]
        try:
            t = int(stem)
        except ValueError:
            continue
        out.append((t, fp))
    out.sort(key=lambda r: r[0])
    return out


def read_regular_snapshot(fp: str) -> dict:
    """
    Read a single ``.txt`` snapshot.

    Returns a dict of numpy arrays (one entry per wet cell). Missing columns
    (for older output files) are returned as ``None``.
    """
    data = np.genfromtxt(fp, delimiter=",", skip_header=1, dtype=float)
    if data.size == 0:
        return {k: (None if k in ("conc_SW", "twetime", "soil_infil_rate")
                    else np.empty(0)) for k in _REG_COLS}
    if data.ndim == 1:
        data = data.reshape(1, -1)
    ncols = data.shape[1]
    out = {}
    for idx, name in enumerate(_REG_COLS):
        if idx < ncols:
            out[name] = (data[:, idx].astype(int)
                         if name in ("irow", "icol") else data[:, idx])
        else:
            out[name] = None
    return out


# --- Triangular mesh readers ----------------------------------------------

def list_triangular_timesteps(results_dir: str) -> list[tuple[int, str]]:
    out = []
    for fp in glob.glob(os.path.join(results_dir, "*.vtu")):
        stem = os.path.basename(fp)[:-4]
        try:
            t = int(stem)
        except ValueError:
            continue
        out.append((t, fp))
    out.sort(key=lambda r: r[0])
    return out


def read_triangular_snapshot(fp: str) -> dict:
    """
    Read a ``.vtu`` timestep.

    Returns dict with cell-centred arrays (h, ux, uy, conc_SW), cell
    coordinates, per-cell areas, and the raw pyvista mesh (for plotting).
    """
    try:
        import pyvista as pv
    except ImportError as e:
        raise ImportError(
            "pyvista is required to read triangular .vtu output. "
            "Install with `pip install pyvista`.") from e

    mesh = pv.read(fp)

    def _arr(name):
        if name in mesh.cell_data.keys():
            return np.asarray(mesh.cell_data[name])
        if name in mesh.point_data.keys():
            return np.asarray(mesh.point_to_cell_data().cell_data[name])
        return None

    centers = mesh.cell_centers().points
    areas = _triangle_areas(mesh)

    return {
        "n_cells": int(mesh.n_cells),
        "cell_x": centers[:, 0],
        "cell_y": centers[:, 1],
        "cell_areas": areas,
        "h":  _arr("h"),
        "ux": _arr("ux"),
        "uy": _arr("uy"),
        "conc_SW": _arr("conc_SW"),
        "pv_mesh": mesh,
    }


def _triangle_areas(mesh) -> np.ndarray:
    """Per-cell triangle area from a pyvista surface/unstructured mesh."""
    pts = np.asarray(mesh.points)
    # Cells stored as (3, i0, i1, i2) tuples in mesh.cells
    try:
        tris = np.asarray(mesh.regular_faces if hasattr(mesh, "regular_faces")
                          else mesh.faces.reshape(-1, 4)[:, 1:])
    except Exception:
        tris = np.asarray(mesh.cells_dict[5]) if 5 in getattr(mesh, "cells_dict", {}) \
            else np.asarray(mesh.cells.reshape(-1, 4)[:, 1:])
    p0 = pts[tris[:, 0]]
    p1 = pts[tris[:, 1]]
    p2 = pts[tris[:, 2]]
    return 0.5 * np.abs(np.cross(p1 - p0, p2 - p0)[..., -1]
                        if (p1 - p0).shape[-1] == 3 else
                        (p1[:, 0] - p0[:, 0]) * (p2[:, 1] - p0[:, 1])
                        - (p1[:, 1] - p0[:, 1]) * (p2[:, 0] - p0[:, 0]))


# --- Unified snapshot iterator --------------------------------------------

def iter_snapshots(results_dir: str, mesh_type: str,
                   stride: int = 1) -> Iterator[tuple[int, dict]]:
    """
    Yield ``(time_s, snapshot_dict)`` in time order.

    ``stride > 1`` samples every Nth timestep for faster preview stats.
    """
    if mesh_type == "triangular":
        tsteps = list_triangular_timesteps(results_dir)
        reader = read_triangular_snapshot
    else:
        tsteps = list_regular_timesteps(results_dir)
        reader = read_regular_snapshot

    for t, fp in tsteps[::max(1, int(stride))]:
        yield t, reader(fp)


def count_timesteps(results_dir: str, mesh_type: str) -> int:
    fn = list_triangular_timesteps if mesh_type == "triangular" \
        else list_regular_timesteps
    return len(fn(results_dir))
