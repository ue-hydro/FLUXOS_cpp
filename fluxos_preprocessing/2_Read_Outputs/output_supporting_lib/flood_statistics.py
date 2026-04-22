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
Flood-simulation statistics computed from a stream of FLUXOS snapshots.

These analyses complement the dynamic WebGL animation with aggregated
numbers and maps that are hard to read off a moving picture:

- per-timestep aggregates (volume, flooded area, peak depth, ...)
- per-cell maxima (the classic "maximum inundation map")
- per-cell first-inundation time (flood-front propagation)
- cell-wise hazard classification (depth × velocity)
- depth histogram across the wet domain
- mass-balance / peak-volume diagnostic

All functions are designed to **stream** the snapshots so we never hold all
timesteps in memory at once — important because a 30-hour river simulation
can easily produce 10+ GB of per-timestep files.
"""

from __future__ import annotations

import numpy as np

# H × V hazard thresholds (metres × m/s = m²/s)
# These follow the Australian Rainfall & Runoff 2019 / UK DEFRA criterion.
_HAZARD_EDGES = [0.3, 0.6, 1.2]       # Low | Medium | High | Extreme
_HAZARD_LABELS = ["Low (H1)", "Medium (H2)", "High (H3)", "Extreme (H4)"]


# --- Small helpers --------------------------------------------------------

def _speed(snap: dict, cellsize: float | None = None) -> np.ndarray | None:
    """
    Velocity magnitude per cell. Prefers ``ux/uy``, but the regular-mesh
    writer sometimes emits zeros there — fall back to ``qx*dxy / qy*dxy``
    (which is the unit discharge times the cell edge) divided by h·cellsize.
    """
    ux, uy = snap.get("ux"), snap.get("uy")
    if ux is not None and uy is not None:
        u = np.asarray(ux); v = np.asarray(uy)
        mag = np.hypot(u, v)
        if mag.size and mag.max() > 0:
            return mag
    # Fall back: recover velocity from discharge columns when available.
    # Only evaluate where h > 5 cm — at the wet/dry interface, very thin films
    # combined with small discharges give unphysical ~100 m/s spikes that
    # dominate the "max velocity" plot and the hazard classification.
    qx = snap.get("qx_dxy"); qy = snap.get("qy_dxy")
    h  = snap.get("h")
    if qx is not None and qy is not None and h is not None and cellsize:
        qx = np.asarray(qx); qy = np.asarray(qy); h = np.asarray(h)
        ok = h > 0.05
        with np.errstate(divide="ignore", invalid="ignore"):
            ux = np.where(ok, qx / (h * cellsize), 0.0)
            uy = np.where(ok, qy / (h * cellsize), 0.0)
        return np.hypot(ux, uy)
    return None


def _cell_areas_regular(snap: dict, cellsize: float) -> np.ndarray:
    n = int(np.asarray(snap["h"]).size)
    return np.full(n, cellsize * cellsize, dtype=float)


# --- Core: single streaming pass computes everything --------------------

def analyse(iterator, mesh_type: str, *,
            cellsize: float | None = None,
            ncols: int | None = None,
            nrows: int | None = None,
            h_threshold: float = 0.01,
            progress=None) -> dict:
    """
    Consume a ``(time_s, snapshot)`` iterator once and compute every statistic
    the report cares about.

    Returns a dict with keys:

    ``time_s``                     : list of timestamps (seconds)
    ``volume_m3``                  : flood volume per timestep
    ``area_flooded_m2``            : area with h > threshold per timestep
    ``max_h_m``                    : max water depth per timestep
    ``mean_h_wet_m``               : mean depth of wet cells per timestep
    ``max_v_ms``                   : max velocity magnitude per timestep
    ``peak_time_s``                : t at which volume peaks
    ``map_cell_x``, ``map_cell_y`` : cell centre coordinates
    ``map_max_h``                  : per-cell max depth (full simulation)
    ``map_max_v``                  : per-cell max speed
    ``map_first_wet_time_s``       : per-cell first time h > threshold (NaN if never)
    ``hazard_class``               : per-cell hazard bin index (0..3)
    ``hazard_counts``              : count per hazard class (np.int)
    ``depth_hist_edges``           : bin edges for max-depth histogram
    ``depth_hist_counts``          : bin counts (over wet cells only)
    ``has_conc``                   : True if conc_SW data was present
    ``max_conc_SW``                : per-cell max concentration (or None)
    ``mean_conc_SW_t``             : per-timestep mean conc over wet cells (or None)
    ``max_conc_SW_t``              : per-timestep max conc (or None)
    """
    total_cells = None
    if mesh_type != "triangular":
        if ncols is None or nrows is None:
            raise ValueError("ncols and nrows are required for regular mesh.")
        total_cells = int(ncols * nrows)
    time_s: list[int] = []
    vol_m3: list[float] = []
    area_m2: list[float] = []
    max_h_t: list[float] = []
    mean_h_wet_t: list[float] = []
    max_v_t: list[float] = []
    max_conc_t: list[float] = []
    mean_conc_t: list[float] = []

    # Per-cell running accumulators (allocated lazily once we know total_cells)
    map_x: np.ndarray | None = None
    map_y: np.ndarray | None = None
    max_h_cell: np.ndarray | None = None
    max_v_cell: np.ndarray | None = None
    first_wet: np.ndarray | None = None
    max_conc_cell: np.ndarray | None = None
    has_conc = False

    n_steps = 0
    for t, snap in iterator:
        n_steps += 1
        _h_raw = snap.get("h")
        h = np.asarray(_h_raw) if _h_raw is not None else np.empty(0)
        v = _speed(snap, cellsize=cellsize)
        if mesh_type == "triangular":
            areas = np.asarray(snap.get("cell_areas"))
            cx = np.asarray(snap.get("cell_x"))
            cy = np.asarray(snap.get("cell_y"))
            idx = None  # triangular files contain all cells each timestep
        else:
            if cellsize is None:
                raise ValueError("cellsize must be provided for regular mesh.")
            areas = _cell_areas_regular(snap, cellsize)
            cx = np.asarray(snap["x"])
            cy = np.asarray(snap["y"])
            # (icol, irow) linearised index into the full grid
            irow = np.asarray(snap["irow"]).astype(int)
            icol = np.asarray(snap["icol"]).astype(int)
            # Guard against out-of-range indices from corrupted rows
            valid = (irow >= 0) & (irow < nrows) & (icol >= 0) & (icol < ncols)
            if not valid.all():
                h = h[valid]; cx = cx[valid]; cy = cy[valid]
                irow = irow[valid]; icol = icol[valid]
                if v is not None:
                    v = v[valid]
                areas = areas[:len(h)]
            idx = irow * ncols + icol

        if h.size == 0:
            time_s.append(int(t)); vol_m3.append(0.0); area_m2.append(0.0)
            max_h_t.append(0.0); mean_h_wet_t.append(0.0); max_v_t.append(0.0)
            continue

        wet = h > h_threshold
        vol = float(np.sum(h * areas))
        area = float(np.sum(areas[wet]))
        max_h = float(h.max())
        mean_h = float(h[wet].mean()) if wet.any() else 0.0
        max_v = float(v.max()) if v is not None and v.size else 0.0

        time_s.append(int(t))
        vol_m3.append(vol); area_m2.append(area)
        max_h_t.append(max_h); mean_h_wet_t.append(mean_h); max_v_t.append(max_v)

        conc = snap.get("conc_SW")
        if conc is not None and np.asarray(conc).size:
            has_conc = True
            conc = np.asarray(conc)
            max_conc_t.append(float(conc.max()))
            mean_conc_t.append(float(conc[wet].mean()) if wet.any() else 0.0)
        else:
            max_conc_t.append(0.0)
            mean_conc_t.append(0.0)

        # Lazy-init per-cell accumulators
        if mesh_type == "triangular":
            if map_x is None:
                map_x = cx.copy(); map_y = cy.copy()
                max_h_cell = np.zeros_like(h, dtype=float)
                max_v_cell = np.zeros_like(h, dtype=float)
                first_wet = np.full_like(h, np.nan, dtype=float)
                if has_conc:
                    max_conc_cell = np.zeros_like(h, dtype=float)
            np.maximum(max_h_cell, h, out=max_h_cell)
            if v is not None and v.size == max_v_cell.size:
                np.maximum(max_v_cell, v, out=max_v_cell)
            just_wet = wet & np.isnan(first_wet)
            first_wet[just_wet] = float(t)
            if has_conc and conc is not None:
                np.maximum(max_conc_cell, conc, out=max_conc_cell)
        else:  # regular
            if map_x is None:
                # Allocate on first pass based on total_cells
                max_h_cell = np.zeros(total_cells, dtype=float)
                max_v_cell = np.zeros(total_cells, dtype=float)
                first_wet = np.full(total_cells, np.nan, dtype=float)
                map_x = np.full(total_cells, np.nan, dtype=float)
                map_y = np.full(total_cells, np.nan, dtype=float)
                if has_conc:
                    max_conc_cell = np.zeros(total_cells, dtype=float)
            # Record coordinates the first time we see each wet cell
            new_cells = np.isnan(map_x[idx])
            map_x[idx[new_cells]] = cx[new_cells]
            map_y[idx[new_cells]] = cy[new_cells]
            np.maximum.at(max_h_cell, idx, h)
            if v is not None:
                np.maximum.at(max_v_cell, idx, v)
            wet_mask = h > h_threshold
            new_wet = wet_mask & np.isnan(first_wet[idx])
            first_wet[idx[new_wet]] = float(t)
            if has_conc and conc is not None:
                np.maximum.at(max_conc_cell, idx, conc)

        if progress:
            progress(n_steps)

    # Compact regular-mesh per-cell arrays: drop cells that never got wet
    if mesh_type != "triangular" and max_h_cell is not None:
        ever_wet = max_h_cell > 0
        map_x = map_x[ever_wet]
        map_y = map_y[ever_wet]
        max_h_cell = max_h_cell[ever_wet]
        max_v_cell = max_v_cell[ever_wet]
        first_wet = first_wet[ever_wet]
        if max_conc_cell is not None:
            max_conc_cell = max_conc_cell[ever_wet]

    # Hazard classification on max_h × max_v
    hazard_cls = np.zeros_like(max_h_cell, dtype=int) if max_h_cell is not None \
        else np.zeros(0, dtype=int)
    if max_h_cell is not None and max_v_cell is not None:
        hv = max_h_cell * (max_v_cell + 0.5)  # ARR-2019 D-factor
        for i, edge in enumerate(_HAZARD_EDGES):
            hazard_cls = np.where(hv > edge, i + 1, hazard_cls)
        hazard_counts = np.bincount(hazard_cls[max_h_cell > 0],
                                    minlength=len(_HAZARD_LABELS))
    else:
        hazard_counts = np.zeros(len(_HAZARD_LABELS), dtype=int)

    # Max-depth histogram (over cells that were ever wet)
    if max_h_cell is not None and (max_h_cell > 0).any():
        wet_max = max_h_cell[max_h_cell > 0]
        depth_hist_counts, depth_hist_edges = np.histogram(wet_max, bins=25)
    else:
        depth_hist_counts = np.array([])
        depth_hist_edges = np.array([])

    peak_time = time_s[int(np.argmax(vol_m3))] if vol_m3 else 0

    return {
        "n_timesteps": n_steps,
        "time_s": time_s,
        "volume_m3": vol_m3,
        "area_flooded_m2": area_m2,
        "max_h_m": max_h_t,
        "mean_h_wet_m": mean_h_wet_t,
        "max_v_ms": max_v_t,
        "peak_time_s": peak_time,
        "map_cell_x": map_x, "map_cell_y": map_y,
        "map_max_h": max_h_cell,
        "map_max_v": max_v_cell,
        "map_first_wet_time_s": first_wet,
        "hazard_class": hazard_cls,
        "hazard_counts": hazard_counts,
        "hazard_labels": _HAZARD_LABELS,
        "hazard_edges": _HAZARD_EDGES,
        "depth_hist_edges": depth_hist_edges.tolist(),
        "depth_hist_counts": depth_hist_counts.tolist(),
        "has_conc": has_conc,
        "max_conc_SW": max_conc_cell,
        "mean_conc_SW_t": mean_conc_t if has_conc else None,
        "max_conc_SW_t": max_conc_t if has_conc else None,
    }


# --- Summary KPIs ---------------------------------------------------------

def summary_kpis(stats: dict) -> dict:
    """Condense the streaming output into single-value KPIs for the report header."""
    vol = np.asarray(stats["volume_m3"])
    area = np.asarray(stats["area_flooded_m2"])
    max_h = np.asarray(stats["max_h_m"])
    max_v = np.asarray(stats["max_v_ms"])
    time_s = np.asarray(stats["time_s"])

    peak_vol = float(vol.max()) if vol.size else 0.0
    peak_area = float(area.max()) if area.size else 0.0
    peak_max_h = float(max_h.max()) if max_h.size else 0.0
    peak_max_v = float(max_v.max()) if max_v.size else 0.0
    peak_t = int(time_s[int(vol.argmax())]) if vol.size else 0
    total_t = int(time_s[-1] - time_s[0]) if time_s.size > 1 else 0
    n_wet_cells = int(np.sum(np.asarray(stats.get("map_max_h", [])) > 0))

    return dict(
        peak_volume_m3=peak_vol,
        peak_area_m2=peak_area,
        peak_depth_m=peak_max_h,
        peak_velocity_ms=peak_max_v,
        peak_time_s=peak_t,
        sim_duration_s=total_t,
        wet_cells=n_wet_cells,
    )
