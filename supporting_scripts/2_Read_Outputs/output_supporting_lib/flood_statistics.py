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
- per-cell total wet-duration (hours each cell spent above threshold)
- cell-wise hazard classification (depth × velocity)
- per-timestep hazard-class areas (evolution over time)
- depth and velocity histograms across the wet domain
- percentiles of max depth and max velocity over wet cells
- supercritical-flow (Froude > 1) extent
- flood-dynamics metrics (rise time, recession time, mean wet depth at peak)
- mass-balance / peak-volume diagnostic

All functions are designed to **stream** the snapshots so we never hold all
timesteps in memory at once — important because a 30-hour river simulation
can easily produce 10+ GB of per-timestep files.
"""

from __future__ import annotations

import numpy as np

# Gravitational constant for Froude number (v / sqrt(g * h))
_G = 9.80665

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
    ``velocity_hist_edges``        : bin edges for max-velocity histogram
    ``velocity_hist_counts``       : velocity bin counts (over wet cells only)
    ``duration_wet_s``             : per-cell total wet duration (seconds)
    ``duration_hist_edges``        : bin edges for wet-duration histogram (hours)
    ``duration_hist_counts``       : duration bin counts (over ever-wet cells)
    ``hazard_area_per_time``       : list of 4-tuples (L, M, H, E) areas (m²) per timestep
    ``max_h_p50`` / ``max_h_p90`` / ``max_h_p95`` / ``max_h_p99`` : depth percentiles
    ``max_v_p50`` / ``max_v_p90`` / ``max_v_p95`` / ``max_v_p99`` : velocity percentiles
    ``froude_max``                 : global max Froude number observed
    ``froude_super_area_m2``       : area (m²) of cells that reached Fr > 1
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
    # New accumulators
    duration_wet_cell: np.ndarray | None = None      # per-cell seconds wet
    max_froude_cell: np.ndarray | None = None        # per-cell max Froude
    cell_area_arr: np.ndarray | None = None          # per-cell area (m^2)
    last_t_val: float | None = None                   # previous timestep seconds
    first_dt_placeholder: float | None = None         # assigned on 2nd step
    pending_first_wet: np.ndarray | None = None       # wet mask at first step, needs dt
    hazard_area_per_time: list[tuple[float, float, float, float]] = []

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

        # Current-snapshot Froude (where wet); NaN otherwise
        froude_now = None
        if v is not None and v.size == h.size:
            with np.errstate(divide="ignore", invalid="ignore"):
                froude_now = np.where(h > h_threshold,
                                      v / np.sqrt(_G * np.maximum(h, 1e-9)),
                                      0.0)

        # Current-snapshot hazard class for per-time area accumulation
        if v is not None and v.size == h.size:
            hv_now = h * (v + 0.5)
            cls_now = np.zeros_like(h, dtype=int)
            for i, edge in enumerate(_HAZARD_EDGES):
                cls_now = np.where(hv_now > edge, i + 1, cls_now)
            cls_now = np.where(h > h_threshold, cls_now, -1)  # dry -> -1
            haz_areas = [0.0, 0.0, 0.0, 0.0]
            for k in range(4):
                mask_k = cls_now == k
                if mask_k.any():
                    haz_areas[k] = float(np.sum(areas[mask_k]))
            hazard_area_per_time.append(tuple(haz_areas))
        else:
            hazard_area_per_time.append((0.0, 0.0, 0.0, 0.0))

        # Compute dt (seconds) for duration-wet accumulation. First step uses
        # the second step's delta as a proxy (fixed up retroactively).
        t_val = float(t)
        if last_t_val is None:
            dt_step = 0.0           # deferred — fixed on step 2
        else:
            dt_step = t_val - last_t_val
            if first_dt_placeholder is None:
                first_dt_placeholder = dt_step

        # Lazy-init per-cell accumulators
        if mesh_type == "triangular":
            if map_x is None:
                map_x = cx.copy(); map_y = cy.copy()
                max_h_cell = np.zeros_like(h, dtype=float)
                max_v_cell = np.zeros_like(h, dtype=float)
                first_wet = np.full_like(h, np.nan, dtype=float)
                duration_wet_cell = np.zeros_like(h, dtype=float)
                max_froude_cell = np.zeros_like(h, dtype=float)
                cell_area_arr = areas.astype(float).copy()
                if has_conc:
                    max_conc_cell = np.zeros_like(h, dtype=float)
            np.maximum(max_h_cell, h, out=max_h_cell)
            if v is not None and v.size == max_v_cell.size:
                np.maximum(max_v_cell, v, out=max_v_cell)
            just_wet = wet & np.isnan(first_wet)
            first_wet[just_wet] = float(t)
            if dt_step > 0 and wet.size == duration_wet_cell.size:
                duration_wet_cell += wet.astype(float) * dt_step
            else:
                # Store wet mask of first step; we'll add its dt once known
                pending_first_wet = wet.copy()
            if froude_now is not None and froude_now.size == max_froude_cell.size:
                np.maximum(max_froude_cell, froude_now, out=max_froude_cell)
            if has_conc and conc is not None:
                np.maximum(max_conc_cell, conc, out=max_conc_cell)
        else:  # regular
            if map_x is None:
                # Allocate on first pass based on total_cells
                max_h_cell = np.zeros(total_cells, dtype=float)
                max_v_cell = np.zeros(total_cells, dtype=float)
                first_wet = np.full(total_cells, np.nan, dtype=float)
                duration_wet_cell = np.zeros(total_cells, dtype=float)
                max_froude_cell = np.zeros(total_cells, dtype=float)
                cell_area_arr = np.full(total_cells, cellsize * cellsize, dtype=float)
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
            if dt_step > 0:
                # Accumulate wet duration
                duration_wet_cell[idx[wet_mask]] += dt_step
            else:
                # First timestep — we only see wet cells on regular mesh;
                # stash their linear indices to charge with the placeholder dt later.
                pending_first_wet = idx[wet_mask].copy()
            if froude_now is not None:
                np.maximum.at(max_froude_cell, idx, froude_now)
            if has_conc and conc is not None:
                np.maximum.at(max_conc_cell, idx, conc)

        last_t_val = t_val

        if progress:
            progress(n_steps)

    # Retroactively charge the first-step wet cells with dt (using 2nd step delta).
    if first_dt_placeholder is not None and pending_first_wet is not None \
            and duration_wet_cell is not None:
        if mesh_type == "triangular":
            if pending_first_wet.size == duration_wet_cell.size:
                duration_wet_cell += pending_first_wet.astype(float) * first_dt_placeholder
        else:
            duration_wet_cell[pending_first_wet] += first_dt_placeholder

    # Compact regular-mesh per-cell arrays: drop cells that never got wet
    if mesh_type != "triangular" and max_h_cell is not None:
        ever_wet = max_h_cell > 0
        map_x = map_x[ever_wet]
        map_y = map_y[ever_wet]
        max_h_cell = max_h_cell[ever_wet]
        max_v_cell = max_v_cell[ever_wet]
        first_wet = first_wet[ever_wet]
        if duration_wet_cell is not None:
            duration_wet_cell = duration_wet_cell[ever_wet]
        if max_froude_cell is not None:
            max_froude_cell = max_froude_cell[ever_wet]
        if cell_area_arr is not None:
            cell_area_arr = cell_area_arr[ever_wet]
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

    # Max-velocity histogram (linear bins with 99th-percentile cutoff to
    # avoid a long empty tail dominated by numerical spikes).
    if max_v_cell is not None and (max_v_cell > 0).any():
        wet_v = max_v_cell[max_v_cell > 0]
        v_hi = float(np.percentile(wet_v, 99))
        v_hi = max(v_hi, 0.1)
        vel_hist_counts, vel_hist_edges = np.histogram(
            np.clip(wet_v, 0, v_hi), bins=25, range=(0, v_hi))
    else:
        vel_hist_counts = np.array([])
        vel_hist_edges = np.array([])

    # Wet-duration histogram (hours)
    if duration_wet_cell is not None and (duration_wet_cell > 0).any():
        dur_h = duration_wet_cell[duration_wet_cell > 0] / 3600.0
        dur_hist_counts, dur_hist_edges = np.histogram(dur_h, bins=25)
    else:
        dur_hist_counts = np.array([])
        dur_hist_edges = np.array([])

    # Percentiles of max_h / max_v over ever-wet cells
    def _pctiles(arr, mask):
        if arr is None or mask is None or not mask.any():
            return (0.0, 0.0, 0.0, 0.0)
        sub = arr[mask]
        if sub.size == 0:
            return (0.0, 0.0, 0.0, 0.0)
        q = np.percentile(sub, [50, 90, 95, 99])
        return tuple(float(x) for x in q)

    ever_wet_mask = (max_h_cell > 0) if max_h_cell is not None else None
    max_h_p50, max_h_p90, max_h_p95, max_h_p99 = _pctiles(max_h_cell, ever_wet_mask)
    max_v_p50, max_v_p90, max_v_p95, max_v_p99 = _pctiles(max_v_cell, ever_wet_mask)

    # Froude stats
    froude_max_val = float(max_froude_cell.max()) if (max_froude_cell is not None
                                                      and max_froude_cell.size) else 0.0
    froude_super_area = 0.0
    if max_froude_cell is not None and cell_area_arr is not None \
            and max_froude_cell.size == cell_area_arr.size:
        sup = max_froude_cell > 1.0
        if sup.any():
            froude_super_area = float(np.sum(cell_area_arr[sup]))

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
        "map_cell_area_m2": cell_area_arr,
        "map_max_h": max_h_cell,
        "map_max_v": max_v_cell,
        "map_max_froude": max_froude_cell,
        "map_first_wet_time_s": first_wet,
        "duration_wet_s": duration_wet_cell,
        "hazard_class": hazard_cls,
        "hazard_counts": hazard_counts,
        "hazard_labels": _HAZARD_LABELS,
        "hazard_edges": _HAZARD_EDGES,
        "hazard_area_per_time": hazard_area_per_time,
        "depth_hist_edges": depth_hist_edges.tolist(),
        "depth_hist_counts": depth_hist_counts.tolist(),
        "velocity_hist_edges": vel_hist_edges.tolist(),
        "velocity_hist_counts": vel_hist_counts.tolist(),
        "duration_hist_edges": dur_hist_edges.tolist(),
        "duration_hist_counts": dur_hist_counts.tolist(),
        "max_h_p50": max_h_p50, "max_h_p90": max_h_p90,
        "max_h_p95": max_h_p95, "max_h_p99": max_h_p99,
        "max_v_p50": max_v_p50, "max_v_p90": max_v_p90,
        "max_v_p95": max_v_p95, "max_v_p99": max_v_p99,
        "froude_max": froude_max_val,
        "froude_super_area_m2": froude_super_area,
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
    mean_h = np.asarray(stats.get("mean_h_wet_m", []))
    time_s = np.asarray(stats["time_s"])

    peak_vol = float(vol.max()) if vol.size else 0.0
    peak_area = float(area.max()) if area.size else 0.0
    peak_max_h = float(max_h.max()) if max_h.size else 0.0
    peak_max_v = float(max_v.max()) if max_v.size else 0.0
    peak_t = int(time_s[int(vol.argmax())]) if vol.size else 0
    total_t = int(time_s[-1] - time_s[0]) if time_s.size > 1 else 0
    n_wet_cells = int(np.sum(np.asarray(stats.get("map_max_h", [])) > 0))

    # --- Flood-dynamics timing (based on flooded-area curve) ---------------
    rise_time_s = None
    recession_time_s = None
    total_inundation_duration_s = None
    if area.size and time_s.size and peak_area > 0:
        thr10 = 0.10 * peak_area
        thr50 = 0.50 * peak_area
        thr90 = 0.90 * peak_area
        above10 = np.where(area >= thr10)[0]
        above90 = np.where(area >= thr90)[0]
        if above10.size and above90.size:
            t_rise_start = float(time_s[above10[0]])
            t_rise_end = float(time_s[above90[0]])
            rise_time_s = max(0.0, t_rise_end - t_rise_start)
        # recession: from peak to first-time area < thr50
        peak_idx = int(area.argmax())
        post = area[peak_idx:]
        rec_rel = np.where(post < thr50)[0]
        if rec_rel.size:
            t_peak = float(time_s[peak_idx])
            t_rec = float(time_s[peak_idx + int(rec_rel[0])])
            recession_time_s = max(0.0, t_rec - t_peak)
        if above10.size:
            total_inundation_duration_s = float(time_s[above10[-1]] - time_s[above10[0]])

    # Mean wet depth at the peak timestep
    mean_wet_depth_at_peak_m = 0.0
    if mean_h.size and vol.size:
        mean_wet_depth_at_peak_m = float(mean_h[int(vol.argmax())])

    # H3+H4 area from the max-field hazard classification
    hazard_h3_h4_area_m2 = 0.0
    hcls = np.asarray(stats.get("hazard_class", []))
    cell_area = stats.get("map_cell_area_m2")
    if hcls.size and cell_area is not None and len(cell_area) == hcls.size:
        mask = (hcls == 2) | (hcls == 3)
        if mask.any():
            hazard_h3_h4_area_m2 = float(np.sum(np.asarray(cell_area)[mask]))

    return dict(
        peak_volume_m3=peak_vol,
        peak_area_m2=peak_area,
        peak_depth_m=peak_max_h,
        peak_velocity_ms=peak_max_v,
        peak_time_s=peak_t,
        sim_duration_s=total_t,
        wet_cells=n_wet_cells,
        # New metrics
        rise_time_s=rise_time_s,
        recession_time_s=recession_time_s,
        total_inundation_duration_s=total_inundation_duration_s,
        mean_wet_depth_at_peak_m=mean_wet_depth_at_peak_m,
        hazard_h3_h4_area_m2=hazard_h3_h4_area_m2,
        max_h_p95=float(stats.get("max_h_p95") or 0.0),
        max_v_p95=float(stats.get("max_v_p95") or 0.0),
        froude_max=float(stats.get("froude_max") or 0.0),
        froude_super_area_m2=float(stats.get("froude_super_area_m2") or 0.0),
    )
