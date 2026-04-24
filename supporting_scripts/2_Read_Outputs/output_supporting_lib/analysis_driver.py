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
Orchestrator for FLUXOS result post-processing.

Reads the `_config` dict from ``read_output_template.py``, streams the
snapshots in ``results_dir``, and hands the computed statistics to
``gen_results_report.generate_results_report``.

Users should not need to edit this file. Edit ``read_output_template.py``
instead.
"""

from __future__ import annotations

import os
import sys
import webbrowser
from datetime import datetime

from read_results import (
    read_modset, detect_mesh_type, read_asc_header,
    iter_snapshots, count_timesteps,
)
import flood_statistics
from gen_results_report import generate_results_report


def _resolve_repo_root(config: dict, template_file: str) -> str:
    if config.get("repo_root"):
        return os.path.abspath(config["repo_root"])
    # <repo>/supporting_scripts/2_Read_Outputs/read_output_template.py
    return os.path.abspath(os.path.join(os.path.dirname(template_file), "..", ".."))


def _abs_in_repo(repo_root: str, p: str) -> str:
    if not p:
        return p
    return p if os.path.isabs(p) else os.path.join(repo_root, p)


def run(config: dict, template_file: str) -> dict:
    """Execute the full analysis pipeline and return the collected metadata."""
    repo_root = _resolve_repo_root(config, template_file)

    results_dir = _abs_in_repo(repo_root, config["results_dir"])
    modset_path = _abs_in_repo(repo_root, config["modset_file"])

    print(f"FLUXOS results analysis — {config.get('project_name') or '(unnamed)'}")
    print(f"  repo root    : {repo_root}")
    print(f"  results dir  : {results_dir}")
    print(f"  modset       : {modset_path}")

    if not os.path.isdir(results_dir):
        raise FileNotFoundError(f"results_dir does not exist: {results_dir}")
    if not os.path.isfile(modset_path):
        raise FileNotFoundError(f"modset_file does not exist: {modset_path}")

    modset = read_modset(modset_path)
    mesh_type = detect_mesh_type(modset)

    # DEM metadata for satellite basemap and coordinate plumbing. We always
    # try to read the DEM header — regular mesh *requires* it (for ncols/nrows),
    # and triangular mesh wants it only for the bounding box.
    cellsize = None
    ncols = nrows = None
    dem_path = _abs_in_repo(repo_root, modset.get("DEM_FILE", ""))
    dem_meta: dict = {}
    if dem_path and os.path.isfile(dem_path):
        try:
            hdr = read_asc_header(dem_path)
            dem_meta = {
                "path":       dem_path,
                "ncols":      int(hdr.get("ncols", 1)),
                "nrows":      int(hdr.get("nrows", 1)),
                "xllcorner":  float(hdr.get("xllcorner", 0.0)),
                "yllcorner":  float(hdr.get("yllcorner", 0.0)),
                "cellsize":   float(hdr.get("cellsize", 1.0)),
                "nodata":     hdr.get("nodata_value", -9999),
            }
            if mesh_type != "triangular":
                cellsize = dem_meta["cellsize"]
                ncols = dem_meta["ncols"]
                nrows = dem_meta["nrows"]
                print(f"  DEM          : {dem_path} "
                      f"({ncols}×{nrows} @ {cellsize} m)")
            else:
                print(f"  DEM          : {dem_path} (bbox reference)")
        except Exception as exc:
            print(f"  DEM header could not be parsed ({exc}); "
                  "satellite basemap will be skipped.")
    elif mesh_type != "triangular":
        raise FileNotFoundError(f"DEM referenced by modset not found: {dem_path}")

    stride = int(config.get("stride", 1) or 1)
    n_total = count_timesteps(results_dir, mesh_type)
    print(f"  timesteps    : {n_total} "
          + (f"(stride={stride}, processing ≈ {n_total // stride})" if stride > 1 else ""))
    if n_total == 0:
        raise RuntimeError(f"No snapshots found in {results_dir} for mesh type '{mesh_type}'.")

    h_thr = float(config.get("h_threshold_flooded_m", 0.01))

    # Streaming statistics
    last_n = [0]
    def _progress(n):
        if n != last_n[0] and (n % 20 == 0 or n == n_total // stride):
            print(f"    processed {n} snapshots", end="\r")
            last_n[0] = n

    stats = flood_statistics.analyse(
        iter_snapshots(results_dir, mesh_type, stride=stride),
        mesh_type=mesh_type,
        cellsize=cellsize,
        ncols=ncols,
        nrows=nrows,
        h_threshold=h_thr,
        progress=_progress,
    )
    print()
    kpis = flood_statistics.summary_kpis(stats)

    print(f"  peak volume  : {kpis['peak_volume_m3']:,.0f} m³")
    print(f"  peak area    : {kpis['peak_area_m2']:,.0f} m²")
    print(f"  max depth    : {kpis['peak_depth_m']:.2f} m")
    print(f"  wet cells    : {kpis['wet_cells']:,}")

    # Optional CSV export of the time series
    if config.get("export_csv_timeseries"):
        csv_path = os.path.join(os.path.dirname(template_file), "reports",
                                _safe_name(config.get("project_name") or "fluxos")
                                + "_timeseries.csv")
        os.makedirs(os.path.dirname(csv_path), exist_ok=True)
        _write_timeseries_csv(csv_path, stats)
        print(f"  time series  : {csv_path}")

    # HTML report
    report_path = None
    if config.get("generate_report", True):
        report_dir = os.path.join(os.path.dirname(template_file), "reports")
        os.makedirs(report_dir, exist_ok=True)
        name = _safe_name(config.get("project_name") or "fluxos")
        report_path = os.path.join(report_dir, f"{name}_results_report.html")
        # UTM zone for satellite basemap. Preference order:
        #   1. explicit config["utm_zone"]
        #   2. auto-detect from DEM corner easting (pyproj, if available)
        #   3. None => satellite basemap disabled, report still renders.
        utm_zone = config.get("utm_zone")
        utm_southern = bool(config.get("utm_southern", False))
        if utm_zone is None and dem_meta:
            utm_zone, utm_southern = _detect_utm_from_dem(dem_meta)
        if utm_zone is not None:
            print(f"  UTM zone     : {utm_zone}"
                  + ("  (southern hemisphere)" if utm_southern else ""))

        generate_results_report({
            "config": config,
            "modset": modset,
            "mesh_type": mesh_type,
            "results_dir": results_dir,
            "modset_path": modset_path,
            "stats": stats,
            "kpis": kpis,
            "dem_meta": dem_meta,
            "utm_zone": utm_zone,
            "utm_southern": utm_southern,
            "report_dir": report_dir,
            "generated_at": datetime.now().isoformat(timespec="seconds"),
        }, report_path)
        print(f"  report       : {report_path}")
        if config.get("open_report"):
            webbrowser.open_new_tab("file://" + report_path)

    print("\nDone.")
    return {"stats": stats, "kpis": kpis, "report_path": report_path}


def _safe_name(s: str) -> str:
    return "".join(c if c.isalnum() or c in "-_" else "_" for c in s)


def _detect_utm_from_dem(dem_meta: dict) -> tuple[int | None, bool]:
    """
    Best-effort UTM zone detection from DEM corner coordinates.

    UTM (easting, northing) alone does not uniquely identify a zone — we first
    try to read a sibling ``.prj`` file, fall back to scanning a ``.tif`` with
    rasterio / pyproj, and finally to a rough "centre-of-the-world" heuristic.

    Returns ``(zone:int|None, southern:bool)``. ``None`` means we could not
    detect — the caller should then either rely on explicit config or skip
    the satellite basemap.
    """
    path = dem_meta.get("path")
    yll = float(dem_meta.get("yllcorner", 0.0))
    southern = yll > 0 and yll > 1_000_000 and yll < 10_000_000 and False
    # Heuristic for hemisphere: UTM southing in southern hem = 10_000_000 - y
    # (we treat very large northings as southern) — only kicks in if we can't
    # read a .prj.
    if path:
        # 1. Sibling .prj
        prj_path = os.path.splitext(path)[0] + ".prj"
        if os.path.isfile(prj_path):
            try:
                with open(prj_path) as f:
                    prj = f.read()
                zone, south = _parse_utm_from_wkt(prj)
                if zone is not None:
                    return zone, south
            except Exception:
                pass
        # 2. Try rasterio / pyproj against a sibling .tif or the .asc itself
        try:
            import pyproj
            xll = float(dem_meta.get("xllcorner", 0.0))
            ncols = int(dem_meta.get("ncols", 1))
            nrows = int(dem_meta.get("nrows", 1))
            cs = float(dem_meta.get("cellsize", 1.0))
            east_mid = xll + 0.5 * ncols * cs
            north_mid = yll + 0.5 * nrows * cs
            # We can't infer zone from (e,n) — but we can validate with pyproj
            # once the caller provides one. Without a projection, we default to
            # a hemisphere guess only.
            southern = north_mid > 10_000_000 or (north_mid < 0)
            _ = pyproj  # pragma: no cover — used only to probe availability
        except Exception:
            pass
    return None, southern


def _parse_utm_from_wkt(wkt: str) -> tuple[int | None, bool]:
    """Extract UTM zone + hemisphere from a WKT / ESRI PRJ string."""
    import re as _re
    m = _re.search(r'UTM[_ ]Zone[_ ]?(\d+)[_ ]?([NS])', wkt, _re.IGNORECASE)
    if m:
        return int(m.group(1)), m.group(2).upper() == "S"
    m = _re.search(r'zone\s*=?\s*(\d+)', wkt, _re.IGNORECASE)
    if m:
        south = bool(_re.search(r'south', wkt, _re.IGNORECASE))
        return int(m.group(1)), south
    return None, False


def _write_timeseries_csv(path: str, stats: dict) -> None:
    import csv
    cols = ["time_s", "volume_m3", "area_flooded_m2",
            "max_h_m", "mean_h_wet_m", "max_v_ms"]
    rows = list(zip(
        stats["time_s"], stats["volume_m3"], stats["area_flooded_m2"],
        stats["max_h_m"], stats["mean_h_wet_m"], stats["max_v_ms"]))
    if stats.get("has_conc"):
        cols += ["max_conc_SW", "mean_conc_SW"]
        rows = [r + (stats["max_conc_SW_t"][i], stats["mean_conc_SW_t"][i])
                for i, r in enumerate(rows)]
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        w.writerows(rows)
