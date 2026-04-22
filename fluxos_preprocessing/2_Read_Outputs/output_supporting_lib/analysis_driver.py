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
    # <repo>/fluxos_preprocessing/2_Read_Outputs/read_output_template.py
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

    # For regular mesh we need DEM dimensions (to linearise (irow,icol) indices)
    cellsize = None
    ncols = nrows = None
    if mesh_type != "triangular":
        dem_path = _abs_in_repo(repo_root, modset.get("DEM_FILE", ""))
        hdr = read_asc_header(dem_path)
        cellsize = float(hdr.get("cellsize", 1.0))
        ncols = int(hdr.get("ncols", 1))
        nrows = int(hdr.get("nrows", 1))
        print(f"  DEM          : {dem_path} "
              f"({ncols}×{nrows} @ {cellsize} m)")

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
        generate_results_report({
            "config": config,
            "modset": modset,
            "mesh_type": mesh_type,
            "results_dir": results_dir,
            "modset_path": modset_path,
            "stats": stats,
            "kpis": kpis,
            "generated_at": datetime.now().isoformat(timespec="seconds"),
        }, report_path)
        print(f"  report       : {report_path}")
        if config.get("open_report"):
            webbrowser.open_new_tab("file://" + report_path)

    print("\nDone.")
    return {"stats": stats, "kpis": kpis, "report_path": report_path}


def _safe_name(s: str) -> str:
    return "".join(c if c.isalnum() or c in "-_" else "_" for c in s)


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
