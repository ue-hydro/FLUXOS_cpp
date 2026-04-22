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
Orchestrator for FLUXOS model configuration.

Reads the `_config` dict from `model_config_template.py`, runs the three
preprocessing steps (DEM → mesh → modset JSON), and hands the collected
metadata off to `gen_report.generate_report()` for the HTML summary.

Users should not edit this file. Edit `model_config_template.py` instead.
"""

from __future__ import annotations

import os
import sys
import json
import traceback
import webbrowser
from datetime import datetime


# --- Helpers ----------------------------------------------------------------

def _resolve_repo_root(config: dict, template_file: str) -> str:
    """
    Resolve the FLUXOS repo root, with `repo_root=None` meaning
    "two directories above the template file" (i.e. the parent of
    `supporting_scripts/`).
    """
    if config.get("repo_root"):
        return os.path.abspath(config["repo_root"])
    # template_file lives at <repo>/supporting_scripts/1_Model_Config/model_config_template.py
    return os.path.abspath(os.path.join(os.path.dirname(template_file),
                                        "..", ".."))


def _abs_in_repo(repo_root: str, maybe_relative: str) -> str:
    """Join `maybe_relative` to `repo_root` only if it's not already absolute."""
    if not maybe_relative:
        return maybe_relative
    return maybe_relative if os.path.isabs(maybe_relative) \
        else os.path.join(repo_root, maybe_relative)


# --- Steps ------------------------------------------------------------------

def _step_dem(config: dict, repo_root: str, bin_dir: str) -> dict:
    """Run DEM preprocessing. Returns a dict of captured metadata.

    Supports two source modes:
    - ``dem_source_mode = "file"`` — use the GeoTIFF at ``dem_source_geotiff``
      (default, keeps the original workflow intact).
    - ``dem_source_mode = "download"`` — fetch a GeoTIFF from an external
      provider (OpenTopography / USGS 3DEP) for ``download_bbox_wgs84``.
    """
    from dem_operations import (
        read_geotiff, downscale_dem, write_esri_ascii, write_dummy_asc,
    )

    src, download_meta = _resolve_dem_source(config, repo_root, bin_dir)
    out_asc = os.path.join(bin_dir, config["dem_output_asc"])

    print(f"\n[1/3] DEM — reading {src}")
    elevation, transform, crs, nodata = read_geotiff(src)
    src_res = abs(transform.a)
    src_rows, src_cols = elevation.shape

    target_res = config.get("dem_target_resolution_m")
    if target_res and abs(target_res - src_res) > 1e-9:
        print(f"      downscaling {src_res:.3f} m → {target_res:.3f} m")
        elevation, transform = downscale_dem(elevation, transform, target_res, nodata)
        cellsize = float(target_res)
    else:
        cellsize = float(src_res)

    nrows, ncols = elevation.shape
    xll = transform.c
    yll = transform.f + transform.e * nrows  # e is negative → lower-left y

    # For triangular meshes the DEM is only used by get_domain_size() on the
    # C++ side — a 10x10 dummy grid is enough and keeps on-disk size tiny.
    if config["mesh_type"] == "triangular":
        print(f"      writing DUMMY ASC (trimesh path): {out_asc}")
        write_dummy_asc(out_asc, xll, yll, cellsize, nodata_value=-9999)
    else:
        print(f"      writing {out_asc}")
        write_esri_ascii(out_asc, elevation, xll, yll, cellsize, nodata)

    valid = elevation[elevation != nodata]
    zmin = float(valid.min()) if valid.size else None
    zmax = float(valid.max()) if valid.size else None

    # Downsample the (post-resample, pre-dummy-ASC) elevation for the preview
    # map in the HTML report. Keeps JSON size bounded regardless of DEM size.
    preview = _dem_preview(elevation, xll, yll, cellsize, nodata,
                           target_max=200)

    return dict(
        source=src,
        output_path=out_asc,
        crs=str(crs) if crs is not None else "unknown",
        source_resolution_m=float(src_res),
        target_resolution_m=float(cellsize),
        source_rows=int(src_rows),
        source_cols=int(src_cols),
        output_rows=int(nrows),
        output_cols=int(ncols),
        xllcorner=float(xll),
        yllcorner=float(yll),
        bbox=[float(xll), float(yll),
              float(xll + ncols * cellsize),
              float(yll + nrows * cellsize)],
        elev_min=zmin,
        elev_max=zmax,
        dummy_asc_for_trimesh=(config["mesh_type"] == "triangular"),
        download=download_meta,
        preview=preview,
    )


def _dem_preview(elevation, xll: float, yll: float, cellsize: float,
                 nodata: float, target_max: int = 200) -> dict:
    """
    Return a downsampled DEM suitable for embedding in the HTML report.

    ``target_max`` bounds the longer side (rows or cols) so the JSON stays
    small for multi-megapixel DEMs. Row 0 of ``z`` is the northernmost row
    (matching ESRI ASCII); ``y_vals`` is correspondingly in descending
    order so Plotly's heatmap renders the map with north up.
    """
    import math
    import numpy as np
    arr = np.asarray(elevation)
    nrows, ncols = arr.shape
    long_side = max(nrows, ncols)
    stride = max(1, int(np.ceil(long_side / float(target_max))))
    sub = arr[::stride, ::stride].astype(float)
    new_rows, new_cols = sub.shape
    dx = cellsize * stride
    y_top = yll + nrows * cellsize        # northernmost Y of the full DEM

    # Build the 2-D list with Python ``None`` for nodata cells — this
    # serialises as JSON ``null``, which Plotly renders as a transparent
    # cell. ``NaN`` is NOT valid JSON and is rejected by JSON.parse in the
    # browser, so we must not let it leak into the encoded figure.
    z_list = []
    for row in sub:
        r = []
        for v in row:
            if math.isnan(v) or math.isclose(v, nodata):
                r.append(None)
            else:
                r.append(float(v))
        z_list.append(r)

    return {
        "z": z_list,
        "x_vals": [float(xll + dx * (i + 0.5)) for i in range(new_cols)],
        "y_vals": [float(y_top - dx * (i + 0.5)) for i in range(new_rows)],
        "rows": int(new_rows),
        "cols": int(new_cols),
        "stride": int(stride),
    }


def _resolve_dem_source(config: dict, repo_root: str, bin_dir: str
                        ) -> tuple[str, dict | None]:
    """
    Return ``(geotiff_path, download_metadata_or_None)`` for the DEM to use.
    When ``dem_source_mode == "download"`` we fetch + reproject before the
    existing pipeline ever sees a file.
    """
    mode = (config.get("dem_source_mode") or "file").lower()
    if mode == "file":
        return _abs_in_repo(repo_root, config["dem_source_geotiff"]), None
    if mode != "download":
        raise ValueError(f"dem_source_mode must be 'file' or 'download', got {mode!r}")

    from dem_download import fetch_dem

    provider = config["download_provider"]
    bbox = tuple(config["download_bbox_wgs84"])
    if len(bbox) != 4:
        raise ValueError("download_bbox_wgs84 must be (lon_min, lat_min, lon_max, lat_max)")

    cache_dir = _abs_in_repo(repo_root, config.get("download_cache_dir",
                                                   "bin/dem_cache"))
    key_env = config.get("download_api_key_env", "OPENTOPOGRAPHY_API_KEY")
    api_key = os.environ.get(key_env)

    target_crs = config.get("download_target_crs", "auto")
    print(f"      DEM download: {provider} bbox={bbox}")
    meta = fetch_dem(provider, bbox, cache_dir,
                     api_key=api_key, target_crs=target_crs,
                     force=bool(config.get("download_force", False)))
    print(f"      downloaded → {meta['output_path']}")
    print(f"      {meta['label']} · native ≈ {meta['native_res_m']} m"
          + (f" · vertical RMSE ≈ {meta['vertical_rmse_m']} m"
             if meta.get("vertical_rmse_m") else ""))
    return meta["output_path"], meta


def _step_mesh(config: dict, repo_root: str, bin_dir: str) -> dict | None:
    """Run mesh generation (triangular only). Returns captured metadata or None."""
    if config["mesh_type"] != "triangular":
        print("\n[2/3] Mesh — regular grid (skipping Gmsh step)")
        return None

    from dem_operations import read_geotiff
    from mesh_generation import (
        compute_slope, build_domain_boundary, create_size_field,
        generate_adaptive_mesh,
    )
    import numpy as np

    src = _abs_in_repo(repo_root, config["dem_source_geotiff"])
    out_msh = os.path.join(bin_dir, config["mesh_output_msh"])

    print(f"\n[2/3] Mesh — {out_msh}")
    elevation, transform, crs, nodata = read_geotiff(src)
    cellsize = abs(transform.a)

    slope = compute_slope(elevation, cellsize, nodata)
    size_field = create_size_field(
        slope,
        config["trimesh_min_size"],
        config["trimesh_max_size"],
        config["trimesh_slope_factor"],
    )
    boundary = build_domain_boundary(elevation, transform, nodata)

    generate_adaptive_mesh(
        boundary, size_field, transform, elevation, nodata,
        config["trimesh_min_size"], config["trimesh_max_size"], out_msh,
    )

    # Parse the produced .msh to capture cell/edge/vertex counts (same stats
    # FLUXOS prints at runtime — useful to surface in the report).
    stats = _parse_gmsh_stats(out_msh)

    # Extract vertex coords + triangle connectivity for the preview map.
    # Gracefully no-op if the mesh is huge — we downsample to ≤ 4k triangles
    # so the embedded JSON doesn't balloon.
    preview = _mesh_preview(out_msh, max_triangles=4000)

    stats.update(dict(
        output_path=out_msh,
        min_size_m=float(config["trimesh_min_size"]),
        max_size_m=float(config["trimesh_max_size"]),
        slope_factor=float(config["trimesh_slope_factor"]),
        max_slope=float(np.nanmax(slope)) if slope.size else 0.0,
        preview=preview,
    ))
    return stats


def _mesh_preview(msh_path: str, max_triangles: int = 4000) -> dict:
    """
    Extract vertex (x, y) arrays + triangle indices from a Gmsh 2.2 ASCII
    .msh file. If there are more than ``max_triangles`` triangles, a stride
    is applied to the triangle list so the embedded preview stays light.

    Returns ``{"x": [...], "y": [...], "tris": [[i0,i1,i2], ...], "stride": N}``
    with vertex indices referring to positions in ``x`` / ``y`` (only the
    vertices touched by the kept triangles are included, so the arrays are
    compact even after subsampling).
    """
    try:
        xs, ys, zs = [], [], []
        tris = []
        with open(msh_path) as f:
            mode = None
            remaining = 0
            for line in f:
                line = line.strip()
                if line == "$Nodes":
                    mode = "nodes_hdr"; continue
                if line == "$EndNodes":
                    mode = None; continue
                if line == "$Elements":
                    mode = "elem_hdr"; continue
                if line == "$EndElements":
                    mode = None; continue
                if mode == "nodes_hdr":
                    remaining = int(line.split()[0])
                    mode = "nodes"; continue
                if mode == "nodes":
                    parts = line.split()
                    if len(parts) >= 4:
                        # id x y z — id is 1-indexed; z carries the DEM elevation
                        xs.append(float(parts[1]))
                        ys.append(float(parts[2]))
                        zs.append(float(parts[3]))
                    continue
                if mode == "elem_hdr":
                    mode = "elements"; continue
                if mode == "elements":
                    parts = line.split()
                    if len(parts) >= 2 and int(parts[1]) == 2:
                        # 3-node triangle: node IDs are the last 3 fields (1-indexed)
                        i0, i1, i2 = (int(parts[-3]) - 1,
                                      int(parts[-2]) - 1,
                                      int(parts[-1]) - 1)
                        tris.append((i0, i1, i2))
                    continue
    except Exception:
        return {}

    if not tris or not xs:
        return {}

    stride = max(1, (len(tris) + max_triangles - 1) // max_triangles)
    kept = tris[::stride]
    # Compact to only the vertices touched by kept triangles
    seen = {}
    comp_x, comp_y, comp_z, comp_tris = [], [], [], []
    for (a, b, c) in kept:
        new_ids = []
        for v in (a, b, c):
            if v not in seen:
                seen[v] = len(comp_x)
                comp_x.append(xs[v])
                comp_y.append(ys[v])
                comp_z.append(zs[v] if v < len(zs) else 0.0)
            new_ids.append(seen[v])
        comp_tris.append(new_ids)
    return {
        "x": comp_x,
        "y": comp_y,
        "z": comp_z,
        "tris": comp_tris,
        "stride": stride,
        "total_triangles": len(tris),
        "shown_triangles": len(kept),
    }


def _parse_gmsh_stats(msh_path: str) -> dict:
    """Quick scan of a Gmsh 2.2 ASCII .msh to count nodes and triangles."""
    n_nodes, n_tri, n_line = 0, 0, 0
    try:
        with open(msh_path) as f:
            mode = None
            for line in f:
                line = line.strip()
                if line == "$Nodes":
                    mode = "nodes_hdr"; continue
                if line == "$EndNodes":
                    mode = None; continue
                if line == "$Elements":
                    mode = "elem_hdr"; continue
                if line == "$EndElements":
                    mode = None; continue
                if mode == "nodes_hdr":
                    n_nodes = int(line.split()[0]); mode = "nodes"; continue
                if mode == "elem_hdr":
                    mode = "elements"; continue
                if mode == "elements":
                    parts = line.split()
                    if len(parts) >= 2:
                        etype = int(parts[1])
                        if etype == 2:  # 3-node triangle
                            n_tri += 1
                        elif etype == 1:  # 2-node line (boundary)
                            n_line += 1
    except Exception:
        pass
    return dict(
        vertex_count=n_nodes,
        cell_count=n_tri,
        boundary_edge_count=n_line,
    )


def _step_config(config: dict, repo_root: str, bin_dir: str,
                 dem_meta: dict, mesh_meta: dict | None) -> dict:
    """Write modset.json. Returns metadata including the output path."""
    modset_name = config["modset_name"]
    out_path = os.path.join(bin_dir, modset_name)

    modset = {
        "COMMENT": f"Generated by model_config_template.py ({config['project_name']})",
        "DEM_FILE": _rel(dem_meta["output_path"], repo_root),
        "SIM_DATETIME_START": config["sim_datetime_start"],
        "RESTART": False,
        # Note: matches the existing typo in read_functions.cpp
        "ROUGNESS_HEIGHT": float(config["roughness_height"]),
        "OUTPUT": {
            "OUTPUT_FOLDER": config["output_folder"],
            "PRINT_STEP": int(config["print_step_s"]),
            "H_MIN_TO_PRINT": float(config["h_min_to_print_m"]),
        },
        "EXTERNAL_MODULES": {
            "ADE_TRANSPORT": {
                "STATUS": bool(config["ade_transport"]["enabled"]),
                "D_COEF": float(config["ade_transport"]["d_coef"]),
            },
        },
    }

    if config["mesh_type"] == "triangular" and mesh_meta is not None:
        modset["MESH_TYPE"] = "triangular"
        modset["MESH_FILE"] = _rel(mesh_meta["output_path"], repo_root)
        modset["MESH_FORMAT"] = "gmsh"
        modset["BOUNDARY_CONDITIONS"] = {"1": {"type": "wall"}}

    if config.get("meteo_file"):
        modset["METEO_FILE"] = config["meteo_file"]

    if config.get("inflow_file"):
        modset["INFLOW_FILE"] = {
            "FILENAME": config["inflow_file"],
            "DISCHARGE_LOCATION": {"X_COORDINATE": 0, "Y_COORDINATE": 0},
        }

    # Soil infiltration (optional — only written if user enables it)
    soil = config.get("soil_infiltration") or {}
    if soil.get("enabled"):
        modset["SOIL_INFILTRATION"] = {
            "STATUS": True,
            "DEFAULT_KS_MM_HR": float(soil.get("default_ks_mm_hr", 10.0)),
        }

    print(f"\n[3/3] Config — writing {out_path}")
    with open(out_path, "w") as f:
        json.dump(modset, f, indent=4)

    return dict(
        output_path=out_path,
        modset=modset,
    )


def _rel(path: str, root: str) -> str:
    """Return `path` relative to `root` with forward slashes."""
    try:
        return os.path.relpath(path, root).replace(os.sep, "/")
    except ValueError:
        return path


# --- Main entry -------------------------------------------------------------

def run(config: dict, template_file: str) -> dict:
    """
    Execute DEM → mesh → modset → report for the given template config.

    Parameters
    ----------
    config : dict
        The `_config` dict from `model_config_template.py`.
    template_file : str
        Absolute path of the template file (used to locate the repo root
        and the directory to write the report into).

    Returns
    -------
    dict
        Collected metadata (the same dict passed to gen_report).
    """
    repo_root = _resolve_repo_root(config, template_file)
    bin_dir = _abs_in_repo(repo_root, config.get("output_bin_dir", "bin"))
    os.makedirs(bin_dir, exist_ok=True)
    os.makedirs(os.path.join(repo_root, "Results"), exist_ok=True)

    print(f"FLUXOS model config — {config['project_name']}")
    print(f"  repo root : {repo_root}")
    print(f"  bin dir   : {bin_dir}")

    report_data = dict(
        config=config,
        repo_root=repo_root,
        bin_dir=bin_dir,
        generated_at=datetime.now().isoformat(timespec="seconds"),
        errors=[],
    )

    # Each step is independently guarded so a failure in one does not block the
    # report — the report surfaces what succeeded and what did not.
    try:
        report_data["dem"] = _step_dem(config, repo_root, bin_dir)
    except Exception as e:
        report_data["dem"] = None
        report_data["errors"].append(("DEM", str(e), traceback.format_exc()))
        print(f"ERROR in DEM step: {e}", file=sys.stderr)

    try:
        report_data["mesh"] = _step_mesh(config, repo_root, bin_dir)
    except Exception as e:
        report_data["mesh"] = None
        report_data["errors"].append(("Mesh", str(e), traceback.format_exc()))
        print(f"ERROR in mesh step: {e}", file=sys.stderr)

    try:
        if report_data["dem"] is not None:
            report_data["config_out"] = _step_config(
                config, repo_root, bin_dir,
                report_data["dem"], report_data.get("mesh"),
            )
        else:
            report_data["config_out"] = None
            report_data["errors"].append(
                ("Config", "DEM step failed — skipping modset generation", ""))
    except Exception as e:
        report_data["config_out"] = None
        report_data["errors"].append(("Config", str(e), traceback.format_exc()))
        print(f"ERROR in config step: {e}", file=sys.stderr)

    # Report
    if config.get("generate_report", True):
        from gen_report import generate_report
        report_dir = os.path.join(os.path.dirname(template_file), "reports")
        os.makedirs(report_dir, exist_ok=True)
        safe_name = "".join(c if c.isalnum() or c in "-_" else "_"
                            for c in config["project_name"])
        report_path = os.path.join(report_dir, f"{safe_name}_report.html")
        generate_report(report_data, report_path)
        report_data["report_path"] = report_path
        print(f"\nReport written to: {report_path}")

        if config.get("open_report", False):
            webbrowser.open_new_tab("file://" + report_path)

    print("\nDone.")
    return report_data
