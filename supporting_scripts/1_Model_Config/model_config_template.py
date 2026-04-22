#!/usr/bin/env python3
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
FLUXOS model configuration template
====================================

HOW TO USE
----------
1. Copy this file (or edit in place for a one-off) and fill in the sections
   below.
2. Run it from any working directory:

       python model_config_template.py

3. The script will:
     - downscale your GeoTIFF DEM to an ESRI-ASCII .asc file in `<repo>/bin/`
     - generate a Gmsh .msh file in `<repo>/bin/` (if `mesh_type = "triangular"`)
     - write a `modset_*.json` config file in `<repo>/bin/`
     - generate an HTML report under `1_Model_Config/reports/` with the
       summary of what was built and copy-paste Docker commands to run
       the simulation.

Everything below is read as a plain Python dict — no YAML, no env vars.
Only this file should need editing for a new project.
"""

_config = dict(

    # ------------------------------------------------------------------
    # 1. Project metadata
    # ------------------------------------------------------------------
    project_name = "Rosa Creek",
    authors      = ["Diogo Costa"],
    date         = "2026-04-21",
    description  = (
        "Snowmelt-driven flood simulation over Rosa Creek DEM."
    ),

    # ------------------------------------------------------------------
    # 2. Paths
    # ------------------------------------------------------------------
    # repo_root: absolute path to the FLUXOS checkout. `None` means
    #   "auto-detect" (assumed to be two directories above this file, i.e.
    #   the parent of `supporting_scripts/`).
    repo_root       = None,

    # Directory (relative to repo_root) where generated .asc / .msh /
    # modset.json files will be written. This is the self-contained
    # "working example" directory — the Docker container mounts it
    # read-only at /work/Working_example, and modsets reference files
    # in here with paths like "Working_example/Rosa_2m.asc".
    output_bin_dir  = "Working_example",

    # Final name of the modset JSON inside output_bin_dir.
    modset_name     = "modset_rosa.json",

    # ------------------------------------------------------------------
    # 3. DEM source
    # ------------------------------------------------------------------
    # Two modes are supported:
    #
    #   "file"      — point `dem_source_geotiff` at an existing GeoTIFF
    #                 you already have on disk (default, traditional path).
    #
    #   "download"  — fetch a DEM automatically for a lat/lon bounding box
    #                 from an external provider, reproject to UTM, then
    #                 run it through the normal pipeline. No manual GIS
    #                 work needed for a first-pass floodplain study.
    dem_source_mode = "file",

    # --- Only used when dem_source_mode == "file" ---------------------
    # Path to your GeoTIFF DEM. Can be absolute or relative to repo_root.
    # Must be in a projected CRS (UTM / equivalent) — geographic (degrees)
    # DEMs will emit a warning and produce incorrect slopes.
    dem_source_geotiff      = "Working_example/Rosa_2m.tif",

    # --- Only used when dem_source_mode == "download" -----------------
    # Bounding box in WGS84 lon/lat: (lon_min, lat_min, lon_max, lat_max).
    # i.e. (upper-left longitude, lower-right latitude,
    #       lower-right longitude, upper-left latitude)
    # Pick a small bbox for a first run — 10×10 km is ≈ 0.1° × 0.1°.
    download_bbox_wgs84 = (-124.10, 52.20, -124.05, 52.23),

    # Which DEM product to fetch.
    # -- OpenTopography (global, free API key required): -----------------
    #    "SRTMGL1"   — SRTM 30 m (NASA, lat ±60°)   vertical RMSE ≈ 16 m
    #    "SRTMGL3"   — SRTM 90 m (NASA, lat ±60°)   vertical RMSE ≈ 16 m
    #    "COP30"     — Copernicus GLO-30 (ESA)      vertical RMSE ≈ 4 m   ★ recommended global default
    #    "COP90"     — Copernicus GLO-90 (ESA)      vertical RMSE ≈ 4 m
    #    "AW3D30"    — ALOS AW3D30 (JAXA)           vertical RMSE ≈ 5 m
    #    "NASADEM"   — NASADEM (reprocessed SRTM)   vertical RMSE ≈ 10 m
    #    "EU_DTM"    — EU-DTM 10 m (Copernicus)     vertical RMSE ≈ 3 m   — Europe (EEA39) only
    #    "GEBCOIceTopo" — GEBCO ice surface         — bathymetry/land combo
    # -- USGS 3DEP (US only, no key, needs `pip install py3dep`): --------
    #    "USGS_10M"  — USGS 3DEP 10 m               vertical RMSE ≈ 1 m   — CONUS
    #    "USGS_3M"   — USGS 3DEP 3 m                vertical RMSE ≈ 0.5 m — partial CONUS
    #    "USGS_1M"   — USGS 3DEP 1 m LiDAR          vertical RMSE ≈ 0.15 m — LiDAR-covered only
    download_provider     = "COP30",

    # OpenTopography requires a free API key:
    #   1. Register at https://portal.opentopography.org/ (free)
    #   2. Put your key in an environment variable and reference its name here:
    #        export OPENTOPOGRAPHY_API_KEY=xxxxxxxxxxxxx
    download_api_key_env  = "OPENTOPOGRAPHY_API_KEY",

    # Downloaded GeoTIFFs are cached in this directory (skips re-downloading
    # on repeated runs). Relative to repo_root.
    download_cache_dir    = "bin/dem_cache",

    # Target CRS for reprojection of the downloaded DEM.
    #   "auto"  — UTM zone inferred from the bbox centre (recommended)
    #   int     — explicit EPSG code, e.g. 32610 for UTM 10N
    #   str     — any PROJ-parseable string, e.g. "EPSG:32629"
    download_target_crs   = "auto",
    download_force        = False,         # True → ignore cache, re-fetch

    # --- Applied to BOTH modes ---------------------------------------
    # Target cell size for the downscaled .asc (metres). Set equal to the
    # GeoTIFF native resolution (or `None`) to skip resampling. Beware: if
    # you set this below the native DEM resolution you are *interpolating*,
    # not gaining information — a 2 m grid built from a 30 m SRTM tile
    # looks smoother than it should.
    dem_target_resolution_m = 2.0,

    # Output .asc filename inside output_bin_dir.
    dem_output_asc          = "Rosa_2m.asc",

    # UTM zone of the DEM — surfaced in the HTML report and used by the
    # post-processing viewer. Not read by the FLUXOS binary itself.
    dem_utm_zone            = 10,

    # ------------------------------------------------------------------
    # 4. Mesh
    # ------------------------------------------------------------------
    # "regular"    -> fluxos uses the .asc directly as a Cartesian grid.
    # "triangular" -> fluxos uses a Gmsh .msh with DEM elevations baked
    #                 into the vertex z-coordinates.
    mesh_type               = "triangular",

    # Only used when mesh_type == "triangular":
    trimesh_min_size        = 2.0,      # finest triangle edge (m)
    trimesh_max_size        = 30.0,     # coarsest triangle edge (m)
    trimesh_slope_factor    = 2.0,      # higher => more refinement in steep terrain
    mesh_output_msh         = "Rosa_trimesh.msh",

    # ------------------------------------------------------------------
    # 5. Forcing files
    #
    # Paths are written into the modset as-is. FLUXOS interprets them
    # relative to its working directory (for the container workflow
    # that's /work, which maps to the repo root — so leave these as
    # `bin/<filename>`).
    # ------------------------------------------------------------------
    meteo_file   = "Working_example/Qmelt_synthetic.fluxos",
    inflow_file  = None,   # e.g. "Working_example/Flow_river.fluxos" — set to None to disable

    # ------------------------------------------------------------------
    # 6. Simulation settings
    # ------------------------------------------------------------------
    sim_datetime_start = "2009-01-01 00:00:00",
    roughness_height   = 0.005,   # meters — uniform Manning-equivalent

    # ------------------------------------------------------------------
    # 7. Output settings
    # ------------------------------------------------------------------
    output_folder    = "Results/",   # relative to FLUXOS working dir
    print_step_s     = 1800,          # seconds between .vtu snapshots
    h_min_to_print_m = 0.001,         # cells with h < this are not written

    # ------------------------------------------------------------------
    # 8. Optional modules
    # ------------------------------------------------------------------
    # ADE_TRANSPORT solves an advection-dispersion-reaction equation alongside
    # the shallow-water solver, tracking a scalar concentration field (e.g.
    # dissolved solute) through the flow. Output variable is `conc_SW`.
    # Set `enabled=False` to skip transport entirely.
    ade_transport      = dict(enabled=True, d_coef=0.5),
    soil_infiltration  = dict(enabled=False, default_ks_mm_hr=10.0),

    # ------------------------------------------------------------------
    # 9. Docker / run settings
    # (these affect the snippets shown in the report, not the files
    # written to disk)
    # ------------------------------------------------------------------
    use_mpi             = False,   # True -> report suggests USE_MPI=ON build
    mpi_np              = 4,       # only used if use_mpi=True
    use_trimesh_build   = True,    # True -> USE_TRIMESH=ON build arg

    # ------------------------------------------------------------------
    # 10. Report options
    # ------------------------------------------------------------------
    generate_report = True,
    open_report     = True,        # False to skip webbrowser.open()
)


# ============================================================================
# You should not need to edit anything below this line.
# ============================================================================

if __name__ == "__main__":
    import os, sys
    _here = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, os.path.join(_here, "config_support_lib"))
    from driver import run          # noqa: E402
    run(_config, template_file=os.path.abspath(__file__))
