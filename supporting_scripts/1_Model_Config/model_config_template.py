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
1. Copy this file to a project-specific name (recommended) or edit in place
   for a one-off:

       cp model_config_template.py model_config_<project>.py

2. Fill in the ``_config = dict(...)`` block below. Every field has inline
   documentation describing the valid values and how to pick one.

3. Run it from any working directory:

       python model_config_<project>.py

4. The script will:
     - fetch or read the DEM (GeoTIFF) and downscale to an ESRI-ASCII ``.asc``
     - generate a Gmsh ``.msh`` (if ``mesh_type = "triangular"``) with DEM
       elevations baked into the vertex z-coordinates
     - write a ``modset_*.json`` FLUXOS configuration file
     - generate an HTML report under ``1_Model_Config/reports/`` summarising
       the build and offering copy-paste Docker / Apptainer commands to
       run the simulation and visualise the results

Everything below is read as a plain Python dict — no YAML, no env vars.
Only this file should need editing for a new project.

COORDINATES: WGS84 IN, UTM INVISIBLE
------------------------------------
Every location you specify in this config is in WGS84 (decimal degrees):

  • ``download_bbox_wgs84``                       — (lon_min, lat_min, lon_max, lat_max)
  • ``inflow_file = dict(..., lon=..., lat=...)`` — point inflow location
  • the HTML report shows lon/lat on hover in every map

Internally FLUXOS is a physics solver: the shallow-water equations need
distances in METRES (a 1° step at lat 0° is 111 km but only ~2 km at lat 89°,
so a solver working in degrees would be dimensionally wrong). So the DEM
you provide — or download — is reprojected to a metric CRS (auto-UTM by
default) before the solver ever sees it. You do not need to know your UTM
zone; the pipeline derives it from the DEM bbox.

To sum up: **users speak WGS84; the solver speaks UTM; the pipeline
handles the conversion**.

EXTERNAL RESOURCES
------------------
  • OpenTopography (free API keys, required for most OpenTopography providers):
      https://portal.opentopography.org/
  • Copernicus GLO-30 on AWS (no key, the default provider below):
      https://registry.opendata.aws/copernicus-dem/
  • USGS 3DEP (US only, no key):
      https://www.usgs.gov/3d-elevation-program
  • UTM-zone finder (plug in your lat/lon):
      https://www.dmap.co.uk/utmworld.htm
  • Manning's n / roughness-height typical values:
      https://www.engineeringtoolbox.com/mannings-roughness-d_799.html
"""

_config = dict(

    # ==================================================================
    # 1. Project metadata — shown in the HTML report header.
    # ==================================================================

    # Short human-readable name for the project. Used in the report title,
    # the generated filenames (slug-cased), and the simulation log.
    # Examples: "Rosa Creek", "Torres Vedras", "Thames Estuary 2026".
    project_name = "Rosa Creek 30h",

    # List of authors / modellers responsible for the run. Free-form.
    # Example: ["Alice Wong", "Bob Ng"]
    authors      = ["Diogo Costa"],

    # ISO-8601 date the configuration was prepared (YYYY-MM-DD).
    # Recorded in the report header as provenance.
    date         = "2026-04-21",

    # One- or two-sentence description — what the simulation is trying to
    # capture (event, purpose, key modelling assumptions). Shown at the
    # top of the report AND used verbatim as the modset's ``COMMENT``
    # field, so keep it short and informative.
    description  = (
        "30-hour multi-pulse flood simulation with sediment transport"
    ),

    # ==================================================================
    # 2. Paths — where generated files land relative to your checkout.
    # ==================================================================

    # Absolute path to the FLUXOS repository root. Two options:
    #   None  — auto-detect (assumed to be two directories above this file,
    #           i.e. the parent of ``supporting_scripts/``). Works for the
    #           in-repo workflow.
    #   str   — an absolute path, e.g. "/home/alice/models/FLUXOS_cpp".
    #           Useful when running this script from a different checkout.
    repo_root       = None,

    # Directory (relative to repo_root) where generated files (.asc, .msh,
    # modset.json) will be written. Common choices:
    #   "Working_example" — the default "canonical inputs" folder shipped in
    #                       the repo. Docker / Apptainer workflows bind-mount
    #                       this at /work/Working_example, so modsets end up
    #                       referencing files as "Working_example/<name>".
    #   "my_project/"     — a separate subfolder if you want project-specific
    #                       inputs alongside the repo.
    output_bin_dir  = "Working_example",

    # Final name of the modset .json inside ``output_bin_dir``.
    # Pattern suggestion: ``modset_<project>.json``.
    modset_name     = "modset_river_30h.json",

    # ==================================================================
    # 3. DEM source — where the elevation model comes from.
    # ==================================================================

    # How to obtain the DEM:
    #   "file"      → point ``dem_source_geotiff`` below at a GeoTIFF on disk.
    #                 Recommended when you have local LiDAR / high-res data,
    #                 or when you are iterating and want fully offline runs.
    #   "download"  → fetch a tile automatically for ``download_bbox_wgs84``
    #                 from one of the providers below (OpenTopography, USGS
    #                 3DEP, or the AWS Copernicus mirror). Results are cached
    #                 so re-runs skip the network. No GIS setup required.
    dem_source_mode = "file",

    # --- Only used when dem_source_mode == "file" ---------------------
    # Absolute or repo-relative path to your GeoTIFF DEM.
    #
    # IMPORTANT: must be in a projected CRS (UTM / State Plane / equivalent)
    # where pixel spacing is in METRES. Geographic (degrees) DEMs will emit
    # a warning and produce wildly wrong slopes because FLUXOS computes
    # gradients as Δz/Δx with Δx in metres.
    #
    # If your DEM is in degrees, reproject first with:
    #   gdalwarp -t_srs EPSG:326XX source.tif utm_source.tif
    # where 326XX is your UTM zone (e.g. 32629 for Portugal, 32610 for
    # British Columbia).
    dem_source_geotiff      = "Working_example/Rosa_2m.tif",

    # --- Only used when dem_source_mode == "download" -----------------
    # Bounding box in WGS84 lon/lat (decimal degrees), four numbers:
    #     (lon_min, lat_min, lon_max, lat_max)
    # i.e. (southwest lon, southwest lat, northeast lon, northeast lat).
    #
    # Tip on sizing: at mid-latitudes 0.01° ≈ 1.1 km (latitude) and
    #                0.01° × cos(lat) km (longitude).
    # A ~10 × 10 km first-pass domain is ≈ 0.1° × 0.1°.
    # Pick a box large enough to contain the channel + floodplain + a buffer.
    #
    # Use Google Maps or the bounding-box picker at
    #   http://bboxfinder.com/
    # to dial in coordinates.
    download_bbox_wgs84 = (-124.10, 52.20, -124.05, 52.23),

    # Which DEM product to fetch. Grouped by whether an API key is required.
    #
    # === NO KEY REQUIRED (recommended starting point) ======================
    # "COP30_AWS"  — Copernicus GLO-30 (ESA) from the public AWS S3 bucket.
    #                30 m native, vertical RMSE ≈ 4 m, GLOBAL coverage.    ★
    #                Identical data to OpenTopography's "COP30" but without
    #                the key hassle.
    #
    # === OpenTopography (free API key required, richer catalogue) ==========
    #                1. Register at https://portal.opentopography.org/ (free)
    #                2. ``export OPENTOPOGRAPHY_API_KEY=xxxxxxxxxxxxx``
    #                3. Set ``download_api_key_env`` below to the env-var name.
    #
    # "SRTMGL1"    — SRTM 30 m (NASA)         RMSE ≈ 16 m   global ±60° lat
    # "SRTMGL3"    — SRTM 90 m (NASA)         RMSE ≈ 16 m   global ±60° lat
    # "COP30"      — Copernicus GLO-30 (ESA)  RMSE ≈ 4 m    GLOBAL
    # "COP90"      — Copernicus GLO-90 (ESA)  RMSE ≈ 4 m    GLOBAL
    # "AW3D30"     — ALOS AW3D30 (JAXA)       RMSE ≈ 5 m    GLOBAL
    # "NASADEM"    — NASADEM (SRTM re-processed)  RMSE ≈ 10 m  ±60° lat
    # "EU_DTM"     — EU-DTM 10 m (Copernicus LS)  RMSE ≈ 3 m   ★ Europe only (EEA39)
    # "GEBCOIceTopo" — GEBCO ice surface                      land + ocean
    #
    # === USGS 3DEP (US only, no key — needs ``pip install py3dep``) ========
    # "USGS_10M"   — 3DEP 10 m             RMSE ≈ 1 m      CONUS
    # "USGS_3M"    — 3DEP 3 m              RMSE ≈ 0.5 m    CONUS (partial)
    # "USGS_1M"    — 3DEP 1 m LiDAR        RMSE ≈ 0.15 m   US LiDAR-covered only ★
    #
    # Rules of thumb — pick the highest-resolution product covering your site:
    #   • USA  → USGS_1M  (if LiDAR),  else USGS_3M  / USGS_10M
    #   • Europe → EU_DTM (10 m)
    #   • elsewhere → COP30_AWS (no key) or COP30 (with key)
    download_provider     = "COP30_AWS",

    # Name of the environment variable that holds your OpenTopography API key.
    # Only consulted when ``download_provider`` is an OpenTopography product.
    # The key itself is NEVER written to disk — only its variable name is.
    # To set the key before running:
    #   bash / zsh:   export OPENTOPOGRAPHY_API_KEY=xxxxxxxxxxxxx
    #   fish:         set -x OPENTOPOGRAPHY_API_KEY xxxxxxxxxxxxx
    #   Windows CMD:  set OPENTOPOGRAPHY_API_KEY=xxxxxxxxxxxxx
    download_api_key_env  = "OPENTOPOGRAPHY_API_KEY",

    # Directory (relative to repo_root) where downloaded GeoTIFFs are cached.
    # Re-runs with the same provider + bbox hit the cache and skip the
    # network. Safe to delete any time to force re-downloads.
    download_cache_dir    = "bin/dem_cache",

    # Target CRS for the reprojection step. Three forms accepted:
    #   "auto"  — pick UTM zone from the bbox centre (recommended; always
    #             gives a local metric grid).
    #   int     — an explicit EPSG code, e.g. 32610 for UTM 10N (BC, Canada),
    #             32629 for UTM 29N (mainland Portugal), 25830 for ETRS89
    #             UTM 30N (Iberian Peninsula).
    #   str     — any PROJ-parseable string, e.g. "EPSG:32629" or a +proj=…
    #             WKT. Useful for national grids (BNG, Lambert-93, …).
    # Find your UTM zone by lat/lon at https://www.dmap.co.uk/utmworld.htm
    download_target_crs   = "auto",

    # Force-refresh the DEM cache:
    #   False — use cached file if present (fast, default).
    #   True  — re-download even if the tile is already on disk.
    download_force        = False,

    # --- Applied to BOTH "file" and "download" modes ------------------
    # Target cell size in METRES for the output .asc that FLUXOS reads.
    #   None / 0      — keep native resolution, no resampling.
    #   <native_res>  — BILINEAR downscaling (lowers resolution). Useful
    #                   for speed when native is overkill.
    #   >native_res   — bilinear up-sampling. WARNING: this is smoothing,
    #                   not added information. Setting e.g. 2 m when your
    #                   source is 30 m SRTM gives you a smooth 2 m grid
    #                   carrying no new detail.
    # Typical choices:
    #   • 1-2 m for urban / street-level flood modelling (LiDAR input)
    #   • 5-10 m for basin-scale floodplain studies
    #   • 30 m for regional / continental screening studies
    dem_target_resolution_m = 2.0,

    # Output .asc filename inside ``output_bin_dir``. Convention:
    #   ``<project>_<res>m.asc``  (e.g. "Rosa_2m.asc", "TorresVedras_30m.asc")
    dem_output_asc          = "Rosa_2m.asc",

    # You DO NOT need a ``dem_utm_zone`` field here any more. UTM is used
    # only internally by the FLUXOS solver (which needs metric distances
    # for physics) and the pipeline auto-derives the zone from the DEM's
    # CRS. Your config can stay fully in WGS84 lon/lat — the bbox, the
    # inflow location, and the report's coordinate tooltips are all in
    # decimal degrees.
    #
    # Only set this to override the auto-derived value — e.g. if you
    # deliberately want a non-UTM metric CRS (ETRS89 / BNG / Lambert-93).
    # Accept an int zone number or leave it out entirely.

    # ==================================================================
    # 4. Mesh — regular Cartesian grid vs. unstructured triangles.
    # ==================================================================
    # Two options:
    #
    #   "regular"    — FLUXOS consumes the .asc directly as a Cartesian
    #                  grid. Cells = ncols × nrows of the .asc. Simplest
    #                  path; works best when the DEM resolution is
    #                  appropriate over the whole domain and the boundary
    #                  is roughly rectangular. No Gmsh step required.
    #
    #   "triangular" — FLUXOS reads a Gmsh .msh with DEM elevations baked
    #                  into vertex z-coordinates. Generated from the DEM
    #                  by this template (slope-adaptive size field).
    #                  Best when:
    #                    • the domain has an irregular boundary
    #                    • terrain mixes steep and flat regions (the
    #                      adaptive grader puts more triangles where
    #                      gradients matter)
    #                    • cell count with a regular grid would be
    #                      wasteful (millions of cells, most of them
    #                      identical in flat floodplains)
    #                  Produces ~5-10× fewer cells than an equivalent
    #                  regular grid with similar accuracy near channels.
    #
    # Same binary handles both — FLUXOS dispatches at runtime based on
    # this field. No recompile needed to switch.
    mesh_type               = "regular",

    # --- Only read when mesh_type == "triangular" ---------------------
    # (These fields are harmless for regular-mesh runs — driver.py
    #  skips the Gmsh step and the modset omits MESH_FILE.)
    #
    # Minimum triangle edge length (metres). Controls the finest detail.
    # Rule: do NOT set below the DEM cell size — that would pretend the
    # mesh has more information than the DEM actually carries. Typical:
    #   1-5 m   — urban / street-level (LiDAR input)
    #   10-30 m — basin-scale flood studies
    #   30 m+   — regional / continental screening
    trimesh_min_size        = 2.0,

    # Maximum triangle edge length (metres). Controls the coarsest element
    # in flat areas. Typical ratio ``max:min`` is 3-10× :
    #   larger  → fewer cells, faster simulation, coarser far from channels
    #   smaller → more uniform grading, cell count closer to regular-mesh
    trimesh_max_size        = 30.0,

    # Slope-based refinement multiplier — more triangles where terrain is
    # steep (where water moves fast and gradients matter).
    #   1.0 — no slope sensitivity; uniform grading between min/max
    #   2.0 — moderate refinement in steep areas  (recommended default)
    #   4.0 — aggressive refinement; for mountainous / stepped terrain
    trimesh_slope_factor    = 2.0,

    # Output Gmsh .msh filename inside ``output_bin_dir``.
    # Convention: ``<project>_trimesh.msh``.
    mesh_output_msh         = "Rosa_trimesh.msh",

    # Default boundary condition for the triangular mesh perimeter.
    #   "outflow" — transmissive (zero-gradient). Water leaves the domain
    #               freely. Matches the regular-mesh default behaviour;
    #               use this for open-ended river reaches or any domain
    #               where water is expected to flow through.
    #   "wall"    — reflective (zero normal flux). Water bounces off the
    #               domain edges and accumulates. Use for closed basins
    #               or when the domain perimeter is a physical barrier
    #               (e.g. a fully-enclosed reservoir).
    # Leaving this unset defaults to "outflow" — the apples-to-apples
    # match with ``mesh_type = "regular"``.
    trimesh_boundary_default = "outflow",

    # ==================================================================
    # 5. Forcing files — precipitation / snowmelt + optional point inflow.
    # ==================================================================
    # Paths are written into the modset exactly as given; FLUXOS resolves
    # them relative to its CWD (for the container workflow that is /work,
    # which maps to the repo root — so relative paths starting with
    # ``Working_example/…`` work unchanged in Docker / Apptainer / native).
    #
    # --- .fluxos forcing file format (both meteo and inflow) --------------
    # Plain CSV, no header, one row per time point. Columns in order:
    #
    #   column 1   time           [s]           seconds from sim start
    #   column 2   rate / flow    [see below]   value at that time
    #   column 3+  concentration  [mg/L]        one column per chemical
    #                                           species (only read when
    #                                           ADE_TRANSPORT is enabled)
    #
    # Units of column 2 depend on the file role:
    #
    #   METEO  (whole-domain precipitation / snowmelt, applied uniformly)
    #       → column 2 in **mm / day**
    #         (solver converts via  rate / 86400000  [m/s])
    #
    #   INFLOW (point source at the X/Y coordinates set in this config)
    #       → column 2 in **m³ / s**
    #         (solver converts via  rate × dt / cellarea  [m of water])
    #
    # The first row is conventionally ``0,0,0`` (t=0, no flow/rate, no
    # concentration) followed by a duplicate ``0,0,0`` — FLUXOS does a
    # linear interpolation between rows, so the duplicate guarantees a
    # flat baseline before the first real value.
    #
    # Example meteo file (``Qmelt_synthetic.fluxos``):
    #     0,0,0            ← t=0 s, rate=0 mm/day, conc=0 mg/L
    #     900,500,0        ← t=900 s (15 min), rate=500 mm/day
    #     3600,0,0         ← rainfall off at t=1 h
    #
    # Example inflow file (``Flow_river.fluxos``):
    #     0,0,0            ← dry start
    #     600,15.0,0       ← 15 m³/s point inflow by t=10 min
    #     21600,15.0,0     ← still 15 m³/s at t=6 h

    # Meteo forcing — whole-domain precipitation / snowmelt time series.
    # Path is resolved relative to FLUXOS's CWD. See the column-format
    # notes above; column 2 must be in mm/day.
    # Set to "" or None to disable (no atmospheric forcing).
    # Rosa Creek default: no atmospheric forcing — the flood is driven
    # by the point inflow below, so METEO_FILE stays empty. Swap in
    # e.g. "Working_example/Qmelt_synthetic.fluxos" to add rainfall.
    meteo_file   = None,

    # OPTIONAL point inflow — e.g. a river or culvert entering the domain.
    # Three accepted forms (pick one):
    #
    #   (1) Dict with WGS84 lon/lat — ★ RECOMMENDED. Same CRS as
    #       `download_bbox_wgs84`, no need to know your UTM zone. The
    #       driver reprojects these to the DEM's projected CRS
    #       automatically using pyproj.
    #           inflow_file = dict(
    #               path = "Working_example/Flow_river.fluxos",
    #               lon  =  -9.25176,   # degrees E (WGS84)
    #               lat  =  39.09674,   # degrees N (WGS84)
    #           )
    #
    #   (2) Dict with projected UTM coords — use when you already know
    #       the easting/northing in the DEM's CRS (e.g. from a QGIS
    #       click-and-query):
    #           inflow_file = dict(
    #               path  = "Working_example/Flow_river.fluxos",
    #               x_utm = 478229.0,   # Easting  [m]  (UTM zone of DEM)
    #               y_utm = 4327542.0,  # Northing [m]  (UTM zone of DEM)
    #           )
    #
    #   (3) Plain string form — legacy. Discharge is applied at (0, 0)
    #       which is almost never what you want. Prefer (1) or (2).
    #           inflow_file = "Working_example/Flow_river.fluxos"
    #
    # The .fluxos file's column 2 must be in **m³/s**. See the format
    # notes above for the full column schema.
    # Set to None to disable.
    #
    # Rosa Creek default: a 30-hour multi-pulse river hydrograph at the
    # valley's upstream inlet. The .fluxos forcing ships with the repo
    # at the path below. We specify the inlet as integer easting/northing
    # in the DEM's UTM zone (10N) so the generated modset matches the
    # hand-tuned reference byte-for-byte; swap in ``lon=..., lat=...``
    # if you prefer WGS84 and don't care about integer precision.
    inflow_file  = dict(
        path  = "Working_example/Flow_river_30h.fluxos",
        x_utm = 426107,   # UTM 10N easting  [m]
        y_utm = 5785362,  # UTM 10N northing [m]
    ),

    # ==================================================================
    # 6. Simulation settings.
    # ==================================================================

    # Start timestamp of the simulation, in ISO-8601 ("YYYY-MM-DD HH:MM:SS").
    # Used (a) in the HTML report header; (b) by the KML / WebGL viewer
    # to label frames with real-world time; (c) to align forcing file time
    # offsets. The simulation advances in seconds from this moment.
    sim_datetime_start = "2009-01-01 00:00:00",

    # Uniform bed-roughness HEIGHT in metres (NOT Manning's n, NOT Chézy).
    # FLUXOS converts this internally to a friction coefficient via a
    # log-law closure. Typical values by land cover:
    #   0.001 m — smooth surfaces (concrete, asphalt, open water)
    #   0.005 m — bare soil, short grass, paved streets with debris
    #   0.010 m — short vegetation, urban mixed surfaces
    #   0.030 m — shrubs, dense grass, gravel
    #   0.100 m — trees, tall vegetation, rough rubble
    # If unsure, 0.005 is a reasonable urban-flood starting point.
    # See the external link at the top of this file for a roughness
    # reference table.
    roughness_height   = 0.005,

    # ==================================================================
    # 7. Output settings — what FLUXOS writes to disk.
    # ==================================================================

    # Directory (relative to FLUXOS's working dir) where per-timestep output
    # files land. Regular mesh → ``<t_sec>.txt`` per print step.
    # Triangular → ``<t_sec>.vtu`` per print step (plus a .pvd time series
    # index for ParaView / VisIt).
    # Convention: ``Results_<project>/`` keeps each project isolated.
    output_folder    = "Results/",

    # Interval (in SECONDS) between snapshot writes. Trade-off:
    #   • smaller → smoother animations, bigger disk footprint
    #   • larger  → fewer files, coarser time resolution
    # Practical values:
    #   60     — 1 min,   fine-grained flash-flood capture
    #   900    — 15 min,  normal urban flood runs
    #   1800   — 30 min,  default
    #   3600   — 1 hr,    long-duration storm runs
    print_step_s     = 300,

    # Minimum water depth (in METRES) to write a cell to the output.
    # Cells with h below this threshold are omitted (they are "dry").
    # Saves huge amounts of disk on regular-mesh outputs for
    # largely-dry domains. Typical:
    #   0.0001 — keep everything ("debug" mode)
    #   0.001  — default; filters numerical trace depths
    #   0.01   — aggressive; only reports cells with visible water
    h_min_to_print_m = 0.001,

    # ==================================================================
    # 8. Optional modules — swappable bolt-ons to the shallow-water solver.
    # ==================================================================

    # ADE transport — Advection-Dispersion-Reaction for a scalar tracer
    # (dissolved solute concentration travelling with the flow). Output
    # variable in the .vtu / .txt is ``conc_SW``.
    #
    #   enabled  — True to activate, False to skip entirely.
    #   d_coef   — dispersion coefficient in m²/s. Controls numerical
    #              smearing vs. sharpness of the front.
    #                0.0   — pure advection, sharpest fronts, some
    #                        numerical oscillations near gradients.
    #                0.1-1 — typical urban-flood range (moderate smoothing).
    #                1-10  — larger basin-scale or regional runs where
    #                        sub-grid mixing is significant.
    ade_transport      = dict(enabled=True, d_coef=0.5),

    # Soil infiltration — simplified constant-rate loss per cell per
    # timestep ( infil = min(Ks·dt, h) ). Useful for urban / pervious
    # mixed surfaces where infiltration is non-negligible on the
    # simulation timescale.
    #
    #   enabled           — True to activate, False to skip (default).
    #   default_ks_mm_hr  — saturated hydraulic conductivity in mm/h.
    #                       Typical values (Rawls et al. 1983, USDA):
    #                         sand       ~ 210  mm/h
    #                         loamy sand ~  60  mm/h
    #                         sandy loam ~  25  mm/h
    #                         loam       ~  13  mm/h
    #                         silt loam  ~   6  mm/h
    #                         clay loam  ~   2  mm/h
    #                         clay       ~   1  mm/h
    #                         asphalt / concrete ≈ 0
    # Horton-decay infiltration. ``soil_types`` is a dict of {id: USDA
    # class name} — every cell that doesn't have a ``SOIL_MAP`` override
    # falls back to type 1. Valid USDA classes are listed in
    # ``src/fluxos/soil_infiltration.cpp::get_usda_soil_table()``.
    # Set ``enabled=False`` to skip infiltration entirely.
    soil_infiltration  = dict(
        enabled    = True,
        soil_map   = "",
        soil_types = {1: "sandy_loam"},
    ),

    # Steady-state convergence check (the sim exits once ``max|Δh|``
    # between two adjacent timesteps drops below ``tolerance`` and
    # ``tim >= min_time_s``). Keep disabled for a fixed-duration run.
    steady_state       = dict(
        enabled      = False,
        tolerance    = 0.001,
        min_time_s   = 1800,
    ),

    # ==================================================================
    # 9. Docker / run settings — only affect the snippets shown in the
    #    HTML report. Do NOT change the files written to disk.
    # ==================================================================

    # Suggest a Hybrid MPI + OpenMP build / run in the report?
    #   False — OpenMP-only (recommended for single machine, < 1000² cells).
    #   True  — the report's build + run snippets include ``-DUSE_MPI=ON``
    #           and ``mpirun -n N`` for HPC clusters / multi-node runs.
    use_mpi             = False,

    # Number of MPI ranks. Only used when ``use_mpi = True``. Practical
    # choice depends on the HPC environment: 4-16 for workstations,
    # 16-64 for cluster nodes, 64+ for large distributed runs. See
    # ``wikipage/source/HPC.rst`` for domain-decomposition guidance.
    mpi_np              = 4,

    # Legacy / compat field. Triangular-mesh support is now ALWAYS
    # compiled in — a single binary handles both regular and triangular
    # meshes, dispatching at runtime based on ``mesh_type`` above. This
    # field is ignored by the report generator; keep it or delete it,
    # your choice.
    use_trimesh_build   = True,

    # Path to the compiled FLUXOS binary — used by the "Run the
    # simulation" snippet in the generated HTML report.
    #
    # Accepts:
    #   • path RELATIVE to ``repo_root`` (default ``"bin/fluxos"``) —
    #     matches the container-compile convention where
    #     ``-DCMAKE_RUNTIME_OUTPUT_DIRECTORY=/work/bin`` drops the
    #     binary at ``<repo>/bin/fluxos`` on the host via the compose
    #     bind mount.
    #   • ABSOLUTE path — for a system-installed binary or an HPC
    #     module path, e.g. ``"/opt/fluxos/bin/fluxos_mpi"``.
    #
    # Typical alternatives:
    #   "bin/fluxos_mpi"       — hybrid MPI + OpenMP build
    #   "bin/fluxos_cuda"      — CUDA GPU build
    #   "build/bin/fluxos"     — native (non-container) debug build tree
    fluxos_executable   = "bin/fluxos",

    # Visualisation time-stride for the interactive WebGL viewer (the
    # "6a. Interactive WebGL viewer" snippet in the report). The viewer
    # keeps every Nth FLUXOS output snapshot as an animation frame:
    #
    #   1  — keep every snapshot (smoothest animation; largest bundle).
    #   2  — every other snapshot (½ bundle, ½ playback duration).
    #   5  — every 5th snapshot (default in ``fluxos_viewer.py``; good
    #        when ``print_step_s`` is short and there are hundreds of
    #        outputs).
    #   N  — every Nth snapshot.
    #
    # Rule of thumb: pick so that
    #   (number_of_outputs / webgl_time_stride) ≈ 30-100 frames.
    # Fewer than ~5 frames yields a jerky animation; more than ~300 makes
    # the bundle unnecessarily large with no perceptual gain.
    #
    # For a short debug run with only 8 snapshots the default of 5 would
    # keep just 2 frames — drop this to 1.
    webgl_time_stride   = 1,

    # ==================================================================
    # 10. Report options — what happens after the template finishes.
    # ==================================================================

    # Produce the HTML configuration report:
    #   True  — write ``reports/<project>_report.html`` with KPI tiles,
    #           DEM + mesh preview map, full modset JSON, module status,
    #           and the copy-paste Next Steps for Docker / Apptainer / run
    #           / visualise. Recommended.
    #   False — skip (useful for CI / automation).
    generate_report = True,

    # After generating the report, open it in the system's default browser
    # (via ``webbrowser.open_new_tab``).
    #   True  — auto-open. Great for interactive use.
    #   False — only print the path (useful for headless runs).
    open_report     = True,
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
