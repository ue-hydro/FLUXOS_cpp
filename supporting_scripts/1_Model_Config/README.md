# FLUXOS — Model Configuration

One-stop template that turns a raw GeoTIFF DEM into the inputs FLUXOS needs
(`Working_example/*.asc`, `Working_example/*.msh`, `Working_example/modset_*.json`)
and produces an HTML summary with copy-paste container-run commands for
building and running the simulation.

The pattern mirrors OpenWQ's `1_Model_Config/` template: **you edit one file**
(`model_config_template.py`), everything else is library code.

## Layout

```
1_Model_Config/
├── model_config_template.py      ← edit this
├── reports/                      ← generated HTML reports land here
└── config_support_lib/           ← internal helpers (don't edit)
    ├── driver.py                 ← orchestrator (writes modset.json inline)
    ├── dem_operations.py         ← GeoTIFF → .asc
    ├── dem_download.py           ← optional: auto-fetch a DEM for a bbox
    ├── mesh_generation.py        ← Gmsh adaptive triangular mesh
    └── gen_report.py             ← HTML report builder
```

## DEM source: bring your own, or let the template fetch one

Set `dem_source_mode` in the template:

- **`"file"`** *(default)* — point `dem_source_geotiff` at an existing
  GeoTIFF on disk. The classic workflow.
- **`"download"`** — fetch a DEM for `download_bbox_wgs84` (a lon/lat
  bounding box) from an external provider, reproject to UTM automatically,
  then continue the normal pipeline. Results are cached under
  `bin/dem_cache/` so repeated runs skip the network.

**Supported providers** (pick one as `download_provider`):

| Provider | Resolution | Coverage | Vertical RMSE | Key? |
|---|---|---|---|---|
| `COP30` ★ | 30 m | Global | 4 m | OpenTopography free key |
| `COP90` | 90 m | Global | 4 m | OpenTopography |
| `SRTMGL1` | 30 m | ±60° latitude | 16 m | OpenTopography |
| `SRTMGL3` | 90 m | ±60° latitude | 16 m | OpenTopography |
| `AW3D30` | 30 m | Global | 5 m | OpenTopography |
| `NASADEM` | 30 m | ±60° latitude | 10 m | OpenTopography |
| `EU_DTM` | 10 m | Europe (EEA39) | 3 m | OpenTopography |
| `GEBCOIceTopo` | 463 m | Global (land + ocean) | — | OpenTopography |
| `USGS_10M` | 10 m | CONUS | 1 m | None (`pip install py3dep`) |
| `USGS_3M` | 3 m | CONUS partial | 0.5 m | None |
| `USGS_1M` | 1 m | US LiDAR areas | 0.15 m | None |

**OpenTopography API key** (free, one-time):
1. Register at <https://portal.opentopography.org/>.
2. `export OPENTOPOGRAPHY_API_KEY=your_key_here` (or set a different env
   var name via `download_api_key_env`).

**Caveat**: don't set `dem_target_resolution_m` below the provider's
native resolution (e.g. 2 m target from a 30 m SRTM source) — that's
interpolation, not added information. The report flags this with a warning.
For urban / street-level flood modelling where you genuinely need sub-10 m
DEMs, you must provide LiDAR-derived data yourself via `"file"` mode.

## How to use

1. **Install Python deps** (once per machine):

   ```
   pip install -r ../requirements.txt
   ```

2. **Copy the template** (optional — you can edit in place for a one-off):

   ```
   cp model_config_template.py my_project.py
   ```

3. **Edit the `_config = dict(...)` block**. The ten sections are commented
   inline; typical edits for a new project are the project name, the DEM
   path, the output modset name, and mesh sizing if you are going triangular.

4. **Run it**:

   ```
   python model_config_template.py     # or: python my_project.py
   ```

   This will:
   - write the `.asc` (and `.msh` if triangular) into `<repo>/bin/`,
   - write the modset into `<repo>/bin/<modset_name>`,
   - write an HTML report into `reports/<project_name>_report.html`,
   - open the report in your default browser (disable via `open_report=False`).

5. **Copy the Docker commands from the report** into your terminal to build
   and run the model. They are generated for your host OS and point directly
   at the modset you just produced.

## What the report contains

| Section | Contents |
|---|---|
| Summary | KPI tiles: mesh type, resolution, cell count, print step, active modules |
| Domain | DEM source / CRS / bbox / elevation range |
| Mesh   | Vertex / cell / edge counts (triangular only) |
| Configuration | All non-default modset values, plus the full JSON in a collapsible |
| Modules | Status of ADE Transport, Soil Infiltration, OpenWQ |
| Errors | (if any step failed) — the rest of the report still generates |
| Next Steps | Copy-paste shell snippets: `cd`, build, run, check outputs, visualise |

## Post-processing

Not part of this template — see the sibling `2_Read_Outputs/` folder, which
has the same template-driven pattern:

- `2_Read_Outputs/read_output_template.py` — generates a flood-statistics
  HTML report (volume / flooded-area time series, max-inundation map,
  hazard classification, depth histogram, first-inundation map, chemical
  transport panel).
- `2_Read_Outputs/output_supporting_lib/fluxos_viewer.py` — KML / MP4 /
  WebGL animated visualisations.

The report's *Next Steps* section already surfaces both with commands
pre-filled with the correct DEM, UTM zone, and modset path for the project.
