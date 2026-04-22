# FLUXOS — Read Outputs

Template-driven analysis of a FLUXOS simulation: reads the per-timestep
output files, computes flood statistics, and produces a self-contained
HTML report that **complements** the WebGL animation with aggregated
numbers and maps that are hard to read off a moving picture.

Mirrors the layout of `1_Model_Config/`: **you edit one file** (`read_output_template.py`),
everything else is library code under `output_supporting_lib/`.

## Layout

```
2_Read_Outputs/
├── read_output_template.py     ← edit this
├── reports/                    ← generated HTML reports + CSV land here
└── output_supporting_lib/      ← internal helpers (don't edit)
    ├── read_results.py         ← regular .txt + triangular .vtu readers
    ├── flood_statistics.py     ← streaming per-timestep + per-cell stats
    ├── gen_results_report.py   ← HTML report (dark theme, Plotly, per-line copy)
    ├── analysis_driver.py      ← orchestrator invoked by the template
    ├── fluxos_viewer.py        ← existing KML / MP4 / WebGL exporter (unchanged)
    └── plot_horton_infiltration.py
```

## Report contents

| # | Section | What it shows | Why |
|---|---|---|---|
| 1 | **Summary KPIs** | Peak volume, peak flooded area, max depth, max velocity, time-to-peak, total wet cells | At-a-glance severity |
| 2 | **Time-series** (5 panels) | Flood volume, flooded area, max depth, mean wet-cell depth, max velocity | Quantifies onset / peak / recession — the animation shows this qualitatively |
| 3 | **Maximum inundation map** | Per-cell max depth over the whole run | The classic floodplain-mapping output |
| 4 | **Hazard classification map** | Cells classified as H1 / H2 / H3 / H4 using the ARR-2019 `D = h·(v + 0.5)` criterion | Standard depth×velocity danger index |
| 5 | **Max-depth histogram** | Depth distribution across wet cells | Shallow spreading vs. deep channelised |
| 6 | **Time-of-first-inundation map** | First t at which each cell became wet | Flood-front propagation |
| 7 | **Chemical transport panel** *(if ADE_TRANSPORT was on)* | `max(conc_SW)` and `mean(conc_SW)` over time | Numeric companion to the WebGL `--variable conc_SW` view |

## How to use

```bash
cd fluxos_preprocessing/2_Read_Outputs
# open read_output_template.py, point it at your results_dir + modset
python read_output_template.py
```

The script:
- streams snapshots (no need to hold all timesteps in memory — a 30-hour sim
  with ~10 GB of `.txt` files analyses fine)
- computes time-series + per-cell aggregates in a single pass
- writes `reports/<project>_results_report.html` and optionally
  `reports/<project>_timeseries.csv`
- opens the HTML in your default browser

## Dependencies

`numpy` for everything; `pyvista` only if `mesh_type = "triangular"` (reads `.vtu`).
Plotly is loaded via CDN from inside the HTML so the Python side doesn't need it.

## Other tools here

- [`output_supporting_lib/fluxos_viewer.py`](output_supporting_lib/fluxos_viewer.py) —
  existing KML / MP4 / WebGL exporter. Run directly:
  `python output_supporting_lib/fluxos_viewer.py --results-dir … --dem …`
- [`output_supporting_lib/plot_horton_infiltration.py`](output_supporting_lib/plot_horton_infiltration.py) —
  infiltration-rate plots for soil module runs.
