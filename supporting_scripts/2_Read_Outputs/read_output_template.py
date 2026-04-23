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
FLUXOS results analysis template
==================================

HOW TO USE
----------
Two equivalent paths — both run the same analysis pipeline:

A) EDIT THIS FILE — point `results_dir` / `modset_file` at your run,
   then just run it:

       python read_output_template.py

B) RUN-TIME OVERRIDES — leave the `_config` dict alone and pass the
   per-project paths as command-line flags (useful for batch / scripted
   runs, and what the HTML configuration report emits):

       python read_output_template.py \
           --results-dir Results_TorresVedras \
           --modset-file Working_example/modset_torres_vedras.json \
           --project-name "Torres Vedras"

Whichever you choose, the script will:
  - stream the per-timestep snapshots from `results_dir` (no need to load
    everything into memory at once)
  - compute flood statistics (volume, flooded area, peak depth, hazard, ...)
  - generate an interactive HTML report under `2_Read_Outputs/reports/`
  - open the report in your default browser (disable with `open_report=False`)
  - optionally export the time-series as CSV for external post-processing

All the computation / plotting code lives in
`2_Read_Outputs/output_supporting_lib/`. You should not need to touch it.
"""

_config = dict(

    # ------------------------------------------------------------------
    # 1. Project metadata (used in the report header)
    # ------------------------------------------------------------------
    project_name = "Rosa Creek",
    authors      = ["Diogo Costa"],

    # ------------------------------------------------------------------
    # 2. Paths
    # ------------------------------------------------------------------
    # `repo_root = None` auto-detects (two directories above this file).
    repo_root    = None,

    # Directory containing the per-timestep output files FLUXOS wrote.
    #   Regular mesh → `<seconds>.txt`
    #   Triangular   → `<seconds>.vtu`
    results_dir  = "Results_river_30h",

    # The modset .json that drove the run — used to detect mesh type, CRS,
    # whether transport was enabled, etc.
    modset_file  = "Working_example/modset_river_30h.json",

    # ------------------------------------------------------------------
    # 3. Analysis settings
    # ------------------------------------------------------------------
    # Cells with h below this threshold (metres) are treated as dry when
    # computing flooded-area, mean-depth, and first-wet-time statistics.
    h_threshold_flooded_m = 0.01,

    # Analyse every Nth timestep for a faster preview. Use 1 for full detail.
    # (A 30-hour river sim with 360 snapshots at stride=5 processes ~72 files.)
    stride = 1,

    # ------------------------------------------------------------------
    # 4. Output options
    # ------------------------------------------------------------------
    generate_report       = True,
    open_report           = True,
    export_csv_timeseries = True,
)


# ============================================================================
# You should not need to edit anything below this line.
# ============================================================================

if __name__ == "__main__":
    import argparse
    import os
    import sys

    # Overlay CLI flags on top of the `_config` dict above. Flags are all
    # optional — if you have edited the file, just run with no args.
    ap = argparse.ArgumentParser(
        description="FLUXOS results analysis (edit _config or pass flags)",
    )
    ap.add_argument(
        "--results-dir", dest="results_dir", default=None,
        help="Folder containing FLUXOS output (.txt for regular mesh, "
             ".vtu for triangular). Overrides _config['results_dir'].",
    )
    ap.add_argument(
        "--modset-file", dest="modset_file", default=None,
        help="Path to the modset .json that produced the results. "
             "Overrides _config['modset_file'].",
    )
    ap.add_argument(
        "--project-name", dest="project_name", default=None,
        help="Human-readable project name used in the HTML report header.",
    )
    ap.add_argument(
        "--stride", type=int, default=None,
        help="Analyse every Nth timestep. Overrides _config['stride'].",
    )
    ap.add_argument(
        "--no-open", dest="no_open", action="store_true",
        help="Do NOT open the HTML report in a browser after generation.",
    )
    args = ap.parse_args()

    # Apply overrides (only when provided)
    if args.results_dir  is not None: _config["results_dir"]  = args.results_dir
    if args.modset_file  is not None: _config["modset_file"]  = args.modset_file
    if args.project_name is not None: _config["project_name"] = args.project_name
    if args.stride       is not None: _config["stride"]       = args.stride
    if args.no_open:                  _config["open_report"]  = False

    _here = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, os.path.join(_here, "output_supporting_lib"))
    from analysis_driver import run      # noqa: E402
    run(_config, template_file=os.path.abspath(__file__))
