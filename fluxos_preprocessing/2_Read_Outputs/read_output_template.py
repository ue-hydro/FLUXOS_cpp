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
1. Edit the `_config` dict below:
    - point `results_dir` at a folder that FLUXOS wrote (e.g. `Results_river_30h`)
    - point `modset_file` at the `.json` that generated those results
    - tweak thresholds / flags if needed

2. Run it from any working directory:

       python read_output_template.py

3. The script will:
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
    modset_file  = "bin/modset_river_30h.json",

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
    import os, sys
    _here = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, os.path.join(_here, "output_supporting_lib"))
    from analysis_driver import run      # noqa: E402
    run(_config, template_file=os.path.abspath(__file__))
