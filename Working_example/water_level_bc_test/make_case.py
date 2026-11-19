#!/usr/bin/env python3
"""Synthetic tidal-channel test for the regular-mesh WATER_LEVEL boundary.

Generates, next to this script:
  channel_10m.asc        30 x 120 grid, 10 m cells. A 60 m wide channel
                         (rows 13-18 counted from the south) runs west-east;
                         bed rises from 3.0 m (west edge) to 4.0 m (east);
                         floodplain at 6.5 m (never overtopped).
  river_inflow.csv       constant 1 m3/s river inflow near the east end.
  tide_level.csv         west-edge stage: 5.0 - 1.0*cos(2*pi*t/10800) m,
                         i.e. low water 4.0 m at t=0, high water 6.0 m at
                         t=5400 s (3 h "tide" so the test runs in seconds).
  modset_wl_west.json    tide imposed on the whole west edge
  modset_no_bc.json      same case without BOUNDARY_CONDITIONS (regression)
  channel_10m_closed.asc + modset_wl_west_closed.json
                         east end walled off -> exact volume-budget check

All elevations are kept > 0 on purpose: read_geo() treats zb == 0 as
NODATA and takes |zb|, so a real estuary DEM must be datum-shifted.
"""
import json, math, os

here = os.path.dirname(os.path.abspath(__file__))
NR, NC, DX = 30, 120, 10
XLL, YLL = 500000.0, 4000000.0
CH = range(13, 19)                     # channel rows (irow, 1 = south)

def zb(irow, icol, closed=False):
    if closed and icol == NC:
        return 9.0                     # high bank: no east-edge weir outflow
    if irow in CH:
        return 3.0 + 1.0 * (icol - 1) / (NC - 1)
    return 6.5

for name, closed in (("channel_10m.asc", False), ("channel_10m_closed.asc", True)):
    with open(os.path.join(here, name), "w") as f:
        f.write(f"ncols {NC}\nnrows {NR}\nxllcorner {XLL}\nyllcorner {YLL}\n"
                f"cellsize {DX}\nNODATA_value -99999\n")
        for irow in range(NR, 0, -1):     # first data line = northernmost row
            f.write(" ".join(f"{zb(irow, c, closed):.4f}" for c in range(1, NC + 1)) + "\n")

T_END = 14400
with open(os.path.join(here, "river_inflow.csv"), "w") as f:
    f.write("0,0,0\n")                # header row (skipped by FLUXOS)
    for t in range(0, T_END + 1, 600):
        f.write(f"{t},1.0,0\n")

with open(os.path.join(here, "tide_level.csv"), "w") as f:
    f.write("t_seconds,level_m\n")
    for t in range(0, T_END + 1, 300):
        f.write(f"{t},{5.0 - 1.0 * math.cos(2 * math.pi * t / 10800):.4f}\n")

base = {
    "COMMENT": "Synthetic tidal channel - WATER_LEVEL boundary test",
    "DEM_FILE": "channel_10m.asc",
    "SIM_DATETIME_START": "2026-01-01 00:00:00",
    "RESTART": False,
    "INFLOW_FILE": {
        "FILENAME": "river_inflow.csv",
        # icol 117, irow 15 (see read_geo's row/col from coordinates)
        "DISCHARGE_LOCATION": {"X_COORDINATE": XLL + 1175, "Y_COORDINATE": YLL + 155},
    },
    "OUTPUT": {"OUTPUT_FOLDER": "Results_no_bc", "PRINT_STEP": 300, "H_MIN_TO_PRINT": 0.001},
    "ROUGNESS_HEIGHT": 0.005,
    "EXTERNAL_MODULES": {"ADE_TRANSPORT": {"STATUS": False, "D_COEF": 0.5}},
}
with open(os.path.join(here, "modset_no_bc.json"), "w") as f:
    json.dump(base, f, indent=4)

wl = json.loads(json.dumps(base))
wl["OUTPUT"]["OUTPUT_FOLDER"] = "Results_wl_west"
wl["BOUNDARY_CONDITIONS"] = {
    "WATER_LEVEL": {"STATUS": True, "FILENAME": "tide_level.csv", "EDGE": "west"}
}
with open(os.path.join(here, "modset_wl_west.json"), "w") as f:
    json.dump(wl, f, indent=4)

# Mass-balance check: east end closed by a 9 m bank, so the only sources and
# sinks are the river inflow and the WATER_LEVEL boundary. Domain storage
# must equal  Q*t + vol_net_m3  (check_propagation.py prints the residual).
cl = json.loads(json.dumps(wl))
cl["DEM_FILE"] = "channel_10m_closed.asc"
cl["OUTPUT"]["OUTPUT_FOLDER"] = "Results_wl_west_closed"
with open(os.path.join(here, "modset_wl_west_closed.json"), "w") as f:
    json.dump(cl, f, indent=4)
print("case written to", here)
