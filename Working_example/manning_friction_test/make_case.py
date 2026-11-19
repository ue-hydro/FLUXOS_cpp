#!/usr/bin/env python3
"""Manning bed-friction test for the regular mesh (MANNING_N).

A straight channel, 1.2 km long, sloping 0.001 from south to north, carries a
constant 3 m3/s and drains over the legacy north-edge weir. With MANNING_N the
flow must settle at Manning's normal depth
    h_n = (q n / sqrt(S))^(3/5),   q = Q / wetted width.
Without it the flow keeps accelerating and thins out (no friction).

Writes slope_channel.asc, inflow.csv, modset_manning.json, modset_no_friction.json.
Run both, then: python3 check_normal_depth.py Results_manning
"""
import json, os
here = os.path.dirname(os.path.abspath(__file__))
NR, NC, DX, S, Q, N = 120, 7, 10.0, 0.001, 3.0, 0.03
XLL, YLL = 500000.0, 4000000.0
with open(os.path.join(here, "slope_channel.asc"), "w") as f:
    f.write(f"ncols {NC}\nnrows {NR}\nxllcorner {XLL}\nyllcorner {YLL}\ncellsize {DX}\nNODATA_value -9999\n")
    for r in range(NR):                       # first line = northernmost row
        f.write(" ".join("50.0000" if (c in (0, NC - 1) or r >= NR - 2)
                         else f"{10.0 + S * r * DX:.4f}" for c in range(NC)) + "\n")
with open(os.path.join(here, "inflow.csv"), "w") as f:
    f.write("0,0,0\n")                        # header row (skipped by FLUXOS)
    f.writelines(f"{t},{Q},0\n" for t in range(0, 7201, 60))
for name, n in (("modset_manning.json", N), ("modset_no_friction.json", None)):
    d = {"COMMENT": "Manning friction test: sloping channel, normal depth",
         "DEM_FILE": "slope_channel.asc", "SIM_DATETIME_START": "2026-01-01 00:00:00",
         "RESTART": False,
         "INFLOW_FILE": {"FILENAME": "inflow.csv",
                         "DISCHARGE_LOCATION": {"X_COORDINATE": XLL + 35, "Y_COORDINATE": YLL + 55}},
         "OUTPUT": {"OUTPUT_FOLDER": "Results_" + name[7:-5], "PRINT_STEP": 600, "H_MIN_TO_PRINT": 0.0001},
         "ROUGNESS_HEIGHT": 0.005,
         "EXTERNAL_MODULES": {"ADE_TRANSPORT": {"STATUS": False, "D_COEF": 0.5}}}
    if n:
        d["MANNING_N"] = n
    json.dump(d, open(os.path.join(here, name), "w"), indent=4)
print("written")
