#!/usr/bin/env python3
"""Mirror-symmetry regression test for the regular-mesh solver.

Closed 31 x 41 basin (10 m cells) with a west-east channel that is symmetric
north/south (after initiate's half-cell zb averaging: bed rows 14-16, 4.0 m
banks at rows 13 and 17, 5.0 m floodplain). 2 m3/s enters the channel centre
line. Overbank water must reach both floodplains alike: run, then
`python3 banks.py outb`. Inflow is given every 10 s so the dry-start step
does not dump a whole forcing interval into one cell (keep < 2000 rows: the
forcing arrays are fixed at 2000 rows).
"""
import json
NR, NC, DX, XLL, YLL = 31, 41, 10, 500000.0, 4000000.0
# channel orig rows 14..17 -> after initiate's half-cell smoothing: rows 14..16 bed, banks at 13 and 17 (symmetric about 15)
def zb(r,c):
    if r in (1,NR) or c in (1,NC): return 9.0
    return 3.0 if 14<=r<=17 else 5.0
with open("bank.asc","w") as f:
    f.write(f"ncols {NC}\nnrows {NR}\nxllcorner {XLL}\nyllcorner {YLL}\ncellsize {DX}\nNODATA_value -99999\n")
    for r in range(NR,0,-1): f.write(" ".join(f"{zb(r,c):.2f}" for c in range(1,NC+1))+"\n")
with open("qb.csv","w") as f:
    f.write("0,0,0\n")
    for t in range(0,7201,10): f.write(f"{t},2.0,0\n")
m={"DEM_FILE":"bank.asc","SIM_DATETIME_START":"2026-01-01 00:00:00","RESTART":False,
   "INFLOW_FILE":{"FILENAME":"qb.csv","DISCHARGE_LOCATION":{"X_COORDINATE":XLL+215,"Y_COORDINATE":YLL+155}},
   "OUTPUT":{"OUTPUT_FOLDER":"outb","PRINT_STEP":1200,"H_MIN_TO_PRINT":0.0001},
   "ROUGNESS_HEIGHT":0.005,"EXTERNAL_MODULES":{"ADE_TRANSPORT":{"STATUS":False,"D_COEF":0.5}}}
json.dump(m,open("mb.json","w"))
