# Manning bed-friction test (regular mesh, `MANNING_N`)

A 1.2 km straight channel with slope S = 0.001 (south → north) carries a constant
3 m³/s and drains over the north-edge weir. With bed friction the flow must settle
at Manning's normal depth `h_n = (q n / √S)^(3/5)`.

```bash
python3 make_case.py
../../build/bin/fluxos modset_manning.json        # MANNING_N = 0.03
../../build/bin/fluxos modset_no_friction.json    # key absent: legacy, frictionless
python3 check_normal_depth.py Results_manning
python3 check_normal_depth.py Results_no_friction
```

Result at t = 2 h (wetted width 40 m: the bed averaging in `initiate()` hides one
of the five channel columns; h_n = 0.205 m):

| irow | with `MANNING_N` 0.03 | without |
|---|---|---|
| 20 | 0.198 m (−3 %) | 0.028 m |
| 40 | 0.206 m (+1 %) | 0.012 m |
| 60 | 0.208 m (+1 %) | 0.021 m |
| 80 | 0.204 m (−0 %) | 0.152 m |
| 100 | 0.200 m (−2 %) | 0.390 m |

Without friction the water accelerates down the slope and thins to a few cm, then
ponds against the weir.
