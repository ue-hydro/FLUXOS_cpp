# WATER_LEVEL boundary — regular-mesh test case

Synthetic 30 × 120 grid (10 m cells): a 60 m wide channel runs west→east,
bed 3.0 m (west) → 4.0 m (east), floodplain 6.5 m. A 3 h "tide"
(4.0 → 6.0 → 4.0 m) is imposed on the whole **west** edge; 1 m³/s river
inflow near the east end.

```bash
python3 make_case.py                      # (re)generate DEMs, CSVs, modsets
mkdir -p Results_wl_west Results_wl_west_closed Results_no_bc
../../build/bin/fluxos modset_wl_west.json         # tide, open east end
../../build/bin/fluxos modset_wl_west_closed.json  # tide, east end walled off
../../build/bin/fluxos modset_no_bc.json           # regression: no BC key
python3 check_propagation.py Results_wl_west_closed
```

Run from this folder (modset paths are relative to the working directory).
Each run takes ~1 min with 4 OpenMP threads.

## What to expect

`expected_stage_along_channel.png`: water surface at x = 15, 95, 395, 795 and
1145 m along the channel centreline. Stage next to the boundary follows the
imposed level; upstream the tide arrives progressively later and damped.

`check_propagation.py` also prints a domain volume budget:

    residual = storage(t) − (Q_in·t + vol_net_m3)

where `vol_net_m3` comes from `water_level_bc_balance.csv`. The residual is the
volume leaving through the *other* (legacy weir) domain edges. In the closed
case it is ≈ −17 200 m³, and an instrumented build showed it is exactly the
north-edge weir outflow (fe_1 at irow = NROWS). So the BC ledger closes the
budget exactly; the loss is a pre-existing solver behaviour (see below).

## Known pre-existing solver behaviour (not caused by the BC)

Without any boundary condition, the unmodified solver (master 65ce73d) moves
overbank water only to the north (+irow) side of this symmetric channel:
1402 wet cells north vs 0 south after 4 h. With the tide, this shows up as a
water surface ≈ 2 m higher on the north bank than on the south bank near the
boundary, and as outflow over the north edge. Worth investigating separately.
