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

The modsets use `"MANNING_N": 0.03` (bed friction) and the default Courant number 0.5.

`expected_stage_along_channel.png`: water surface at x = 15, 95, 395, 795 and
1145 m along the channel centreline (closed case). The stage next to the boundary
follows the imposed level; upstream the tide arrives later and slightly damped
(e.g. at t = 9000 s the tide is 4.50 m and the closed upstream end 4.58 m).

`check_propagation.py` also prints a domain volume budget:

    residual = storage(t) − (Q_in·t + vol_net_m3)

where `vol_net_m3` comes from `water_level_bc_balance.csv`. In the closed case
the residual is within ±1 m³ at every print step: the boundary ledger closes the
budget. No water reaches the 6.5 m floodplain.

## History

Earlier versions of this test showed water on the floodplain, a north-edge weir
loss of ≈ 17 000 m³ and a north/south asymmetry. These came from the solver, not
from the boundary:
- a north-biased solver (fixed on `fix/solver-symmetry`);
- the legacy Courant number 1, which drove a grid-scale oscillation between the
  channel banks and the channel. The banks spilled onto the floodplain, which
  then drained over the north-edge weir. The default is now 0.5 (`COURANT_NUMBER`).
