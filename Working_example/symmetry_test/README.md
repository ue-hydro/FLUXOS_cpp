# Solver mirror-symmetry test

```bash
python3 make_bank.py && mkdir -p outb
../../build/bin/fluxos mb.json
python3 banks.py outb      # S(r<=12) and N(r>=18) floodplain volumes [m3]
```

Expected (6 m3/s, CFL 0.5): no overbank flow until ~4800 s, then both
floodplains wet. At 7200 s S = 4572 m3 and N = 6114 m3, banks 4215/4222 m3, and the
mirrored floodplain rows are within a few percent near the channel. The
north side is larger (13 vs 11 rows) and the inflow cell sits half a row north
of the centre line, so the totals are not identical. Before the fix (master 65ce73d + clock init): all overbank water
went north and the south floodplain stayed dry.
