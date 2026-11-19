# Solver mirror-symmetry test

```bash
python3 make_bank.py && mkdir -p outb
../../build/bin/fluxos mb.json
python3 banks.py outb      # S(r<=12) and N(r>=18) floodplain volumes [m3]
```

Expected after fix/solver-symmetry: both floodplains wet, similar volumes
(the outer walls are not exactly symmetric about the channel, so not
identical). Before the fix (master 65ce73d + clock init): all overbank water
went north and the south floodplain stayed dry.
