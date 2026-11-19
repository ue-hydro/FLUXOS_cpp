#!/usr/bin/env python3
"""Post-process the WATER_LEVEL test: stage along the channel centreline vs
time, plus a domain volume budget. Usage: python3 check_propagation.py [Results_dir]
Writes stage_along_channel.csv (and a PNG if matplotlib is available)."""
import csv, glob, math, os, sys

d = sys.argv[1] if len(sys.argv) > 1 else "Results_wl_west"
CELL_A, ROW, STATIONS, Q_IN = 100.0, 15, [1, 2, 10, 40, 80, 115], 1.0

def tide(t):
    return 5.0 - 1.0 * math.cos(2 * math.pi * t / 10800)

rows = []
for fn in sorted(glob.glob(os.path.join(d, "*.txt")), key=lambda p: int(os.path.basename(p)[:-4])):
    t = int(os.path.basename(fn)[:-4])
    stage, vol = {}, 0.0
    with open(fn) as f:
        for line in f:
            p = line.replace(",", " ").split()
            try:
                irow, icol, z, h = int(float(p[0])), int(float(p[1])), float(p[4]), float(p[5])
            except (ValueError, IndexError):
                continue
            vol += h * CELL_A
            if irow == ROW and icol in STATIONS:
                stage[icol] = z
    rows.append((t, stage, vol))

bal = {}
bf = os.path.join(d, "water_level_bc_balance.csv")
if os.path.exists(bf):
    for r in csv.DictReader(open(bf)):
        bal[int(float(r["t_s"]))] = float(r["vol_net_m3"])

out = os.path.join(d, "stage_along_channel.csv")
with open(out, "w") as f:
    f.write("t_s,tide_m," + ",".join(f"z_icol{c}" for c in STATIONS) + ",storage_m3,bc_net_m3,inflow_m3,residual_m3\n")
    print(f"{'t':>6} {'tide':>6} " + " ".join(f"{'z@'+str(c):>7}" for c in STATIONS) + f" {'storage':>9} {'BCnet+Qt':>9} {'resid':>8}")
    for t, st, vol in rows:
        zs = [st.get(c, float("nan")) for c in STATIONS]
        bn = bal.get(t, 0.0)
        res = vol - (bn + Q_IN * t)   # = -(edge-weir outflow) if the budget closes
        f.write(f"{t},{tide(t):.4f}," + ",".join(f"{z:.4f}" for z in zs) + f",{vol:.2f},{bn:.2f},{Q_IN*t:.1f},{res:.2f}\n")
        if t % 900 == 0:
            print(f"{t:6d} {tide(t):6.3f} " + " ".join(f"{z:7.3f}" for z in zs) + f" {vol:9.0f} {bn+Q_IN*t:9.0f} {res:8.0f}")
print("written", out)

try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    ts = [r[0] / 3600 for r in rows]
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(ts, [tide(r[0]) for r in rows], "k--", lw=2, label="imposed level (west edge)")
    for c in STATIONS[1:]:
        ax.plot(ts, [r[1].get(c, float("nan")) for r in rows], label=f"stage, x = {c*10-5} m")
    ax.set_xlabel("time [h]"); ax.set_ylabel("water surface [m]"); ax.legend(fontsize=8)
    ax.set_title("Tide propagating up the channel (FLUXOS regular mesh)")
    fig.tight_layout(); fig.savefig(os.path.join(d, "stage_along_channel.png"), dpi=120)
except ImportError:
    pass
