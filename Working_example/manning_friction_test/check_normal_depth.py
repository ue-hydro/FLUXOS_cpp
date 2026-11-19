#!/usr/bin/env python3
"""Compare the last snapshot's along-channel depth with Manning's normal depth.
Usage: python3 check_normal_depth.py Results_manning [n=0.03]"""
import glob, os, sys
d = sys.argv[1] if len(sys.argv) > 1 else "Results_manning"
n = float(sys.argv[2]) if len(sys.argv) > 2 else 0.03
S, Q = 0.001, 3.0
fn = max(glob.glob(os.path.join(d, "*.txt")), key=lambda p: int(os.path.basename(p)[:-4]))
H = {}
for line in open(fn):
    p = line.split(",")
    try:
        r, h = int(p[0]), float(p[5])
    except (ValueError, IndexError):
        continue
    if h > 0.005:
        H.setdefault(r, []).append(h)
width = 10.0 * max(len(v) for v in H.values())          # wetted width (zb averaging may hide a column)
hn = (Q / width * n / S ** 0.5) ** 0.6
print(f"{fn}: wetted width {width:.0f} m, Manning normal depth {hn:.3f} m")
for r in (20, 40, 60, 80, 100):
    if r in H:
        h = sum(H[r]) / len(H[r])
        print(f"  irow {r:3d}: h = {h:.3f} m  ({100 * (h - hn) / hn:+.1f} %)")
