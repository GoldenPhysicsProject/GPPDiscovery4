#!/usr/bin/env python3
"""
Weil window scan, pass 2: refine the F_theta sign crossing and audit the
small-L end. Turbo, 2026-10-13. NUMERICAL.
"""
import sys, time
import numpy as np

sys.path.insert(0, "/home/hatch/workspace")
sys.path.insert(0, "/home/hatch/workspace/user/files")
from ritz_analysis import setup, neutral_basis
from theta_control_ritz_4_spk4 import trigbasis

FUNCS = trigbasis(30, 0.75, (0, 1, 2))

def m(L, th, h=0.01):
    d = setup(th, FUNCS, L=L, h=h, ycut=6.0)
    ev = neutral_basis(d["E"], d["Wm"], d["Cm"])[1]
    return float(ev[0]) - 1.0

print("crossing refinement, h=0.01", flush=True)
for L in [0.5, 0.55, 0.6, 0.65, 0.6931, 0.72, 0.75]:
    vals = [m(L, th) for th in (0.0, 0.1, 0.25)]
    print(f"  L={L:6.4f} lam={np.exp(L):6.3f}  th0 {vals[0]: .6e}  th.1 {vals[1]: .6e}  th.25 {vals[2]: .6e}", flush=True)

print("grid audit at L=0.5 (is the th=0 dip a grid artifact?)", flush=True)
for h in (0.01, 0.005, 0.0025):
    print(f"  h={h}: th=0 {m(0.5, 0.0, h): .6e}   th=0.25 {m(0.5, 0.25, h): .6e}", flush=True)
print("done", flush=True)
