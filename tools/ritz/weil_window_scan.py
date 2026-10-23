#!/usr/bin/env python3
"""
Weil small-support window scan (Turbo, 2026-10-13).

Claude Code's PR #11 queue item 1: reproduce small-support Weil positivity
numerically. Take the pole-neutral Weil/Ritz form on test functions supported
in [-L, L] (log scale; support ratio lambda = e^L), compute the reduced
neutral minimum eigenvalue lam_min - 1 for zeta (theta=0) and for the
shifted product F_theta, and find where each first goes negative as the
window grows. The gap between the two crossings is the invariant.

Machinery: the project's verified Ritz rig (theta_control_ritz_4_spk4 via
ritz_analysis), the same setup whose theta-scan reproduces the banked
-0.052*theta^2 law and matched the 2026-10-08 baseline digit-for-digit
from the seeded Discovery4 workbench.

Label: NUMERICAL (finite basis, finite grid).
"""
import sys, time, os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ritz_analysis import setup, neutral_basis
from theta_control_ritz_4_spk4 import trigbasis

FUNCS = trigbasis(30, 0.75, (0, 1, 2))
Ls = [0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
THETAS = [0.0, 0.25]

print(f"basis dim {len(FUNCS)}; h=0.01 ycut=6.0", flush=True)
print("     L    lambda      lam_min-1 (th=0)   lam_min-1 (th=0.25)", flush=True)
for L in Ls:
    row = []
    for th in THETAS:
        t0 = time.time()
        d = setup(th, FUNCS, L=L, h=0.01, ycut=6.0)
        evals = neutral_basis(d["E"], d["Wm"], d["Cm"])[1]
        row.append(float(evals[0]) - 1.0)
        print(f"  [L={L} th={th}: {time.time()-t0:.0f}s]", flush=True)
    print(f"  {L:4.2f}  {np.exp(L):8.2f}   {row[0]: .6e}      {row[1]: .6e}", flush=True)
print("done", flush=True)
