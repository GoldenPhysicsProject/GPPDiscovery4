#!/usr/bin/env python3
"""
Supplement to ground_vector.py: quantify the near-null bundle at theta=0
(how many reduced eigenvalues sit within tol of 1) and check whether the
jump fingerprint is stable when averaged over the bundle. NUMERICAL.
"""
import sys
import numpy as np
sys.path.insert(0, "/home/hatch/workspace")
sys.path.insert(0, "/home/hatch/workspace/user/files")
sys.path.insert(0, "/home/hatch/workspace/gppdiscovery4/tools/ritz")
from ritz_analysis import setup, neutral_basis
from theta_control_ritz_4_spk4 import trigbasis

FUNCS = trigbasis(30, 0.75, (0, 1, 2))

for L in (1.0, 2.0, 4.0):
    d = setup(0.0, FUNCS, L=L, h=0.01, ycut=6.0)
    B, ev, U = neutral_basis(d["E"], d["Wm"], d["Cm"])
    dev = np.abs(ev - 1.0)
    print(f"L={L}: neutral dim={B.shape[1]}")
    for tol in (1e-7, 1e-6, 1e-5, 1e-4, 1e-3):
        print(f"   #ev with |lam-1|<{tol:.0e}: {(dev < tol).sum()}")
    print(f"   smallest 8 (lam-1): {np.array2string(ev[:8]-1, precision=2, suppress_small=False)}")
    if L == 4.0:
        k = int((dev < 1e-6).sum())
        print(f"   bundle-averaged fingerprint over {k} near-null dirs (L=4):")
        Ea = Ep = 0.0
        pn = {}
        for j in range(k):
            c = B @ U[:, j]
            for kind, y, n, M in d["per_jump"]:
                val = float(c @ M @ c) / float(ev[j])
                if kind == "arch":
                    Ea += val
                else:
                    Ep += val
                    if n is not None:
                        pn[int(n)] = pn.get(int(n), 0.0) + val
        print(f"   arch avg {Ea/k:.6f}, prime avg {Ep/k:.6f}")
        for n, v in sorted(pn.items(), key=lambda kv: -kv[1])[:8]:
            print(f"      n={n:>4d}  avg {v/k:.6e}")
print("done")
