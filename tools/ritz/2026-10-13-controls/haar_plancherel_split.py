#!/usr/bin/env python3
"""
Haar vs Plancherel positivity in the finite Ritz model (Turbo, 2026-10-13).

E is assembled as a sum of per-jump matrices M_j = D_j^T diag(Pu*Ps*h*wt_j) D_j.
If every weight profile Pu*Ps*wt_j >= 0 pointwise, each M_j is PSD: E is a
TERMWISE Haar square (Haar positivity holds by construction). Then any
failure of the reduced form is purely Plancherel-side: the spectrum of the
reduced form (the discrete spectral measure) vs the constraint/W side.

Checks: (1) sign of rho(y) on the archimedean panel; (2) min eigenvalue of
each jump-class sum (arch-only, prime-only) at theta=0 and theta=0.25;
(3) min eigenvalue of total E (unreduced) for contrast with the reduced
neutral value. NUMERICAL.
"""
import sys, os
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import theta_control_ritz_4_spk4 as R
from ritz_analysis import setup, neutral_basis
from theta_control_ritz_4_spk4 import trigbasis

ys, wy = R.ypan()
rho = np.array([R.rho(y) for y in ys])
print(f"arch panel: {len(ys)} pts, rho min {rho.min():.6e}, rho max {rho.max():.6e}, quad w min {wy.min():.6e}")

FUNCS = trigbasis(30, 0.75, (0, 1, 2))
for th in (0.0, 0.25):
    d = setup(th, FUNCS, L=4.0, h=0.01, ycut=6.0)
    E, Wm, Cm = d["E"], d["Wm"], d["Cm"]
    kinds = {}
    for kind, y, n, M in d["per_jump"]:
        kinds.setdefault(kind, []).append(M)
    print(f"theta={th}: jump classes " + ", ".join(f"{k}x{len(v)}" for k, v in kinds.items()))
    for k, Ms in kinds.items():
        S = np.sum(Ms, axis=0)
        ev = np.linalg.eigvalsh((S + S.T) / 2)
        print(f"   {k:6s} sum: min eig {ev[0]: .6e}   (termwise PSD iff >= 0)")
    evE = np.linalg.eigvalsh((E + E.T) / 2)
    ev_red = neutral_basis(E, Wm, Cm)[1]
    print(f"   TOTAL E unreduced min eig {evE[0]: .6e};  reduced neutral min {ev_red[0]: .6e} (min-1 = {ev_red[0]-1: .3e})")
print("done")
