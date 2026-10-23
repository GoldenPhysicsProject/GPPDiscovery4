#!/usr/bin/env python3
"""
Extract the ground (bottom) neutral eigenvector of the Ritz/Weil reduced form
for zeta (theta=0): the direction on which the balance closes exactly
(reduced lam_min = 1). Also contrast with theta=0.25 and compute the
jump-class fingerprint. Label: NUMERICAL.

Method: setup via ritz_analysis (243-dim trig basis trigbasis(30,0.75,(0,1,2)),
h=0.01, ycut=6.0). neutral_basis whitens Wm, imposes the pole-neutral
constraints, and diagonalizes the reduced form. c = B @ U[:,0] is the
full-space coefficient vector; v(u) = sum_j c_j f_j(u) is the test function.
"""
import sys, time, math, os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import theta_control_ritz_4_spk4 as R
from ritz_analysis import setup, neutral_basis
from theta_control_ritz_4_spk4 import trigbasis

FUNCS = trigbasis(30, 0.75, (0, 1, 2))
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "out")

def ground(theta, L):
    t0 = time.time()
    d = setup(theta, FUNCS, L=L, h=0.01, ycut=6.0)
    B, ev, U = neutral_basis(d["E"], d["Wm"], d["Cm"])
    c = B @ U[:, 0]
    v = d["G"].T @ c
    print(f"theta={theta} L={L}: lam_min-1={ev[0]-1: .3e} "
          f"lam_2-1={ev[1]-1: .3e} dim={B.shape[1]} [{time.time()-t0:.0f}s]",
          flush=True)
    return d, B, ev, U, c, v

def mass_bins(u, v, h):
    w = v*v*h
    w = w/w.sum()
    tot = lambda cond: w[cond].sum()
    return {
        "center(|u|<0.5)": tot(np.abs(u) < 0.5),
        "log2-band(0.55-0.85)": tot((np.abs(u) >= 0.55) & (np.abs(u) <= 0.85)),
        "mid(0.85-1.5)": tot((np.abs(u) > 0.85) & (np.abs(u) <= 1.5)),
        "tail(>1.5)": tot(np.abs(u) > 1.5),
    }

def even_odd(u, v):
    idx = {round(x, 8): i for i, x in enumerate(u)}
    ve = np.zeros_like(v); vo = np.zeros_like(v)
    for i, x in enumerate(u):
        j = idx.get(round(-x, 8))
        if j is not None:
            ve[i] = 0.5*(v[i]+v[j]); vo[i] = 0.5*(v[i]-v[j])
        else:
            ve[i] = v[i]
    n = lambda a: float(np.sqrt((a*a).sum()))
    return n(ve)/n(v), n(vo)/n(v)

def corr(a, b):
    a = a - a.mean(); b = b - b.mean()
    return float((a @ b)/math.sqrt((a @ a)*(b @ b)))

def fit_gauss(u, v):
    # grid search sigma; a = least squares scale
    best = (None, None, -1.0)
    for sig in np.linspace(0.1, 3.0, 300):
        g = np.exp(-u*u/(2*sig*sig))
        a = float((v @ g)/(g @ g))
        r = abs(float((v @ (a*g))/math.sqrt((v @ v)*(a*g @ (a*g)))))
        if r > best[2]:
            best = (a, sig, r)
    return best

def common_corr(d1, v1, d2, v2):
    # correlate on intersection of grids by interpolation
    lo = max(d1["u"][0], d2["u"][0]); hi = min(d1["u"][-1], d2["u"][-1])
    uc = np.arange(lo, hi + 1e-9, 0.01)
    a = np.interp(uc, d1["u"], v1); b = np.interp(uc, d2["u"], v2)
    return corr(a, b)

print("=== STEP 1: ground vectors theta=0, L=1,2,4 ===", flush=True)
D = {}
for L in (1.0, 2.0, 4.0):
    d, B, ev, U, c, v = ground(0.0, L)
    D[L] = (d, c, v, ev)
    np.savetxt(f"{OUT}/ground_vector_L{L:.0f}_th0.csv", np.column_stack([d["u"], v]),
               delimiter=",", header="u,v", comments="")
    print(f"  mass bins L={L}: " +
          ", ".join(f"{k}={x:.3f}" for k, x in mass_bins(d["u"], v, 0.01).items()))
    ee, oo = even_odd(d["u"], v)
    print(f"  L={L}: ||even||/||v||={ee:.4f}  ||odd||/||v||={oo:.4f}")

for a, b in ((1.0, 2.0), (2.0, 4.0), (1.0, 4.0)):
    r = common_corr(D[a][0], D[a][2], D[b][0], D[b][2])
    print(f"shape corr L={a} vs L={b} on common support: {r:.6f}")

print("=== STEP 2: closed-form candidate overlaps (L=4, th=0) ===", flush=True)
d, c, v, ev = D[4.0]
u = d["u"]; PhiF = d["PhiF"]
phi = PhiF(u)
print(f"corr(v, Phi(u))            = {corr(v, phi):+.6f}")
print(f"corr(v, Phi(u)*cosh(u/2))  = {corr(v, phi*np.cosh(u/2)):+.6f}")
phip = (PhiF(u+0.01)-PhiF(u-0.01))/0.02
print(f"corr(v, Phi'(u))           = {corr(v, phip):+.6f}")
a, sig, r = fit_gauss(u, v)
print(f"best Gaussian: a={a:+.6f} sigma={sig:.4f} |corr|={r:.6f}")
print(f"v(0)={v[len(v)//2]:+.6f} (sign check vs a={a:+.6f})")

print("=== STEP 3: theta=0.25, L=4 contrast ===", flush=True)
d25, B25, ev25, U25, c25, v25 = ground(0.25, 4.0)
np.savetxt(f"{OUT}/ground_vector_L4_th025.csv", np.column_stack([d25["u"], v25]),
           delimiter=",", header="u,v", comments="")
r = corr(np.interp(u, d25["u"], v25), v)
print(f"corr(v_th0, v_th025) on [-4,4] = {r:+.6f}")
print("mass th=0.25: " + ", ".join(f"{k}={x:.3f}" for k, x in mass_bins(d25["u"], v25, 0.01).items()))
ee, oo = even_odd(d25["u"], v25)
print(f"th=0.25: ||even||/||v||={ee:.4f} ||odd||/||v||={oo:.4f}")
print(f"v_th025(0)={v25[len(v25)//2]:+.6f}; E-exp={c25@d25['E']@c25:.6f} Wm-exp={c25@d25['Wm']@c25:.6f}")

print("=== STEP 4: jump fingerprint of the th=0, L=4 ground vector ===", flush=True)
Eexp = c @ d["E"] @ c; Wexp = c @ d["Wm"] @ c
print(f"ground: E-exp={Eexp:.10f}  Wm-exp={Wexp:.10f}  (E-W)exp={Eexp-Wexp: .3e}")
arch = [(y, n, M) for kind, y, n, M in d["per_jump"] if kind == "arch"]
prime = [(y, n, M) for kind, y, n, M in d["per_jump"] if kind == "prime"]
Ea = sum(float(c @ M @ c) for _, _, M in arch)
Ep = sum(float(c @ M @ c) for _, _, M in prime)
print(f"arch class total  {Ea:.10f} ({Ea/Eexp*100:.2f}% of E)")
print(f"prime class total {Ep:.10f} ({Ep/Eexp*100:.2f}% of E)")
pj = sorted([(float(c @ M @ c), int(n)) for _, n, M in prime if n is not None], reverse=True)
print("top prime jumps by E-value (n, value):")
for val, n in pj[:10]:
    print(f"   n={n:>4d}  {val:.6e}")
print("done")
