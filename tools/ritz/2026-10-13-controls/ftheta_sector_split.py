#!/usr/bin/env python3
"""
F_theta ARCHIMEDEAN-SEWING SECTOR TEST (Turbo's numerics, 2026-10-13).
Assigned by Codex on GPP-bridge PR #11 (comment 6152631994): "test explicit
theta branches in your paired entropy/zero-count rig; isolate where completed
Arch sewing breaks."   LABEL: NUMERICAL throughout.

Question: for the fake control F_th(s) = xi(s-th) xi(s+th), WHERE does the
negative mass of the completed Weil form live, in sector terms?
  - a single sector (archimedean-only or prime-only restriction already
    violates E >= W on the neutral space)?
  - only in cross-sector combinations / the fully sewn form E_a + E_p vs W?

Rig: copy-and-adapt of haar_plancherel_split.py / ground_vector.py (originals
untouched). Same basis trigbasis(30, 0.75, (0,1,2)) (243 functions), same
window L = 4.0, h = 0.01, ycut = 6.0, same pole-neutral reduction
(ritz_analysis.setup + neutral_basis). The reduced form tests E >= W on the
pole-neutral space; lam_min - 1 is the defect (zeta = 0, F_th ~ -0.052 th^2
banked). Sector split: E = E_a + E_p from the labeled per-jump list
(arch = Gamma-factor jumps rho(y) * 2cosh(th y);  prime = jumps at y = log n
with weight Lambda(n)/sqrt(n) * 2cosh(th log n)); W is the pole weight.

Diagnostics per theta:
  (0) full-space sector PSD check: min eig of E_a, E_p (termwise Haar-square
      construction predicts >= 0 for BOTH zeta and F_th);
  (1) neutral-sector bottoms: min eig of B' E_a B and B' E_p B (Rayleigh
      quotients against the W-metric, since B whitens Wm). NOTE these do NOT
      sum to lam_min; each measures how small that sector alone can go on the
      neutral space. < 1 in a sector = that sector alone carries defect mass;
  (2) sewn directions: bottom eigenvector c of the full reduced form, the
      pair-restricted bottom eigenvectors c_a, c_p; sector expectations
      (c' E_a c, c' E_p c, c' W c) for each direction -- the sector shares of
      each defect;
  (3) sewing cross term: for the normalized pair c_a, c_p, the cross-sector
      quadratic cross_a,p = c_a' E_p c_a + c_p' E_a c_p (see sector_report).

Theta sweep th in {0, 0.05, 0.1, 0.15, 0.2, 0.25} for the th^2 coefficient of
lam_min - 1 and the sector data as a function of th.

DH control (cheap reuse): setup_DH is copied verbatim from dh_control.py
(b-coefficients recomputed, Mellin spot check kept as a gate). DH is ENTIRE:
no pole weight W and no neutral space, so the matched diagnostics are the
full-space sector PSD check and the bottom generalized eigenvalue of
(E_sector, A_arch) >= 0 test of dh_control. DH's prime-jump weights b(n) have
MIXED signs, so its prime sector is not a termwise Haar square.
"""
import sys, time, os
import numpy as np
from math import log, sqrt, pi

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import theta_control_ritz_4_spk4 as R
from ritz_analysis import setup, neutral_basis
from theta_control_ritz_4_spk4 import trigbasis

FUNCS = trigbasis(30, 0.75, (0, 1, 2))
L, H, YCUT = 4.0, 0.01, 6.0
sym = lambda M: (M + M.T) / 2

def sector_report(th):
    t0 = time.time()
    d = setup(th, FUNCS, L=L, h=H, ycut=YCUT)
    E, Wm = d["E"], d["Wm"]
    Ea = sum(M for kind, y, n, M in d["per_jump"] if kind == "arch")
    Ep = sum(M for kind, y, n, M in d["per_jump"] if kind == "prime")
    resid = np.linalg.norm(E - (Ea + Ep)) / np.linalg.norm(E)
    B, ev, U = neutral_basis(E, Wm, d["Cm"])
    lam = ev[0]
    # noise scale: dense-eigvalsh backward error bound for the reduced form
    noise = 1e-13 * max(1.0, np.linalg.norm(B.T @ E @ B, 2))
    out = dict(theta=th, resid=resid, lam_min=lam, dim=B.shape[1],
               noise=noise, sectors={}, dirs={}, wres=None)
    out["eigmin_full_arch"] = float(np.linalg.eigvalsh(sym(Ea))[0])
    out["eigmin_full_prime"] = float(np.linalg.eigvalsh(sym(Ep))[0])
    out["eigmin_full_E"] = float(np.linalg.eigvalsh(sym(E))[0])
    out["wres"] = float(np.linalg.norm(B.T @ Wm @ B - np.eye(B.shape[1])))
    Ared = sym(B.T @ Ea @ B); Pred = sym(B.T @ Ep @ B)
    ea_a, Ua = np.linalg.eigh(Ared)
    ea_p, Up = np.linalg.eigh(Pred)
    out["sectors"] = dict(arch=float(ea_a[0]), prime=float(ea_p[0]))
    for tag, c in (("full", B @ U[:, 0]),
                  ("archdir", B @ Ua[:, 0]),
                  ("primedir", B @ Up[:, 0])):
        xa = float(c @ Ea @ c); xp = float(c @ Ep @ c); w = float(c @ Wm @ c)
        out["dirs"][tag] = dict(Ea=xa, Ep=xp, W=w, lam=(xa + xp) / w,
                                share_a=xa / (xa + xp), share_p=xp / (xa + xp))
    ca, cp = out["dirs"]["archdir"], out["dirs"]["primedir"]
    out["cross"] = ca["Ep"] + cp["Ea"]        # cross-sector sewing term
    out["sew_sum"] = ca["Ea"] + cp["Ep"] + out["cross"]
    out["time"] = time.time() - t0
    return out

print("=== F_theta sector split, matched zeta / F_theta rig "
      "(trigbasis(30,0.75,(0,1,2)), L=4.0, h=0.01, ycut=6.0) ===", flush=True)
res = {}
for th in (0.0, 0.25, 0.05, 0.1, 0.15, 0.2):
    r = sector_report(th)
    res[th] = r
    s, dr = r["sectors"], r["dirs"]
    print(f"\n--- theta = {th:.2f}   (neutral dim {r['dim']}, "
          f"sector-sum resid {r['resid']:.1e}, W-whitening resid {r['wres']:.1e}, "
          f"[{r['time']:.0f}s]) ---", flush=True)
    print(f"  full-space PSD: min eig E_arch {r['eigmin_full_arch']: .3e}   "
          f"E_prime {r['eigmin_full_prime']: .3e}   E {r['eigmin_full_E']: .3e}")
    print(f"  neutral-sector bottoms (vs W-metric): arch {s['arch']: .6e}   "
          f"prime {s['prime']: .6e}   [1 = threshold]")
    print(f"  FULL form: lam_min = {r['lam_min']:.8f}  (lam-1 = {r['lam_min']-1: .3e}; "
          f"noise floor ~{r['noise']:.0e})")
    for tag, q in dr.items():
        print(f"    dir {tag:9s}: E_arch {q['Ea']: .6e}  E_prime {q['Ep']: .6e}  "
              f"W {q['W']: .6e}  lam {(q['lam']): .6f}  "
              f"shares arch {q['share_a']*100:6.2f}% / prime {q['share_p']*100:6.2f}%")
    print(f"  sewing: c_a'E_p c_a + c_p'E_a c_p = {r['cross']: .6e};  "
          f"E_a[c_a]+cross+E_p[c_p] = {r['sew_sum']: .6e} < 2 iff a sewn defect exists")

ths = np.array(sorted(res))
defect = np.array([res[t]["lam_min"] - 1.0 for t in ths])
c2 = np.sum(ths**2 * defect) / np.sum(ths**4)
A = np.vstack([ths**2, ths**4]).T
(c2f, c4f), *_ = np.linalg.lstsq(A, defect, rcond=None)
print(f"\n=== theta^2 coefficient of lam_min - 1 (full sewn form) ===")
print(f"  one-param fit defect ~ c2 th^2:        c2 = {c2: .4f} "
      f"(banked pole-neutral law: -0.052)")
print(f"  two-param fit defect ~ c2 th^2 + c4 th^4: c2 = {c2f: .4f}, c4 = {c4f: .4f}")
for t, dfc in zip(ths, defect):
    print(f"    th={t:.2f}: lam-1 = {dfc: .3e}   (lam-1)/th^2 = "
          f"{(dfc/t**2 if t > 0 else float('nan')): .4f}")

# ================================================================ DH control
print("\n=== DH control (copied from dh_control.py; entire: no pole weight W, "
      "no neutral space) ===", flush=True)
from mpmath import mp, mpc, zeta as mpzeta, gamma as mpgamma, pi as mpi, e as me

mp.dps = 50
_chi = {1: 1, 2: 1j, 3: -1j, 4: -1}
_tau = sum(_chi[a] * me**(2j*mpi*a/5) for a in range(1, 5))
_eps = _tau / (1j*mp.sqrt(5))
XVAL = float((1j*(1-_eps)/(1+_eps)).real)

def dh_acoef(n):
    return [0.0, 1.0, XVAL, -XVAL, -1.0][n % 5]

def dh_F_mp(s):
    s = mpc(s)
    f = sum(dh_acoef(r) * mpzeta(s, mpc(r)/5) for r in range(1, 5)) / (5**s)
    return (5**(s/2)) * (mpi**(-(s+1)/2)) * mpgamma((s+1)/2) * f

NMAX = 500
_a = np.array([dh_acoef(n) for n in range(NMAX + 1)])
_ainv = np.zeros(NMAX + 1); _ainv[1] = 1.0
for n in range(2, NMAX + 1):
    s_ = 0.0
    for dd in range(2, n + 1):
        if n % dd == 0:
            s_ += _a[dd] * _ainv[n // dd]
    _ainv[n] = -s_
_alog = np.array([_a[n] * log(n) if n >= 1 else 0.0 for n in range(NMAX + 1)])
BCOEF = np.zeros(NMAX + 1)
for n in range(1, NMAX + 1):
    s_ = 0.0
    for dd in range(1, n + 1):
        if n % dd == 0:
            s_ += _alog[dd] * _ainv[n // dd]
    BCOEF[n] = s_

def Phi_DH(u):
    u = np.asarray(u, float); au = np.abs(u)
    y = np.exp(2*au)
    out = np.zeros_like(au)
    for n in range(1, 80):
        an = dh_acoef(n)
        if an == 0.0:
            continue
        out += an * (n/sqrt(5.0)) * np.exp(-pi*n*n*y/5.0)
    return 2.0*np.exp(1.5*au)*out

# Mellin gate (verbatim check from dh_control.py)
ug = np.arange(-9.0, 9.0 + 1e-9, 0.002)
Pg = Phi_DH(ug)
mel = [abs(np.trapz(Pg*np.exp((s-0.5)*ug), ug) - complex(dh_F_mp(s))) /
       abs(complex(dh_F_mp(s))) for s in [1.5+5j, 0.7+3j, 2.0+10j, 0.3+8j]]
print(f"  Phi_DH Mellin spot checks: max rel err {max(mel):.3e}", flush=True)
assert max(mel) <= 1e-8, "DH MELLIN GATE FAILED"

def rho_DH(y):
    return np.exp(-1.5*y) / (1.0 - np.exp(-2.0*y))

def setup_DH(funcs, L=4.0, h=0.01, ycut=6.0):
    """Verbatim copy of dh_control.setup_DH, extended to also return E_prime."""
    u = np.arange(-L, L + h/2, h)
    G = np.array([f(u) for f in funcs]); Pu = Phi_DH(u)
    m = len(funcs)
    E = np.zeros((m, m)); A = np.zeros((m, m)); P = np.zeros((m, m))
    ys, wy = R.ypan()
    jumps = [("arch", y, w*rho_DH(y), None) for y, w in zip(ys, wy)]
    jumps += [("prime", log(n), BCOEF[n]/sqrt(n), n) for n in range(2, NMAX + 1)
              if log(n) <= ycut and BCOEF[n] != 0.0]
    for kind, y, wt, _n in jumps:
        sh = u - y; Ps = Phi_DH(sh)
        if abs(wt)*np.max(np.abs(Pu*Ps)) < 1e-40:
            continue
        D = G - np.array([f(sh) for f in funcs])
        Mj = (D*(Pu*Ps*h*wt)) @ D.T
        E += Mj
        if kind == "arch":
            A += Mj
        else:
            P += Mj
    return dict(u=u, G=G, E=E, A=A, P=P)

def gen_min(M, S, ridge=1e-11):
    """verbatim from dh_control.py"""
    w, V = np.linalg.eigh(sym(S))
    keep = w > ridge*w.max()
    T = V[:, keep]/np.sqrt(w[keep])
    Mm = sym(T.T @ M @ T)
    ev = np.linalg.eigvalsh(Mm)
    return float(ev[0]), int(keep.sum()), float(w[keep].max()/w[keep].min())

ddh = setup_DH(FUNCS, L=L, h=H, ycut=YCUT)
mu_full, keptA, condA = gen_min(ddh["E"], ddh["A"])
mu_prime, _, _ = gen_min(ddh["P"], ddh["A"])
print(f"  DH full-space PSD: min eig E_arch {np.linalg.eigvalsh(sym(ddh['A']))[0]: .3e}   "
      f"E_prime {np.linalg.eigvalsh(sym(ddh['P']))[0]: .3e}   "
      f"E {np.linalg.eigvalsh(sym(ddh['E']))[0]: .3e}")
print(f"  DH bottoms vs arch scale A (arch kept {keptA}, cond {condA:.1e}): "
      f"mu(E) = {mu_full: .4e}   mu(E_prime only) = {mu_prime: .4e}   "
      f"[mu(E_arch only) = 1 by construction]")
print("done")
