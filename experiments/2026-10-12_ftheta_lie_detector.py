#!/usr/bin/env python3
"""F_theta lie-detector demo (Turbo seat, GPP).

Daniel's falsifier: F_theta(s) = xi(s-theta)*xi(s+theta), 0<theta<1/2.
It satisfies xi's functional equation EXACTLY and is built from xi's Euler
factors, yet has explicit OFF-LINE zeros at rho +- theta for every xi-zero rho.
Hence no argument using only {functional equation + Euler-product structure}
can prove RH: any such argument would also "prove" RH for F_theta, which is
false. This script exhibits that numerically.

Checks:
 1. FE: |F(s) - F(1-s)| ~ 0 at test points (relative).
 2. Off-line zeros: with rho1 = 1/2 + i*gamma1 the first xi-zero (refined by
    findroot), |F(1/2-theta+i*gamma1)| and |F(1/2+theta+i*gamma1)| are ~0,
    while |F(1/2+i*gamma1)| is NOT ~0 (xi itself vanishes only on the line).
 3. Line blindness: a sign-change scan of Re F(1/2+it) over t in [10,20] is
    compared against a 2D findroot sweep of the box [0.3,0.7]x[10,20]; the
    off-line pair is invisible to any line-only procedure.

Prints parameters + results to stdout and appends one JSON line to the
private bridge results log.
"""
import mpmath as mp
import json
import time
import os

t0 = time.time()
mp.mp.dps = 30
THETA = float(os.environ.get("THETA", "0.1"))
assert 0 < THETA < 0.5, "need 0<theta<1/2"

print(f"params: dps={mp.mp.dps} theta={THETA}", flush=True)


def xi(s):
    s = mp.mpc(s)
    return mp.mpc(1) / 2 * s * (s - 1) * mp.power(mp.pi, -s / 2) * mp.gamma(s / 2) * mp.zeta(s)


def F(s):
    return xi(s - THETA) * xi(s + THETA)


# ---- 1. functional equation -------------------------------------------
fe_pts = [0.3 + 10j, 0.7 + 14.13j, 2 + 5j, -1 + 7j, 1.5 + 0.5j]
fe_worst = 0.0
for s in fe_pts:
    a, b = F(s), F(1 - s)
    rel = abs(a - b) / max(abs(a), mp.mpf("1e-300"))
    fe_worst = max(fe_worst, float(rel))
print(f"FE check: max |F(s)-F(1-s)|/|F(s)| over {len(fe_pts)} pts = {fe_worst:.3e}", flush=True)

# ---- 2. refine first xi zero, exhibit off-line F zeros -----------------
gamma1 = mp.findroot(lambda t: mp.re(xi(mp.mpc("0.5") + mp.j * t)), 14.13)
gamma1 = float(gamma1)
print(f"refined gamma1 = {gamma1:.12f}", flush=True)

z_lo = 0.5 - THETA + 1j * gamma1
z_hi = 0.5 + THETA + 1j * gamma1
z_on = 0.5 + 1j * gamma1
v_lo, v_hi = abs(F(z_lo)), abs(F(z_hi))
v_xi_on = abs(xi(z_on))
v_F_on = abs(F(z_on))
print(f"|F(1/2-theta+i g1)| = {mp.nstr(v_lo, 4)}   (off-line zero)", flush=True)
print(f"|F(1/2+theta+i g1)| = {mp.nstr(v_hi, 4)}   (off-line zero)", flush=True)
print(f"|xi(1/2+i g1)|      = {mp.nstr(v_xi_on, 4)}   (xi on-line zero)", flush=True)
print(f"|F(1/2+i g1)|       = {mp.nstr(v_F_on, 4)}   (NOT a zero: xi has none there)", flush=True)

# ---- 3. 2D findroot sweep of the box -----------------------------------
found = []
for rs in (0.35, 0.45, 0.55, 0.65):
    for ts in (11.0, 13.0, 15.0, 17.0, 19.0):
        try:
            r = mp.findroot(lambda s: F(s), rs + 1j * ts, tol=mp.mpf("1e-25"), maxsteps=60)
            if 0.3 <= mp.re(r) <= 0.7 and 10 <= mp.im(r) <= 20 and abs(F(r)) < 1e-8:
                key = (round(float(mp.re(r)), 6), round(float(mp.im(r)), 6))
                if key not in [k for k, _ in found]:
                    found.append((key, float(abs(F(r)))))
        except Exception:
            pass
print(f"2D sweep: {len(found)} distinct zeros of F in [0.3,0.7]x[10,20]:", flush=True)
for (rs_, ts_), vv in sorted(found):
    print(f"   s = {rs_:.6f} + {ts_:.6f} i   |F| = {vv:.2e}"
          + ("   <-- OFF the critical line" if abs(rs_ - 0.5) > 1e-4 else ""), flush=True)

# ---- 4. line-only scan --------------------------------------------------
N = 4001
ts = [10 + 10 * i / (N - 1) for i in range(N)]
prev = mp.re(F(mp.mpc("0.5") + mp.j * ts[0]))
sign_changes = 0
for t in ts[1:]:
    cur = mp.re(F(mp.mpc("0.5") + mp.j * t))
    if prev == 0:
        sign_changes += 1
    elif cur * prev < 0:
        sign_changes += 1
    prev = cur
print(f"line scan: Re F(1/2+it) sign changes on t in [10,20]: {sign_changes}", flush=True)
print("moral: F obeys xi's FE exactly and is built from xi's Euler factors, yet has",
      flush=True)
print("off-line zeros at rho+-theta. FE + Euler-product structure alone CANNOT prove RH.",
      flush=True)
print("Any argument that would also 'work' for F_theta is broken. (Daniel's lie detector.)",
      flush=True)

dt = time.time() - t0
print(f"runtime_s = {dt:.1f}", flush=True)

rec = {
    "date": "2026-10-12", "seat": "turbo", "kind": "experiment",
    "label": "NUMERICAL", "title": "F_theta lie-detector demo",
    "ref": "experiments/2026-10-12_ftheta_lie_detector.py",
    "detail": f"theta={THETA} dps=30: FE rel err {fe_worst:.2e}; "
              f"|F(1/2+-theta+i g1)|={float(v_lo):.2e},{float(v_hi):.2e}; "
              f"|F(1/2+i g1)|={float(v_F_on):.2e} (nonzero); "
              f"2D box zeros found={len(found)} (all off-line); "
              f"line-scan sign changes [10,20]={sign_changes}; runtime {dt:.0f}s",
}
logp = os.path.expanduser("~/workspace/gpp-bridge/MuseResearch/turbo-results-log.md")
with open(logp, "a") as fh:
    fh.write(json.dumps(rec) + "\n")
print(f"logged -> {logp}", flush=True)
