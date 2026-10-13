"""
test23.py -- Tests 2 and 3: loss-of-positivity prime profiles and the
degree-two residual prediction.
"""
import numpy as np, time, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # staged: import sibling scripts locally
from math import log, sqrt
import theta_control_ritz_4_spk4 as R
import ritz_analysis as A

def profiles(thetas=(0.0, 0.1, 0.2, 0.3), kmax=30, L=4.0, h=0.01, ycut=6.0):
    """Per-jump contributions of the lowest Ritz direction at each theta."""
    out = {}
    for th in thetas:
        out[th] = A.loss_profile(th, L=L, h=h, ycut=ycut, kmax=kmax)
        print(f"theta={th}: lam_min={out[th]['lam_min']:.6f}", flush=True)
    return out

def summarize_loss(profs):
    """Rank prime-power jumps by loss of contribution vs theta=0."""
    base = {n: c for kind, y, n, c in profs[0.0]['per_jump'] if kind == 'prime'}
    arch0 = sum(c for kind, y, n, c in profs[0.0]['per_jump'] if kind == 'arch')
    rows = []
    for th in sorted(profs):
        if th == 0.0:
            continue
        cur = {n: c for kind, y, n, c in profs[th]['per_jump'] if kind == 'prime'}
        losses = []
        for n in set(base) | set(cur):
            losses.append((n, cur.get(n, 0.0) - base.get(n, 0.0), cur.get(n, 0.0)))
        arch = sum(c for kind, y, n, c in profs[th]['per_jump'] if kind == 'arch')
        losses.sort(key=lambda t: t[1])   # most negative first
        rows.append((th, losses, arch - arch0, profs[th]['lam_min']))
    return rows

def residual_profile(ns):
    """Predicted per-prime loss shape: (1-p^-2)/p for prime n, else None."""
    out = {}
    for n in ns:
        # factor n to check primality
        p = None
        m = n
        for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
            if m % q == 0:
                c = 0
                while m % q == 0:
                    m //= q; c += 1
                if m == 1:
                    p = q; pw = c
                break
        if p is not None and pw == 1:
            out[n] = (1 - p**-2)/p
    return out

if __name__ == "__main__":
    profs = profiles()
    rows = summarize_loss(profs)
    print("\n=== TEST 2: top loss channels per theta ===")
    for th, losses, darch, lm in rows:
        print(f"\ntheta={th} lam_min-1={lm-1:.3e}  d_arch={darch:.3e}")
        print("  top 12 prime losses (n, delta, contrib):")
        for n, d, c in losses[:12]:
            print(f"    n={n:4d}  delta={d:+.3e}  contrib={c:.3e}")
    # TEST 3b: correlation of loss shape with residual (1-p^-2)/p
    print("\n=== TEST 3b: loss-shape vs residual correlation ===")
    resp = residual_profile(
        [n for kind, y, n, c in profs[0.0]['per_jump'] if kind == 'prime'])
    for th, losses, darch, lm in rows:
        d = {n: dl for n, dl, c in losses}
        ns = [n for n in d if n in resp]
        x = np.array([d[n] for n in ns]); y = np.array([resp[n] for n in ns])
        # loss should be NEGATIVE multiple of residual: correlate -x with y
        cc = np.corrcoef(-x, y)[0, 1]
        print(f"theta={th}: corr(-loss, residual) = {cc:.4f}  (n_primes={len(ns)})")
    # shape stability across theta
    print("\n=== shape stability: corr of loss profiles across theta ===")
    dmap = {}
    for th, losses, darch, lm in rows:
        dmap[th] = {n: dl for n, dl, c in losses}
    ths = sorted(dmap)
    for i in range(len(ths)):
        for j in range(i+1, len(ths)):
            ns = [n for n in dmap[ths[i]] if n in dmap[ths[j]] and n in resp]
            x = np.array([dmap[ths[i]][n] for n in ns])
            y = np.array([dmap[ths[j]][n] for n in ns])
            print(f"  corr(loss@{ths[i]}, loss@{ths[j]}) = {np.corrcoef(x,y)[0,1]:.4f}")
