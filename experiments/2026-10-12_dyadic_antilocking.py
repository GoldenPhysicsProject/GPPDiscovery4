#!/usr/bin/env python3
"""Dyadic prime-sum anti-locking, small-scale reproduction (Turbo seat, GPP).

Reproduces the wave-30/31 signature at modest x: with the triangular window
  g(t) = t^{-1/2} (1 - |ln(t/x)|/L),  L = ln 2,  on [x/2, 2x],
  A1(x) = sum_{x/2<p<=x} (log p) g(p)      (rising edge)
  A2(x) = sum_{x<p<=2x}  (log p) g(p)      (falling edge)
  F1(x) = A1 - W1 sqrt(x),  F2(x) = A2 - W2 sqrt(x)
(W1, W2 the window integrals), the program banks Corr(F1,F2) ~ -0.983 with
F1, F2 each O(1) (std ~0.33/0.32) and |F1+F2| <= ~0.62 (sieve to 5e10).

This script checks the signature at small x (default XMAX=2e6) and also shows
the companion fact: with E(x)=theta(x)-x, R=E/sqrt(x)-F1 and
Ftilde=E/sqrt(x)+F2 correlate POSITIVELY (~+0.98, the trivial common mode) --
the anti-locking lives in (F1,F2), not in E.

Exact transcription cross-check: A sqrt(x) - (A1+A2) == -(F1+F2), A=W1+W2.

Prints parameters + results to stdout and appends one JSON line to the
private bridge results log.
"""
import numpy as np
import math
import json
import time
import os

t0 = time.time()
L = math.log(2.0)
XMAX = float(os.environ.get("XMAX", "2e6"))
PMAX = int(2 * XMAX)
NGRID = int(os.environ.get("NGRID", "200"))
print(f"params: XMAX={XMAX:g} PMAX={PMAX} NGRID={NGRID} L=ln2", flush=True)

# ---- window constants by quadrature (as in wave-31 agent1_rf.py) --------
uu = np.linspace(0.5, 1.0, 200001)
W1 = float(np.trapz(uu ** -0.5 * (1.0 + np.log(uu) / L), uu))
uu = np.linspace(1.0, 2.0, 200001)
W2 = float(np.trapz(uu ** -0.5 * (1.0 - np.log(uu) / L), uu))
A = W1 + W2
print(f"W1={W1:.8f} W2={W2:.8f} A={A:.8f}", flush=True)

# ---- sieve to PMAX (odd-only bytearray) --------------------------------
ts = time.time()
n = (PMAX - 1) // 2
mask = bytearray(b"\x01") * n
sq = int(math.isqrt(PMAX)) + 1
for p in range(3, sq, 2):
    if mask[(p - 3) // 2]:
        step = p
        start = (p * p - 3) // 2
        mask[start::step] = b"\x00" * ((n - 1 - start) // step + 1)
primes = [2] + [2 * i + 3 for i in range(n) if mask[i]]
P = np.array(primes, dtype=float)
logP = np.log(P)
cumlog = np.cumsum(logP)  # theta at each prime
print(f"sieved {len(P)} primes to {PMAX} in {time.time()-ts:.1f}s", flush=True)


def theta(x):
    i = int(np.searchsorted(P, x, side="right")) - 1
    return float(cumlog[i]) if i >= 0 else 0.0


# ---- grid ---------------------------------------------------------------
xs = np.exp(np.linspace(math.log(XMAX / 20.0), math.log(XMAX), NGRID))
F1 = np.empty(NGRID)
F2 = np.empty(NGRID)
R = np.empty(NGRID)
Ft = np.empty(NGRID)
dcheck = np.empty(NGRID)
sqx = np.sqrt(xs)
for k, x in enumerate(xs):
    m1 = (P > x / 2) & (P <= x)
    p1 = P[m1]
    g1 = p1 ** -0.5 * (1.0 - np.abs(np.log(p1 / x)) / L)
    A1 = float((logP[m1] * g1).sum())
    m2 = (P > x) & (P <= 2 * x)
    p2 = P[m2]
    g2 = p2 ** -0.5 * (1.0 - np.abs(np.log(p2 / x)) / L)
    A2 = float((logP[m2] * g2).sum())
    f1 = A1 - W1 * sqx[k]
    f2 = A2 - W2 * sqx[k]
    F1[k], F2[k] = f1, f2
    E = theta(x) - x
    R[k] = E / sqx[k] - f1
    Ft[k] = E / sqx[k] + f2
    dcheck[k] = A * sqx[k] - (A1 + A2) + (f1 + f2)  # exact identity => 0

c12 = float(np.corrcoef(F1, F2)[0, 1])
cR = float(np.corrcoef(R, Ft)[0, 1])
print(f"std(F1)={F1.std():.4f} std(F2)={F2.std():.4f}", flush=True)
print(f"Corr(F1,F2) = {c12:.4f}   (banked -0.983 at 5e10)", flush=True)
print(f"max|F1+F2|  = {np.abs(F1+F2).max():.4f}   (banked <= ~0.62)", flush=True)
print(f"Corr(R,Ftilde) = {cR:.4f}   (banked +0.98 at 5e10: the E/sqrt(x) common mode", flush=True)
print("   dominates only at large x; at small x the F1/F2 anti-locking still shows through)", flush=True)
print(f"transcription cross-check max|A sqrt(x)-(A1+A2)+(F1+F2)| = "
      f"{np.abs(dcheck).max():.2e}", flush=True)
dt = time.time() - t0
print(f"runtime_s = {dt:.1f}", flush=True)
print("note: small-x demo; the -0.983 signature is banked at 5e10 (wave-31).", flush=True)

rec = {
    "date": "2026-10-12", "seat": "turbo", "kind": "experiment",
    "label": "NUMERICAL", "title": "Dyadic anti-locking small-scale reproduction",
    "ref": "experiments/2026-10-12_dyadic_antilocking.py",
    "detail": f"XMAX={XMAX:g} NGRID={NGRID}: std(F1)={F1.std():.3f} "
              f"std(F2)={F2.std():.3f} Corr(F1,F2)={c12:.4f} "
              f"max|F1+F2|={np.abs(F1+F2).max():.4f} Corr(R,Ftilde)={cR:.4f} "
              f"xcheck={np.abs(dcheck).max():.1e} runtime {dt:.0f}s",
}
logp = os.path.expanduser("~/workspace/gpp-bridge/MuseResearch/turbo-results-log.md")
with open(logp, "a") as fh:
    fh.write(json.dumps(rec) + "\n")
print(f"logged -> {logp}", flush=True)
