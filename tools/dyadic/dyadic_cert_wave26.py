"""Wave-26 Agent 2: interval-certified windowed dyadic remainder bound.

Target: B_w(x) = sum_{x/2<=n<=2x} (Lambda(n)/n^{3/2}) * w(n/x),
        w(t) = (1 - |log t|/log 2)_+,  x >= 2.
C_1 = 1 + sup_{x>=2} B_w(x)  (the "1" absorbs D_l/C_l(0) terms, wave-19).

PROOF STRUCTURE (all steps rigorous):
 (P1) For y = log x, each summand a_n*(1-|log n - y|/log2)_+ is a tent
      function of y: continuous, piecewise linear with kinks at
      y in {log(n/2), log n, log(2n)}. Finite sum (n<=400 prime powers
      suffice on [2,200]) => B_w is continuous piecewise-linear in y on
      [log2, log200]; attains its max at a kink or an endpoint.
      => sup_{[2,200]} B_w = max over kinks {n/2, n, 2n} cap [2,200]
         plus endpoints {2, 200}. Pure math, no computation.
 (P2) At each kink x=A/B (exact Fraction), activity n in [x/2,2x] is decided
      by EXACT integer arithmetic (2Bn>=A, Bn<=2A). Each term's interval
      enclosure uses the rigorous log table (exact-rational Taylor) and the
      wave-20 outward-rounded interval engine. ipos is rigorous because the
      true triangular weight is >= 0. Sum via proved isum_hi.
      => max of upper ends is a rigorous upper bound of the true sup.
 (P3) Tail x>200: B_w(x) <= B_unw(x) <= 3.3*sqrt(2/x) < 0.33 via partial
      summation with psi(t) <= 1.1t (Rosser-Schoenfeld 1962; algebra
      re-verified in-report; validated computationally on [100,1e6]).
      Robust: even K=1.11 (Chebyshev) gives < 0.334.

Also re-certifies the wave-19 unwindowed sup (<=0.822791) as machinery
validation (kinks + midpoints, since the unwindowed B is a step function).
"""
import json
import math
from fractions import Fraction

import numpy as np

from agent2_iv import (iadd, isub, imul, imul_nn, idiv_nn, iabs, ipos,
                       isqrt, isum_hi)

LOG = {int(k): (v[0], v[1]) for k, v in json.load(open("agent2_logtable.json")).items()}
LOG2 = LOG[2]

X0 = 200

# ---- exact prime-power catalog <= 400 -------------------------------------
def prime_powers(N):
    sieve = bytearray(b'\x01') * (N + 1)
    sieve[0:2] = b'\x00\x00'
    for i in range(2, int(N ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i:N + 1:i] = b'\x00' * ((N - i * i) // i + 1)
    primes = [i for i in range(2, N + 1) if sieve[i]]
    out = []
    for p in primes:
        pk = p
        while pk <= N:
            out.append((pk, p))
            if pk > N // p:
                break
            pk *= p
    return sorted(out)

PP = prime_powers(400)
print(f"prime powers <= 400: {len(PP)}")
NS = [n for n, p in PP]

# ---- a_n = Lambda(n)/n^{3/2} intervals ------------------------------------
A_IV = {}
for n, p in PP:
    nf = float(n)
    nn = (nf, nf)
    n32 = imul_nn(nn, isqrt(nn))          # n*sqrt(n) = n^{3/2}, encloses true
    A_IV[n] = idiv_nn(LOG[p], n32)        # (log p)/n^{3/2}

# ---- kink enumeration (exact Fractions) ------------------------------------
kinks = set()
for n, p in PP:
    for x in (Fraction(n, 2), Fraction(n, 1), Fraction(2 * n, 1)):
        if Fraction(2, 1) <= x <= Fraction(X0, 1):
            kinks.add(x)
kinks = sorted(kinks)
print(f"kinks in [2,{X0}]: {len(kinks)}")

def logx_iv(fr):
    A, B = fr.numerator, fr.denominator
    return isub(LOG[A], LOG[B])

def eval_windowed(fr):
    """Rigorous upper bound of B_w at x = fr (exact rational)."""
    A, B = fr.numerator, fr.denominator
    lx = logx_iv(fr)
    his = []
    for n, p in PP:
        if 2 * B * n >= A and B * n <= 2 * A:      # x/2 <= n <= 2x, exact
            d = iabs(isub(LOG[n], lx))            # |log(n/x)|
            w = ipos(isub((1.0, 1.0), idiv_nn(d, LOG2)))
            his.append(imul(A_IV[n], w)[1])
    return isum_hi(np.array(his))

results = {}
best = 0.0
bestx = None
for xk in kinks + [Fraction(2, 1), Fraction(X0, 1)]:
    hi = eval_windowed(xk)
    results[str(xk)] = hi
    if hi > best:
        best, bestx = hi, xk
print(f"windowed certified sup over [2,{X0}]: <= {best:.12f} at x={bestx}")

# ---- unwindowed validation (step function: kinks + midpoints) --------------
def eval_unwindowed_kink(fr):
    A, B = fr.numerator, fr.denominator
    his = [A_IV[n][1] for n, p in PP if 2 * B * n >= A and B * n <= 2 * A]
    return isum_hi(np.array(his))

def eval_unwindowed_mid(m):
    A, B = m.numerator, m.denominator
    his = [A_IV[n][1] for n, p in PP if 2 * B * n > A and B * n < 2 * A]
    return isum_hi(np.array(his))

bestu = 0.0
for xk in kinks:
    bestu = max(bestu, eval_unwindowed_kink(xk))
for i in range(len(kinks) - 1):
    m = (kinks[i] + kinks[i + 1]) / 2
    bestu = max(bestu, eval_unwindowed_mid(m))
bestu = max(bestu, eval_unwindowed_kink(Fraction(2, 1)),
            eval_unwindowed_kink(Fraction(X0, 1)))
print(f"unwindowed certified sup over [2,{X0}]: <= {bestu:.12f} (wave-19: 0.822791)")

# ---- tail ------------------------------------------------------------------
tail = 3.3 * math.sqrt(2.0 / X0)
print(f"analytic tail for x>{X0}: < {tail:.12f} (exact real value 0.33)")

C1 = 1.0 + max(best, 0.33)
print(f"CERTIFIED: sup B_w <= {max(best,0.33):.10f}  =>  C_1 <= {C1:.10f}")

json.dump({"sup_windowed": best, "sup_windowed_at": str(bestx),
           "sup_unwindowed": bestu, "tail": tail, "C1": C1,
           "n_kinks": len(kinks), "n_pp": len(PP)},
          open("agent2_dyadic_cert.json", "w"), indent=1)

# ---- mpmath cross-check: true kink values (50-digit) must lie below cert ---
import mpmath
mpmath.mp.dps = 50
mplog2 = mpmath.log(2)
viol = 0
for xk in kinks:
    x = mpmath.mpf(xk.numerator) / mpmath.mpf(xk.denominator)
    s = mpmath.mpf(0)
    for n, p in PP:
        if xk / 2 <= n <= 2 * xk:
            a = mpmath.log(p) / mpmath.mpf(n) ** mpmath.mpf('1.5')
            w = 1 - abs(mpmath.log(mpmath.mpf(n) / x)) / mplog2
            if w > 0:
                s += a * w
    if float(s) > results[str(xk)]:
        viol += 1
        print("VIOLATION at", xk, float(s), results[str(xk)])
print(f"mpmath cross-check: {len(kinks)} kinks, {viol} violations")

# ---- dense-grid sanity (NUMERICAL): no value above kink max ----------------
xs = np.concatenate([np.linspace(2, 20, 20000), np.linspace(20, 200, 20000)])
ns = np.array(NS, dtype=float)
avs = np.array([float(mpmath.log(p)) / n ** 1.5 for n, p in PP])
gmax, gx = 0.0, 0.0
lw2 = float(mplog2)
for x in xs:
    m = (ns >= x / 2) & (ns <= 2 * x)
    nn = ns[m]
    w = 1.0 - np.abs(np.log(nn / x)) / lw2
    w = np.clip(w, 0, 1)
    s = float((avs[m] * w).sum())
    if s > gmax:
        gmax, gx = s, x
print(f"dense-grid (NUMERICAL) max: {gmax:.10f} at x={gx:.4f} "
      f"(cert bound {best:.10f})")
assert gmax <= best + 1e-9, "grid exceeds certified bound!"
print("OK: grid max below certified bound")
