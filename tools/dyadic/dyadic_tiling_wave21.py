"""Wave 21 Agent 2 — dyadic induction barrier: numerical verification.

Verifies the EXACT dyadic tiling identities and quantifies why the induction
P+(x) ==> P+(2x) cannot close with allowed (F_theta-blind classical) inputs.

S(x)   = sum_{x/2<=p<=2x} (log p)/sqrt(p) * w_x(p),  w_x(t)=(1-|log(t/x)|/log2)_+
U(x)   = A*sqrt(x) - S(x),  A = (6*sqrt2-8)/log2
Tiling on [x,2x]: w_x + w_{2x} = 1  ==>  S(2x) = T1+T2+T3 - S(x), Ti >= 0
  T1 = sum_{(2x,4x]} a_p w_{2x}(p)   (new territory, windowed)
  T2 = sum_{[x,2x]}  a_p             (old territory, UNWINDOWED)
  T3 = sum_{[x/2,x]} a_p w_x(p)      (old territory, windowed)
Deficit identity: U(2x) + U(x) = A*sqrt(2x)+A*sqrt(x) - (T1+T2+T3) =: -D(x)
i.e. U(2x) = -U(x) - D(x).   <-- the IH quantity enters with MINUS sign.
"""
import numpy as np, math, sys, json

N = 30_000_000
print("sieving to", N, file=sys.stderr)
sieve = bytearray(b"\x01") * (N + 1)
sieve[0:2] = b"\x00\x00"
for i in range(2, int(N ** 0.5) + 1):
    if sieve[i]:
        sieve[i * i:N + 1:i] = b"\x00" * ((N - i * i) // i + 1)
primes = np.array([i for i in range(2, N + 1) if sieve[i]], dtype=np.float64)
logp = np.log(primes)
ap = logp / np.sqrt(primes)          # a_p = (log p)/sqrt(p), m=1
pre_ap = np.cumsum(ap)               # prefix sums for unwindowed blocks
pre_logp = np.cumsum(logp)           # prefix sums for theta
print(f"{len(primes)} primes", file=sys.stderr)

LOG2 = math.log(2.0)
A = (6 * math.sqrt(2.0) - 8) / LOG2
print(f"A = {A:.12f}", file=sys.stderr)

def w(x, t):
    return max(0.0, 1.0 - abs(math.log(t / x)) / LOG2)

def S_of(x):
    lo, hi = x / 2.0, 2.0 * x
    i0 = int(np.searchsorted(primes, lo, side="left"))
    i1 = int(np.searchsorted(primes, hi, side="right"))
    s = 0.0
    for p, a in zip(primes[i0:i1], ap[i0:i1]):
        s += a * (1.0 - abs(math.log(p / x)) / LOG2)
    return s

def T1_of(x):  # (2x,4x] with w_{2x}
    z = 2.0 * x
    i0 = int(np.searchsorted(primes, z, side="right"))
    i1 = int(np.searchsorted(primes, 2.0 * z, side="right"))
    s = 0.0
    for p, a in zip(primes[i0:i1], ap[i0:i1]):
        s += a * (1.0 - abs(math.log(p / z)) / LOG2)
    return s

def T3_of(x):  # [x/2,x] with w_x
    i0 = int(np.searchsorted(primes, x / 2.0, side="left"))
    i1 = int(np.searchsorted(primes, x, side="right"))
    s = 0.0
    for p, a in zip(primes[i0:i1], ap[i0:i1]):
        s += a * (1.0 - abs(math.log(p / x)) / LOG2)
    return s

def T2_of(x):  # [x,2x] unwindowed
    i0 = int(np.searchsorted(primes, x, side="left"))
    i1 = int(np.searchsorted(primes, 2.0 * x, side="right"))
    return float(pre_ap[i1 - 1] - (pre_ap[i0 - 1] if i0 > 0 else 0.0))

def theta(t):
    i1 = int(np.searchsorted(primes, t, side="right"))
    return float(pre_logp[i1 - 1]) if i1 > 0 else 0.0

scales = [100, 1000, 10_000, 100_000, 1_000_000, 5_000_000]
rows = []
print(f"{'x':>10} {'S(x)':>12} {'S(2x)':>12} {'tiling_err':>11} {'U(x)':>10} {'U(2x)':>10} "
      f"{'D(x)':>10} {'ident_err':>10} {'e2=T2main-T2':>13} {'E(2x)/sqrt(2x)':>14}")
for x in scales:
    Sx, S2x = S_of(x), S_of(2.0 * x)
    t1, t2, t3 = T1_of(x), T2_of(x), T3_of(x)
    tiling_err = abs(S2x - (t1 + t2 + t3 - Sx))
    Ux, U2x = A * math.sqrt(x) - Sx, A * math.sqrt(2 * x) - S2x
    Dx = (t1 + t2 + t3) - (A * math.sqrt(2 * x) + A * math.sqrt(x))
    ident_err = abs(U2x + Ux + Dx)
    m2 = 2.0 * (math.sqrt(2 * x) - math.sqrt(x))
    e2 = m2 - t2
    E2x = (theta(2 * x) - 2 * x) / math.sqrt(2 * x)
    rows.append(dict(x=x, Sx=Sx, S2x=S2x, tiling_err=tiling_err, Ux=Ux, U2x=U2x,
                     Dx=Dx, ident_err=ident_err, e2=e2, m2=m2, E2x=E2x))
    print(f"{x:>10} {Sx:>12.6f} {S2x:>12.6f} {tiling_err:>11.2e} {Ux:>10.4f} {U2x:>10.4f} "
          f"{Dx:>10.4f} {ident_err:>10.2e} {e2:>13.4f} {E2x:>14.4f}")

# Wrong-sign demonstration: the IH upper bound on U(x) becomes a LOWER bound on U(2x).
# Best O(1)-free upper bound derivable: -U(x) <= 0.9*sqrt(x) [trivial w<=1 + theta<=1.1t]
# and -D(x) <= K*E0/sqrt(x) [PNT]; show this "derivable bound" >> needed C.
print("\n--- induction gap: best derivable upper bound on U(2x) vs needed O(1) ---")
print(f"{'x':>10} {'needed C~':>10} {'0.9*sqrt(x)':>12} {'ratio':>10}")
for r in rows:
    x = r["x"]
    print(f"{x:>10} {0.68:>10.2f} {0.9*math.sqrt(x):>12.1f} {0.9*math.sqrt(x)/0.68:>10.0f}")

with open("agent2_dyadic.json", "w") as f:
    json.dump(dict(A=A, rows=rows), f)
print("\nwrote agent2_dyadic.json", file=sys.stderr)
