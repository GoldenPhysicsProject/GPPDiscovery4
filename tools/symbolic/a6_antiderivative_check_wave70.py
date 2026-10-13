#!/usr/bin/env python3
"""
Wave 70, Agent 3 -- Check 3: symbolic (sympy) verification of A6's
antiderivatives Psi_r, Psi_f and the closed forms A_rise, A_fall.

Claims verified by computer algebra (label: PROVED, computer-algebra
verification):
  Psi_r'(u) = g_rise(u)/u,  Psi_f'(u) = g_fall(u)/u   (symbolic diff, simplify to 0)
  A_rise = Psi_r(1) - Psi_r(1/2) = 2 + (2 sqrt2 - 4)/log 2
  A_fall = Psi_f(2) - Psi_f(1)   = -2 + (4 sqrt2 - 4)/log 2
plus independent 30-digit numeric evaluation of int g(u)/u du on each half
(mpmath quadrature of the smooth integrand g(u)/u) against the closed forms.
"""
import sympy as sp
from mpmath import mp, mpf, log, sqrt

u = sp.symbols("u", positive=True)
L = sp.log(2)

g_rise = sp.sqrt(u) * (1 + sp.log(u) / L)
g_fall = sp.sqrt(u) * (1 - sp.log(u) / L)

Psi_r = 2 * sp.sqrt(u) + (2 * sp.sqrt(u) * sp.log(u) - 4 * sp.sqrt(u)) / L
Psi_f = 2 * sp.sqrt(u) - (2 * sp.sqrt(u) * sp.log(u) - 4 * sp.sqrt(u)) / L

for name, Psi, g in [("rise", Psi_r, g_rise), ("fall", Psi_f, g_fall)]:
    diff = sp.simplify(sp.diff(Psi, u) - g / u)
    print(f"Psi_{name}' - g/u  simplifies to: {diff}")
    assert diff == 0

Ar_sym = sp.simplify(Psi_r.subs(u, 1) - Psi_r.subs(u, sp.Rational(1, 2)))
Af_sym = sp.simplify(Psi_f.subs(u, 2) - Psi_f.subs(u, 1))
Ar_closed = 2 + (2 * sp.sqrt(2) - 4) / L
Af_closed = -2 + (4 * sp.sqrt(2) - 4) / L
print("Psi_r(1)-Psi_r(1/2) - [2+(2sqrt2-4)/log2] =",
      sp.simplify(Ar_sym - Ar_closed))
print("Psi_f(2)-Psi_f(1)   - [-2+(4sqrt2-4)/log2] =",
      sp.simplify(Af_sym - Af_closed))
assert sp.simplify(Ar_sym - Ar_closed) == 0
assert sp.simplify(Af_sym - Af_closed) == 0

print("\nNumeric values (sympy, 30 digits):")
print("A_rise =", sp.N(Ar_closed, 30))
print("A_fall =", sp.N(Af_closed, 30))
print("A_rise+A_fall =", sp.N(Ar_closed + Af_closed, 30))
print("A = (6sqrt2-8)/log2 =", sp.N((6 * sp.sqrt(2) - 8) / L, 30))

# Independent numeric quadrature of int g(u)/u du (smooth integrand, mpmath)
mp.dps = 30
Lm = log(2)

def w_rise(uu):
    return 1 + log(uu) / Lm

def w_fall(uu):
    return 1 - log(uu) / Lm

Ar_num = mp.quad(lambda uu: sqrt(uu) * w_rise(uu) / uu,
                 [mpf("0.5"), mpf(1)])
Af_num = mp.quad(lambda uu: sqrt(uu) * w_fall(uu) / uu,
                 [mpf(1), mpf(2)])
Ar_cf = 2 + (2 * sqrt(2) - 4) / Lm
Af_cf = -2 + (4 * sqrt(2) - 4) / Lm
print("\nmpmath quadrature int g(u)/u du vs closed form:")
print(f"rise: quad={mp.nstr(Ar_num, 25)} closed={mp.nstr(Ar_cf, 25)} "
      f"diff={mp.nstr(Ar_num - Ar_cf, 5)}")
print(f"fall: quad={mp.nstr(Af_num, 25)} closed={mp.nstr(Af_cf, 25)} "
      f"diff={mp.nstr(Af_num - Af_cf, 5)}")
print("\nALL SYMBOLIC CHECKS PASSED")
