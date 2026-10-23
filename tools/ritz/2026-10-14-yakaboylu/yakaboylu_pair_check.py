#!/usr/bin/env python3
"""
Yakaboylu reflection-pairing check (Turbo, 2026-10-14). LABEL: NUMERICAL
throughout. This script makes NO claim about RH.

Context: Yakaboylu (arXiv:2408.15135v17) reduces RH to positivity of a
reflection pairing over the zeros. Turbo's audit of v17 (bridge
MuseResearch/turbo-note-2026-10-14-yakaboylu-v17-audit.md) observes that
in the biorthonormal zero basis his intertwiner W acts as the reflection
rho -> 1 - conj(rho), and predicts that for the shifted product
F_theta(s) = xi(s-theta) xi(s+theta) the Weil zero-side pairing splits,
pair by mirror pair, into a symmetric (positive) part and an
antisymmetric (defect) part. This script verifies that decomposition
numerically with the existing rigs and quantifies it.

Conventions (identical to weil_window_scan.py / ftheta_sector_split.py /
ritz_analysis.py, which see):
  - Basis: trigbasis(30, 0.75, (0,1,2)) (243 real functions of u).
  - Window L = 4.0, grid h = 0.01 on u in [-L, L], ycut = 6.0.
  - Scalar transform of a real test function f = sum_j c_j g_j at a point z
    (this is ritz_analysis.F_matrix, per basis function):
        ghat_theta(z) = int Phi_theta(u) f(u) e^{(z-1/2) u} du
    evaluated as the rig's rectangle sum on the grid. Phi_theta is the
    F_theta archimedean profile (make_labeled; Phi_0 = zeta's Phi).
  - F_theta zero pairs from zeta_ordinates_140.npy (gamma_k, k = 1..140):
        a_k = 1/2 + theta + i gamma_k,   b_k = 1/2 - theta + i gamma_k,
    with b_k = 1 - conj(a_k).
  - For real f, ghat(1-a) = conj(ghat(b)) and ghat(1-b) = conj(ghat(a))
    (Phi_theta and f are real), so the Weil pair contribution is
        pair(a,b) = ghat(a) ghat(1-a) + ghat(b) ghat(1-b)
                  = 2 Re[ ghat(a) conj(ghat(b)) ].
    With m = (ghat(a)+ghat(b))/2, d = (ghat(a)-ghat(b))/2, the
    parallelogram identity gives the decomposition under test:
        pair(a,b) = 2 ( |m|^2 - |d|^2 ).
    S = sum_pairs 2|m|^2 (symmetric mass), A = sum_pairs 2|d|^2
    (antisymmetric/defect mass); zero-side pair sum = S - A.
  - Relation to the rig: ritz_analysis.zero_side_matrix builds the full
    zero-side matrix Z (all zeros, both half-planes) from exactly these
    cross terms; for a coefficient vector c, c' Z c = 2 (S - A) with the
    sums truncated at the same zero set. Checked below at theta = 0.25.

Test functions (coefficient vectors in the trig basis, Wm-normalized
neutral Ritz directions of the theta = 0.25 rig, c' Wm c = 1):
  - CARRIER: bottom neutral direction at theta = 0.25, L = 4
    (completed form value lam_min - 1 = -3.254e-3, banked: the direction
    where the F_theta completed Weil form is negative).
  - CONTROL: top neutral direction at theta = 0.25 (largest reduced
    eigenvalue; completed form strongly positive).
Both coefficient vectors are then held FIXED while theta varies; only
Phi_theta and the pair locations move, exactly as in the rig's zero side.
Profile-family note: for theta > 0, Phi_theta is the archimedean profile
of the product xi(s-theta) xi(s+theta) (max |Phi_theta| ~ 0.320 on the
grid, nearly theta-independent), while the theta = 0 row uses zeta's own
profile Phi_0 (max ~ 0.894). The theta = 0 row is therefore on a
different overall scale and is NOT the theta -> 0 limit of the F_theta
family's S; it is included for the d = 0 sanity check, which is
profile-independent (a = b pointwise at theta = 0).

CAVEAT (read before quoting any number): S - A is the ZERO-SIDE pair sum
only, truncated at 140 pairs (gamma <= 302.70). It is NOT the full Weil
form. The archimedean/prime-side balance (the explicit formula, the
E - W completed form of the rigs) is a separate, existing computation;
nothing here re-derives or replaces it.
"""
import sys, time, os
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from ritz_analysis import setup, neutral_basis, make_labeled, F_matrix, \
    zero_side_matrix
from theta_control_ritz_4_spk4 import trigbasis

HERE = os.path.dirname(os.path.abspath(__file__))
FUNCS = trigbasis(30, 0.75, (0, 1, 2))
L, H, YCUT = 4.0, 0.01, 6.0
GAMMA = np.load(f"{HERE}/zeta_ordinates_140.npy")
NPAIR = len(GAMMA)

print("=== Yakaboylu reflection-pairing check (NUMERICAL) ===", flush=True)
print(f"basis dim {len(FUNCS)}, L={L}, h={H}, ycut={YCUT}; "
      f"{NPAIR} zero pairs, gamma in [{GAMMA[0]:.4f}, {GAMMA[-1]:.4f}]",
      flush=True)

# ---------------------------------------------------------------- rig at 0.25
t0 = time.time()
d25 = setup(0.25, FUNCS, L=L, h=H, ycut=YCUT)
B, ev, U = neutral_basis(d25["E"], d25["Wm"], d25["Cm"])
u = d25["u"]; G = d25["G"]
FUN = {
    "carrier (bottom neutral dir, th=0.25)": (B @ U[:, 0], float(ev[0]) - 1.0),
    "control (top neutral dir, th=0.25)":    (B @ U[:, -1], float(ev[-1]) - 1.0),
}
print(f"\nneutral dim {B.shape[1]}; rig setup [{time.time()-t0:.0f}s]", flush=True)
for name, (c, lam) in FUN.items():
    print(f"  {name}: completed-form value lam-1 = {lam: .6e} "
          f"(Wm-normalization c'Wm c = {float(c @ d25['Wm'] @ c):.6f})",
          flush=True)

_Pu_cache = {}
def Pu_of(theta):
    if theta not in _Pu_cache:
        if theta == 0.25:
            _Pu_cache[theta] = d25["Pu"]
        else:
            PhiF = make_labeled(theta, YCUT)[0]
            _Pu_cache[theta] = PhiF(u)
    return _Pu_cache[theta]

def ghat(theta, c, zpts):
    """Rig transform: int Phi_theta(u) f(u) e^{(z-1/2)u} du (rectangle sum)."""
    w = Pu_of(theta) * (c @ G) * H
    return w @ np.exp(np.outer(u, np.asarray(zpts, complex) - 0.5))

def pair_terms(theta, c, gammas):
    """Per-pair (2|m|^2, 2|d|^2) arrays for the given ordinates."""
    a = 0.5 + theta + 1j * gammas
    b = 0.5 - theta + 1j * gammas
    ga, gb = ghat(theta, c, a), ghat(theta, c, b)
    m = 0.5 * (ga + gb); dd = 0.5 * (ga - gb)
    return 2.0 * np.abs(m) ** 2, 2.0 * np.abs(dd) ** 2, ga, gb

# ------------------------------------------------- Task 2: identity, th=0.25
print("\n--- pair identity check, theta = 0.25, first 50 pairs ---", flush=True)
g50 = GAMMA[:50]
a50 = 0.5 + 0.25 + 1j * g50
b50 = 0.5 - 0.25 + 1j * g50
for name, (c, _lam) in FUN.items():
    ga = ghat(0.25, c, a50); gb = ghat(0.25, c, b50)
    g1a = ghat(0.25, c, 1.0 - a50)          # evaluated directly, no shortcut
    g1b = ghat(0.25, c, 1.0 - b50)
    pair_direct = ga * g1a + gb * g1b       # definition
    pair_cross = 2.0 * np.real(ga * np.conj(gb))
    m = 0.5 * (ga + gb); dd = 0.5 * (ga - gb)
    pair_md = 2.0 * (np.abs(m) ** 2 - np.abs(dd) ** 2)
    scale = np.max(np.abs(ga) ** 2 + np.abs(gb) ** 2)
    e1 = np.max(np.abs(pair_direct - pair_cross))
    e2 = np.max(np.abs(pair_direct - pair_md))
    e3 = np.max(np.abs(g1a - np.conj(gb)))  # reality relation ghat(1-a)=conj(ghat(b))
    print(f"  {name}:", flush=True)
    print(f"    max |pair_direct - 2Re[ga conj(gb)]|  = {e1:.3e} "
          f"(rel. to max|ga|^2+|gb|^2: {e1/scale:.3e})", flush=True)
    print(f"    max |pair_direct - 2(|m|^2-|d|^2)|    = {e2:.3e} "
          f"(rel. {e2/scale:.3e})", flush=True)
    print(f"    max |ghat(1-a) - conj(ghat(b))|       = {e3:.3e} "
          f"(rel. to max|ghat|: {e3/np.max(np.abs(gb)):.3e})", flush=True)

# --------------------------------------- Task 4 (first): theta = 0 sanity
print("\n--- theta = 0 sanity: a = b on the line, d must vanish ---", flush=True)
for name, (c, _lam) in FUN.items():
    _s, _a, ga, gb = pair_terms(0.0, c, GAMMA)
    dd = 0.5 * (ga - gb)
    print(f"  {name}: max |d| over {NPAIR} pairs = {np.max(np.abs(dd)):.3e} "
          f"(max |ghat| = {np.max(np.abs(ga)):.3e})", flush=True)

# ------------------------------------------------- Task 3: S / A theta table
print(f"\n--- S (symmetric mass), A (antisymmetric mass), S-A "
      f"({NPAIR} pairs) ---", flush=True)
THETAS = [0.01, 0.02, 0.05, 0.1, 0.25]   # fit set (all in the F_theta family)
table = {}
for name, (c, _lam) in FUN.items():
    print(f"  {name}:", flush=True)
    print(f"    {'theta':>6} {'S':>14} {'A':>14} {'S-A':>14} {'A/S':>10} "
          f"{'A/theta^2':>12}", flush=True)
    rows = {}
    for th in [0.0] + THETAS:
        s_arr, a_arr, _ga, _gb = pair_terms(th, c, GAMMA)
        S, A = float(s_arr.sum()), float(a_arr.sum())
        rows[th] = (S, A)
        print(f"    {th:6.2f} {S:14.8e} {A:14.8e} {S-A:14.8e} "
              f"{(A/S if S != 0 else float('nan')):10.3e} "
              f"{(A/th**2 if th > 0 else 0.0):12.6e}", flush=True)
    table[name] = rows
    ths = np.array(THETAS); As = np.array([rows[t][1] for t in THETAS])
    alpha = float(np.sum(ths**2 * As) / np.sum(ths**4))     # A ~ alpha th^2
    pred = alpha * ths**2
    slopes = [float(np.log(As[i+1]/As[i]) / np.log(ths[i+1]/ths[i]))
              for i in range(len(ths) - 1)]
    print(f"    fit A ~ alpha*theta^2 (through origin): alpha = {alpha:.6e}; "
          f"pred/actual = {np.array2string(pred/As, precision=4)}", flush=True)
    print(f"    successive log-log slopes dlogA/dlogth: "
          f"{['%.3f' % s for s in slopes]}  (2.000 = pure theta^2)", flush=True)

# --------------------------------------- rig cross-check: c'Zc vs 2(S-A)
print("\n--- rig cross-check at theta = 0.25: zero_side_matrix vs pairs ---",
      flush=True)
zeros = 0.5 + 1j * GAMMA
Fa, Fb = F_matrix(d25["PhiF"], FUNCS, u, H, zeros + 0.25, zeros - 0.25)
Z = zero_side_matrix(Fa, Fb)
for name, (c, _lam) in FUN.items():
    zc = float(c @ Z @ c)
    S, A = table[name][0.25]
    print(f"  {name}: c'Zc = {zc: .8e};  2(S-A) = {2*(S-A): .8e}; "
          f"rel diff = {abs(zc - 2*(S-A))/abs(zc):.3e}", flush=True)
print("  (Z sums over ALL zeros, both half-planes; the pair sum covers the",
      flush=True)
print("   upper half-plane, hence the factor 2. Same 140-zero truncation.)",
      flush=True)

# partial-sum stability at 50 pairs (carrier, theta = 0.25)
s50, a50m, _ga, _gb = pair_terms(0.25, FUN["carrier (bottom neutral dir, th=0.25)"][0], g50)
S, A = table["carrier (bottom neutral dir, th=0.25)"][0.25]
print(f"\ncarrier theta=0.25, first 50 pairs only: S50 = {s50.sum():.8e}, "
      f"A50 = {a50m.sum():.8e}  (vs 140 pairs: S = {S:.8e}, A = {A:.8e})",
      flush=True)
print("done", flush=True)
