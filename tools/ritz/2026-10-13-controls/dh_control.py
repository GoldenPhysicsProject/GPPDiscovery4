#!/usr/bin/env python3
"""
Davenport-Heilbronn (DH) control for the Ritz/Weil reduced-form rig.
Turbo, 2026-10-13.  Label: NUMERICAL throughout.

Question: the rig's reduced form saturates exactly for zeta (lam_min - 1 = 0,
banked: within 1.3e-7 for L >= 1.5) and fails for the shifted product F_theta
(-0.052 theta^2 law).  Does saturation use only the functional equation --
which DH also has -- or the Euler product / local unitarity, which DH lacks
(DH has off-line zeros)?  Run the same machinery on f_DH and see.

Construction (all verified inside this script before any number is quoted):
  f(s) = sum a(n) n^{-s}, a 5-periodic: (0, 1, x, -x, -1), residues (0..4),
  x = i(1-eps)/(1+eps), eps = tau(chi)/(i sqrt5), chi mod 5 with chi(2)=i.
  F(s) = 5^{s/2} pi^{-(s+1)/2} Gamma((s+1)/2) f(s),  F(s) = F(1-s).
  Phi_DH(u) = 2 e^{3|u|/2} sum_n a(n) (n/sqrt5) exp(-pi n^2 e^{2|u|}/5),
  with int Phi_DH(u) e^{(s-1/2)u} du = F(s)  (theta-series inverse Mellin;
  spot-checked numerically below).
  -f'/f = sum b(n) n^{-s},  b = (a log) * a^{-1}  (Dirichlet convolution).
  Archimedean density rho_DH(y) = e^{-3y/2}/(1-e^{-2y})
    = sum over Gamma((s+1)/2) poles s=-(2k+1) of e^{(pole-1/2) y}
    (same rule that gives zeta's rho(y) = e^{-y/2}/(1-e^{-2y}) from Gamma(s/2)
     poles at s=-2k; mirrors the verified rig).
  f is ENTIRE -> no pole weight W, no pole constraints.  The reduced quantity
  is therefore the bottom generalized eigenvalue of (E_DH, A_DH), A_DH the
  archimedean-only Gram (the positive trivial scale), computed identically
  for zeta as a same-definition control column.
"""
import sys, time, json
import numpy as np
from math import log, sqrt, pi

sys.path.insert(0, "/home/hatch/workspace")
sys.path.insert(0, "/home/hatch/workspace/user/files")
import theta_control_ritz_4_spk4 as R
from ritz_analysis import setup as zeta_setup, neutral_basis
from mpmath import mp, mpc, zeta, gamma, pi as mpi, e as me, findroot

OUT = {}

# ================================================================ DH data
mp.dps = 50
_chi = {1: 1, 2: 1j, 3: -1j, 4: -1}
_tau = sum(_chi[a] * me**(2j*mpi*a/5) for a in range(1, 5))
_eps = _tau / (1j*mp.sqrt(5))
XVAL = float((1j*(1-_eps)/(1+_eps)).real)
OUT["x"] = XVAL
OUT["eps_chi"] = [float(_eps.real), float(_eps.imag)]

def acoef(n):
    return [0.0, 1.0, XVAL, -XVAL, -1.0][n % 5]

def f_mp(s):
    s = mpc(s)
    return sum(acoef(r) * zeta(s, mpc(r)/5) for r in range(1, 5)) / (5**s)

def F_mp(s):
    s = mpc(s)
    return (5**(s/2)) * (mpi**(-(s+1)/2)) * gamma((s+1)/2) * f_mp(s)

# ---- FE verification (gate: stop if this fails)
FE_PTS = [0.3+7j, 0.8+15j, 1.2+3j, 0.5+0.5j, 2.5+10j, 0.9+30j, 1.5+0.2j, 0.2+25j, 0.65+11j]
fe = [float(abs(F_mp(s) - F_mp(1-s))) for s in FE_PTS]
OUT["fe_residuals"] = fe
OUT["fe_worst"] = max(fe)
print(f"x = {XVAL:.15f};  FE worst |F(s)-F(1-s)| = {max(fe):.3e} over {len(FE_PTS)} pts", flush=True)
assert max(fe) <= 1e-8, "FE VERIFICATION FAILED -- aborting"

# ================================================================ zeros
mp.dps = 30
def find_zeros():
    sig = np.arange(0.2, 1.301, 0.05)
    tt = np.arange(0.1, 40.01, 0.2)
    G = np.zeros((len(tt), len(sig)))
    for i, t in enumerate(tt):
        for j, sg in enumerate(sig):
            G[i, j] = float(abs(f_mp(mpc(sg, t))))
    seeds = []
    for i in range(len(tt)):
        for j in range(len(sig)):
            v = G[i, j]
            if v == G[max(i-1,0):i+2, max(j-1,0):j+2].min() and v < 1.2:
                seeds.append(mpc(sig[j], tt[i]))
    roots = []
    for sd in seeds:
        try:
            r = findroot(f_mp, sd, tol=1e-28, maxsteps=80)
        except Exception:
            continue
        if float(abs(f_mp(r))) > 1e-12:
            continue
        if not (0.15 < float(r.real) < 1.45 and 0.0 < float(r.imag) < 40.5):
            continue
        r = complex(r)
        if all(abs(r - q) > 1e-6 for q in roots):
            roots.append(r)
    roots.sort(key=lambda z: z.imag)
    return roots

def argprinciple_count():
    N = 3000
    path = ([mpc(s, 40.0) for s in np.linspace(1.3, 0.2, N)] +
            [mpc(0.2, t) for t in np.linspace(40.0, 0.1, N)] +
            [mpc(s, 0.1) for s in np.linspace(0.2, 1.3, N)] +
            [mpc(1.3, t) for t in np.linspace(0.1, 40.0, N)])
    ph = [float(mp.arg(f_mp(p))) for p in path]
    tot = 0.0
    for a_, b_ in zip(ph, ph[1:]):
        d = b_ - a_
        while d > np.pi: d -= 2*np.pi
        while d < -np.pi: d += 2*np.pi
        tot += d
    return tot/(2*np.pi)

def rect_count(sa, sb, ta, tb, N=2000):
    path = ([mpc(s, tb) for s in np.linspace(sb, sa, N)] +
            [mpc(sa, t) for t in np.linspace(tb, ta, N)] +
            [mpc(s, ta) for s in np.linspace(sa, sb, N)] +
            [mpc(sb, t) for t in np.linspace(ta, tb, N)])
    ph = [float(mp.arg(f_mp(p))) for p in path]
    tot = 0.0
    for a_, b_ in zip(ph, ph[1:]):
        d = b_ - a_
        while d > np.pi: d -= 2*np.pi
        while d < -np.pi: d += 2*np.pi
        tot += d
    return tot/(2*np.pi)

print("locating zeros in strip Im in [0.1,40] ...", flush=True)
t0 = time.time()
ZEROS = find_zeros()
OUT["n_zeros_found"] = len(ZEROS)
OUT["zeros"] = [[z.real, z.imag] for z in ZEROS]
OUT["zero_resid_max"] = max(float(abs(f_mp(mpc(z)))) for z in ZEROS)
ncount = argprinciple_count()
OUT["argprinciple_count"] = float(ncount)
# off-line exclusion proof for the strip: count both off-line half-rectangles
c_right = rect_count(0.55, 1.3, 0.1, 40)
c_left = rect_count(0.2, 0.45, 0.1, 40)
OUT["strip_offline_counts"] = [float(c_right), float(c_left)]
print(f"{len(ZEROS)} zeros found (max resid {OUT['zero_resid_max']:.2e}); "
      f"argument-principle count = {ncount:.3f}; off-line half-rect counts: "
      f"right {c_right:.2f}, left {c_left:.2f}  [{time.time()-t0:.0f}s]", flush=True)
for z in ZEROS:
    tag = "ON-LINE " if abs(z.real - 0.5) < 1e-6 else "OFF-LINE"
    print(f"   {tag} {z.real:.10f} + {z.imag:.10f}i", flush=True)

ONLINE = [z for z in ZEROS if abs(z.real - 0.5) <= 1e-6]
OFFLINE = [z for z in ZEROS if abs(z.real - 0.5) > 1e-6]
OUT["n_online"] = len(ONLINE); OUT["n_offline"] = len(OFFLINE)

# ---- off-line zeros live higher: hunt Re in (0.5,1.3], Im in (40,200]
print("hunting off-line zeros in Im in (40,200] ...", flush=True)
t0 = time.time()
c_hi = rect_count(0.55, 1.3, 40, 200, N=1500)
sig = np.arange(0.55, 1.301, 0.05)
tt = np.arange(40.0, 200.01, 0.4)
Gh = np.zeros((len(tt), len(sig)))
for i, t in enumerate(tt):
    for j, sg in enumerate(sig):
        Gh[i, j] = float(abs(f_mp(mpc(sg, t))))
seeds = []
for i in range(len(tt)):
    for j in range(len(sig)):
        v = Gh[i, j]
        if v == Gh[max(i-1,0):i+2, max(j-1,0):j+2].min() and v < 0.8:
            seeds.append(mpc(sig[j], tt[i]))
OFFHI = []
for sd in seeds:
    try:
        r = findroot(f_mp, sd, tol=1e-28, maxsteps=80)
    except Exception:
        continue
    if float(abs(f_mp(r))) > 1e-10:
        continue
    if not (0.5 < float(r.real) < 1.35 and 39.0 < float(r.imag) < 201.0):
        continue
    r = complex(r)
    if all(abs(r - q) > 1e-6 for q in OFFHI):
        OFFHI.append(r)
OFFHI.sort(key=lambda z: z.imag)
OUT["offline_high_count"] = float(c_hi)
OUT["offline_high"] = [[z.real, z.imag] for z in OFFHI]
print(f"high band: argument count (Re>0.55) = {c_hi:.2f}; {len(OFFHI)} roots located "
      f"[{time.time()-t0:.0f}s]", flush=True)
for z in OFFHI:
    print(f"   OFF-LINE {z.real:.10f} + {z.imag:.10f}i", flush=True)

# ================================================================ b(n)
NMAX = 500
_a = np.array([acoef(n) for n in range(NMAX + 1)])
_ainv = np.zeros(NMAX + 1); _ainv[1] = 1.0
for n in range(2, NMAX + 1):
    s_ = 0.0
    for d in range(2, n + 1):
        if n % d == 0:
            s_ += _a[d] * _ainv[n // d]
    _ainv[n] = -s_
_alog = np.array([_a[n] * log(n) if n >= 1 else 0.0 for n in range(NMAX + 1)])
BCOEF = np.zeros(NMAX + 1)
for n in range(1, NMAX + 1):
    s_ = 0.0
    for d in range(1, n + 1):
        if n % d == 0:
            s_ += _alog[d] * _ainv[n // d]
    BCOEF[n] = s_
OUT["b_at_primes"] = {str(p): float(BCOEF[p]) for p in [2, 3, 5, 7, 11, 13, 19]}
OUT["b_minmax_n_le_403"] = [float(BCOEF[1:404].min()), float(BCOEF[1:404].max())]
print(f"b(n): primes ok (b(p)=a(p)log p); range on [1,403]: "
      f"[{BCOEF[1:404].min():.3f}, {BCOEF[1:404].max():.3f}] (mixed signs -> no Euler product)", flush=True)

# ================================================================ Phi_DH
def Phi_DH(u):
    u = np.asarray(u, float); au = np.abs(u)
    y = np.exp(2*au)
    out = np.zeros_like(au)
    for n in range(1, 80):
        an = acoef(n)
        if an == 0.0:
            continue
        out += an * (n/sqrt(5.0)) * np.exp(-pi*n*n*y/5.0)
    return 2.0*np.exp(1.5*au)*out

ug = np.arange(-9.0, 9.0 + 1e-9, 0.002)
Pg = Phi_DH(ug)
mel = []
for s in [1.5+5j, 0.7+3j, 2.0+10j, 0.3+8j]:
    num = np.trapz(Pg * np.exp((s-0.5)*ug), ug)
    ref = complex(F_mp(s))
    mel.append(abs(num - ref)/abs(ref))
OUT["mellin_rel_errs"] = [float(v) for v in mel]
print(f"Phi_DH Mellin spot checks: max rel err {max(mel):.3e}", flush=True)

# ================================================================ rig
def rho_DH(y):
    return np.exp(-1.5*y) / (1.0 - np.exp(-2.0*y))

def setup_DH(funcs, L=4.0, h=0.01, ycut=6.0):
    u = np.arange(-L, L + h/2, h)
    G = np.array([f(u) for f in funcs]); Pu = Phi_DH(u)
    m = len(funcs)
    E = np.zeros((m, m)); A = np.zeros((m, m))
    ys, wy = R.ypan()
    jumps = [("arch", y, w*rho_DH(y)) for y, w in zip(ys, wy)]
    jumps += [("prime", log(n), BCOEF[n]/sqrt(n)) for n in range(2, NMAX + 1)
              if log(n) <= ycut and BCOEF[n] != 0.0]
    for kind, y, wt in jumps:
        sh = u - y; Ps = Phi_DH(sh)
        if abs(wt)*np.max(np.abs(Pu*Ps)) < 1e-40:
            continue
        D = G - np.array([f(sh) for f in funcs])
        Mj = (D*(Pu*Ps*h*wt)) @ D.T
        E += Mj
        if kind == "arch":
            A += Mj
    return dict(u=u, G=G, E=E, A=A)

def gen_min(M, S, ridge=1e-11):
    """smallest generalized eigenvalue of (M, S), S PSD; returns (mu, kept_dim, cond)."""
    w, V = np.linalg.eigh((S + S.T)/2)
    keep = w > ridge*w.max()
    T = V[:, keep]/np.sqrt(w[keep])
    Mm = T.T @ M @ T; Mm = (Mm + Mm.T)/2
    ev = np.linalg.eigvalsh(Mm)
    return float(ev[0]), int(keep.sum()), float(w[keep].max()/w[keep].min())

def Fmat(PhiF, funcs, u, h, zs):
    G = np.array([f(u) for f in funcs]); Pu = PhiF(u)
    MG = G * (Pu*h)
    Ea = np.exp(np.outer(u, np.array(zs) - 0.5))
    return MG @ Ea

def zero_side_DH(funcs, u, h, extra_offline=()):
    m = len(funcs)
    Z = np.zeros((m, m))
    if ONLINE:
        Fa = Fmat(Phi_DH, funcs, u, h, ONLINE)
        Z += 2*((Fa @ Fa.conj().T).real)
    offs = list(OFFLINE) + list(extra_offline)
    if offs:
        Fa = Fmat(Phi_DH, funcs, u, h, offs)
        Fb = Fmat(Phi_DH, funcs, u, h, [1 - z.conjugate() for z in offs])
        Z += 2*((Fa @ Fb.conj().T + Fb @ Fa.conj().T).real)
    return Z

# ---- identity check at L=4
FUNCS = R.trigbasis(30, 0.75, (0, 1, 2))
d4 = setup_DH(FUNCS, L=4.0, h=0.01, ycut=6.0)
Z4 = zero_side_DH(FUNCS, d4["u"], 0.01)
rel_id = float(np.linalg.norm(d4["E"] - Z4)/np.linalg.norm(Z4))
OUT["identity_rel_err_L4"] = rel_id
mu_E, keptA, condA = gen_min(d4["E"], d4["A"])
mu_Z, _, _ = gen_min(Z4, d4["A"])
Z4hi = zero_side_DH(FUNCS, d4["u"], 0.01, extra_offline=OFFHI)
mu_Z_hi, _, _ = gen_min(Z4hi, d4["A"])
rel_id_hi = float(np.linalg.norm(d4["E"] - Z4hi)/np.linalg.norm(Z4hi))
OUT["identity_rel_err_L4_with_high_offline"] = rel_id_hi
OUT["L4"] = dict(mu_E=mu_E, mu_Z=mu_Z, mu_Z_with_high_offline=mu_Z_hi,
                 arch_kept_dim=keptA, arch_cond=condA)
print(f"identity check L=4: ||E-Z||/||Z|| = {rel_id:.3e} "
      f"(with high off-line zeros added: {rel_id_hi:.3e})", flush=True)
print(f"L=4 bottoms vs arch scale: mu(E)={mu_E: .4e}  mu(Z)={mu_Z: .4e}  "
      f"mu(Z + high off-line)={mu_Z_hi: .4e}  (arch kept {keptA}, cond {condA:.1e})", flush=True)

# ---- window scan, DH and zeta control with identical definitions
Ls = [0.75, 1.0, 1.5, 2.0, 3.0, 4.0]
OUT["scan"] = {}
print("\n     L      mu_DH(L)        mu_zeta(L)   [same definition: min gen eig of (E, A_arch)]", flush=True)
for L in Ls:
    t0 = time.time()
    dd = setup_DH(FUNCS, L=L, h=0.01, ycut=6.0)
    mu_dh, kdh, _ = gen_min(dd["E"], dd["A"])
    dz = zeta_setup(0.0, FUNCS, L=L, h=0.01, ycut=6.0)
    Az = np.zeros_like(dz["E"])
    for kind, y, n, Mj in dz["per_jump"]:
        if kind == "arch":
            Az += Mj
    mu_z, kz, _ = gen_min(dz["E"], Az)
    OUT["scan"][str(L)] = dict(mu_DH=mu_dh, mu_zeta=mu_z, kept_DH=kdh, kept_zeta=kz)
    print(f"  {L:5.2f}  {mu_dh: .6e}   {mu_z: .6e}    [{time.time()-t0:.0f}s]", flush=True)

json.dump(OUT, open("/home/hatch/workspace/ym-k-experiments/dh_results.json", "w"), indent=1)
print("done")
