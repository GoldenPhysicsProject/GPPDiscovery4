#!/usr/bin/env python3
"""
bundle_mechanism.py -- adjudicate two readings of the near-null bundle in the
Ritz/Weil reduced form at theta=0. NUMERICAL throughout.

Readings:
  (P) Claude Code / CCM prolate: the bundle is the finite shadow of the
      prolate near-null space (functions whose Phi-weighted Mellin transform
      ~vanishes at the zeta zeros) + construction kernel of the quotiented
      pole direction. Bottom eigenvalue positive but exponentially small.
  (I) Turbo identity: shadow of an exact operator identity E == W on a
      subspace of the neutral space.

Machinery identical to ground_vector.py / weil_window_scan.py:
trigbasis(30, 0.75, (0,1,2)) (243-dim), h=0.01, theta=0.

Tests:
 1. Bundle dimension vs L (ycut=6) and vs ycut (L=4), tolerances 1e-5/1e-6/1e-7.
 2. Counting: bundle dim vs N(T) zero-counting laws; sampling-map singular
    values (the mechanism check); Fourier-zero signature of bundle vectors
    and bundle-averaged spectral density S(gamma) dips at zeta ordinates.
 3. Per-jump bundle action ||P_B M P_B||_F/||M||_F (+ neutral-restricted
    version and bundle trace share) vs y at L=4.
 4. Perturbation split: ycut 4/8, cosh-only constraint (sinh dropped),
    alt basis polys (0,1). Principal cosines vs the reference bundle and
    average-projector core dimension.
"""
import sys, time
import numpy as np

sys.path.insert(0, "/home/hatch/workspace")
sys.path.insert(0, "/home/hatch/workspace/user/files")
sys.path.insert(0, "/home/hatch/workspace/gppdiscovery4/tools/ritz")
import theta_control_ritz_4_spk4 as R
from ritz_analysis import setup, neutral_basis, F_matrix
from theta_control_ritz_4_spk4 import trigbasis

OUT = "/home/hatch/workspace/ym-k-experiments"
FUNCS = trigbasis(30, 0.75, (0, 1, 2))
FUNCS_ALT = trigbasis(30, 0.75, (0, 1))
TOLS = (1e-5, 1e-6, 1e-7)
H = 0.01

def log(s=""):
    print(s, flush=True)

# ------------------------------------------------------------------ zeros
def zero_ordinates(K=140):
    path = f"{OUT}/zeta_ordinates_{K}.npy"
    try:
        return np.load(path)
    except Exception:
        from mpmath import mp, zetazero
        mp.dps = 30
        g = np.array([float(zetazero(k).imag) for k in range(1, K + 1)])
        np.save(path, g)
        return g

t0 = time.time()
GAM = zero_ordinates()
log(f"zeta ordinates cached: {len(GAM)} zeros, gamma_max={GAM[-1]:.1f} "
    f"[{time.time()-t0:.0f}s]")
def NofT(T):
    return int((GAM <= T).sum())

# ------------------------------------------------------------------ helpers
def bundle_idx(ev, tol):
    return np.where(np.abs(ev - 1.0) < tol)[0]

def run(L, ycut, funcs=FUNCS, cm_rows=None):
    d = setup(0.0, funcs, L=L, h=H, ycut=ycut)
    Cm = d["Cm"] if cm_rows is None else d["Cm"][cm_rows]
    B, ev, U = neutral_basis(d["E"], d["Wm"], Cm)
    na = sum(1 for j in d["per_jump"] if j[0] == "arch")
    np_ = sum(1 for j in d["per_jump"] if j[0] == "prime")
    return d, B, ev, U, na, np_

def counts_str(ev):
    return "  ".join(f"tol {t:g}: {len(bundle_idx(ev, t))}" for t in TOLS)

# ======================================================================
log("\n================ TEST 1: bundle dimension scaling ================")
log("--- L scan at ycut=6 (theta=0) ---")
log(f"{'L':>4} {'ndim':>5} {'#arch':>6} {'#prime':>6} {'lam_min-1':>11}  bundle counts")
scan = {}
for L in (1.0, 1.5, 2.0, 3.0, 4.0):
    d, B, ev, U, na, np_ = run(L, 6.0)
    scan[L] = dict(d=d, B=B, ev=ev, U=U, ndim=B.shape[1])
    log(f"{L:>4.1f} {B.shape[1]:>5d} {na:>6d} {np_:>6d} {ev[0]-1:>11.2e}  {counts_str(ev)}")
    if L not in (2.0, 4.0):
        del scan[L]["d"]  # free the per-jump matrices; keep spectrum
log("--- ycut scan at L=4 ---")
yc_runs = {}
for yc in (4.0, 8.0):
    d, B, ev, U, na, np_ = run(4.0, yc)
    yc_runs[yc] = dict(d=d, B=B, ev=ev, U=U, ndim=B.shape[1])
    log(f"ycut={yc:g} ndim={B.shape[1]} #arch={na} #prime={np_} "
        f"lam_min-1={ev[0]-1: .2e}  {counts_str(ev)}")
ref = scan[4.0]
log(f"reference (L=4, ycut=6): ndim={ref['ndim']}  {counts_str(ref['ev'])}")
log("bottom 20 reduced eigenvalues (L=4, ycut=6), ev-1:")
log("  " + " ".join(f"{x-1:.1e}" for x in ref["ev"][:20]))
# prime jump n-values actually present per run (exposes the overlap filter)
for name, rr in (("L=4 yc=6", ref), ("L=4 yc=4", yc_runs[4.0]), ("L=4 yc=8", yc_runs[8.0])):
    ns = sorted(int(n) for k, y, n, M in rr["d"]["per_jump"] if k == "prime")
    log(f"prime jumps present [{name}]: n = {ns}")

# ======================================================================
log("\n================ TEST 2a: counting laws ================")
Ls = [1.0, 1.5, 2.0, 3.0, 4.0]
bL = np.array([len(bundle_idx(scan[L]["ev"], 1e-6)) for L in Ls], float)
nd = np.array([scan[L]["ndim"] for L in Ls], float)
LsA = np.array(Ls)
log(f"{'L':>4} {'ndim':>5} {'bundle':>6} {'implied rank ndim-bundle':>24}")
for i, L in enumerate(Ls):
    log(f"{L:>4.1f} {int(nd[i]):>5d} {int(bL[i]):>6d} {int(nd[i]-bL[i]):>24d}")
def fit_law(Tfun, label, complement=False):
    cs = np.logspace(-2, 2.5, 900)
    best = None
    for c in cs:
        pred = np.array([NofT(Tfun(c, L)) for L in Ls], float)
        if complement:
            pred = nd - pred
        err = float(((pred - bL) ** 2).sum())
        if best is None or err < best[0]:
            best = (err, c, pred)
    err, c, pred = best
    log(f"{label}: best c={c:.3f}  pred={list(pred.astype(int))} "
        f"obs={list(bL.astype(int))}  SSE={err:.1f}")
    return c
log("-- direct: bundle(L) vs N(T(L)) --")
fit_law(lambda c, L: c * np.exp(L), "T = c*e^L      ")
fit_law(lambda c, L: c * L * np.exp(L), "T = c*L*e^L    ")
log("-- complement: bundle(L) vs ndim - N(T(L)) --")
fit_law(lambda c, L: c * np.exp(L), "T = c*e^L      ", complement=True)
fit_law(lambda c, L: c * L * np.exp(L), "T = c*L*e^L    ", complement=True)
fit_law(lambda c, L: c * L, "T = c*L        ", complement=True)
# implied rank vs L: linear (phase-space with fixed bandlimit) fit
r = nd - bL
A = np.vstack([LsA, np.ones_like(LsA)]).T
slope, icept = np.linalg.lstsq(A, r, rcond=None)[0]
pred = slope * LsA + icept
log(f"implied rank r(L)=ndim-bundle: linear fit r = {slope:.2f}*L + {icept:.2f}; "
    f"pred={np.round(pred,1)} obs={r.astype(int)}; "
    f"phase-space bandlimit T_eff = slope*pi/2 = {slope*np.pi/2:.1f}")

# ======================================================================
log("\n================ TEST 2b: sampling-map mechanism (L=2, L=4) ================")
NZ = 60
zeros_c = 0.5 + 1j * GAM[:NZ]
for L in (2.0, 4.0):
    rr = scan[L]; d, B, ev, U = rr["d"], rr["B"], rr["ev"], rr["U"]
    Fa, _ = F_matrix(d["PhiF"], FUNCS, d["u"], H, zeros_c)     # (243, NZ)
    Ms = Fa.T @ B                                              # (NZ, ndim) samples
    sv = np.linalg.svd(Ms, compute_uv=False)
    log(f"L={L}: sampling map {Ms.shape}; singular values log10:")
    log("  " + " ".join(f"{np.log10(max(s,1e-300)):.1f}" for s in sv))
    for rt in (1e-2, 1e-3, 1e-4, 1e-6):
        rk = int((sv > rt * sv[0]).sum())
        log(f"   rank(sv > {rt:g}*sv0) = {rk:2d}  -> ndim-rank = {B.shape[1]-rk}")
    ib = bundle_idx(ev, 1e-6)
    mask = np.ones(B.shape[1], bool); mask[ib] = False
    sb = Ms @ U[:, ib]; sc = Ms @ U[:, mask]
    log(f"   bundle dirs ({len(ib)}): max|sample|={np.abs(sb).max():.2e} "
        f"rms={np.sqrt((np.abs(sb)**2).mean()):.2e}")
    log(f"   complement ({int(mask.sum())}): max|sample|={np.abs(sc).max():.2e} "
        f"rms={np.sqrt((np.abs(sc)**2).mean()):.2e}")
    # sanity: E - Wm vs zero-side Gram with NZ zeros
    Z = 2 * ((Fa @ Fa.conj().T).real)
    Pside = B.T @ (d["E"] - d["Wm"]) @ B; Zside = B.T @ Z @ B
    rel = np.linalg.norm(Pside - Zside) / np.linalg.norm(Zside)
    log(f"   ||(E-Wm) - Z_{NZ}||/||Z_{NZ}|| on neutral space = {rel:.2e}")

# ======================================================================
log("\n================ TEST 2c: Fourier-zero signature (L=4) ================")
def orth_funcs(d, B, U, idx):
    """L2-orthonormal bundle functions on the grid."""
    V = d["G"].T @ (B @ U[:, idx])                    # (nu, db)
    Q, _ = np.linalg.qr(V * np.sqrt(H))
    return Q / np.sqrt(H)                              # columns: int f_i f_j = delta

def ft_profiles(d, F, gam):
    """F: (nu, k) functions on grid -> |int Phi f e^{i gamma u} du| and plain."""
    u = d["u"]; Pu = d["PhiF"](u)
    E = np.exp(1j * np.outer(gam, u))
    Fw = (E * (Pu * H)[None, :]) @ F
    Fp = (E * H) @ F
    return Fw, Fp

d = ref["d"]; B, ev, U = ref["B"], ref["ev"], ref["U"]
ib = bundle_idx(ev, 1e-6)
mask = np.ones(B.shape[1], bool); mask[ib] = False
Fb_orth = orth_funcs(d, B, U, ib)
Fc_orth = orth_funcs(d, B, U, np.where(mask)[0])
gam = np.arange(0.0, 90.001, 0.02)
Fw_b, _ = ft_profiles(d, Fb_orth, gam)
Fw_c, _ = ft_profiles(d, Fc_orth, gam)
S_b = (np.abs(Fw_b) ** 2).sum(1)
S_c = (np.abs(Fw_c) ** 2).sum(1)
def at(garr, S, g0):
    return float(S[np.argmin(np.abs(garr - g0))])
med_b = np.median(S_b); med_c = np.median(S_c)
log(f"bundle S(gamma): median over [0,90] = {med_b:.3e}; complement median = {med_c:.3e}")
log(f"{'k':>3} {'gamma_k':>8} {'S_b(gamma_k)/med':>17} {'S_c(gamma_k)/med':>17} "
    f"{'S_b(gamma_k+3.7)/med':>21}")
for k in range(15):
    g0 = GAM[k]
    log(f"{k+1:>3} {g0:>8.2f} {at(gam,S_b,g0)/med_b:>17.2e} {at(gam,S_c,g0)/med_c:>17.2e} "
        f"{at(gam,S_b,g0+3.7)/med_b:>21.2e}")
# individual bundle vectors: FT minima locations (Phi-weighted and plain)
Vb = d["G"].T @ (B @ U[:, ib])
for j in (0, 1):
    v = Vb[:, j]
    Fw1, Fp1 = ft_profiles(d, v[:, None], gam)
    for name, F1 in (("Phi-weighted", Fw1), ("plain", Fp1)):
        a = np.abs(F1[:, 0])
        mins = [i for i in range(5, len(a) - 5)
                if a[i] == a[i-5:i+6].min()]
        locs = [float(gam[i]) for i in mins[:8]]
        log(f"bundle vec {j} ({name} FT): first |F| minima at gamma = "
            f"{[round(x,2) for x in locs]}")
log(f"   zeta ordinates for comparison: {[round(float(x),2) for x in GAM[:8]]}")

# ======================================================================
log("\n================ TEST 3: per-jump bundle action (L=4, ycut=6) ================")
E = d["E"]
Cb = B @ U[:, ib]                                     # (243, db) coeff bundle
P = Cb @ np.linalg.pinv(Cb)                           # coeff-space projector
Bt = B.T
trE_bundle = float(np.trace(Cb.T @ E @ Cb))
rows = []
for kind, y, n, M in d["per_jump"]:
    fro = np.linalg.norm(P @ M @ P) / max(np.linalg.norm(M), 1e-300)
    MB = Bt @ M @ B
    fro_n = np.linalg.norm((U[:, ib].T @ MB @ U[:, ib])) / max(np.linalg.norm(MB), 1e-300)
    share = float(np.trace(Cb.T @ M @ Cb)) / trE_bundle
    rows.append((kind, float(y), n, fro, fro_n, share))
arch = sorted([r for r in rows if r[0] == "arch"], key=lambda r: r[1])
prime = sorted([r for r in rows if r[0] == "prime"], key=lambda r: r[2])
log("arch jumps, binned by y (mean Frobenius ratio full / neutral-restricted):")
bins = [(0, 0.5), (0.5, 1.0), (1.0, 1.5), (1.5, 2.0), (2.0, 3.0), (3.0, 6.0), (6.0, 13.0)]
for lo_, hi_ in bins:
    sel = [r for r in arch if lo_ <= r[1] < hi_]
    if sel:
        log(f"  y in [{lo_:>4},{hi_:>4}): n_jumps={len(sel):3d}  "
            f"fro={np.mean([r[3] for r in sel]):.3e}  fro_neutral={np.mean([r[4] for r in sel]):.3e}  "
            f"share_sum={sum(r[5] for r in sel): .3e}")
log("prime jumps by n:")
log(f"{'n':>5} {'y=log n':>8} {'fro_full':>10} {'fro_neutral':>12} {'bundle trace share':>19}")
for kind, y, n, fro, fro_n, share in prime:
    log(f"{int(n):>5} {y:>8.3f} {fro:>10.3e} {fro_n:>12.3e} {share:>19.3e}")
tot_share_prime = sum(r[5] for r in prime); tot_share_arch = sum(r[5] for r in arch)
log(f"trace shares: arch total {tot_share_arch:.4f}, prime total {tot_share_prime:.4f}")

# ======================================================================
log("\n================ TEST 4: perturbation split (L=4) ================")
def bundle_funcs_on_grid(dd, BB, UU, idx):
    V = dd["G"].T @ (BB @ UU[:, idx])
    Q, _ = np.linalg.qr(V * np.sqrt(H))
    return Q  # columns orthonormal in plain Euclidean on grid (weight folded in)

perts = {}
perts["reference ycut=6"] = (ref["d"], ref["B"], ref["U"], ref["ev"])
perts["ycut=4"] = (yc_runs[4.0]["d"], yc_runs[4.0]["B"], yc_runs[4.0]["U"], yc_runs[4.0]["ev"])
perts["ycut=8"] = (yc_runs[8.0]["d"], yc_runs[8.0]["B"], yc_runs[8.0]["U"], yc_runs[8.0]["ev"])
# cosh-only constraint on the reference setup (drop sinh quotient)
Bc, evc, Uc = neutral_basis(ref["d"]["E"], ref["d"]["Wm"], ref["d"]["Cm"][:1])
perts["cosh-only constraint"] = (ref["d"], Bc, Uc, evc)
# alt basis
da, Ba, eva, Ua, naa, npa = run(4.0, 6.0, funcs=FUNCS_ALT)
perts["alt basis polys(0,1)"] = (da, Ba, Ua, eva)
Q = {}
for name, (dd, BB, UU, ee) in perts.items():
    idx = bundle_idx(ee, 1e-6)
    Q[name] = bundle_funcs_on_grid(dd, BB, UU, idx)
    log(f"{name:>24}: ndim={BB.shape[1]:3d}  bundle(1e-6)={len(idx):2d}  "
        f"bundle(1e-5)={len(bundle_idx(ee,1e-5)):2d}")
Qref = Q["reference ycut=6"]
log("principal cosines of each perturbed bundle vs reference bundle:")
for name, Qp in Q.items():
    if name == "reference ycut=6":
        continue
    sv = np.linalg.svd(Qref.T @ Qp, compute_uv=False)
    sv = np.sort(sv)[::-1]
    log(f"  {name:>24}: cosines={np.round(sv,3)}  #(>0.9)={int((sv>0.9).sum())} "
        f"#(>0.7)={int((sv>0.7).sum())}")
# average projector core
Pavg = sum(Qp @ Qp.T for Qp in Q.values()) / len(Q)
pw = np.linalg.eigvalsh(Pavg)[::-1]
log(f"average-projector eigenvalues (top 20): {np.round(pw[:20],3)}")
log(f"core dim: #(eig>0.8)={int((pw>0.8).sum())}  #(eig>0.5)={int((pw>0.5).sum())} "
    f"(bundles have dims {[Q[n].shape[1] for n in Q]})")
# per-reference-vector robustness: min over perturbations of projection fraction
log("reference bundle vectors: min-over-perturbations projection fraction:")
fr = []
for i in range(Qref.shape[1]):
    q = Qref[:, i]
    fr.append(min(float(np.sum((Qp.T @ q) ** 2)) for Qp in Q.values()
                  if Qp is not Qref))
fr = np.array(fr)
log(f"  fractions={np.round(np.sort(fr)[::-1],3)}")
log(f"  robust (>=0.8 in ALL perturbations): {int((fr>=0.8).sum())} of {len(fr)}")
log(f"\n[total elapsed {time.time()-t0:.0f}s] done")
