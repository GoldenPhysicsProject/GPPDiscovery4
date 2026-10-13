#!/usr/bin/env python3
"""
bundle_counting.py -- test Claude Code's counting law for the near-null
bundle of the Ritz/Weil reduced form:  bundle_dim ~ D - 2*N(Omega),
D = neutral dim, N(Omega) = #{zeta zeros with ordinate <= Omega},
as the test-basis BANDWIDTH Omega varies at fixed window L=4, theta=0.

Machinery identical to bundle_mechanism.py / ground_vector.py:
ritz_analysis.setup + neutral_basis, h=0.01, ycut=6.
Effective Omega is MEASURED per rung from the zero-sampling map on the
neutral space (row norms of samples at zeta ordinates), not asserted
from the basis parameters.  NUMERICAL throughout.
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
H = 0.01
L = 4.0
YCUT = 6.0
TOLS = (1e-5, 1e-6, 1e-7)
NZ = 80          # zeros used for the sampling map / Omega measurement

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

GAM = zero_ordinates()
def NofT(T):
    return int((GAM <= T).sum())

# ------------------------------------------------- neutral basis, cut param
def neutral_basis_cut(E, Wm, Cm, cut=1e-12):
    """neutral_basis with the whitening/SVD cut exposed (default = rig's)."""
    w, V = np.linalg.eigh(Wm)
    keep = w > cut * w.max()
    T = V[:, keep] / np.sqrt(w[keep])
    Ct = Cm @ T
    _, s, Vt = np.linalg.svd(Ct)
    r = (s > cut * s.max()).sum()
    Z = Vt[r:].T
    B = T @ Z
    Eb = B.T @ E @ B
    Eb = (Eb + Eb.T) / 2
    ev, U = np.linalg.eigh(Eb)
    return B, ev, U

def bundle_count(ev, tol):
    return int((np.abs(ev - 1.0) < tol).sum())

# ------------------------------------------------------------- one rung
def run_rung(label, funcs, cut=1e-12):
    t0 = time.time()
    d = setup(0.0, funcs, L=L, h=H, ycut=YCUT)
    B, ev, U = neutral_basis_cut(d["E"], d["Wm"], d["Cm"], cut=cut)
    D = B.shape[1]
    # sampling map on the neutral space, Euclidean-orthonormalized columns
    Q, _ = np.linalg.qr(B)
    zeros_c = 0.5 + 1j * GAM[:NZ]
    Fa, _ = F_matrix(d["PhiF"], funcs, d["u"], H, zeros_c)   # (m, NZ)
    Ms = Fa.T @ Q                                            # (NZ, D) samples
    rown = np.linalg.norm(Ms, axis=1)                        # per-zero reach
    sv = np.linalg.svd(Ms, compute_uv=False)
    out = dict(label=label, m=len(funcs), D=D, ev=ev,
               lam_min=float(ev[0]), rown=rown, sv=sv,
               bundle={t: bundle_count(ev, t) for t in TOLS},
               elapsed=time.time() - t0)
    del d
    return out

def omega_meas(rown, eps):
    """Largest ordinate whose sampling row norm >= eps * max row norm."""
    thr = eps * rown.max()
    idx = np.where(rown >= thr)[0]
    return float(GAM[idx.max()]) if len(idx) else 0.0

def report_rung(r):
    line = (f"{r['label']:>22} m={r['m']:>3} D={r['D']:>3} "
            f"lam_min-1={r['lam_min']-1: .2e}  bundle: " +
            " ".join(f"{t:g}:{r['bundle'][t]}" for t in TOLS))
    log(line)
    for eps in (1e-2, 1e-3):
        Om = omega_meas(r["rown"], eps)
        N = NofT(Om)
        log(f"{'':>22}  Omega(eps={eps:g}) = {Om:6.2f}  N(Omega)={N:>2}   "
            f"D-2N = {r['D']-2*N:>3}   D-N = {r['D']-N:>3}   "
            f"bundle@1e-6 = {r['bundle'][1e-6]}")
    return r

# ======================================================================
log("================ LADDER A: kmax at fixed dk=0.75, polys (0,1,2) ================")
ladderA = []
for kmax in (16, 20, 24, 28, 30, 34, 38, 44, 50):
    funcs = trigbasis(kmax, 0.75, (0, 1, 2))
    r = run_rung(f"kmax={kmax} dk=0.75", funcs)
    ladderA.append(r)
    report_rung(r)
    log(f"{'':>22}  [{r['elapsed']:.0f}s]")

log("\n================ LADDER B: fixed basis size 243, dk varies (Omega knob only) ================")
ladderB = []
for kmax, dk in ((20, 0.5), (30, 0.75), (44, 1.1)):
    funcs = trigbasis(kmax, dk, (0, 1, 2))
    r = run_rung(f"kmax={kmax} dk={dk}", funcs)
    ladderB.append(r)
    report_rung(r)
    log(f"{'':>22}  [{r['elapsed']:.0f}s]")

log("\n================ SENSITIVITY: whitening cut at kmax=30, dk=0.75 ================")
funcs30 = trigbasis(30, 0.75, (0, 1, 2))
for cut in (1e-10, 1e-12, 1e-14):
    r = run_rung(f"cut={cut:g}", funcs30, cut=cut)
    log(f"cut={cut:g}: D={r['D']}  bundle: " +
        " ".join(f"{t:g}:{r['bundle'][t]}" for t in TOLS) +
        f"  Omega(1e-3)={omega_meas(r['rown'],1e-3):.2f}  [{r['elapsed']:.0f}s]")

# ======================================================================
log("\n================ CROSSING ANALYSIS (ladder A, sorted by Omega) ================")
for eps in (1e-3,):
    rows = []
    for r in ladderA:
        Om = omega_meas(r["rown"], eps)
        rows.append((Om, r))
    rows.sort(key=lambda x: x[0])
    log(f"eps={eps:g}:  rung-to-rung changes as Omega grows")
    log(f"{'Omega':>7} {'D':>4} {'bundle@1e-5':>11} {'bundle@1e-6':>11} "
        f"{'D-2N':>5} {'D-N':>5} {'crossed gammas':>30}")
    prev = None
    for Om, r in rows:
        N = NofT(Om)
        crossed = ""
        if prev is not None:
            g0, g1 = prev[0], Om
            cg = [float(g) for g in GAM if g0 < g <= g1]
            db6 = r["bundle"][1e-6] - prev[1]["bundle"][1e-6]
            crossed = (f"crossed {[round(x,1) for x in cg]}: "
                       f"dbundle@1e-6={db6:+d} (pred {-2*len(cg):+d} for -2N, "
                       f"{-len(cg):+d} for -N)")
        log(f"{Om:>7.2f} {r['D']:>4} {r['bundle'][1e-5]:>11} {r['bundle'][1e-6]:>11} "
            f"{r['D']-2*N:>5} {r['D']-N:>5}  {crossed}")
        prev = (Om, r)

log("\n================ SAMPLING SINGULAR VALUES (ladder A) ================")
for r in ladderA:
    sv = r["sv"]
    log(f"{r['label']:>22}: log10 sv (first 24): " +
        " ".join(f"{np.log10(max(s,1e-300)):.1f}" for s in sv[:24]))

np.savez(f"{OUT}/bundle_counting_data.npz",
         **{f"A{i}_ev": r["ev"] for i, r in enumerate(ladderA)},
         **{f"A{i}_rown": r["rown"] for i, r in enumerate(ladderA)},
         **{f"B{i}_ev": r["ev"] for i, r in enumerate(ladderB)},
         **{f"B{i}_rown": r["rown"] for i, r in enumerate(ladderB)},
         gam=GAM)
log(f"\n[done] data saved to {OUT}/bundle_counting_data.npz")
