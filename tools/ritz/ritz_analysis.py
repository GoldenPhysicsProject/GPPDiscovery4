"""
ritz_analysis.py -- three numerical tests for the fixed-window / theta-control program.

Test 1: prime-side matrix vs zero-side matrix agreement (normalization check).
Test 2: lowest Ritz eigenvector vs theta; per-prime loss-of-positivity profile.
Test 3: degree-two residual vs the negative direction (correlation test).

Imports the discretization machinery from theta_control_ritz_4_spk4.py so the
setup is IDENTICAL to the baseline script.
"""
import numpy as np, time, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # staged: import sibling scripts locally
from math import log, sqrt, pi, cosh
import theta_control_ritz_4_spk4 as R

# ---------------------------------------------------------------- labeled setup
def make_labeled(theta, ycut=12.0):
    """Same as R.make_problem but with labeled jumps."""
    if theta == 0:
        PhiF = lambda u: R.Phi(u)
        coshw = lambda y: 1.0
        W = lambda u: 0.5*(np.exp(u/2)+np.exp(-u/2))*R.Phi(u)
        cons = [lambda u: np.cosh(u/2), lambda u: np.sinh(u/2)]
    else:
        sp, _, _ = R.build(theta)
        PhiF = lambda u: np.nan_to_num(sp(np.clip(u, -9, 9)))*(np.abs(u) <= 9)
        coshw = lambda y: 2*np.cosh(theta*y)
        a = 0.5-theta; b = 0.5+theta
        xa = R.xi(2*theta).real; xb = R.xi(1+2*theta).real
        W = lambda u: PhiF(u)*(0.5*xa*2*np.cosh(a*u)+0.5*xb*2*np.cosh(b*u))
        cons = [lambda u: np.cosh(a*u), lambda u: np.sinh(a*u),
                lambda u: np.cosh(b*u), lambda u: np.sinh(b*u)]
    ys, wy = R.ypan()
    jumps = [('arch', y, w*R.rho(y)*coshw(y), None) for y, w in zip(ys, wy)]
    jumps += [('prime', log(n), Lm/sqrt(n)*coshw(log(n)), n)
              for n, Lm in R.LAM.items() if log(n) <= ycut]
    return PhiF, jumps, W, cons

def setup(theta, funcs, L=8.0, h=0.01, ycut=12.0):
    PhiF, jumps, W, cons = make_labeled(theta, ycut)
    u = np.arange(-L, L+h/2, h)
    G = np.array([f(u) for f in funcs]); Pu = PhiF(u)
    m = len(funcs)
    E = np.zeros((m, m))
    per_jump = []          # (kind, y, n, M_jump)
    for kind, y, wt, n in jumps:
        sh = u - y; Ps = PhiF(sh)
        if wt*np.max(Pu*Ps) < 1e-40:
            continue
        D = G - np.array([f(sh) for f in funcs])
        M = (D*(Pu*Ps*h*wt)) @ D.T
        E += M
        per_jump.append((kind, y, n, M))
    Wm = (G*(W(u)*h)) @ G.T
    Cm = np.array([(G*(PhiF(u)*c(u)*h)).sum(1) for c in cons])
    return dict(u=u, G=G, Pu=Pu, PhiF=PhiF, E=E, Wm=Wm, Cm=Cm,
                per_jump=per_jump, theta=theta, funcs=funcs, L=L, h=h)

def neutral_basis(E, Wm, Cm):
    """Replicates reduce_eig but also returns the neutral-space basis B and eigenpairs."""
    w, V = np.linalg.eigh(Wm); keep = w > 1e-12*w.max(); T = V[:, keep]/np.sqrt(w[keep])
    Ct = Cm @ T; _, s, Vt = np.linalg.svd(Ct); r = (s > 1e-12*s.max()).sum(); Z = Vt[r:].T
    B = T @ Z
    Eb = B.T @ E @ B; Eb = (Eb+Eb.T)/2
    ev, U = np.linalg.eigh(Eb)
    return B, ev, U

# ---------------------------------------------------------------- zeros
_zcache = {}
def zeta_zeros(N):
    if N not in _zcache:
        from mpmath import mp, zetazero
        mp.dps = 30
        _zcache[N] = np.array([complex(zetazero(k)) for k in range(1, N+1)])
    return _zcache[N]

def F_matrix(PhiF, funcs, u, h, zeros_a, zeros_b=None):
    """F[j,k] = int Phi(u) g_j(u) e^{(z_k-1/2) u} du.  zeros_b=None -> single set."""
    G = np.array([f(u) for f in funcs]); Pu = PhiF(u)
    MG = G * (Pu*h)                       # (nbasis, nu)
    def block(zz):
        Ea = np.exp(np.outer(u, zz-0.5))  # (nu, nz)
        return (MG @ Ea)                  # (nbasis, nz)
    Fa = block(zeros_a)
    Fb = block(zeros_b) if zeros_b is not None else None
    return Fa, Fb

def zero_side_matrix(Fa, Fb=None):
    """Z[i,j] = sum_k F_i(z_k^+) conj(F_j(1-conj(z_k^+))) + (paired term).

    Fa covers z_k^+ for k=1..N (upper half-plane).  The Weil sum runs over ALL
    zeros, so conjugates must be included: for real test functions this doubles
    the upper half-plane sum.  With Fb given (the z_k^- = 1-conj(z_k^+) set),
    Z = sum_k [F_i(z^+)conj(F_j(z^-)) + F_i(z^-)conj(F_j(z^+))] + conjugates.
    """
    if Fb is None:                        # on-line zeros: 1-conj(z)=z
        return 2*((Fa @ Fa.conj().T).real)
    # paired zeros; conjugate zeros contribute conj(...) -> 2*Re (real symmetric)
    return 2*((Fa @ Fb.conj().T + Fb @ Fa.conj().T).real)

# ---------------------------------------------------------------- TEST 1
def test1(theta, Nzero=300, L=4.0, h=0.01, ycut=6.0):
    t0 = time.time()
    funcs = R.trigbasis(30, 0.75, (0, 1, 2))
    S = setup(theta, funcs, L=L, h=h, ycut=ycut)
    B, ev, U = neutral_basis(S['E'], S['Wm'], S['Cm'])
    zeros = zeta_zeros(Nzero)
    if theta == 0:
        Fa, _ = F_matrix(S['PhiF'], funcs, S['u'], h, zeros)
        Z = zero_side_matrix(Fa)
    else:
        zp = zeros + theta; zm = zeros - theta
        Fa, Fb = F_matrix(S['PhiF'], funcs, S['u'], h, zp, zm)
        Z = zero_side_matrix(Fa, Fb)
    prime_side = B.T @ (S['E'] - S['Wm']) @ B
    zero_side = B.T @ Z @ B
    num = np.linalg.norm(prime_side - zero_side)
    den = np.linalg.norm(zero_side)
    # convergence check: half the zeros
    return dict(theta=theta, rel_err=num/den,
                lam_min=ev[0], dim=B.shape[1], time=time.time()-t0)

# ---------------------------------------------------------------- TEST 2
def loss_profile(theta, L=4.0, h=0.01, ycut=6.0, kmax=30):
    """Per-jump contributions of the lowest Ritz direction; returns dicts."""
    funcs = R.trigbasis(kmax, 0.75, (0, 1, 2))
    S = setup(theta, funcs, L=L, h=h, ycut=ycut)
    B, ev, U = neutral_basis(S['E'], S['Wm'], S['Cm'])
    c = B @ U[:, 0]                       # full-space coefficient vector, lowest dir
    out = []
    for kind, y, n, M in S['per_jump']:
        out.append((kind, y, n, float(c @ M @ c)))
    return dict(theta=theta, lam_min=float(ev[0]), per_jump=out,
                Wexp=float(c @ S['Wm'] @ c), Eexp=float(c @ S['E'] @ c))

# ---------------------------------------------------------------- TEST 3a: exact residual algebra
def test3a(primes=(2, 3, 5, 7, 11, 13), thetas=(0.1, 0.2, 0.3)):
    rng = np.random.default_rng(0)
    rep = []
    for p in primes:
        for th in thetas:
            r_plus = p**(-0.5+th); r_minus = p**(-0.5-th); q = p**-1
            for _ in range(5):
                u = rng.uniform(0.01, 0.9)
                two = (1-r_plus*u)*(1-r_minus*u)
                comp = (1-q*r_plus*u)*(1-q*r_minus*u) - two
                pred = (1-q)*(r_plus+r_minus)*u - (1-q*q)*r_plus*r_minus*u*u
                # one-radius must have NO u^2 term
                one = (1-q*(p**-0.5)*u) - (1-(p**-0.5)*u)
                one_pred = (1-q)*(p**-0.5)*u
                rep.append((abs(comp-pred), abs(one-one_pred),
                            abs((1-q*q)*r_plus*r_minus - (1-p**-2)*p**-1)))
    rep = np.array(rep)
    return dict(max_resid_err=float(rep[:, 0].max()),
                max_one_err=float(rep[:, 1].max()),
                max_coef_err=float(rep[:, 2].max()))

if __name__ == "__main__":
    print("== TEST 3a: exact residual algebra ==")
    r = test3a(); print(r)
    print("== TEST 1: theta=0 ==")
    r = test1(0.0); print({k: (f"{v:.3e}" if isinstance(v, float) else v) for k, v in r.items()})
    print("== TEST 1: theta=0.2 ==")
    r = test1(0.2); print({k: (f"{v:.3e}" if isinstance(v, float) else v) for k, v in r.items()})
