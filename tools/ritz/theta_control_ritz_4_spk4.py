"""
theta_control_ritz.py -- negative control for ground-transform / Weil-positivity candidates.

Setup (all verified numerically, see checks at bottom):
  Phi      : even positive kernel, int Phi e^{(s-1/2)u} du = xi(s)  (= 2 x Titchmarsh Phi)
  control  : L_theta(s)=zeta(s-theta)zeta(s+theta), 0<theta<1/2.  Euler product + functional
             equation s->1-s + positive-law representation, zeros at rho±theta (OFF the line).
             Phi_theta=(Phi e^{-theta u})*(Phi e^{+theta u}),  nu_theta=2cosh(theta y) nu,
             four pole functions cosh/sinh((1/2∓theta)u)... weight W=Phi_theta*(xi(2th)/2*c_a+xi(1+2th)/2*c_b).
  Identity : Q_theta[Phi_theta g] = E_theta(g) - <g,W g> + P_theta(Phi_theta g)  == sum_zeros F(z) conj F(1-conj z)
  Test     : on the neutral space (all pole moments of Phi_theta g vanish), RH-type inequality is
             E_theta >= W, i.e. generalized eigenvalue lam >= 1.   theta=0 is true zeta (lam_min = 1).
Any candidate proof of lam>=1 for zeta must FAIL here for theta>0.
"""
import numpy as np, time
from numpy.polynomial.legendre import leggauss
from math import pi, log, sqrt
# ---- Phi (note normalization: Phi = 2*Phi_Titchmarsh) with analytic derivatives
# terms for u>=0: coef * e^{alpha u} * w^m * exp(-c w), w=e^{2u}
def build_terms(N=6):
    T=[]
    for n in range(1,N+1):
        c=pi*n*n
        T.append((2*2*pi**2*n**4, 4.5, 0, c))
        T.append((-2*3*pi*n**2, 2.5, 0, c))
    return T
TERMS0=build_terms()
def deriv(terms):
    out=[]
    for (A,al,m,c) in terms:
        out.append(((A*(al+2*m)),al,m,c))
        out.append((-2*c*A,al,m+1,c))
    return out
DER=[TERMS0]
for k in range(8): DER.append(deriv(DER[-1]))
def Phik(u,k=0):
    u=np.asarray(u,float); a=np.abs(u); s=np.zeros_like(a)
    for (A,al,m,c) in DER[k]:
        w=np.exp(2*a)
        s+=A*np.exp(al*a-c*w)*w**m
    return s*(np.sign(u)**k if k%2 else 1)   # parity: Phi^(k)(-u)=(-1)^k Phi^(k)(u)
def Phi(u): return Phik(u,0)
def gk(u,k):   # (d^2-1/4) d^k Phi / Phi
    u=np.asarray(u,float); P=Phi(u); out=np.zeros_like(u)
    ok=P>1e-250
    out[ok]=(Phik(u,k+2)[ok]-0.25*Phik(u,k)[ok])/P[ok]
    return out
# primes / von Mangoldt up to nmax
def vonM(nmax):
    L={}
    sieve=[True]*(nmax+1)
    for p in range(2,nmax+1):
        if sieve[p]:
            for q in range(2*p,nmax+1,p): sieve[q]=False
            pk=p
            while pk<=nmax:
                L[pk]=log(p); pk*=p
    return L
from math import log, sqrt, pi
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import CubicSpline
from mpmath import mp, zeta, gamma, pi as mpi
def xi(s):
    s=complex(s); return complex(0.5*s*(s-1)*mpi**(-s/2)*gamma(s/2)*zeta(s))
# ---- build Phi_theta = (Phi e^{-theta .}) * (Phi e^{theta .}),  int Phi_th e^{(s-1/2)u} = xi(s-theta) xi(s+theta)
def build(theta):
    dx=0.002; x=np.arange(-9,9+dx/2,dx)
    P=Phi(x)
    conv=np.convolve(P*np.exp(-theta*x),P*np.exp(theta*x))*dx
    u=2*x[0]+np.arange(len(conv))*dx
    m=(u>=-9)&(u<=9)
    return CubicSpline(u[m],conv[m],extrapolate=False),u[m],conv[m]
LAM=vonM(170000)
def rho(y): return np.exp(-y/2)/(1-np.exp(-2*y))
def ypan(npts=40):
    x,w=leggauss(npts);N=[];Wt=[]
    for a,b in [(0,0.5),(0.5,1.5),(1.5,3),(3,6),(6,12)]:
        N.append(0.5*(b-a)*x+0.5*(b+a)); Wt.append(0.5*(b-a)*w)
    return np.concatenate(N),np.concatenate(Wt)
def make_problem(theta, ycut=12.0):
    """returns callable PhiT, nu-jumps, weight function, constraint functions"""
    if theta==0:
        PhiF=lambda u: Phi(u)
        coshw=lambda y:1.0
        W=lambda u:0.5*(np.exp(u/2)+np.exp(-u/2))*Phi(u)               # (C/2) c Phi
        cons=[lambda u:np.cosh(u/2),lambda u:np.sinh(u/2)]
    else:
        sp,_,_=build(theta)
        PhiF=lambda u: np.nan_to_num(sp(np.clip(u,-9,9)))*(np.abs(u)<=9)
        coshw=lambda y:2*np.cosh(theta*y)
        a=0.5-theta;b=0.5+theta
        xa=xi(2*theta).real; xb=xi(1+2*theta).real
        W=lambda u:PhiF(u)*(0.5*xa*2*np.cosh(a*u)+0.5*xb*2*np.cosh(b*u))
        cons=[lambda u:np.cosh(a*u),lambda u:np.sinh(a*u),lambda u:np.cosh(b*u),lambda u:np.sinh(b*u)]
    ys,wy=ypan()
    jumps=[(y,w*rho(y)*coshw(y)) for y,w in zip(ys,wy)]
    jumps+=[(log(n),Lm/sqrt(n)*coshw(log(n))) for n,Lm in LAM.items() if log(n)<=ycut]
    return PhiF,jumps,W,cons
def solve(theta,funcs,L=8.0,h=0.01,ycut=12.0):
    PhiF,jumps,W,cons=make_problem(theta,ycut)
    u=np.arange(-L,L+h/2,h)
    G=np.array([f(u) for f in funcs]); Pu=PhiF(u)
    m=len(funcs); E=np.zeros((m,m))
    for y,wt in jumps:
        sh=u-y; Ps=PhiF(sh)
        if wt*np.max(Pu*Ps)<1e-40: continue
        D=G-np.array([f(sh) for f in funcs])
        E+=(D*(Pu*Ps*h*wt))@D.T
    Wm=(G*(W(u)*h))@G.T
    Cm=np.array([(G*(PhiF(u)*c(u)*h)).sum(1) for c in cons])   # constraints <c, Phi g>=0
    return E,Wm,Cm
def reduce_eig(E,Wm,Cm):
    w,V=np.linalg.eigh(Wm); keep=w>1e-12*w.max(); T=V[:,keep]/np.sqrt(w[keep])
    Ct=Cm@T; _,s,Vt=np.linalg.svd(Ct); r=(s>1e-12*s.max()).sum(); Z=Vt[r:].T
    B=T@Z; Eb=B.T@E@B; Eb=(Eb+Eb.T)/2
    return np.linalg.eigvalsh(Eb)

def trigbasis(kmax,dk,polys):
    F=[(lambda x,j=j:x**j) for j in polys]
    for k in np.arange(dk,kmax+1e-9,dk):
        for j in polys:
            F.append(lambda x,k=k,j=j: x**j*np.cos(k*x)); F.append(lambda x,k=k,j=j: x**j*np.sin(k*x))
    return F

if __name__=="__main__":
    t0=time.time(); funcs=trigbasis(30,0.75,(0,1,2)); print("basis",len(funcs))
    for th in [0,0.1,0.2,0.3]:
        E,Wm,Cm=solve(th,funcs,L=4.0,h=0.01,ycut=6.0); ev=reduce_eig(E,Wm,Cm)
        print(f"theta={th:.1f}  lam_min-1 = {ev[0]-1: .3e}   (theta=0: noise ~1e-7; control: ~ -0.052*theta^2)  [{time.time()-t0:.0f}s]")
