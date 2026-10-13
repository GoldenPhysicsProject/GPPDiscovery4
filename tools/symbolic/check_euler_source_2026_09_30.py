import json, math
import numpy as np
from scipy.integrate import quad
from scipy.special import digamma, bernoulli

def cq(f,a,b):
    return quad(lambda x: complex(f(x)).real,a,b,epsabs=2e-11,epsrel=2e-11,limit=300)[0]+1j*quad(lambda x: complex(f(x)).imag,a,b,epsabs=2e-11,epsrel=2e-11,limit=300)[0]

def zeta_pair(s,N=100,K=12):
    n=np.arange(1,N,dtype=float); l=np.log(n)
    val=np.sum(np.exp(-s*l)); der=-np.sum(l*np.exp(-s*l))
    ln=math.log(N)
    q=N**(1-s)/(s-1); val+=q
    der+=q*(-ln-1/(s-1))
    q=.5*N**(-s); val+=q; der-=ln*q
    bs=bernoulli(2*K)
    for k in range(1,K+1):
        a=np.array([s+j for j in range(2*k-1)])
        q=bs[2*k]/math.factorial(2*k)*np.prod(a)*N**(-s-2*k+1)
        val+=q;der+=q*(np.sum(1/a)-ln)
    return val,der

def prime_powers(X):
    sieve=np.ones(X+1,dtype=bool);sieve[:2]=False
    for p in range(2,math.isqrt(X)+1):
        if sieve[p]:sieve[p*p::p]=False
    result=[]
    for p in np.flatnonzero(sieve):
        q=int(p)
        while q<=X:
            result.append((q,math.log(p)))
            q*=int(p)
    return sorted(result)

def rho(y): return np.exp(-y/2)/(-np.expm1(-2*y))
def r(y): return np.exp(-y)/(-np.expm1(-2*y))

out=[]
for L in [2.,4.,6.,8.,10.]:
    ns,lam=np.array(prime_powers(math.floor(math.exp(L)))).T
    ln=np.log(ns);weights=lam/np.sqrt(ns)
    for w in [1j, .7+.8j, 2+1.3j]:
        s=.5-1j*w;E=np.exp(2j*L*w)
        D=lambda y:E/(1-E)*(np.exp(1j*w*y)-np.exp(-1j*w*y))
        pole=2*w/(w*w+.25)*(np.cosh(L/2)/np.sin(L*w)-1/np.tan(L*w))
        arch=cq(lambda y:rho(y)*np.sin((L-y)*w)/np.sin(L*w)-r(y),0,L)
        arch+=(math.log(4*math.pi)+np.euler_gamma)/2-.5*math.log(1/math.tanh(L/2))
        prime=np.sum(weights*np.sin((L-ln)*w)/np.sin(L*w))
        H=pole-arch-prime
        z,zp=zeta_pair(s)
        target=1/s+1/(s-1)-.5*math.log(math.pi)+.5*digamma(s/2)+zp/z
        # Independent forward-sum/tail decomposition; zeta evaluated by EM.
        pole_error=4j*w/(w*w+.25)*(E-np.cosh(L/2)*np.exp(1j*L*w))/(1-E)
        arch_error=-cq(lambda y:rho(y)*D(y),0,L)+cq(lambda y:rho(y)*np.exp(1j*w*y),L,np.inf)
        prime_error=-np.sum(weights*D(ln))+(-zp/z-np.sum(lam*np.exp(-s*ln)))
        err=pole_error+arch_error+prime_error
        out.append(dict(L=L,w=str(w),source_error=float(abs(H-target)),error_formula_residual=float(abs(H-target-err))))
print(json.dumps(out,indent=2))
