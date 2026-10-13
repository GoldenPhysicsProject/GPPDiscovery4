import math,json
import numpy as np
from scipy.integrate import quad
from check_euler_source_2026_09_30 import prime_powers
pi=math.pi

def phi_score(u):
    a=abs(u);v=math.exp(a);p=pi*v*v
    if p>740:
        return 0., math.copysign(1.,u)*(4.5-2*p)
    n=np.arange(1.,6.); t=pi*(n*v)**2
    e=np.exp(-t)
    val=math.exp(a/2)*np.sum((4*t*t-6*t)*e)
    der=math.exp(a/2)*np.sum((-8*t**3+30*t*t-15*t)*e)
    return val,(der/val if u>=0 else -der/val)

def mu(u):
    p,_=phi_score(u)
    return 2*math.cosh(u/2)*p

def integral(f): return quad(f,-4,0,epsabs=1e-11,limit=200)[0]+quad(f,0,4,epsabs=1e-11,limit=200)[0]
def jump(y):
    return integral(lambda u:phi_score(u)[0]*phi_score(u-y)[0]*(phi_score(u)[1]-phi_score(u-y)[1])**2)
norm=integral(mu)
var=integral(lambda u:mu(u)*phi_score(u)[1]**2)
rg=integral(lambda u:mu(u)*math.tanh(u/2)*phi_score(u)[1])
arch=quad(lambda y: math.exp(-y/2)/(-math.expm1(-2*y))*jump(y),0,4,epsabs=2e-10,limit=100)[0]
prime=sum(l/math.sqrt(n)*jump(math.log(n)) for n,l in prime_powers(54))
target=.5*(var+rg*rg)
print(json.dumps(dict(mu_mass=norm,score_variance=var,boundary_score_mean=rg,arch_jump=arch,prime_jump=prime,jump_energy=arch+prime,variance_target=target,difference=arch+prime-target,cutoffs=dict(u=4,y=4,max_prime_power=54),status='finite numerical consistency check; no interval error enclosure'),indent=2))
