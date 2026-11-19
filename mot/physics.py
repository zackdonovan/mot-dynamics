"""Incoherent F=2 -> F'=3 population model; all quantities SI except detuning/Gamma."""
from dataclasses import dataclass
import numpy as np
HBAR=1.054571817e-34
MASS=1.44316060e-25
GAMMA=2*np.pi*6.065e6
K=2*np.pi/780.241209e-9
MUB_HBAR=9.2740100783e-24/HBAR
DIRECTIONS=np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1.]])
MG=np.arange(-2,3)
Q=np.arange(-1,2)
ME=MG[:,None]+Q
CG=np.stack([(3-MG)*(4-MG)/30,(9-MG**2)/15,(3+MG)*(4+MG)/30],axis=-1)
EINDEX=ME+3
@dataclass
class Controls:
    detuning: tuple=(-2.,-2.,-2.)
    purity: tuple=(1.,1.,1.)
    power: tuple=(1.,1.,1.)
    gradient: float=.12 # T/m radial; axial gradient -2*radial
    waist: float=.012
    geometry: str="axial6"
    def __post_init__(self):
        if self.geometry not in ('axial6','oblique8'): raise ValueError('unknown geometry')
        n=3 if self.geometry=='axial6' else 4
        if any(len(x)!=n for x in (self.detuning,self.purity,self.power)): raise ValueError('one control per beam pair required')
        if not all(np.isfinite(x).all() for x in (self.detuning,self.purity,self.power)): raise ValueError('finite controls required')
        if np.any(np.abs(self.purity)>1): raise ValueError('helicity mixture must lie in [-1,1]')
        if not np.isfinite([self.waist,self.gradient]).all() or np.any(np.array(self.power)<0) or self.waist<=0 or self.gradient<0: raise ValueError('invalid intensity, waist or gradient')

def beam_directions(c):
    if c.geometry=="axial6": return DIRECTIONS
    if c.geometry=="oblique8":
        d=np.array([1.,1.,1.])/np.sqrt(3)
        return np.vstack([DIRECTIONS,d,-d])
    raise ValueError("unknown geometry")

def polarization_weights(cosine,helicity):
    """Incoherent mixture of opposite circular helicities. q order -1,0,+1."""
    c=np.asarray(cosine); h=np.asarray(helicity)
    return np.stack(((1+c*c-2*h*c)/4,(1-c*c)/2,(1+c*c+2*h*c)/4),axis=-1)

def rates(r,v,c):
    r=np.atleast_2d(r); v=np.atleast_2d(v)
    b=r*np.array([c.gradient,c.gradient,-2*c.gradient]); bm=np.linalg.norm(b,axis=1)
    axis=np.divide(b,bm[:,None],out=np.tile([0.,0.,1.],(len(r),1)),where=bm[:,None]>1e-15)
    directions=beam_directions(c)
    cosine=axis@directions.T
    signs=np.array([1,1,-1] if c.geometry=="axial6" else [1,1,-1,1])
    helicity=np.repeat(np.array(c.purity)*signs,2)
    weights=polarization_weights(cosine,helicity)
    transverse=np.maximum(0,(r*r).sum(1)[:,None]-(r@directions.T)**2)
    sat=np.repeat(c.power,2)[None,:]*np.exp(-2*transverse/c.waist**2)
    zeeman=MUB_HBAR*bm[:,None,None]*((2/3)*ME-.5*MG[:,None])
    delta=GAMMA*np.repeat(c.detuning,2)[None,:,None,None]-K*(v@directions.T)[:,:,None,None]-zeeman[:,None,:,:]
    return GAMMA/2*sat[:,:,None,None]*weights[:,:,None,:]*CG[None,None,:,:]/(1+4*(delta/GAMMA)**2)

def rate_generator(w):
    """Column generator: dp/dt = A p, prior to normalization replacement."""
    n=len(w); total=w.sum(axis=1); a=np.zeros((n,12,12))
    for g in range(5):
        for q in range(3):
            e=5+EINDEX[g,q]; t=total[:,g,q]; decay=GAMMA*CG[g,q]
            a[:,e,g]+=t; a[:,g,g]-=t
            a[:,g,e]+=t+decay; a[:,e,e]-=t+decay
    return a

def populations(w):
    a=rate_generator(w); n=len(w)
    a[:,-1,:]=1.; rhs=np.zeros((n,12)); rhs[:,-1]=1.
    return np.linalg.solve(a,rhs[...,None])[...,0]

def force_and_diffusion(r,v,c,gravity=True):
    w=rates(r,v,c); active=w.sum(axis=(1,2,3))>GAMMA*1e-8
    p=np.zeros((len(w),12)); p[:,:5]=.2
    if active.any(): p[active]=populations(w[active])
    difference=p[:,:5,None]-p[:,5+EINDEX]
    scattering=(w*difference[:,None,:,:]).sum(axis=(2,3))
    f=HBAR*K*(scattering@beam_directions(c))
    if gravity: f[:,2]-=MASS*9.80665
    # Independent-event diffusion approximation: stimulated absorption/emission
    # both add variance; spontaneous emission is approximated isotropic.
    event_rate=(w*(p[:,:5,None]+p[:,5+EINDEX])[:,None,:,:]).sum(axis=(2,3))
    d=beam_directions(c)
    covariance=(HBAR*K)**2*(np.einsum('nb,bi,bj->nij',event_rate,d,d)
                + GAMMA*p[:,5:].sum(1)[:,None,None]*np.eye(3)[None,:,:]/3)
    return f,covariance

def force(r,v,c,gravity=True):
    return force_and_diffusion(r,v,c,gravity)[0]

def sample_atoms(n,seed):
    """Declared synthetic incoming cloud, not a thermal vapor capture rate."""
    rng=np.random.default_rng(seed)
    return rng.normal(0,.006,(n,3)),rng.normal(0,8.,(n,3))

def simulate(c,r,v,dt=1e-4,duration=.025,history=False,recoil=False,seed=None):
    """Midpoint mean force; optional Gaussian recoil kicks at step end.

    Recoil is an independent-event, isotropic-spontaneous approximation, not an
    exact quantum-jump model. Mean-force mode is retained for v1 reproducibility.
    Last step is shortened so requested duration is never rounded away.
    """
    if not np.isfinite([dt,duration]).all() or dt<=0 or duration<0:
        raise ValueError('positive finite dt and nonnegative finite duration required')
    r=np.array(r,dtype=float,copy=True);v=np.array(v,dtype=float,copy=True)
    if r.ndim!=2 or r.shape[1]!=3 or r.shape!=v.shape or not np.isfinite(r).all() or not np.isfinite(v).all():
        raise ValueError('finite matching (n,3) positions and velocities required')
    records=[];rng=np.random.default_rng(seed); elapsed=0.
    for step in range(int(np.ceil(duration/dt))):
        h=min(dt,duration-elapsed)
        if h<=0: break
        a=force(r,v,c)/MASS
        vm=v+a*h/2; rm=r+v*h/2
        f,cov=force_and_diffusion(rm,vm,c)
        v+=f/MASS*h; r+=vm*h
        if recoil:
            eig,vec=np.linalg.eigh(cov)
            kick=np.einsum('nij,nj->ni',vec,np.sqrt(np.maximum(eig,0)*h)*rng.normal(size=v.shape))/MASS
            v+=kick
        elapsed+=h
        if history and step%10==0: records.append(np.concatenate([r,v],axis=1).copy())
    retained=(np.linalg.norm(r,axis=1)<.002)&(np.linalg.norm(v,axis=1)<.5)
    score=np.exp(-np.sum((r/.003)**2,axis=1)-np.sum((v/1.)**2,axis=1))
    return {'r':r,'v':v,'retained':retained,'score':score,'history':np.array(records),'elapsed':elapsed}
