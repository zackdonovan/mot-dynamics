"""Independent PyLCP comparison of actual 12-state production force (no coherence)."""
import json,sys,time
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np,pylcp
from mot.physics import *

def pylcp_rate_force(r,v,c):
    dirs=beam_directions(c);transverse=np.maximum(0,np.dot(r,r)-(dirs@r)**2)
    sat=np.repeat(c.power,2)*np.exp(-2*transverse/c.waist**2)
    h=np.repeat(np.array(c.purity)*([1,1,-1] if c.geometry=='axial6' else [1,1,-1,1]),2)
    beams=[]
    # Splitting each mixture into two beams is appropriate for rate equations
    # only. It would imply coherent fields in OBE and is not used there.
    for d,s,helicity,det in zip(dirs,sat,h,np.repeat(c.detuning,2)):
        for pol in [-1,1]:
            power=s*(1+pol*helicity)/2
            if power>0:beams.append(dict(kvec=d,pol=pol,delta=float(det),s=float(power)))
    lasers=pylcp.laserBeams(beams,beam_type=pylcp.infinitePlaneWaveBeam)
    Hg,mg=pylcp.hamiltonians.singleF(F=2,gF=.5,muB=1)
    He,me=pylcp.hamiltonians.singleF(F=3,gF=2/3,muB=1)
    ham=pylcp.hamiltonian(Hg,He,mg,me,pylcp.hamiltonians.dqij_two_bare_hyperfine(2,3))
    b=r*np.array([c.gradient,c.gradient,-2*c.gradient])*MUB_HBAR/GAMMA
    re=pylcp.rateeq(lasers,pylcp.magField(b),ham,include_mag_forces=False,v0=v*K/GAMMA)
    return re.find_equilibrium_force()*HBAR*K*GAMMA

if __name__=='__main__':
    out=[];rng=np.random.default_rng(741902);start=time.time()
    for geom,n in [('axial6',3),('oblique8',4)]:
        for i in range(20):
            c=Controls(detuning=tuple(rng.uniform(-4,-.3,n)),purity=tuple(rng.uniform(-1,1,n)),power=tuple(rng.uniform(.1,1.3,n)),geometry=geom)
            r=rng.normal(0,.004,3);v=rng.normal(0,5,3)
            a=force(r,v,c,gravity=False)[0];b=pylcp_rate_force(r,v,c)
            out.append(dict(geometry=geom,r=r.tolist(),v=v.tolist(),controls=c.__dict__,production=a.tolist(),pylcp=b.tolist(),relative_error=float(np.linalg.norm(a-b)/max(np.linalg.norm(b),1e-35))))
    dest=Path('results/pylcp');dest.mkdir(exist_ok=True)
    (dest/'production_rate_comparison.json').write_text(json.dumps(out,indent=2))
    print('n',len(out),'max_relative_error',max(x['relative_error'] for x in out),'runtime_s',time.time()-start)
