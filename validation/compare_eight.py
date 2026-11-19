from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from mot.physics import *
from reference import steady
rng=np.random.default_rng(183);worst=0.
for i in range(20):
    c=Controls(detuning=tuple(rng.uniform(-4,-.8,4)),purity=tuple(rng.uniform(-1,1,4)),power=(.75,)*4,geometry='oblique8')
    r=rng.uniform(-.008,.008,3);v=rng.uniform(-10,10,3);d=beam_directions(c)
    s=np.repeat(c.power,2)*np.exp(-2*((r*r).sum()-(r@d.T)**2)/c.waist**2)
    b=r*np.array([c.gradient,c.gradient,-2*c.gradient])
    f,p,L=steady(d*K,np.repeat(c.detuning,2)*GAMMA,s,np.repeat(np.array(c.purity)*[1,1,-1,1],2),b,v)
    actual=force(r,v,c,False)[0];np.testing.assert_allclose(actual,f,atol=1e-33,rtol=2e-9)
    worst=max(worst,np.linalg.norm(actual-f)/np.linalg.norm(f))
print('20 randomized eight-beam independent force comparisons passed; maximum relative error:',worst)
