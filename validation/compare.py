import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from mot.physics import Controls,force,rates,populations,DIRECTIONS,K,GAMMA,CG
from reference import steady
import numpy as np
rng=np.random.default_rng(923)
max_relative_error=0
for i in range(30):
    c=Controls(detuning=tuple(rng.uniform(-4,-.5,3)),purity=tuple(rng.uniform(.3,1,3)),power=tuple(rng.uniform(.3,2,3)))
    r=rng.uniform(-.008,.008,3);v=rng.uniform(-10,10,3)
    transverse=(r*r).sum()-(r@DIRECTIONS.T)**2
    s=np.repeat(c.power,2)*np.exp(-2*transverse/c.waist**2)
    b=r*np.array([c.gradient,c.gradient,-2*c.gradient])
    f,p,L=steady(DIRECTIONS*K,np.repeat(c.detuning,2)*GAMMA,s,np.repeat(np.array(c.purity)*[1,1,-1],2),b,v)
    candidate=force(r,v,c,gravity=False)[0]
    np.testing.assert_allclose(candidate,f,atol=1e-33,rtol=2e-9)
    np.testing.assert_allclose(populations(rates(r,v,c))[0],p,atol=2e-12)
    max_relative_error=max(max_relative_error,np.linalg.norm(candidate-f)/np.linalg.norm(f))
print('30 randomized independent force/population comparisons passed')
print('max force relative error:',max_relative_error)
