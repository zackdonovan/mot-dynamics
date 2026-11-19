"""Distributional, not pathwise, half-step check of three predeclared cases."""
import json
from pathlib import Path
import numpy as np
from mot.physics import *
p=Path('results/v2');data=json.loads(Path('results/heldout.json').read_text());rng=np.random.default_rng(20261118);r=rng.normal(0,.006,(512,3));unitv=rng.normal(size=(512,3));out=[]
for name,sigma in [('standard',8.),('gp_ucb_11',8.),('gp_ucb_11',4.)]:
 c=Controls(**data[name]['controls']);a=simulate(c,r,unitv*sigma,dt=1.25e-5,duration=.025,recoil=True,seed=17890);b=simulate(c,a['r'],a['v'],dt=1.25e-5,duration=.05,recoil=True,seed=17891)
 out.append({'geometry':'axial6','controls_name':name,'velocity_sigma':sigma,'dt':1.25e-5,'count25':int(a['retained'].sum()),'count75':int(b['retained'].sum()),'fraction75':float(b['retained'].mean()),'note':'not Brownian path-coupled: different increments at half step; distributional comparison only, single seed'})
 (p/'recoil_halfstep.json').write_text(json.dumps(out,indent=2));print(out[-1],flush=True)
