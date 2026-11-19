import json
from pathlib import Path
import numpy as np
from spatial_motion import run,system
out=Path('results/spatial');rows=[]
for geom in ['axial6','oblique8']:
 for name,c in [('origin',[0.,0,0]),('off_axis',[.006,.003,.002])]:
  row,*_=run(geom,c,rtol=2e-10);row.update(check='tight_tolerance',position_name=name);rows.append(row)
  # Absolute-position PyLCP reference; compare physical displacement.
  n=6 if geom=='axial6' else 8;p=np.random.default_rng(721).uniform(0,2*np.pi,n);df=np.repeat([0.,.12,.29,.47][:n//2],2)
  ob,*rest=system(geom,p,df,c);center=rest[2]
  ob.set_initial_position(center);ob.set_initial_velocity(np.array([.05,0,0]));ob.set_initial_rho_from_populations(np.r_[np.ones(5)/5,np.zeros(7)])
  ob.evolve_motion([0,100.],random_recoil=False,t_eval=[100.],method='DOP853',rtol=2e-8,atol=2e-10)
  short,*_=run(geom,c,T=100.)
  from mot.physics import K
  rows.append(dict(check='public_motion_reference',geometry=geom,position_name=name,velocity_error=float(np.linalg.norm(ob.sol.v[:,-1]-short['final_v'])),displacement_error_m=float(np.linalg.norm((ob.sol.r[:,-1]-center)/K-short['final_displacement_m']))))
  (out/'checks.json').write_text(json.dumps(rows,indent=2));print(geom,name,flush=True)
