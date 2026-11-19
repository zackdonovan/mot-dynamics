import json
from pathlib import Path
import numpy as np
from moving_coherence import run
from pylcp_multibeam import solve
out=Path('results/moving');rows=[]
for geom,n in [('axial6',6),('oblique8',8)]:
 p=np.random.default_rng(721).uniform(0,2*np.pi,n)
 for architecture,offsets in [('common',np.zeros(n)),('offset_pairs',np.repeat([0.,.12,.29,.47][:n//2],2))]:
  r,t,f=run(geom,p,offsets,[.05,0,0],T=8000.)
  r.update(check='longer_window',architecture=architecture);rows.append(r)
  r,t,f=run(geom,p,offsets,[.05,0,0],T=2000.,rtol=2e-10)
  r.update(check='tighter_tolerance',architecture=architecture);rows.append(r)
  (out/'checks.json').write_text(json.dumps(rows,indent=2))
for geom,n in [('axial6',6),('oblique8',8)]:
 p=np.random.default_rng(721).uniform(0,2*np.pi,n)
 r=solve(geom,[.1,0,0],p);rows.append(dict(check='stationary_direct_steady',geometry=geom,force=r['force']))
(out/'checks.json').write_text(json.dumps(rows,indent=2))
