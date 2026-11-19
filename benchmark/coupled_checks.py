import json
from pathlib import Path
import numpy as np
from coupled_motion import run
out=Path('results/coupled');rows=[]
for geom,n in [('axial6',6),('oblique8',8)]:
 p=np.random.default_rng(721).uniform(0,2*np.pi,n)
 for arch,df in [('common',np.zeros(n)),('offset_pairs',np.repeat([0.,.12,.29,.47][:n//2],2))]:
  for check,T,rtol,ref in [('tight_tolerance',2000.,2e-10,False),('public_reference',200.,2e-8,True)]:
   row,*_=run(geom,p,df,T=T,rtol=rtol,reference=ref);row.update(check=check,architecture=arch);rows.append(row)
   (out/'checks.json').write_text(json.dumps(rows,indent=2));print(geom,arch,check,flush=True)
