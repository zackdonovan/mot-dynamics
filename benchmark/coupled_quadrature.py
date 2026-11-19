import json
from pathlib import Path
import numpy as np
from coupled_motion import run
out=Path('results/coupled');rows=[]
for geom,n in [('axial6',6),('oblique8',8)]:
 p=np.random.default_rng(721).uniform(0,2*np.pi,n)
 for arch,df in [('common',np.zeros(n)),('offset_pairs',np.repeat([0.,.12,.29,.47][:n//2],2))]:
  row,*_=run(geom,p,df,sample_dt=.25);row.update(architecture=arch);rows.append(row)
(out/'quadrature.json').write_text(json.dumps(rows,indent=2))
