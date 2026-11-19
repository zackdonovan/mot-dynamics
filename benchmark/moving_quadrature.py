import json
from pathlib import Path
import numpy as np
from moving_coherence import run
rows=[]
for geom,n in [('axial6',6),('oblique8',8)]:
 p=np.random.default_rng(721).uniform(0,2*np.pi,n)
 for arch,offsets in [('common',np.zeros(n)),('offset_pairs',np.repeat([0.,.12,.29,.47][:n//2],2))]:
  r,_,_=run(geom,p,offsets,[.05,0,0],sample_dt=.5)
  r['architecture']=arch;rows.append(r)
Path('results/moving/quadrature.json').write_text(json.dumps(rows,indent=2))
