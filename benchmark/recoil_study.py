import json,time
from pathlib import Path
import numpy as np
from spatial_motion import run
out=Path('results/recoil_validity');out.mkdir(exist_ok=True);rows=[]
for geom in ['axial6','oblique8']:
 for name,c in [('origin',[0.,0,0]),('off_axis',[.006,.003,.002])]:
  for T in [2000.,8000.]:
   start=time.time();r,t,x,v,f=run(geom,c,T=T);r.update(position_name=name,wall_time_s=time.time()-start);rows.append(r)
   np.savez_compressed(out/f'{geom}_{name}_{int(T)}.npz',time=t,displacement=x,v=v,force=f)
   (out/'study.json').write_text(json.dumps(rows,indent=2));print(geom,name,T,r['spontaneous_recoil_diagnostic']['expected_spontaneous_photons'],r['wall_time_s'],flush=True)
