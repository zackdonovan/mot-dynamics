import json
from pathlib import Path
from spatial_motion import run
rows=[]
for name,c in [('origin',[0.,0,0]),('off_axis',[.006,.003,.002])]:
 r,*_=run('axial6',c,sample_dt=.25);r['position_name']=name;rows.append(r)
Path('results/recoil_validity/quadrature.json').write_text(json.dumps(rows,indent=2))
