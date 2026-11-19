import json,os
from pathlib import Path
import numpy as np
os.environ.setdefault('MPLCONFIGDIR','/tmp/mot-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path('results/spatial');rows=json.loads((out/'study.json').read_text());checks=json.loads((out/'checks.json').read_text());pairs=[];conv=[]
for a in rows:
 if not a['frozen']:
  b=next(x for x in rows if x['geometry']==a['geometry'] and x['position_name']==a['position_name'] and x['frozen'])
  pairs.append(dict(geometry=a['geometry'],position_name=a['position_name'],final_velocity_difference=float(np.linalg.norm(np.array(a['final_v'])-b['final_v'])),max_field_change=a['max_field_change'],max_relative_intensity_change=a['max_relative_intensity_change']))
for a in checks:
 if a['check']=='tight_tolerance':
  b=next(x for x in rows if x['geometry']==a['geometry'] and x['position_name']==a['position_name'] and not x['frozen'])
  conv.append(dict(geometry=a['geometry'],position_name=a['position_name'],velocity_error=float(np.linalg.norm(np.array(a['final_v'])-b['final_v'])),displacement_error_m=float(np.linalg.norm(np.array(a['final_displacement_m'])-b['final_displacement_m']))))
full=[x for x in rows if not x['frozen']];refs=[x for x in checks if x['check']=='public_motion_reference']
s=dict(paired_frozen_comparisons=pairs,tolerance_checks=conv,max_public_force_error=max(x['max_public_force_error'] for x in full),max_public_hamiltonian_error=max(x['max_public_hamiltonian_error'] for x in full),max_trace_error=max(x['max_trace_error'] for x in rows),min_eigenvalue=min(x['min_eigenvalue'] for x in rows),max_public_velocity_error=max(x['velocity_error'] for x in refs),max_public_displacement_error_m=max(x['displacement_error_m'] for x in refs))
(out/'summary.json').write_text(json.dumps(s,indent=2))
fig,axs=plt.subplots(2,2,figsize=(10,7),constrained_layout=True)
for col,geom in enumerate(['axial6','oblique8']):
 for row,name in enumerate(['origin','off_axis']):
  ax=axs[row,col]
  for frozen in [False,True]:
   a=np.load(out/f'{geom}_{name}_{frozen}.npz');ax.plot(a['time'],a['v'][:,0],ls='--' if frozen else '-',label='locally frozen' if frozen else 'full spatial')
  ax.set(title=geom+' / '+name,xlabel='time × Γ',ylabel='vx / (Γ/k)');ax.legend(fontsize=8)
fig.savefig(out/'spatial_vs_local.png',dpi=160)
print(json.dumps(s,indent=2))
