import json,os
from pathlib import Path
import numpy as np
os.environ.setdefault('MPLCONFIGDIR','/tmp/mot-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path('results/coupled');rows=json.loads((out/'study.json').read_text());checks=json.loads((out/'checks.json').read_text());quad=json.loads((out/'quadrature.json').read_text());diff=[]
for c in checks:
 if c['check']=='tight_tolerance':
  b=next(r for r in rows if r['geometry']==c['geometry'] and r['architecture']==c['architecture'] and r['seed']==721)
  diff.append(dict(geometry=c['geometry'],architecture=c['architecture'],velocity_error=float(np.linalg.norm(np.array(c['final_v'])-b['final_v'])),position_error=float(np.linalg.norm(np.array(c['final_r'])-b['final_r']))))
refs=[x for x in checks if x['check']=='public_reference']
s=dict(n_main_runs=len(rows),max_public_force_error=max(r['max_public_force_error'] for r in rows),max_trace_error=max(r['max_trace_error'] for r in rows),min_density_eigenvalue=min(r['min_eigenvalue'] for r in rows),max_public_density_error=max(r['public_density_error'] for r in refs),max_public_velocity_error=max(r['public_velocity_error'] for r in refs),max_public_position_error=max(r['public_position_error'] for r in refs),tolerance_checks=diff,coarse_max_momentum_error=max(r['momentum_quadrature_error'] for r in rows),fine_max_momentum_error=max(r['momentum_quadrature_error'] for r in quad))
(out/'summary.json').write_text(json.dumps(s,indent=2))
fig,axs=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
for ax,geom in zip(axs,['axial6','oblique8']):
 for seed in [721,722]:
  for arch,color in [('common','tab:blue'),('offset_pairs','tab:orange')]:
   a=np.load(out/f'{geom}_{seed}_{arch}.npz');ax.plot(a['time'],a['v'][:,0],color=color,alpha=.7,label=arch if seed==721 else None)
 ax.axhline(.05,color='black',ls='--',label='initial / previously imposed vx');ax.set(title=geom,xlabel='time × Γ',ylabel='vx / (Γ/k)');ax.legend(fontsize=7)
fig.savefig(out/'coupled_velocity.png',dpi=170)
print(json.dumps(s,indent=2))
