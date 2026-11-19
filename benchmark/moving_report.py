import json,os
from pathlib import Path
import numpy as np
os.environ.setdefault('MPLCONFIGDIR','/tmp/mot-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path('results/moving');rows=json.loads((out/'study.json').read_text());checks=json.loads((out/'checks.json').read_text());quad=json.loads((out/'quadrature.json').read_text());comparisons=[]
for c in checks+quad:
 arch=c.get('architecture','common');b=next(r for r in rows if r['geometry']==c['geometry'] and r['seed']==721 and r['architecture']==arch and r['velocity'][0]==(0 if c.get('check')=='stationary_direct_steady' else .05))
 f=c['force'] if 'force' in c else c['mean_tail_force']
 comparisons.append(dict(geometry=c['geometry'],architecture=arch,check=c.get('check','quadrature_dt_0.5'),force_difference_norm=float(np.linalg.norm(np.array(f)-b['mean_tail_force'])),reference_force=b['mean_tail_force'],checked_force=f))
summary=dict(n_main_runs=len(rows),max_public_pylcp_density_error=max(r.get('pylcp_time_density_error',0) for r in rows),max_hamiltonian_error=max(r['hamiltonian_reconstruction_error'] for r in rows),max_trace_error=max(r['max_trace_error'] for r in rows),min_density_eigenvalue=min(r['min_density_eigenvalue'] for r in rows),max_half_tail_window_change=max(np.linalg.norm(np.array(r['last_quarter_force'])-r['previous_quarter_force']) for r in rows),checks=comparisons)
(out/'summary.json').write_text(json.dumps(summary,indent=2))
fig,axs=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
for ax,geom in zip(axs,['axial6','oblique8']):
 for arch,color in [('common','tab:blue'),('offset_pairs','tab:orange')]:
  for seed in [721,722]:
   group=[r for r in rows if r['geometry']==geom and r['architecture']==arch and r['seed']==seed]
   ax.plot([r['velocity'][0] for r in group],[r['mean_tail_force'][0] for r in group],'-o',color=color,alpha=.7,label=arch+' coherent' if seed==721 else None)
   if seed==721:ax.plot([r['velocity'][0] for r in group],[r['rate_force'][0] for r in group],'--',color=color,label=arch+' rates')
 ax.set(title=geom,xlabel='prescribed vx / (Γ/k)',ylabel='tail mean x force / ℏkΓ');ax.legend(fontsize=7)
fig.savefig(out/'moving_force_comparison.png',dpi=170)
print(json.dumps(summary,indent=2))
