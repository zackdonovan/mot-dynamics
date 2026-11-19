import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path('results/multibeam');data=json.loads((out/'phase_ensemble.json').read_text());checks=json.loads((out/'checks.json').read_text())
summary=[];paired=[]
for geom in ['axial6','oblique8']:
 rows=[x for x in data if x['geometry']==geom];zero=np.array([x['force'] for x in rows[0]['rows']]);near=np.array([x['force'] for x in rows[1]['rows']])
 for a in rows[2:]:
  d=np.array([x['force'] for x in a['rows']])-zero
  paired.append(dict(geometry=geom,field=a['field'],mean_delta_force=d.mean(0).tolist(),standard_error=(d.std(0,ddof=1)/np.sqrt(len(d))).tolist(),rate_delta_force=(np.array(a['rate_force'])-rows[0]['rate_force']).tolist()))
 summary.append(dict(geometry=geom,max_paired_zero_field_force_change=float(np.linalg.norm(near-zero,axis=1).max())))
flat=[r for x in data for r in x['rows']]
stats=dict(n_steady_states=len(flat),max_residual=max(r['residual'] for r in flat),max_trace_error=max(r['trace_error'] for r in flat),max_hermiticity_error=max(r['hermiticity_error'] for r in flat),minimum_density_eigenvalue=min(r['min_eigenvalue'] for r in flat),near_zero_checks=summary,paired_field_response=paired,checks=checks)
(out/'summary.json').write_text(json.dumps(stats,indent=2))
fig,axes=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
for ax,geom in zip(axes,['axial6','oblique8']):
 d=next(x for x in data if x['geometry']==geom and x['field']==[.1,0,0]);f=np.array([r['force'] for r in d['rows']])
 ax.hist(f[:,0],bins=18,alpha=.65,label='128 static phase realizations')
 ax.axvline(d['rate_force'][0],color='black',ls='--',label='population rate model')
 ax.axvline(f[:,0].mean(),color='tab:red',label='coherent phase mean')
 ax.set(title=geom+'; B = (0.1, 0, 0)',xlabel='x force / ℏkΓ',ylabel='count');ax.legend(fontsize=7)
fig.savefig(out/'phase_force_distribution.png',dpi=180)
print(json.dumps(stats,indent=2))
