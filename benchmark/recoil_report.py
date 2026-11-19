import json,os
from pathlib import Path
import numpy as np
os.environ.setdefault('MPLCONFIGDIR','/tmp/mot-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mot.physics import GAMMA,K
out=Path('results/recoil_validity');rows=json.loads((out/'study.json').read_text());q=json.loads((out/'quadrature.json').read_text());qc=[]
for a in q:
 b=next(x for x in rows if x['geometry']==a['geometry'] and x['position_name']==a['position_name'] and x['T']==2000)
 aa=a['spontaneous_recoil_diagnostic'];bb=b['spontaneous_recoil_diagnostic']
 qc.append(dict(position_name=a['position_name'],relative_count_change=abs(aa['expected_spontaneous_photons']/bb['expected_spontaneous_photons']-1),covariance_change_norm=float(np.linalg.norm(np.array(aa['free_recoil_velocity_covariance'])-bb['free_recoil_velocity_covariance']))))
long=[r for r in rows if r['T']==8000];seconds=[r['wall_time_s'] for r in long];scale=.075/(8000/GAMMA)
s=dict(quadrature_checks=qc,max_rate_closure_error=max(r['spontaneous_recoil_diagnostic']['max_emission_rate_error'] for r in rows),min_tensor_eigenvalue=min(r['spontaneous_recoil_diagnostic']['min_instantaneous_tensor_eigenvalue'] for r in rows),long_run_wall_seconds=seconds,estimated_75ms_single_trajectory_hours=[x*scale/3600 for x in seconds],estimated_512_trajectory_serial_days=[x*scale*512/86400 for x in seconds])
(out/'summary.json').write_text(json.dumps(s,indent=2));print(json.dumps(s,indent=2))
fig,ax=plt.subplots(figsize=(9,4),constrained_layout=True)
labels=[];values=[]
for r in rows:
 labels.append(r['geometry']+'\n'+r['position_name']+'\n'+str(round(r['T']/GAMMA*1e6))+' μs');values.append(r['spontaneous_recoil_diagnostic']['free_recoil_rms_norm']*GAMMA/K)
ax.bar(range(len(values)),values,color=['tab:blue' if r['T']==2000 else 'tab:orange' for r in rows]);ax.axhline(.05*GAMMA/K,color='black',ls='--',label='initial deterministic speed')
ax.set_xticks(range(len(values)),labels,fontsize=8);ax.set(ylabel='free spontaneous-recoil RMS norm (m/s)',title='Noise budget, not a simulated final velocity spread');ax.legend(fontsize=8)
fig.savefig(out/'recoil_budget.png',dpi=170)
