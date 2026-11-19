import json,os
from pathlib import Path
import numpy as np
from scipy.stats import spearmanr
os.environ.setdefault('MPLCONFIGDIR','/tmp/mot-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
out=Path('results/multifidelity');rows=json.loads((out/'coherent.json').read_text());rate=json.loads((out/'rate.json').read_text());ds=[-1.,-2.,-3.,-4.]
assert len(rows)==64,'wait for full predeclared study'
summary={};rng=np.random.default_rng(721993);means=[];bands=[]
for batch in ['discovery','heldout']:
 seeds=sorted({r['phase_seed'] for r in rows if r['batch']==batch});a=np.array([[next(r['energy_removal_fraction'] for r in rows if r['batch']==batch and r['phase_seed']==seed and r['base_detuning']==d) for d in ds] for seed in seeds])
 boot=a[rng.integers(0,len(seeds),(20000,len(seeds)))].mean(1);mean=a.mean(0);ci=np.quantile(boot,[.025,.975],axis=0)
 summary[batch]=dict(mean=mean.tolist(),bootstrap95=ci.tolist(),phase_sample_std=a.std(0,ddof=1).tolist(),ranking=[ds[i] for i in np.argsort(-mean)],spearman_with_rate=float(spearmanr(mean,[r['energy_removal_fraction'] for r in rate]).statistic),bootstrap_winner_frequency=[float(np.mean(np.argmax(boot,axis=1)==i)) for i in range(4)])
 means.append(mean);bands.append(ci)
 if batch=='heldout':
  best=ds.index(summary['discovery']['ranking'][0]);competitors=np.delete(np.arange(4),best)
  summary['heldout']['selected_vs_best_alternative_bootstrap95']=np.quantile(boot[:,best]-np.max(boot[:,competitors],axis=1),[.025,.975]).tolist()
summary['rate_ranking']=[r['detuning'] for r in sorted(rate,key=lambda x:-x['energy_removal_fraction'])]
summary['coherent_total_wall_s']=sum(r['wall_s'] for r in rows);summary['rate_total_wall_s']=sum(r['wall_s'] for r in rate)
summary['median_single_evaluation_speed_ratio']=float(np.median([r['wall_s'] for r in rows])/np.median([r['wall_s'] for r in rate]))
summary['max_free_recoil_rms_norm']=max(r['spontaneous_recoil_diagnostic']['free_recoil_rms_norm'] for r in rows)
stress=json.loads((out/'stress_coherent.json').read_text());sr=json.loads((out/'stress_rate.json').read_text())
assert len(stress)==16
summary['origin_stress']={str(d):dict(coherent_mean=float(np.mean([r['energy_removal_fraction'] for r in stress if r['base_detuning']==d])),samples=[r['energy_removal_fraction'] for r in stress if r['base_detuning']==d],rate=next(r['energy_removal_fraction'] for r in sr if r['detuning']==d)) for d in ds}
(out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
fig,ax=plt.subplots(figsize=(8,4),constrained_layout=True)
x=np.arange(4);ax.plot(x,[r['energy_removal_fraction'] for r in rate],'s--',label='instantaneous population model')
for j,batch in enumerate(['discovery','heldout']):
 m=means[j];ci=bands[j];ax.errorbar(x+(j-.5)*.08,m,yerr=[m-ci[0],ci[1]-m],fmt='o',capsize=4,label='coherent '+batch+' (8 phases)')
ax.set_xticks(x,[str(d) for d in ds]);ax.set(xlabel='base detuning / Γ',ylabel='fractional kinetic-energy removal',title='Conditional model comparison; not capture or experimental cooling');ax.legend(fontsize=8)
fig.savefig(out/'candidate_comparison.png',dpi=180)
