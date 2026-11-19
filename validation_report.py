import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=Path('results/v2');rows=json.loads((p/'followup.json').read_text())
assert len(rows)==36,'validation matrix incomplete'
lines=['# validation results (v2)','', 'Independent, synthetic follow-up. No laboratory validation; no new optimizer selection.','', '| geometry | setting | input velocity sigma (m/s) | dynamics | retained 25 ms | retained 75 ms | lost from 25 ms cohort |','|---|---|---:|---|---:|---:|---:|']
for x in rows:
 lines.append(f"| {x['geometry']} | {x['controls_name']} | {x['velocity_sigma']:g} | {'recoil approximation' if x['recoil'] else 'mean force'} | {x['25ms']['count']}/512 ({x['25ms']['fraction']:.1%}) | {x['75ms']['count']}/512 ({x['75ms']['fraction']:.1%}) | {x['75ms']['lost_from_25ms']} |")
lines+=['','## paired comparisons at 75 ms','', '95% percentile bootstrap intervals resample paired atoms (5,000 draws); conditional on these controls and approximate physics. No optimization-seed or physical-model uncertainty.','', '| geometry | source sigma | dynamics | GP minus standard (percentage points) | GP minus Sobol (percentage points) |','|---|---:|---|---|---|']
comp=[]
for geom in ['axial6','oblique8']:
 for sigma in [4.,8.,12.]:
  for recoil in [False,True]:
   masks={name:np.load(p/f'{geom}_{name}_s{sigma}_recoil{recoil}.npz')['mask75'] for name in ['standard','sobol_11','gp_ucb_11']}
   rng=np.random.default_rng(582);idx=rng.integers(0,512,(5000,512));res={}
   for name in ['standard','sobol_11']:
    d=masks['gp_ucb_11'].astype(float)-masks[name];ci=np.percentile(d[idx].mean(1)*100,[2.5,97.5]);res[name]={'difference_pp':float(d.mean()*100),'ci95_pp':ci.tolist()}
   comp.append({'geometry':geom,'velocity_sigma':sigma,'recoil':recoil,'comparisons':res})
   def fmt(name):
    z=res[name];return f"{z['difference_pp']:.1f} [{z['ci95_pp'][0]:.1f}, {z['ci95_pp'][1]:.1f}]"
   lines.append(f"| {geom} | {sigma:g} | {'recoil approximation' if recoil else 'mean force'} | {fmt('standard')} | {fmt('sobol_11')} |")
(p/'paired_comparisons.json').write_text(json.dumps(comp,indent=2))
lines+=['','## model limitation exposed','', 'At a fixed velocity (0.1, 0.2, 0.3) m/s, moving the position only 1 nm away from the origin along different directions changes the computed force by up to 33.7% for the six-beam finalists and 15.7% for the eight-beam finalists, relative to the origin fallback. This is a local-axis/incoherent-population approximation failure at vanishing field. It prevents claiming high-fidelity near-center physics; more samples or better optimization do not resolve it.','', 'Initial-cloud population-relaxation medians span 2.3–7.1 microseconds; the slowest sampled case is 425 microseconds. In the declared initial cloud, the 95th-percentile product of relaxation time and magnetic-axis turning rate is below 0.023. That supports timescale separation for many initial states, not all trajectories or near-zero-field coherences.','', 'See V2_METHODS.md for the independent-event recoil approximation and uncertainty scope. Snapshot retention is not proof of continuous or long-term confinement.']
summary=Path('V2_SUMMARY.md').read_text() if Path('V2_SUMMARY.md').exists() else ''
Path('VALIDATION_RESULTS.md').write_text(summary+'\n## complete matrix\n\n'+'\n'.join(lines[4:])+'\n')
fig,axes=plt.subplots(1,2,figsize=(11,4),sharey=True)
for ax,geom in zip(axes,['axial6','oblique8']):
 for name,col in [('standard','gray'),('sobol_11','#dd9128'),('gp_ucb_11','#2878b5')]:
  for recoil in [False,True]:
   ss=sorted([x for x in rows if x['geometry']==geom and x['controls_name']==name and x['recoil']==recoil],key=lambda x:x['velocity_sigma'])
   ax.plot([x['velocity_sigma'] for x in ss],[100*x['75ms']['fraction'] for x in ss],color=col,linestyle='-' if recoil else '--',marker='o',label=name+(' + recoil' if recoil else ' mean'))
 ax.set_title(geom);ax.set_xlabel('initial velocity-component sigma (m/s)');ax.grid(alpha=.2);ax.set_ylim(0,100)
axes[0].set_ylabel('75-ms synthetic retention (%)');axes[1].legend(fontsize=8)
fig.suptitle('Fixed controls: source distribution matters; solid lines use approximate recoil')
fig.tight_layout();fig.savefig(p/'source_sensitivity.png',dpi=160)
