import json,importlib.metadata,platform
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path('results/pylcp')
toy=json.loads((P/'toy.json').read_text());prod=json.loads((P/'production_rate_comparison.json').read_text());zero=json.loads((P/'zero_field_coherent.json').read_text());metrics=json.loads((P/'toy_metrics.json').read_text())
rotation=max(np.ptp([r['pylcp_obe_force'][r['axis']] for r in toy if r['pair'] and r['s']==s and r['beta']==b]) for s in [.05,1.] for b in [-2.,-1.,-.1,0.,.1,1.,2.])
summary={'production_rate_cases':len(prod),'max_production_relative_error':max(r['relative_error']for r in prod),'toy_metrics':metrics,'max_toy_rotation_absolute_spread':float(rotation),'max_zero_steady_residual':max(r['steady_residual']for r in zero),'max_zero_time_evolution_density_error':max(r.get('time_evolution_vs_direct_rho_norm',0)for r in zero),'nearzero_rate_z_range':[min(r['pylcp_rate_force'][2]for r in zero[1:6]),max(r['pylcp_rate_force'][2]for r in zero[1:6])],'nearzero_obe_z_range':[min(r['fixed_basis_obe_force'][2]for r in zero[1:6]),max(r['fixed_basis_obe_force'][2]for r in zero[1:6])]}
assert summary['max_production_relative_error']<1e-9
assert metrics['max_single_obe_exact_error']<1e-8
assert metrics['max_steady_residual']<1e-7
assert summary['max_zero_time_evolution_density_error']<1e-7
assert summary['max_toy_rotation_absolute_spread']<1e-8
(P/'summary.json').write_text(json.dumps(summary,indent=2))
(P/'environment.json').write_text(json.dumps({'python':platform.python_version(),'packages':{name:importlib.metadata.version(name)for name in ['pylcp','numpy','scipy','sympy','numba','llvmlite']}},indent=2))
fig,ax=plt.subplots(1,2,figsize=(10.5,4))
for s,col in [(.05,'#397fb0'),(1.,'#c07a20')]:
 rr=[x for x in toy if x['pair']and x['axis']==2 and x['s']==s]
 ax[0].plot([x['beta'] for x in rr],[x['pylcp_rate_force'][2]for x in rr],'--',color=col,label=f'rate, s={s}')
 ax[0].plot([x['beta'] for x in rr],[x['pylcp_obe_force'][2]for x in rr],'o-',color=col,label=f'OBE, s={s}',markersize=4)
ax[0].set_title("F=0→1 opposite-beam toy, atom at rest")
ax[0].set_xlabel('dimensionless magnetic field');ax[0].set_ylabel('force / (ℏkΓ)');ax[0].legend(fontsize=8);ax[0].grid(alpha=.2)
labels=['x','y','z','(1,1,1)','(1,−1,2)'];x=np.arange(5)
ax[1].bar(x-.18,[r['pylcp_rate_force'][2] for r in zero[1:6]],width=.36,label='rate equation',color='#c07a20')
ax[1].bar(x+.18,[r['fixed_basis_obe_force'][2]for r in zero[1:6]],width=.36,label='fixed-basis OBE',color='#397fb0')
ax[1].set_xticks(x,labels,rotation=20);ax[1].set_ylabel('force / (ℏkΓ)');ax[1].set_xlabel('direction of near-zero magnetic field')
ax[1].set_title('F=2→3 single-beam toy, |B|=10⁻⁸');ax[1].legend(fontsize=8,loc='lower right')
fig.tight_layout();fig.savefig(P/'coherent_comparison.png',dpi=170)
print(json.dumps(summary,indent=2))
