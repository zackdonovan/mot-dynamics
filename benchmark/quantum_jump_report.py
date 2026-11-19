import json,os
from pathlib import Path
import numpy as np
os.environ.setdefault('MPLCONFIGDIR','/tmp/mot-mpl')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=Path('results/quantum_jumps');a=np.load(p/'trajectories.npz');ref=np.load(p/'reference.npz');mom=json.loads((p/'count_moments.json').read_text());rows=json.loads((p/'convergence.json').read_text())
assert len(a['counts'])==16384
count=a['counts'];re=a['recoils'];st=a['terminal_states'];n=len(count)
cm=count-count.mean();var=count.var(ddof=1);var_se=np.sqrt((np.mean(cm**4)-(n-3)/(n-1)*var**2)/n)
center=re-re.mean(0);variance=re.var(axis=0,ddof=1);vse=np.sqrt((np.mean(center**4,axis=0)-(n-3)/(n-1)*variance**2)/n)
blocks=[]
for i in range(4):
 block=st[i*4096:(i+1)*4096];rho=block.T@block.conj()/4096;blocks.append(float(np.linalg.norm(rho-ref['rho'])))
s=dict(n=n,count_mean=float(count.mean()),reference_count_mean=mom['mean'],count_mean_z=float((count.mean()-mom['mean'])/(count.std(ddof=1)/np.sqrt(n))),count_variance=float(var),reference_count_variance=mom['variance'],count_variance_z=float((var-mom['variance'])/var_se),sample_fano=float(var/count.mean()),reference_fano=mom['fano'],recoil_mean_z=(re.mean(0)/(re.std(axis=0,ddof=1)/np.sqrt(n))).tolist(),recoil_variance_z=((variance-np.diag(ref['recoil_covariance']))/vse).tolist(),independent_4096_block_density_errors=blocks,final_density_error=rows[-1]['density_frobenius_error'],final_covariance_relative_error=rows[-1]['covariance_relative_frobenius_error'])
(p/'statistical_checks.json').write_text(json.dumps(s,indent=2));print(json.dumps(s,indent=2))
fig,axs=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
N=np.array([r['n'] for r in rows]);err=np.array([r['density_frobenius_error'] for r in rows]);axs[0].loglog(N,err,'o-',label='sampled density error');axs[0].loglog(N,err[1]*np.sqrt(N[1]/N),'--',label='N^−1/2 guide');axs[0].set(xlabel='quantum trajectories',ylabel='density Frobenius error');axs[0].legend(fontsize=8)
axs[1].hist(count,bins=np.arange(count.max()+2)-.5,density=True,alpha=.7);axs[1].set(xlabel='emitted photons in 200 Γ⁻¹',ylabel='frequency',title=f'count Fano factor {s["sample_fano"]:.2f}; reference {mom["fano"]:.2f}')
fig.savefig(p/'quantum_jump_validation.png',dpi=170)
