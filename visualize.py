"""Deterministic illustrative plots, not used for model selection."""
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from mot.physics import *
from pathlib import Path
p=Path('results');c=Controls(**json.loads((p/'heldout.json').read_text())['gp_ucb_11']['controls'])
cloud=sample_atoms(32,9997)
fig,axes=plt.subplots(1,2,figsize=(11,4))
for ax,label,cfg in zip(axes,['standard','GP seed11'],[Controls(),c]):
    out=simulate(cfg,*cloud,history=True);history=out['history'];
    for j in range(32):ax.plot(history[:,j,0]*1000,history[:,j,2]*1000,alpha=.4,lw=.8)
    ax.set(title=f'{label}, same 32 new atoms',xlabel='x (mm)',ylabel='z (mm)',xlim=(-30,30),ylim=(-30,30));ax.add_patch(plt.Circle((0,0),2,fill=False,color='black'))
fig.suptitle('Trajectory projections, clipped view; circle is spatial proxy only');fig.tight_layout();fig.savefig(p/'trajectories.png',dpi=160)
fig,axes=plt.subplots(1,2,figsize=(11,4))
for i,label in enumerate(['x','y','z']):
    grid=np.linspace(-.006,.006,101);r=np.zeros((101,3));r[:,i]=grid
    axes[0].plot(grid*1e3,force(r,r*0,Controls(),False)[:,i]/MASS,label=label)
    grid=np.linspace(-12,12,101);v=np.zeros((101,3));v[:,i]=grid
    axes[1].plot(grid,force(v*0,v,Controls(),False)[:,i]/MASS,label=label)
axes[0].set(xlabel='displacement (mm)',ylabel='acceleration (m/s²)',title='Restoring force, v=0')
axes[1].set(xlabel='velocity (m/s)',ylabel='acceleration (m/s²)',title='Damping force, r=0')
for ax in axes:ax.legend();ax.axhline(0,color='black',lw=.5)
fig.tight_layout();fig.savefig(p/'force_profiles.png',dpi=160)
