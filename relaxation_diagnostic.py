"""Frozen-field population relaxation diagnostics; not a coherent-state audit."""
import json
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from mot.physics import *
OUT=Path('results/v2');OUT.mkdir(exist_ok=True)
records=[]
for geometry,folder in [('axial6',Path('results')),('oblique8',Path('results/eight_beam'))]:
    data=json.loads((folder/'heldout.json').read_text())
    for name in ['standard','sobol_11','gp_ucb_11']:
        c=Controls(**data[name]['controls']);r,v=sample_atoms(512,20261119)
        w=rates(r,v,c);a=rate_generator(w);eig=np.linalg.eigvals(a)
        real=np.sort(eig.real,axis=1);gap=-real[:,-2];tau=1/np.maximum(gap,1e-300)
        b=r*np.array([c.gradient,c.gradient,-2*c.gradient]);bdot=v*np.array([c.gradient,c.gradient,-2*c.gradient]);norm=np.linalg.norm(b,axis=1)
        turn=np.linalg.norm(np.cross(b,bdot),axis=1)/np.maximum(norm*norm,1e-300)
        p=populations(w);pinit=np.zeros(12);pinit[:5]=.2
        residual=[]
        for j in range(16):
            residual.append([float(np.linalg.norm(expm(a[j]*(tau[j]*mult))@pinit-p[j],ord=1)) for mult in [.1,1,10]])
        def quant(x):return {str(q):float(np.percentile(x,q)) for q in [50,95,99,100]}
        records.append({'geometry':geometry,'controls_name':name,'n':len(r),'relaxation_us_percentiles':quant(tau*1e6),'tau_times_axis_turn_rate_percentiles':quant(tau*turn),'tau_over_waist_transit_percentiles':quant(tau*np.linalg.norm(v,axis=1)/c.waist),'frozen_16point_L1_errors_at_0p1_1_10_tau':residual, 'note':'spectral gap at initial synthetic cloud only; L1 decay with frozen rates does not validate dynamics along trajectories or near B=0 coherences'})
(OUT/'relaxation.json').write_text(json.dumps(records,indent=2))
for row in records:print(row['geometry'],row['controls_name'],row['relaxation_us_percentiles'],row['tau_times_axis_turn_rate_percentiles'])
