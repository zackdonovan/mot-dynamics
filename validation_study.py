"""Locked v2 follow-up: no optimization or finalist selection using these results."""
import json,time
from pathlib import Path
import numpy as np
from mot.physics import *
OUT=Path('results/v2');OUT.mkdir(exist_ok=True)
protocol={'n':512,'cloud_seed':20261118,'velocity_sigma':[4.,8.,12.], 'position_sigma':.006,'snapshots_ms':[25,75], 'dt':2.5e-5,'controls':['standard','sobol_11','gp_ucb_11'],'geometries':['axial6','oblique8'],'recoil_seed':17890,'note':'fixed prior seed11 finalists; no reoptimization; 3 sources x 6 finalists x 2 dynamics modes; paired atoms and recoil RNG stream across comparisons; single recoil seed per condition, not replicated-seed error bars'}
(OUT/'protocol.json').write_text(json.dumps(protocol,indent=2))
rows=[];start=time.time();rng=np.random.default_rng(protocol['cloud_seed']);r0=rng.normal(0,.006,(512,3));unitv=rng.normal(size=(512,3))
def summarize(o,first=None):
    ret=o['retained'];rad=np.linalg.norm(o['r'],axis=1);speed=np.linalg.norm(o['v'],axis=1)
    k=int(ret.sum());n=len(ret);p=k/n;z=1.959963984540054;d=1+z*z/n
    ci=[(p+z*z/(2*n)-z*np.sqrt(p*(1-p)/n+z*z/(4*n*n)))/d,(p+z*z/(2*n)+z*np.sqrt(p*(1-p)/n+z*z/(4*n*n)))/d]
    out={'count':k,'fraction':p,'wilson95':ci,'radius_median_mm':float(np.median(rad)*1000),'speed_median_mps':float(np.median(speed)), 'threshold_counts':{f'r{radius}_v{vel}':int(((rad<radius/1000)&(speed<vel)).sum()) for radius in [1,2,4] for vel in [.25,.5,1.]}}
    if first is not None:out.update({'lost_from_25ms':int((first&~ret).sum()),'retained_at_both':int((first&ret).sum()),'entered_since_25ms':int((~first&ret).sum())})
    return out
for geometry,folder in [('axial6',Path('results')),('oblique8',Path('results/eight_beam'))]:
    data=json.loads((folder/'heldout.json').read_text())
    for name in protocol['controls']:
        c=Controls(**data[name]['controls'])
        for sigma in protocol['velocity_sigma']:
            for recoil in [False,True]:
                a=simulate(c,r0,unitv*sigma,dt=protocol['dt'],duration=.025,recoil=recoil,seed=protocol['recoil_seed'])
                # Fresh RNG seed for the subsequent independent increments.
                b=simulate(c,a['r'],a['v'],dt=protocol['dt'],duration=.050,recoil=recoil,seed=protocol['recoil_seed']+1)
                row={'geometry':geometry,'controls_name':name,'velocity_sigma':sigma,'recoil':recoil,'25ms':summarize(a),'75ms':summarize(b,a['retained'])}
                np.savez_compressed(OUT/f'{geometry}_{name}_s{sigma}_recoil{recoil}.npz',r25=a['r'],v25=a['v'],r75=b['r'],v75=b['v'],mask25=a['retained'],mask75=b['retained'])
                rows.append(row);(OUT/'followup.json').write_text(json.dumps(rows,indent=2));print(len(rows),geometry,name,sigma,recoil,row['25ms']['fraction'],row['75ms']['fraction'],round(time.time()-start,1),flush=True)
(OUT/'runtime.json').write_text(json.dumps({'seconds':time.time()-start}))
