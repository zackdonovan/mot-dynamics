"""Refine final comparisons; fixed seed11, not chosen using holdout performance."""
import json
from pathlib import Path
import numpy as np
from mot.physics import *
cloud=sample_atoms(4096,88291)
for folder in [Path('results'),Path('results/eight_beam')]:
    original=json.loads((folder/'heldout.json').read_text());out={}
    for name in ['standard','gp_ucb_11']:
        cfg=Controls(**original[name]['controls']);coarse=simulate(cfg,*cloud);fine=simulate(cfg,*cloud,dt=2.5e-5)
        d=np.linalg.norm(coarse['r']-fine['r'],axis=1)
        retained=fine['retained'];out[name]={'coarse_fraction':float(coarse['retained'].mean()),'fine_fraction':float(retained.mean()),'changed_labels':int((coarse['retained']!=retained).sum()),'position_error_p95_all_mm':float(np.percentile(d,95)*1e3),'position_error_p95_fine_retained_mm':float(np.percentile(d[retained],95)*1e3)}
        print(folder,name,out[name],flush=True)
    (folder/'timestep_confirm.json').write_text(json.dumps(out,indent=2))
