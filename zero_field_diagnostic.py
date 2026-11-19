"""Reproduce basis sensitivity in local-axis incoherent model at vanishing B."""
import json
from pathlib import Path
import numpy as np
from mot.physics import Controls,force
rows=[]
dirs=np.array([[1,0,0],[0,1,0],[0,0,1],[1,1,1],[1,-1,2]],float);dirs/=np.linalg.norm(dirs,axis=1)[:,None]
for geometry,folder in [('axial6',Path('results')),('oblique8',Path('results/eight_beam'))]:
    data=json.loads((folder/'heldout.json').read_text())
    for name in ['standard','gp_ucb_11']:
        c=Controls(**data[name]['controls']);v=np.tile([.1,.2,.3],(len(dirs),1))
        f=force(dirs*1e-9,v,c,False);f0=force(np.zeros((1,3)),v[:1],c,False)[0]
        rows.append({'geometry':geometry,'name':name,'radius_m':1e-9,'velocity':[.1,.2,.3],'approach_directions':dirs.tolist(),'force':f.tolist(),'origin_force':f0.tolist(),'max_relative_difference_from_origin':float(np.max(np.linalg.norm(f-f0,axis=1))/np.linalg.norm(f0))})
Path('results/v2/zero_field_axis_sensitivity.json').write_text(json.dumps(rows,indent=2))
