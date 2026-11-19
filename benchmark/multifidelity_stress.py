"""Separate condition-shift probe, not used to retune discovery selection."""
import json,time
from pathlib import Path
import numpy as np
from multifidelity_study import rate_model,CANDIDATES,T,V0
from spatial_motion import run
out=Path('results/multifidelity');c=[0.,0,0];rows=[]
(out/'stress_protocol.json').write_text(json.dumps(dict(center_m=c,seeds=list(range(1010,1014)),candidates=CANDIDATES,T=T,purpose='origin condition shift; not used for model selection'),indent=2))
(out/'stress_rate.json').write_text(json.dumps([rate_model(d,c) for d in CANDIDATES],indent=2))
for seed in range(1010,1014):
 for d in CANDIDATES:
  start=time.time();r,*_=run('axial6',c,frozen=True,T=T,base_detuning=d,phase_seed=seed);v=np.array(r['final_v']);r.update(energy_removal_fraction=float(1-v@v/(V0@V0)),wall_s=time.time()-start);rows.append(r)
  (out/'stress_coherent.json').write_text(json.dumps(rows,indent=2));print(seed,d,flush=True)
