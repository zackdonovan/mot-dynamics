"""Predeclared four-candidate pilot; rank conditional deterministic models.
Independent discovery/held-out optical-phase batches; no capture objective.
"""
import json,time,sys
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from spatial_motion import run,ETA
from mot.physics import Controls,force,HBAR,K,GAMMA
OUT=Path('results/multifidelity');OUT.mkdir(exist_ok=True)
CANDIDATES=[-1.,-2.,-3.,-4.];CENTER=[.006,.003,.002];T=1000.;V0=np.array([.05,0,0])

def rate_model(detuning,center=CENTER):
 c=Controls(detuning=tuple(detuning+np.array([0.,.12,.29])))
 def rhs(t,v):return ETA*force(np.array(center),v*GAMMA/K,c,gravity=False)[0]/(HBAR*K*GAMMA)
 start=time.time();sol=solve_ivp(rhs,(0,T),V0,method='DOP853',rtol=2e-9,atol=2e-11)
 if not sol.success:raise RuntimeError(sol.message)
 v=sol.y[:,-1];return dict(detuning=detuning,final_v=v.tolist(),energy_removal_fraction=float(1-v@v/(V0@V0)),wall_s=time.time()-start)

if __name__=='__main__':
 # The batches and candidates are fixed before inspecting results.
 (OUT/'protocol.json').write_text(json.dumps(dict(candidates=CANDIDATES,center_m=CENTER,T=T,initial_v=V0.tolist(),discovery_seeds=list(range(810,818)),heldout_seeds=list(range(910,918)),objective='1 - final squared speed / initial squared speed; larger is better',scope='frozen local profiles, deterministic mean force, six beams, fixed power and pure helicities'),indent=2))
 rates=[rate_model(d) for d in CANDIDATES];(OUT/'rate.json').write_text(json.dumps(rates,indent=2));rows=[]
 for batch,seeds in [('discovery',range(810,818)),('heldout',range(910,918))]:
  for seed in seeds:
   for d in CANDIDATES:
    start=time.time();r,*_=run('axial6',CENTER,frozen=True,T=T,base_detuning=d,phase_seed=seed)
    v=np.array(r['final_v']);r.update(batch=batch,energy_removal_fraction=float(1-v@v/(V0@V0)),wall_s=time.time()-start);rows.append(r)
    (OUT/'coherent.json').write_text(json.dumps(rows,indent=2));print(batch,seed,d,r['energy_removal_fraction'],flush=True)
