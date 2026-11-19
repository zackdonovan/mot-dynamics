"""Matched F=0->F'=1 frozen-atom benchmark. Not the production Rb F=2->3 model."""
import json,time,sys
from pathlib import Path
import numpy as np
import pylcp
from pylcp.common import cart2spherical
OUT=Path('results/pylcp');OUT.mkdir(exist_ok=True)

def analytical_rate(s,delta,beta,pair=True):
    # For +beam sigma+ couples m=+1; -beam sigma+ couples m=-1.
    wplus=s/2/(1+4*(delta-beta)**2)
    wminus=s/2/(1+4*(delta+beta)**2) if pair else 0.
    pg=1/(1+wplus/(1+wplus)+wminus/(1+wminus))
    peplus=pg*wplus/(1+wplus);peminus=pg*wminus/(1+wminus)
    return wplus*(pg-peplus)-wminus*(pg-peminus)

def case(s,delta,beta,axis=2,pair=True,T=100.):
    direction=np.eye(3)[axis]
    beams=[dict(kvec=direction,pol=1,delta=0.,s=s)]
    if pair:beams.append(dict(kvec=-direction,pol=1,delta=0.,s=s))
    lasers=pylcp.laserBeams(beams,beam_type=pylcp.infinitePlaneWaveBeam)
    Hg,mg=pylcp.hamiltonians.singleF(F=0,gF=0,muB=1)
    He,me=pylcp.hamiltonians.singleF(F=1,gF=1,muB=1)
    dq=pylcp.hamiltonians.dqij_two_bare_hyperfine(0,1)
    ham=pylcp.hamiltonian(Hg,He-delta*np.eye(3),mg,me,dq)
    field=pylcp.magField(direction*beta)
    re=pylcp.rateeq(lasers,field,ham,include_mag_forces=False)
    fre=re.find_equilibrium_force()
    ob=pylcp.obe(lasers,field,ham,include_mag_forces=False,transform_into_re_im=False)
    ob.set_initial_rho_from_populations(np.array([1.,0,0,0]))
    ob.evolve_density([0,T],t_eval=np.array([T-10,T]),rtol=1e-9,atol=1e-11)
    fo=ob.force(ob.sol.r,ob.sol.t,ob.sol.rho)
    rho=ob.sol.rho[:,:,-1]
    Eq={'g->e':lasers.total_electric_field(np.zeros(3),T)}
    # PyLCP's OBE uses gamma*E/4 (gamma=1), whereas return_full_H
    # uses E/2. Scale this independently assembled Hamiltonian accordingly.
    H=ham.return_full_H({key:value/2 for key,value in Eq.items()},cart2spherical(direction*beta))
    # PyLCP indexes vectors column-major but sol.rho is reshaped row-major.
    # Transpose for the conventional density-matrix commutator.
    physical_rho=rho.T
    residual=-1j*(H@physical_rho-physical_rho@H)+(ob.ev_mat['decay']@rho.ravel()).reshape((4,4)).T
    row=dict(s=s,delta=delta,beta=beta,axis=axis,pair=pair,T=T,rate_analytic=float(analytical_rate(s,delta,beta,pair)),pylcp_rate_force=fre.tolist(),pylcp_obe_force=fo[:,-1].tolist(),tail_force_difference=float(np.linalg.norm(fo[:,-1]-fo[:,-2])),density_trace_error=float(abs(np.trace(rho)-1)),density_min_eigenvalue=float(np.linalg.eigvalsh(rho).min()),steady_residual_norm=float(np.linalg.norm(residual)))
    if not pair:row['two_level_exact']=float((s/2)/(1+s+4*(delta-beta)**2))
    return row

if __name__=='__main__':
    start=time.time();rows=[]
    for s in [.05,1.]:
        for beta in [-2.,-1.,-.1,0.,.1,1.,2.]:
            for axis in [0,1,2]:
                rows.append(case(s,-2.,beta,axis=axis));print(len(rows),'pair',s,beta,axis,flush=True)
                (OUT/'toy.json').write_text(json.dumps(rows,indent=2))
    for s in [.01,.1,1.]:
        for beta in [-3.,-2.,0.,1.,3.]:
            rows.append(case(s,-2.,beta,pair=False));print(len(rows),'single',s,beta,flush=True)
            (OUT/'toy.json').write_text(json.dumps(rows,indent=2))
    # Repeat selected high-saturation, near-zero points at twice integration time.
    checks=[case(1.,-2.,b,axis=a,T=200.) for b,a in [(.1,2),(1.,0),(0.,1)]]
    (OUT/'toy_doubletime.json').write_text(json.dumps(checks,indent=2))
    metrics={'n_cases':len(rows),'runtime_s':time.time()-start,'max_rate_analytic_error':max(abs(r['rate_analytic']-r['pylcp_rate_force'][r['axis']]) for r in rows),'max_single_obe_exact_error':max(abs(r['two_level_exact']-r['pylcp_obe_force'][r['axis']]) for r in rows if not r['pair']),'max_steady_residual':max(r['steady_residual_norm'] for r in rows),'max_tail_force_difference':max(r['tail_force_difference'] for r in rows),'max_trace_error':max(r['density_trace_error'] for r in rows),'min_density_eigenvalue':min(r['density_min_eigenvalue'] for r in rows)}
    (OUT/'toy_metrics.json').write_text(json.dumps(metrics,indent=2));print(metrics,flush=True)
