"""Frozen-atom, monochromatic six/eight-beam coherent phase ensemble.

Pure circular helicities, equal detunings, no motion, no phase diffusion.
Phase ensemble is quenched disorder: average AFTER solving each steady state.
It is not an incoherent laser model and not a capture calculation.
"""
import json,time,sys
from pathlib import Path
import numpy as np,pylcp
from pylcp.common import cart2spherical
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from mot.physics import Controls,beam_directions

def directions(geometry):
    n=3 if geometry=='axial6' else 4
    return beam_directions(Controls(geometry=geometry,detuning=(-2.,)*n,purity=(1.,)*n,power=(1.,)*n))

def solve(geometry,field,phases,total_s=6.,evolve=False):
    dirs=directions(geometry);nbeam=len(dirs)
    if len(phases)!=nbeam: raise ValueError('one optical phase per beam required')
    h=np.repeat([1,1,-1] if nbeam==6 else [1,1,-1,1],2)
    lasers=pylcp.laserBeams([dict(kvec=d,pol=int(p),s=total_s/nbeam,delta=0.,phase=float(phase)) for d,p,phase in zip(dirs,h,phases)],beam_type=pylcp.infinitePlaneWaveBeam)
    Hg,mg=pylcp.hamiltonians.singleF(F=2,gF=.5,muB=1)
    He,me=pylcp.hamiltonians.singleF(F=3,gF=2/3,muB=1)
    ham=pylcp.hamiltonian(Hg,He+2*np.eye(7),mg,me,pylcp.hamiltonians.dqij_two_bare_hyperfine(2,3))
    ob=pylcp.obe(lasers,pylcp.magField(np.array(field)),ham,include_mag_forces=False,transform_into_re_im=False,use_sparse_matrices=False)
    E=lasers.total_electric_field(np.zeros(3),0.)
    H=ham.return_full_H({'g->e':E/2},cart2spherical(np.array(field)))
    L=-1j*(np.kron(np.eye(12),H)-np.kron(H.T,np.eye(12)))+ob.ev_mat['decay']
    A=np.vstack([L,np.eye(12).ravel()]);rhs=np.zeros(145,complex);rhs[-1]=1
    vec=np.linalg.lstsq(A,rhs,rcond=None)[0];rho=vec.reshape(12,12).T
    row=dict(phases=np.asarray(phases).tolist(),force=ob.force(np.zeros(3),0.,rho.T).tolist(),excited_population=float(np.trace(rho[5:,5:]).real),residual=float(np.linalg.norm(L@vec)),trace_error=float(abs(np.trace(rho)-1)),hermiticity_error=float(np.linalg.norm(rho-rho.conj().T)),min_eigenvalue=float(np.linalg.eigvalsh(rho).min()))
    if evolve:
        ob.set_initial_rho_from_populations(np.r_[np.ones(5)/5,np.zeros(7)])
        ob.evolve_density([0,12000.],t_eval=[12000.],method='BDF',rtol=1e-9,atol=1e-11)
        row.update(time_success=bool(ob.sol.success),time_density_error=float(np.linalg.norm(ob.sol.rho[:,:,-1]-rho.T)))
    return row

def rate(geometry,field,total_s=6.):
    dirs=directions(geometry);h=np.repeat([1,1,-1] if len(dirs)==6 else [1,1,-1,1],2)
    lasers=pylcp.laserBeams([dict(kvec=d,pol=int(p),s=total_s/len(dirs),delta=-2.) for d,p in zip(dirs,h)],beam_type=pylcp.infinitePlaneWaveBeam)
    Hg,mg=pylcp.hamiltonians.singleF(F=2,gF=.5,muB=1);He,me=pylcp.hamiltonians.singleF(F=3,gF=2/3,muB=1)
    ham=pylcp.hamiltonian(Hg,He,mg,me,pylcp.hamiltonians.dqij_two_bare_hyperfine(2,3))
    re=pylcp.rateeq(lasers,pylcp.magField(np.array(field)),ham,include_mag_forces=False)
    f=re.find_equilibrium_force();return f.tolist()

if __name__=='__main__':
    out=Path('results/multibeam');out.mkdir(exist_ok=True)
    rng=np.random.default_rng(271828);records=json.loads((out/'phase_ensemble.json').read_text()) if '--resume' in sys.argv else [];checks=[];start=time.time()
    for geom,n in [('axial6',6),('oblique8',8)]:
        phases=rng.uniform(0,2*np.pi,(128,n));phases[:,0]=0 # fix unobservable common phase
        for field in [[0.,0,0],[1e-8,0,0],[.1,0,0],[0,0,.1]]:
            if any(x['geometry']==geom and x['field']==field for x in records): continue
            rows=[solve(geom,field,p) for p in phases]
            f=np.array([r['force'] for r in rows]);pop=np.array([r['excited_population'] for r in rows])
            records.append(dict(geometry=geom,field=field,rate_force=rate(geom,field),mean_force=f.mean(0).tolist(),force_standard_error=(f.std(0,ddof=1)/np.sqrt(len(f))).tolist(),force_std=f.std(0,ddof=1).tolist(),mean_excited_population=float(pop.mean()),force_mean_first64=f[:64].mean(0).tolist(),rows=rows))
            print(geom,field,records[-1]['mean_force'],flush=True)
            (out/'phase_ensemble.json').write_text(json.dumps(records,indent=2))
        base=solve(geom,[.1,0,0],phases[0],evolve=True);shifted=solve(geom,[.1,0,0],phases[0]+.731)
        checks.append(dict(geometry=geom,global_phase_force_error=float(np.linalg.norm(np.array(base['force'])-shifted['force'])),global_phase_population_error=abs(base['excited_population']-shifted['excited_population']),time_density_error=base['time_density_error'],time_success=base['time_success']))
    (out/'checks.json').write_text(json.dumps(dict(checks=checks,runtime_s=time.time()-start),indent=2))
