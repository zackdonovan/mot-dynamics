"""Fixed-basis coherent check on an F=2->3 single circular beam, stationary atom.

This isolates the basis pathology. It is not the full six/eight-beam geometry,
not a validation of capture, and contains no incoherent helicity mixtures.
"""
import json,time
from pathlib import Path
import numpy as np,pylcp
from pylcp.common import cart2spherical
OUT=Path('results/pylcp');OUT.mkdir(exist_ok=True)

def case(field,evolve=False):
    s=1.;delta=-2.;field=np.array(field,float)
    lasers=pylcp.laserBeams([dict(kvec=np.array([0.,0,1.]),pol=1,delta=0.,s=s)],beam_type=pylcp.infinitePlaneWaveBeam)
    Hg,mg=pylcp.hamiltonians.singleF(F=2,gF=.5,muB=1);He,me=pylcp.hamiltonians.singleF(F=3,gF=2/3,muB=1)
    ham=pylcp.hamiltonian(Hg,He-delta*np.eye(7),mg,me,pylcp.hamiltonians.dqij_two_bare_hyperfine(2,3))
    re=pylcp.rateeq(lasers,pylcp.magField(field),ham,include_mag_forces=False)
    fre=re.find_equilibrium_force()
    ob=pylcp.obe(lasers,pylcp.magField(field),ham,include_mag_forces=False,transform_into_re_im=False,use_sparse_matrices=False)
    E=lasers.total_electric_field(np.zeros(3),0.)
    # Units gamma=1; OBE coupling gamma E/4; public H builder uses E/2.
    H=ham.return_full_H({'g->e':E/2},cart2spherical(field))
    n=12;L=-1j*(np.kron(np.eye(n),H)-np.kron(H.T,np.eye(n)))+ob.ev_mat['decay']
    # Column-major density vector consistent with PyLCP's density_index.
    trace=np.eye(n).ravel();A=np.vstack([L,trace]);rhs=np.zeros(n*n+1,complex);rhs[-1]=1.
    vec=np.linalg.lstsq(A,rhs,rcond=None)[0];rho=vec.reshape((n,n)).T
    fobe=ob.force(np.zeros(3),0.,rho.T)
    row={'field':field.tolist(),'pylcp_rate_force':fre.tolist(),'fixed_basis_obe_force':fobe.tolist(),'steady_residual':float(np.linalg.norm(L@vec)),'trace_error':float(abs(np.trace(rho)-1)),'hermiticity_error':float(np.linalg.norm(rho-rho.conj().T)),'min_density_eigenvalue':float(np.linalg.eigvalsh(rho).min())}
    if evolve:
        ob.set_initial_rho_from_populations(np.r_[np.ones(5)/5,np.zeros(7)])
        ob.evolve_density([0,12000.],t_eval=[12000.],method='BDF',rtol=1e-9,atol=1e-11)
        time_rho=ob.sol.rho[:,:,-1]
        row['time_evolution_force']=ob.force(np.zeros(3),12000.,time_rho).tolist()
        row['time_evolution_vs_direct_rho_norm']=float(np.linalg.norm(time_rho-rho.T))
        row['time_evolution_success']=bool(ob.sol.success)
    return row
if __name__=='__main__':
    start=time.time();rows=[case([0.,0,0.],evolve=True)]
    dirs=[np.array(d,float)/np.linalg.norm(d) for d in [[1,0,0],[0,1,0],[0,0,1],[1,1,1],[1,-1,2]]]
    for magnitude in [1e-8,1e-4,.01,.1,1.]:
        for i,d in enumerate(dirs):
            rows.append(case(d*magnitude,evolve=(magnitude==.1 and i==0)))
            (OUT/'zero_field_coherent.json').write_text(json.dumps(rows,indent=2));print(len(rows),magnitude,i,rows[-1]['pylcp_rate_force'][2],rows[-1]['fixed_basis_obe_force'][2],flush=True)
    print('runtime_s',time.time()-start,flush=True)
