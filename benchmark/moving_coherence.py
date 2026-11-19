"""Prescribed straight paths through phase-coherent monochromatic/offset beams.
Density dynamics is integrated, not locally equilibrated. No force feedback.
"""
import sys,json,time
from pathlib import Path
import numpy as np,pylcp
from scipy.integrate import solve_ivp
from scipy.sparse import csr_matrix
from pylcp.common import cart2spherical
sys.path.insert(0,str(Path(__file__).resolve().parent))
from pylcp_multibeam import directions

def setup(geometry,phases,offsets,velocity,field=(.1,0,0)):
    dirs=directions(geometry);h=np.repeat([1,1,-1] if len(dirs)==6 else [1,1,-1,1],2)
    lasers=pylcp.laserBeams([dict(kvec=d,pol=int(p),s=6/len(dirs),delta=float(df),phase=float(phi)) for d,p,df,phi in zip(dirs,h,offsets,phases)],beam_type=pylcp.infinitePlaneWaveBeam)
    Hg,mg=pylcp.hamiltonians.singleF(F=2,gF=.5,muB=1);He,me=pylcp.hamiltonians.singleF(F=3,gF=2/3,muB=1)
    ham=pylcp.hamiltonian(Hg,He+2*np.eye(7),mg,me,pylcp.hamiltonians.dqij_two_bare_hyperfine(2,3))
    ob=pylcp.obe(lasers,pylcp.magField(np.array(field)),ham,v0=np.asarray(velocity),include_mag_forces=False,transform_into_re_im=False,use_sparse_matrices=False)
    H0=ham.return_full_H({'g->e':np.zeros(3)},cart2spherical(np.asarray(field)))
    # Cache only ground/excited upper blocks; lower blocks are Hermitian adjoints.
    blocks=[]
    for b in lasers.beam_vector:
        H=ham.return_full_H({'g->e':b.electric_field(np.zeros(3),0.)/2},cart2spherical(np.asarray(field)))
        blocks.append(H[:5,5:])
    # The H upper block is proportional to the complex field conjugate.
    # Determine its time sign by direct comparison rather than assuming it.
    omega=np.asarray(offsets)-dirs@velocity
    test_t=.731
    exact=ham.return_full_H({'g->e':lasers.total_electric_field(np.asarray(velocity)*test_t,test_t)/2},cart2spherical(np.asarray(field)))[:5,5:]
    errors=[np.linalg.norm(np.einsum('b,bij->ij',np.exp(sign*1j*omega*test_t),blocks)-exact) for sign in [-1,1]]
    sign=[-1,1][int(np.argmin(errors))]
    decay=csr_matrix(ob.ev_mat['decay'])
    def rhs(t,vec):
        rho=vec.reshape(12,12).T
        H=H0.copy();upper=np.einsum('b,bij->ij',np.exp(sign*1j*omega*t),blocks)
        H[:5,5:]=upper;H[5:,:5]=upper.conj().T
        return (-1j*(H@rho-rho@H)).T.ravel()+decay@vec
    return ob,lasers,ham,rhs,float(min(errors))

def run(geom,phases,offsets,velocity,T=2000.,rtol=2e-8,validate=False,sample_dt=2.):
    velocity=np.asarray(velocity,float);ob,lasers,ham,rhs,herr=setup(geom,phases,offsets,velocity)
    y0=np.diag(np.r_[np.ones(5)/5,np.zeros(7)]).astype(complex).T.ravel()
    ts=np.linspace(0,T,max(1001,int(np.ceil(T/sample_dt))+1))
    sol=solve_ivp(rhs,(0,T),y0,method='DOP853',rtol=rtol,atol=rtol/100,t_eval=ts)
    if not sol.success:raise RuntimeError(sol.message)
    densities=sol.y.T.reshape(-1,12,12)
    force=np.array([ob.force(velocity*t,float(t),rho) for t,rho in zip(ts,densities)])
    # Density arrays above use PyLCP's transposed convention.
    physical=densities.transpose(0,2,1)
    tail=ts>=T/2;last=ts>=3*T/4;prev=(ts>=T/2)&(ts<=3*T/4)
    average=lambda mask:np.trapezoid(force[mask],ts[mask],axis=0)/(ts[mask][-1]-ts[mask][0])
    # Rate comparator uses actual detuning -2+offset and matching Doppler shift.
    dirs=directions(geom);h=np.repeat([1,1,-1] if len(dirs)==6 else [1,1,-1,1],2)
    lb=pylcp.laserBeams([dict(kvec=d,pol=int(p),s=6/len(dirs),delta=-2+float(df)) for d,p,df in zip(dirs,h,offsets)],beam_type=pylcp.infinitePlaneWaveBeam)
    Hg,mg=pylcp.hamiltonians.singleF(F=2,gF=.5,muB=1);He,me=pylcp.hamiltonians.singleF(F=3,gF=2/3,muB=1)
    hr=pylcp.hamiltonian(Hg,He,mg,me,pylcp.hamiltonians.dqij_two_bare_hyperfine(2,3))
    re=pylcp.rateeq(lb,pylcp.magField(np.array([.1,0,0])),hr,v0=velocity,include_mag_forces=False)
    row=dict(geometry=geom,velocity=velocity.tolist(),offsets=list(offsets),phases=list(phases),T=T,rtol=rtol,sample_dt=sample_dt,nfev=sol.nfev,mean_tail_force=average(tail).tolist(),previous_quarter_force=average(prev).tolist(),last_quarter_force=average(last).tolist(),rate_force=re.find_equilibrium_force().tolist(),hamiltonian_reconstruction_error=herr,max_trace_error=float(np.max(np.abs(np.trace(physical,axis1=1,axis2=2)-1))),max_hermiticity_error=float(np.linalg.norm(physical-physical.conj().transpose(0,2,1),axis=(1,2)).max()),min_density_eigenvalue=float(np.linalg.eigvalsh(physical).min()))
    if validate:
        ob.set_initial_rho_from_populations(np.r_[np.ones(5)/5,np.zeros(7)])
        ob.evolve_density([0,min(T,30.)],t_eval=[min(T,30.)],method='DOP853',rtol=rtol,atol=rtol/100)
        direct=solve_ivp(rhs,(0,min(T,30.)),y0,method='DOP853',rtol=rtol,atol=rtol/100)
        row['pylcp_time_density_error']=float(np.linalg.norm(ob.sol.rho[:,:,-1]-direct.y[:,-1].reshape(12,12)))
    return row,ts,force

if __name__=='__main__':
    out=Path('results/moving');out.mkdir(exist_ok=True);rows=[];start=time.time()
    for geom,n in [('axial6',6),('oblique8',8)]:
      for seed in [721,722]:
        phases=np.random.default_rng(seed).uniform(0,2*np.pi,n)
        for architecture,offsets in [('common',np.zeros(n)),('offset_pairs',np.repeat([0.,.12,.29,.47][:n//2],2))]:
          for speed in [0.,.05,.2]:
            row,ts,f=run(geom,phases,offsets,[speed,0,0],validate=(seed==721 and speed==.05))
            row.update(seed=seed,architecture=architecture);rows.append(row)
            np.savez_compressed(out/f'{geom}_{seed}_{architecture}_{speed}.npz',time=ts,force=f)
            (out/'study.json').write_text(json.dumps(rows,indent=2));print(geom,seed,architecture,speed,row['mean_tail_force'],flush=True)
    (out/'runtime.json').write_text(json.dumps({'runtime_s':time.time()-start}))
