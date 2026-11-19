"""Semiclassical mean-force density+Newton motion in a uniform local field.
No jumps/recoil diffusion, trapping envelope, or capture classification.
"""
import sys,json,time
from pathlib import Path
import numpy as np,pylcp
from scipy.integrate import solve_ivp
from scipy.sparse import csr_matrix
from pylcp.common import cart2spherical
sys.path.insert(0,str(Path(__file__).resolve().parent));sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from pylcp_multibeam import directions
from mot.physics import HBAR,K,MASS,GAMMA
ETA=HBAR*K*K/(MASS*GAMMA)

def system(geometry,phases,offsets,total_s=6.,eta=ETA):
 dirs=directions(geometry);h=np.repeat([1,1,-1] if len(dirs)==6 else [1,1,-1,1],2)
 lasers=pylcp.laserBeams([dict(kvec=d,pol=int(p),s=total_s/len(dirs),delta=float(df),phase=float(phi)) for d,p,df,phi in zip(dirs,h,offsets,phases)],beam_type=pylcp.infinitePlaneWaveBeam)
 Hg,mg=pylcp.hamiltonians.singleF(F=2,gF=.5,muB=1);He,me=pylcp.hamiltonians.singleF(F=3,gF=2/3,muB=1)
 ham=pylcp.hamiltonian(Hg,He+2*np.eye(7),mg,me,pylcp.hamiltonians.dqij_two_bare_hyperfine(2,3),mass=1/eta if eta else 1e100)
 field=np.array([.1,0,0]);ob=pylcp.obe(lasers,pylcp.magField(field),ham,include_mag_forces=False,transform_into_re_im=False,use_sparse_matrices=False)
 H0=ham.return_full_H({'g->e':np.zeros(3)},cart2spherical(field))
 blocks=np.array([ham.return_full_H({'g->e':b.electric_field(np.zeros(3),0.)/2},cart2spherical(field))[:5,5:] for b in lasers.beam_vector])
 # This upper Hamiltonian block is linear in E*: determine convention once.
 rr=np.array([.13,.27,-.19]);tt=.731;arg=np.asarray(offsets)*tt-dirs@rr
 exact=ham.return_full_H({'g->e':lasers.total_electric_field(rr,tt)/2},cart2spherical(field))[:5,5:]
 errors=[np.linalg.norm(np.einsum('b,bij->ij',np.exp(s*1j*arg),blocks)-exact) for s in [-1,1]];sign=[-1,1][int(np.argmin(errors))]
 decay=csr_matrix(ob.ev_mat['decay'])
 def hforce(r,t,rho):
  individual=blocks*np.exp(sign*1j*(np.asarray(offsets)*t-dirs@r))[:,None,None]
  H=H0.copy();upper=individual.sum(0);H[:5,5:]=upper;H[5:,:5]=upper.conj().T
  deriv=np.einsum('ba,bij->aij',-sign*1j*dirs,individual)
  force=-2*np.real(np.einsum('ji,aij->a',rho[5:,:5],deriv))
  return H,force
 def rhs(t,y):
  vec=y[:144];rho=vec.reshape(12,12).T;r=y[147:150].real
  H,f=hforce(r,t,rho);dv=eta*f
  return np.r_[(-1j*(H@rho-rho@H)).T.ravel()+decay@vec,dv,y[144:147]]
 return ob,hforce,rhs,min(errors)

def run(geom,phases,offsets,v0=(.05,0,0),T=2000.,rtol=2e-8,eta=ETA,total_s=6.,reference=False,sample_dt=2.):
 ob,hforce,rhs,herr=system(geom,phases,offsets,total_s,eta)
 y0=np.r_[np.diag(np.r_[np.ones(5)/5,np.zeros(7)]).astype(complex).T.ravel(),v0,[0,0,0]]
 ts=np.linspace(0,T,max(501,int(np.ceil(T/sample_dt))+1))
 sol=solve_ivp(rhs,(0,T),y0,method='DOP853',rtol=rtol,atol=rtol/100,t_eval=ts)
 if not sol.success:raise RuntimeError(sol.message)
 rho=sol.y[:144].T.reshape(-1,12,12).transpose(0,2,1);v=sol.y[144:147].real.T;r=sol.y[147:150].real.T
 forces=np.array([hforce(ri,float(t),rh)[1] for ri,t,rh in zip(r,ts,rho)])
 public_errors=[np.linalg.norm(hforce(r[j],float(ts[j]),rho[j])[1]-ob.force(r[j],float(ts[j]),rho[j].T)) for j in np.linspace(0,len(ts)-1,10,dtype=int)]
 # momentum balance uses independent force quadrature, not solver velocity.
 impulse=eta*np.trapezoid(forces,ts,axis=0)
 row=dict(geometry=geom,phases=list(phases),offsets=list(offsets),v0=list(v0),T=T,rtol=rtol,sample_dt=sample_dt,eta=eta,total_s=total_s,final_v=v[-1].tolist(),final_r=r[-1].tolist(),max_speed_change=float(np.linalg.norm(v-np.asarray(v0),axis=1).max()),max_hamiltonian_error=float(herr),max_public_force_error=float(max(public_errors)),momentum_quadrature_error=float(np.linalg.norm(v[-1]-v[0]-impulse)),max_trace_error=float(np.abs(np.trace(rho,axis1=1,axis2=2)-1).max()),min_eigenvalue=float(np.linalg.eigvalsh(rho).min()),max_position_imaginary=float(np.abs(sol.y[144:].imag).max()))
 if reference:
  ob.set_initial_rho_from_populations(np.r_[np.ones(5)/5,np.zeros(7)]);ob.set_initial_velocity(np.asarray(v0));ob.set_initial_position(np.zeros(3))
  ob.evolve_motion([0,T],random_recoil=False,t_eval=[T],method='DOP853',rtol=rtol,atol=rtol/100)
  row.update(public_density_error=float(np.linalg.norm(ob.sol.rho[:,:,-1]-rho[-1].T)),public_velocity_error=float(np.linalg.norm(ob.sol.v[:,-1]-v[-1])),public_position_error=float(np.linalg.norm(ob.sol.r[:,-1]-r[-1])))
 return row,ts,r,v,forces

if __name__=='__main__':
 out=Path('results/coupled');out.mkdir(exist_ok=True);rows=[];start=time.time()
 for geom,n in [('axial6',6),('oblique8',8)]:
  for seed in [721,722]:
   p=np.random.default_rng(seed).uniform(0,2*np.pi,n)
   for arch,df in [('common',np.zeros(n)),('offset_pairs',np.repeat([0.,.12,.29,.47][:n//2],2))]:
    row,t,r,v,f=run(geom,p,df);row.update(seed=seed,architecture=arch);rows.append(row)
    np.savez_compressed(out/f'{geom}_{seed}_{arch}.npz',time=t,r=r,v=v,force=f)
    (out/'study.json').write_text(json.dumps(rows,indent=2));print(geom,seed,arch,row['final_v'],row['max_public_force_error'],flush=True)
 (out/'runtime.json').write_text(json.dumps({'runtime_s':time.time()-start}))
