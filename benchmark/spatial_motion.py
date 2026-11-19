"""Fixed-basis coupled motion with collimated Gaussian beams and quadrupole B.
Deterministic mean-force local trajectories; no recoil or capture claim.
"""
import sys,json
from pathlib import Path
import numpy as np,pylcp
from scipy.integrate import solve_ivp
from scipy.sparse import csr_matrix
from pylcp.common import cart2spherical
sys.path.insert(0,str(Path(__file__).resolve().parent));sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from pylcp_multibeam import directions
from mot.physics import HBAR,K,MASS,GAMMA,MUB_HBAR
ETA=HBAR*K*K/(MASS*GAMMA)

def system(geom,phases,offsets,center,waist=.012,gradient=.12,frozen=False,base_detuning=-2.):
 dirs=directions(geom);h=np.repeat([1,1,-1] if len(dirs)==6 else [1,1,-1,1],2)
 center=np.asarray(center)*K;w=waist*K;G=gradient*MUB_HBAR/(GAMMA*K)*np.diag([1.,1.,-2.])
 # Dynamics uses local displacement from center; optical phases are defined
 # at that center so changing macroscopic position doesn't conflate phases.
 def field(x):return G@(center+np.asarray(x)) if not frozen else G@center
 # Public beams take absolute coordinates; offset initial phases to match local ones.
 lasers=pylcp.laserBeams([dict(kvec=d,pol=int(p),s=6/len(dirs),delta=float(df),phase=float(phi-d@center),wb=w) for d,p,df,phi in zip(dirs,h,offsets,phases)],beam_type=pylcp.gaussianBeam)
 Hg,mg=pylcp.hamiltonians.singleF(F=2,gF=.5,muB=1);He,me=pylcp.hamiltonians.singleF(F=3,gF=2/3,muB=1)
 ham=pylcp.hamiltonian(Hg,He-base_detuning*np.eye(7),mg,me,pylcp.hamiltonians.dqij_two_bare_hyperfine(2,3),mass=1/ETA)
 # Reference object uses absolute positions and exact quadrupole field.
 ob=pylcp.obe(lasers,pylcp.magField(lambda r:G@r),ham,include_mag_forces=True,transform_into_re_im=False,use_sparse_matrices=False)
 zero=np.zeros(3);Hbase=ham.return_full_H({'g->e':zero},cart2spherical(zero))
 Hunit=np.array([ham.return_full_H({'g->e':zero},cart2spherical(e))-Hbase for e in np.eye(3)])
 dHB=np.einsum('ja,jmn->amn',G,Hunit)
 # Unit-intensity plane-wave blocks with specified local phase.
 plane=pylcp.laserBeams([dict(kvec=d,pol=int(p),s=6/len(dirs),delta=float(df),phase=float(phi)) for d,p,df,phi in zip(dirs,h,offsets,phases)],beam_type=pylcp.infinitePlaneWaveBeam)
 blocks=np.array([ham.return_full_H({'g->e':b.electric_field(zero,0.)/2},cart2spherical(zero))[:5,5:] for b in plane.beam_vector])
 # Verified PyLCP convention from prior benchmarks: upper block linear in E.
 rr=np.array([.13,.27,-.19]);tt=.731;arg=np.asarray(offsets)*tt-dirs@rr
 exact=ham.return_full_H({'g->e':plane.total_electric_field(rr,tt)/2},cart2spherical(zero))[:5,5:]
 errors=[np.linalg.norm(np.einsum('b,bij->ij',np.exp(s*1j*arg),blocks)-exact) for s in [-1,1]];sign=[-1,1][int(np.argmin(errors))]
 decay=csr_matrix(ob.ev_mat['decay'])
 def hforce(x,t,rho):
  r=center+(np.zeros(3) if frozen else x)
  perp=r[None,:]-(dirs@r)[:,None]*dirs
  amplitude=np.exp(-np.sum(perp*perp,axis=1)/w**2)
  individual=blocks*(amplitude*np.exp(sign*1j*(np.asarray(offsets)*t-dirs@x)))[:,None,None]
  H=Hbase+np.einsum('j,jmn->mn',field(x),Hunit);upper=individual.sum(0);H[:5,5:]=upper;H[5:,:5]=upper.conj().T
  envelope_gradient=np.zeros_like(perp) if frozen else -2*perp/w**2
  derivative=np.einsum('ba,bij->aij',envelope_gradient-sign*1j*dirs,individual)
  optical=-2*np.real(np.einsum('ji,aij->a',rho[5:,:5],derivative))
  magnetic=np.zeros(3) if frozen else -np.real(np.einsum('ji,aij->a',rho,dHB))
  return H,optical+magnetic,optical,magnetic
 def rhs(t,y):
  rho=y[:144].reshape(12,12).T;x=y[147:].real;H,f,*_=hforce(x,t,rho)
  return np.r_[(-1j*(H@rho-rho@H)).T.ravel()+decay@y[:144],ETA*f,y[144:147]]
 return ob,hforce,rhs,center,G,w

def run(geom,center,frozen=False,T=2000.,rtol=2e-8,gradient=.12,waist=.012,sample_dt=2.,base_detuning=-2.,phase_seed=721):
 n=6 if geom=='axial6' else 8;p=np.random.default_rng(phase_seed).uniform(0,2*np.pi,n);df=np.repeat([0.,.12,.29,.47][:n//2],2)
 ob,hf,rhs,c,G,w=system(geom,p,df,center,waist,gradient,frozen,base_detuning)
 y0=np.r_[np.diag(np.r_[np.ones(5)/5,np.zeros(7)]).astype(complex).T.ravel(),[.05,0,0],[0,0,0]]
 ts=np.linspace(0,T,max(501,int(np.ceil(T/sample_dt))+1));sol=solve_ivp(rhs,(0,T),y0,t_eval=ts,method='DOP853',rtol=rtol,atol=rtol/100)
 if not sol.success:raise RuntimeError(sol.message)
 rho=sol.y[:144].T.reshape(-1,12,12).transpose(0,2,1);x=sol.y[147:].real.T;v=sol.y[144:147].real.T
 f=np.array([hf(a,float(t),r)[1] for a,t,r in zip(x,ts,rho)]);forceerr=[];herr=[]
 if not frozen:
  for j in np.linspace(0,len(ts)-1,10,dtype=int):
   forceerr.append(np.linalg.norm(f[j]-ob.force(c+x[j],float(ts[j]),rho[j].T)))
   E=ob.laserBeams['g->e'].total_electric_field(c+x[j],float(ts[j]))
   H=ob.hamiltonian.return_full_H({'g->e':E/2},cart2spherical(G@(c+x[j])))
   herr.append(np.linalg.norm(H-hf(x[j],float(ts[j]),rho[j])[0]))
 dirs=directions(geom);absr=c+x;transverse=np.sum(absr*absr,axis=1)[:,None]-(absr@dirs.T)**2
 attenuation=np.exp(-2*transverse/w**2)
 row=dict(geometry=geom,base_detuning=base_detuning,phase_seed=phase_seed,center_m=list(center),frozen=frozen,T=T,rtol=rtol,gradient_T_m=gradient,waist_m=waist,final_v=v[-1].tolist(),final_displacement_m=(x[-1]/K).tolist(),max_field_change=float(np.linalg.norm(x@G.T,axis=1).max()),max_relative_intensity_change=float(np.max(np.abs(attenuation/attenuation[0]-1))),initial_beam_intensity_factors=attenuation[0].tolist(),max_public_force_error=max(forceerr,default=None),max_public_hamiltonian_error=max(herr,default=None),max_trace_error=float(np.abs(np.trace(rho,axis1=1,axis2=2)-1).max()),min_eigenvalue=float(np.linalg.eigvalsh(rho).min()))
 from recoil_validity import diagnostic
 row['spontaneous_recoil_diagnostic']=diagnostic(ts,rho,ob.hamiltonian.d_q_bare['g->e'],ETA)
 return row,ts,x,v,f

if __name__=='__main__':
 out=Path('results/spatial');out.mkdir(exist_ok=True);rows=[]
 for geom in ['axial6','oblique8']:
  for name,c in [('origin',[0.,0,0]),('off_axis',[.006,.003,.002])]:
   for frozen in [False,True]:
    row,t,x,v,f=run(geom,c,frozen);row['position_name']=name;rows.append(row)
    np.savez_compressed(out/f'{geom}_{name}_{frozen}.npz',time=t,displacement=x,v=v,force=f)
    (out/'study.json').write_text(json.dumps(rows,indent=2));print(geom,name,frozen,row['final_v'],row['max_public_force_error'],flush=True)
