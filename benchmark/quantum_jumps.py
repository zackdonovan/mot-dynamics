"""Photon-resolved fixed-atom quantum trajectories, exact event root timing.
No external force feedback or absorption/force-noise model.
"""
import sys,json,time
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
sys.path.insert(0,str(Path(__file__).resolve().parent))
from coupled_motion import system,ETA
from recoil_validity import cartesian_dipoles,emission_tensor

def setup():
 phases=np.random.default_rng(721).uniform(0,2*np.pi,6)
 ob,hf,*_=system('axial6',phases,np.zeros(6))
 H=hf(np.zeros(3),0.,np.eye(12)/12)[0]
 dq=ob.hamiltonian.d_q_bare['g->e'];d=cartesian_dipoles(dq)
 closure=sum(a.conj().T@a for a in d)
 A=-1j*H-closure/2
 ev,U=np.linalg.eig(A);Ui=np.linalg.inv(U)
 errors=[np.linalg.norm((U*np.exp(ev*t))@Ui-expm(A*t)) for t in [0.,.1,1.,10.,100.,200.]]
 def evolve(psi,t):return U@(np.exp(ev*t)*(Ui@psi))
 return ob,H,d,closure,evolve,dict(eigenvector_condition=float(np.linalg.cond(U)),max_propagator_error=float(max(errors)),closure_error=float(np.linalg.norm(closure-np.diag(np.r_[np.zeros(5),np.ones(7)]))))

def emit(psi,d,rng):
 """Sample direction and unresolved transverse polarization via jump operators."""
 amp=d@psi;rate=float(np.linalg.norm(amp)**2)
 if rate<=0:raise ValueError('cannot emit from a zero-excited state')
 while True:
  n=rng.normal(size=3);n/=np.linalg.norm(n)
  weight=rate-np.linalg.norm(n@amp)**2
  if rng.random()*rate<weight:break
 axis=np.array([0.,0,1.]) if abs(n[2])<.9 else np.array([1.,0,0])
 e1=np.cross(n,axis);e1/=np.linalg.norm(e1);e2=np.cross(n,e1)
 a=e1@amp;b=e2@amp;pa=float(np.linalg.norm(a)**2);pb=float(np.linalg.norm(b)**2)
 state=a if rng.random()*(pa+pb)<pa else b
 return state/np.linalg.norm(state),n

def trajectory(evolve,d,T,rng,initial_state=None):
 if initial_state is None:
  psi=np.zeros(12,complex);psi[rng.integers(0,5)]=1
 else:psi=np.array(initial_state,complex)/np.linalg.norm(initial_state)
 t=0.;count=0;recoil=np.zeros(3);events=[]
 while t<T:
  threshold=rng.random();remaining=T-t;end=evolve(psi,remaining)
  if np.vdot(end,end).real>threshold:
   psi=end/np.linalg.norm(end);break
  def loss(dt):
   q=evolve(psi,dt);return np.vdot(q,q).real-threshold
  wait=brentq(loss,0.,remaining,xtol=1e-10,rtol=1e-12)
  psi=evolve(psi,wait);psi/=np.linalg.norm(psi);psi,n=emit(psi,d,rng)
  t+=wait;count+=1;recoil-=ETA*n;events.append(t)
 return psi,count,recoil,events

def reference(ob,H,d,T):
 n=12;decay=ob.ev_mat['decay'];rho0=np.diag(np.r_[np.ones(5)/5,np.zeros(7)]).astype(complex)
 def rhs(t,y):
  rho=y[:144].reshape(12,12).T
  dr=(-1j*(H@rho-rho@H)).T.ravel()+decay@y[:144]
  rate,Q=emission_tensor(rho,ob.hamiltonian.d_q_bare['g->e'])
  return np.r_[dr,rate,ETA**2*Q.ravel()]
 sol=solve_ivp(rhs,(0,T),np.r_[rho0.T.ravel(),0.,np.zeros(9)],method='DOP853',rtol=2e-10,atol=2e-12)
 if not sol.success:raise RuntimeError(sol.message)
 y=sol.y[:,-1]
 return y[:144].reshape(12,12).T,float(y[144].real),y[145:].real.reshape(3,3)

if __name__=='__main__':
 out=Path('results/quantum_jumps');out.mkdir(exist_ok=True);ob,H,d,closure,evolve,checks=setup();T=200.
 ref,N,C=reference(ob,H,d,T);np.savez_compressed(out/'reference.npz',rho=ref,expected_count=N,recoil_covariance=C,H=H)
 (out/'propagator_checks.json').write_text(json.dumps(checks,indent=2))
 states=[];counts=[];recoils=[];rows=[];start=time.time()
 for i in range(16384):
  psi,count,recoil,_=trajectory(evolve,d,T,np.random.default_rng(700000+i));states.append(psi);counts.append(count);recoils.append(recoil)
  if i+1 in [128,512,1024,2048,4096,8192,16384]:
   st=np.array(states);co=np.array(counts);re=np.array(recoils);rho=st.T@st.conj()/len(st)
   row=dict(n=len(st),density_frobenius_error=float(np.linalg.norm(rho-ref)),mean_count=float(co.mean()),count_standard_error=float(co.std(ddof=1)/np.sqrt(len(co))),expected_count=N,mean_recoil=re.mean(0).tolist(),sample_recoil_covariance=np.cov(re,rowvar=False).tolist(),reference_recoil_covariance=C.tolist(),covariance_relative_frobenius_error=float(np.linalg.norm(np.cov(re,rowvar=False)-C)/np.linalg.norm(C)),wall_s=time.time()-start);rows.append(row)
   (out/'convergence.json').write_text(json.dumps(rows,indent=2));print(row['n'],row['density_frobenius_error'],row['mean_count'],row['covariance_relative_frobenius_error'],flush=True)
 np.savez_compressed(out/'trajectories.npz',terminal_states=states,counts=counts,recoils=recoils)
