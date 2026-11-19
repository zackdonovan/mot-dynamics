"""Independent counting-field derivatives of the master-equation generator."""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse import csr_matrix
from quantum_jumps import setup
ob,H,d,closure,*_=setup();n=12
L=csr_matrix(-1j*(np.kron(np.eye(n),H)-np.kron(H.T,np.eye(n)))+ob.ev_mat['decay'])
J=csr_matrix(sum(np.kron(a.conj(),a) for a in d))
rho0=np.diag(np.r_[np.ones(5)/5,np.zeros(7)]).astype(complex).T.ravel()
def rhs(t,y):
 r,a,b=y.reshape(3,144)
 return np.r_[L@r,L@a+J@r,L@b+2*(J@a)]
s=solve_ivp(rhs,(0,200.),np.r_[rho0,np.zeros(288)],method='DOP853',rtol=2e-10,atol=2e-12)
r,a,b=s.y[:,-1].reshape(3,12,12);mean=float(np.trace(a).real);factorial=float(np.trace(b).real);var=factorial+mean-mean**2
Path('results/quantum_jumps/count_moments.json').write_text(json.dumps(dict(mean=mean,variance=var,fano=var/mean,integration_success=bool(s.success)),indent=2))
