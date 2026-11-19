import unittest,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
from spatial_motion import system
from recoil_validity import cartesian_dipoles,emission_tensor,diagnostic
class RecoilValidityTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):
  ob,*_=system('axial6',np.zeros(6),np.zeros(6),[0,0,0]);cls.dq=ob.hamiltonian.d_q_bare['g->e']
 def test_isotropic_excited_and_ground_limits(self):
  rho=np.diag(np.r_[np.zeros(5),np.ones(7)/7]).astype(complex);rate,Q=emission_tensor(rho,self.dq)
  self.assertAlmostEqual(float(rate),1);np.testing.assert_allclose(Q,np.eye(3)/3,atol=1e-13)
  rho=np.diag(np.r_[np.ones(5)/5,np.zeros(7)]).astype(complex);rate,Q=emission_tensor(rho,self.dq);np.testing.assert_allclose(Q,0,atol=1e-13)
 def test_angular_quadrature_and_bounds(self):
  rng=np.random.default_rng(822);z=rng.normal(size=12)+1j*rng.normal(size=12);z/=np.linalg.norm(z);rho=np.outer(z,z.conj());rate,Q=emission_tensor(rho,self.dq)
  d=cartesian_dipoles(self.dq);M=np.array([[np.trace(rho@a.conj().T@b) for b in d] for a in d]).real
  n=rng.normal(size=(500000,3));n/=np.linalg.norm(n,axis=1)[:,None]
  w=np.trace(M)-np.einsum('ni,ij,nj->n',n,M,n)
  numeric=np.einsum('n,ni,nj->ij',w,n,n)/w.sum()
  np.testing.assert_allclose(numeric,Q/rate,atol=.0015)
  eig=np.linalg.eigvalsh(Q);self.assertGreaterEqual(eig.min(),.2*rate-1e-13);self.assertLessEqual(eig.max(),.4*rate+1e-13)
 def test_integrated_constant_rate(self):
  rho=np.diag(np.r_[np.zeros(5),np.ones(7)/7]).astype(complex)
  a=diagnostic(np.array([0.,10.]),np.array([rho,rho]),self.dq,.01)
  self.assertAlmostEqual(a['expected_spontaneous_photons'],10)
  np.testing.assert_allclose(a['free_recoil_velocity_covariance'],np.eye(3)*.001/3,atol=1e-13)
if __name__=='__main__':unittest.main()
