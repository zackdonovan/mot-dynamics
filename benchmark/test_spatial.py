import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
from spatial_motion import system
class SpatialTests(unittest.TestCase):
 def test_full_hamiltonian_gradient(self):
  p=np.arange(6)*.3;df=np.repeat([0.,.12,.29],2)
  ob,hf,*_=system('axial6',p,df,[.006,.003,.002])
  rng=np.random.default_rng(14);z=rng.normal(size=(12,12))+1j*rng.normal(size=(12,12));rho=z@z.conj().T;rho/=np.trace(rho)
  x=np.array([1.,2.,3.]);t=2.;_,f,*_=hf(x,t,rho);dx=1e-4;numeric=[]
  for e in np.eye(3):numeric.append(-np.trace(rho@(hf(x+dx*e,t,rho)[0]-hf(x-dx*e,t,rho)[0])/(2*dx)).real)
  np.testing.assert_allclose(f,numeric,rtol=1e-6,atol=1e-9)
 def test_base_detuning_parameter(self):
  p=np.arange(6)*.3;df=np.repeat([0.,.12,.29],2);rho=np.eye(12)/12
  _,a,*_=system('axial6',p,df,[.006,.003,.002],frozen=True,base_detuning=-1.)
  _,b,*_=system('axial6',p,df,[.006,.003,.002],frozen=True,base_detuning=-3.)
  ha=a(np.zeros(3),0.,rho)[0];hb=b(np.zeros(3),0.,rho)[0]
  np.testing.assert_allclose(hb-ha,np.diag(np.r_[np.zeros(5),np.ones(7)*2]),atol=1e-12)
 def test_frozen_uniform_plane_limit(self):
  p=np.arange(8)*.3;df=np.repeat([0.,.12,.29,.47],2);rho=np.eye(12)/12
  _,a,*_=system('oblique8',p,df,[0,0,0],waist=np.inf,gradient=0,frozen=False)
  _,b,*_=system('oblique8',p,df,[0,0,0],waist=np.inf,gradient=0,frozen=True)
  for x in [np.zeros(3),np.array([5.,3.,-2.])]:
   ha,fa,*_=a(x,2.,rho);hb,fb,*_=b(x,2.,rho)
   np.testing.assert_allclose(ha,hb,atol=1e-12);np.testing.assert_allclose(fa,fb,atol=1e-12)
if __name__=='__main__':unittest.main()
