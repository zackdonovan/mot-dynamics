import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
from coupled_motion import run
from moving_coherence import run as fixed
class CoupledTests(unittest.TestCase):
 def test_ballistic_dark_limit(self):
  a,t,r,v,f=run('axial6',np.zeros(6),np.zeros(6),v0=(.05,.02,-.01),T=20,total_s=0)
  np.testing.assert_allclose(r[-1],[1,.4,-.2],atol=1e-10);np.testing.assert_allclose(v[-1],[.05,.02,-.01],atol=1e-12)
  np.testing.assert_allclose(f,0,atol=1e-12)
 def test_infinite_mass_prescribed_path(self):
  p=np.arange(6)*.41;df=np.repeat([0.,.12,.29],2)
  a,t,r,v,f=run('axial6',p,df,T=20,eta=0)
  b,_,ff=fixed('axial6',p,df,[.05,0,0],T=20)
  np.testing.assert_allclose(r[-1],[1,0,0],atol=1e-10)
  np.testing.assert_allclose(f[-1],ff[-1],atol=1e-8)
 def test_public_coupled_reference(self):
  a,*_=run('oblique8',np.arange(8)*.31,np.repeat([0.,.12,.29,.47],2),T=20,reference=True)
  self.assertLess(a['public_density_error'],1e-7);self.assertLess(a['public_velocity_error'],1e-9);self.assertLess(a['public_position_error'],1e-8)
  self.assertLess(a['max_public_force_error'],1e-10)
if __name__=='__main__':unittest.main()
