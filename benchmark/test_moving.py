import unittest,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
from moving_coherence import run,directions
class MovingTests(unittest.TestCase):
 def test_doppler_equivalence(self):
  p=np.arange(6)*.37;v=np.array([.05,.02,-.03]);df=np.repeat([0.,.12,.29],2)
  a,t,fa=run('axial6',p,df,v,T=20)
  b,t,fb=run('axial6',p,df-directions('axial6')@v,[0,0,0],T=20)
  np.testing.assert_allclose(fa,fb,atol=2e-9)
 def test_public_pylcp_and_density(self):
  a,t,f=run('oblique8',np.arange(8)*.41,np.repeat([0.,.12,.29,.47],2),[.05,0,0],T=30,validate=True)
  self.assertLess(a['pylcp_time_density_error'],1e-7)
  self.assertLess(a['hamiltonian_reconstruction_error'],1e-12)
  self.assertLess(a['max_trace_error'],1e-10)
  self.assertGreater(a['min_density_eigenvalue'],-1e-7)
if __name__=='__main__':unittest.main()
