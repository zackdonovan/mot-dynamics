import unittest,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import numpy as np
from pylcp_multibeam import solve
class MultibeamTests(unittest.TestCase):
 def test_global_phase_and_density(self):
  for geom,n in [('axial6',6),('oblique8',8)]:
   p=np.random.default_rng(8).uniform(0,6.28,n)
   a=solve(geom,[.1,.2,.3],p);b=solve(geom,[.1,.2,.3],p+1.13)
   np.testing.assert_allclose(a['force'],b['force'],atol=1e-10)
   self.assertLess(a['residual'],1e-10);self.assertLess(a['trace_error'],1e-10)
   self.assertLess(a['hermiticity_error'],1e-10);self.assertGreater(a['min_eigenvalue'],-1e-10)
 def test_continuity_at_zero(self):
  p=np.random.default_rng(18).uniform(0,6.28,6)
  a=solve('axial6',[0.,0,0],p);b=solve('axial6',[1e-8,0,0],p)
  self.assertLess(np.linalg.norm(np.array(a['force'])-b['force']),1e-7)
if __name__=='__main__':unittest.main()
