import unittest,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from pylcp_toy import case,analytical_rate
from pylcp_production_rates import pylcp_rate_force
from pylcp_zero_field import case as zero_case
from mot.physics import *
class ExternalBenchmarks(unittest.TestCase):
 def test_single_beam_exact(self):
  x=case(.1,-2.,.7,pair=False)
  self.assertLess(abs(x['pylcp_obe_force'][2]-x['two_level_exact']),1e-8)
  self.assertLess(x['steady_residual_norm'],1e-7)
 def test_production_rate_equations(self):
  r,v=sample_atoms(2,80)
  for geom,n in [('axial6',3),('oblique8',4)]:
   c=Controls(detuning=(-2.,)*n,purity=(.4,)*n,power=(.6,)*n,geometry=geom)
   a=force(r,v,c,False)
   for j in range(2):np.testing.assert_allclose(a[j],pylcp_rate_force(r[j],v[j],c),rtol=1e-9,atol=1e-32)
 def test_coherent_zero_field_limit(self):
  z=zero_case([0.,0,0]);x=zero_case([1e-8,0,0])
  self.assertLess(abs(z['fixed_basis_obe_force'][2]-x['fixed_basis_obe_force'][2]),1e-10)
  self.assertGreater(abs(x['pylcp_rate_force'][2]-x['fixed_basis_obe_force'][2]),.01)
  self.assertLess(x['steady_residual'],1e-10)
if __name__=='__main__':unittest.main()
