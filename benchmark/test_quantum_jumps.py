import unittest,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
from quantum_jumps import setup,trajectory,emit
class QuantumJumpTests(unittest.TestCase):
 @classmethod
 def setUpClass(cls):cls.ob,cls.H,cls.d,cls.closure,cls.evolve,cls.checks=setup()
 def test_propagator_and_closure(self):
  self.assertLess(self.checks['max_propagator_error'],1e-11);self.assertLess(self.checks['closure_error'],1e-12)
 def test_exact_single_decay_wait(self):
  psi=np.zeros(12,complex);psi[5]=1
  def free(s,t):return s*np.r_[np.ones(5),np.ones(7)*np.exp(-t/2)]
  for seed in range(10):
   expected=-np.log(np.random.default_rng(seed).random())
   final,count,recoil,events=trajectory(free,self.d,20.,np.random.default_rng(seed),initial_state=psi)
   self.assertEqual(count,1);self.assertAlmostEqual(events[0],expected,places=8);self.assertLess(np.linalg.norm(final[5:]),1e-12)
 def test_jump_maps_to_normalized_ground(self):
  rng=np.random.default_rng(143);p=rng.normal(size=12)+1j*rng.normal(size=12);p/=np.linalg.norm(p)
  for _ in range(30):
   state,n=emit(p,self.d,rng);self.assertAlmostEqual(float(np.linalg.norm(state)),1);self.assertAlmostEqual(float(np.linalg.norm(n)),1);self.assertLess(np.linalg.norm(state[5:]),1e-12)
if __name__=='__main__':unittest.main()
