import unittest
import numpy as np
from mot.physics import *
class PhysicsTests(unittest.TestCase):
    def test_branches(self):
        for e in range(7):self.assertAlmostEqual(CG[EINDEX==e].sum(),1)
    def test_polarization(self):
        c=np.linspace(-1,1,51)
        for h in [-1,-.3,0,.8,1]:
            p=polarization_weights(c,h);self.assertTrue((p>=0).all());np.testing.assert_allclose(p.sum(-1),1)
    def test_restoring_damping(self):
        c=Controls();r=np.eye(3)*1e-4
        self.assertTrue((np.diag(force(r,r*0,c,False))<0).all())
        self.assertTrue((np.diag(force(r*0,np.eye(3)*.1,c,False))<0).all())
    def test_population(self):
        r,v=sample_atoms(40,1);p=populations(rates(r,v,Controls()))
        self.assertTrue((p>=-1e-12).all());np.testing.assert_allclose(p.sum(1),1)
    def test_dark(self):
        f=force(np.zeros((1,3)),np.zeros((1,3)),Controls(power=(0,0,0)))
        np.testing.assert_allclose(f,[[0,0,-MASS*9.80665]])
    def test_inversion(self):
        r,v=sample_atoms(40,3);c=Controls()
        np.testing.assert_allclose(force(r,v,c,False),-force(-r,-v,c,False),atol=1e-33)
    def test_invalid(self):
        with self.assertRaises(ValueError):Controls(purity=(2,1,1))
    def test_eight(self):
        c=Controls(detuning=(-2,)*4,purity=(1,1,1,0),power=(.75,)*4,geometry='oblique8')
        self.assertEqual(beam_directions(c).shape,(8,3));self.assertTrue(np.isfinite(force(*sample_atoms(10,2),c)).all())
if __name__=='__main__':unittest.main()
