import unittest
import numpy as np
from mot.physics import *

class ValidationTests(unittest.TestCase):
    def test_generator(self):
        w=rates(*sample_atoms(20,123),Controls());a=rate_generator(w);p=populations(w)
        np.testing.assert_allclose(a.sum(1),0,atol=2e-8)
        np.testing.assert_allclose(np.einsum('nij,nj->ni',a,p),0,atol=2e-8)
        off=a.copy();off[:,np.arange(12),np.arange(12)]=0
        self.assertTrue(np.all(off>=0))
    def test_ballistic_exact_time_and_integer_inputs(self):
        r=np.array([[1,2,3]]);v=np.array([[1,0,2]])
        t=.00123;a=np.array([0,0,-9.80665])
        for dt in [.0002,.1]:
            out=simulate(Controls(power=(0,0,0)),r,v,dt=dt,duration=t,recoil=True,seed=1)
            np.testing.assert_allclose(out['r'],r+v*t+.5*a*t*t,atol=1e-14)
            np.testing.assert_allclose(out['v'],v+a*t,atol=1e-14)
            self.assertAlmostEqual(out['elapsed'],t)
    def test_invalid_time_and_controls(self):
        for value in [0,-1,np.nan,np.inf]:
            with self.assertRaises(ValueError):simulate(Controls(),np.zeros((1,3)),np.zeros((1,3)),dt=value)
        for kw in [dict(waist=np.nan),dict(gradient=np.inf)]:
            with self.assertRaises(ValueError):Controls(**kw)
    def test_diffusion_zero_dark_and_positive(self):
        r,v=sample_atoms(50,31)
        _,d=force_and_diffusion(r,v,Controls());self.assertTrue(np.min(np.linalg.eigvalsh(d))>=0)
        _,dark=force_and_diffusion(r,v,Controls(power=(0,0,0)))
        np.testing.assert_array_equal(dark,np.zeros((50,3,3)))
    def test_monte_carlo_kick_covariance(self):
        n=24000;h=1e-6;r=np.zeros((n,3));v=r.copy();c=Controls()
        out=simulate(c,r,v,dt=h,duration=h,recoil=True,seed=843)
        deterministic=simulate(c,r[:1],v[:1],dt=h,duration=h)
        residual=out['v']-deterministic['v']
        _,d=force_and_diffusion(r[:1],np.array([[0,0,-9.80665*h/2]]),c)
        expected=d[0]*h/MASS**2; empirical=np.cov(residual.T,bias=True)
        np.testing.assert_allclose(np.diag(empirical),np.diag(expected),rtol=.035)
        self.assertLess(np.max(abs(empirical-np.diag(np.diag(empirical))))/np.max(np.diag(expected)),.025)
    def test_reproducibility(self):
        r,v=sample_atoms(5,32); c=Controls()
        a=simulate(c,r,v,duration=.0003,recoil=True,seed=61)
        b=simulate(c,r,v,duration=.0003,recoil=True,seed=61)
        np.testing.assert_array_equal(a['v'],b['v'])
if __name__=='__main__':unittest.main()
