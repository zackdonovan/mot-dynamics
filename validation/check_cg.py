"""Independent sympy angular-momentum cross-check, not a runtime dependency."""
import numpy as np
from sympy.physics.wigner import clebsch_gordan
F=2
branches={me:0.0 for me in range(-3,4)}
for m in range(-F,F+1):
    for q in [-1,0,1]:
        exact=float(clebsch_gordan(F,1,F+1,m,q,m+q))**2
        if q==0: formula=((F+1)**2-m*m)/((F+1)*(2*F+1))
        else: formula=(F+q*m+1)*(F+q*m+2)/(2*(F+1)*(2*F+1))
        np.testing.assert_allclose(formula,exact,atol=1e-15)
        branches[m+q]+=exact
np.testing.assert_allclose(list(branches.values()),1,atol=1e-15)
print('15 CG coefficients and 7 spontaneous branch sums independently verified')
