"""Slow independent scalar reference for 87Rb F=2→F'=3 rate model.
Uses Wigner CG, eigenspace null solve; meant for validation, not production.
All omega quantities angular frequency in SI. s is cycling Isat ratio.
"""
import numpy as np
from sympy.physics.wigner import clebsch_gordan

def steady(kvecs, detunings, saturations, helicities, B, v,
           gamma=2*np.pi*6.065e6, mu_over_hbar=9.2740100783e-24/1.054571817e-34):
    from scipy.linalg import null_space
    kvecs=np.asarray(kvecs); B=np.asarray(B); v=np.asarray(v)
    Bmag=np.linalg.norm(B)
    axis=B/Bmag if Bmag>1e-15 else np.array([0.,0.,1.])
    L=np.zeros((12,12)); beams=[]
    for k,delta,s,h in zip(kvecs,detunings,saturations,helicities):
        cos=np.dot(k/np.linalg.norm(k),axis)
        weights={0:(1-cos*cos)/2,1:(1+cos*cos+2*h*cos)/4,-1:(1+cos*cos-2*h*cos)/4}
        rates=np.zeros((5,7))
        for m in range(-2,3):
            for q in [-1,0,1]:
                me=m+q; strength=float(clebsch_gordan(2,1,3,m,q,me))**2
                d=delta-np.dot(k,v)-mu_over_hbar*Bmag*((2/3)*me-.5*m)
                w=gamma*s*strength*weights[q]/(2*(1+4*d*d/gamma**2))
                rates[m+2,me+3]=w
                g=m+2;e=me+8
                L[e,g]+=w;L[g,e]+=w;L[g,g]-=w;L[e,e]-=w
        beams.append(rates)
    for m in range(-2,3):
        for q in [-1,0,1]:
            me=m+q; strength=float(clebsch_gordan(2,1,3,m,q,me))**2
            g=m+2;e=me+8
            L[g,e]+=gamma*strength;L[e,e]-=gamma*strength
    null=null_space(L)
    if null.shape[1]!=1: raise ValueError('nonunique stationary distribution')
    pop=null[:,0]/null[:,0].sum()
    force=np.zeros(3)
    for k,rates in zip(kvecs,beams):
        force+=1.054571817e-34*k*np.sum(rates*(pop[:5,None]-pop[None,5:]))
    return force,pop,L

if __name__=='__main__':
    k=2*np.pi/780.241209686e-9
    dirs=np.vstack((np.eye(3),-np.eye(3)))
    hs=[1,1,-1,1,1,-1]
    args=(dirs*k,[-2*2*np.pi*6.065e6]*6,[1]*6,hs)
    for dim in range(3):
        r=np.eye(3)[dim]*1e-3; B=.1*r*np.array([1,1,-2])
        f,p,L=steady(*args,B,[0,0,0]); assert f[dim]<0
        np.testing.assert_allclose(L.sum(axis=0),0,atol=1e-7)
        assert p.min()>-1e-12
        print('restoring',dim, f)
        f,p,L=steady(*args,[0,0,0],np.eye(3)[dim]); assert f[dim]<0
        print('damping',dim, f)
