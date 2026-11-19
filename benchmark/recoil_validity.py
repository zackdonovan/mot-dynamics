"""Spontaneous recoil budget along deterministic trajectories.
Not total momentum diffusion: excludes absorption/force noise and feedback.
"""
import numpy as np

def cartesian_dipoles(dq):
 return np.array([(dq[0]-dq[2])/np.sqrt(2),1j*(dq[0]+dq[2])/np.sqrt(2),dq[1]])

def emission_tensor(rho,dq):
 """Angular second moment times emission rate, with Gamma=1.
 M_ij=Tr(rho d_i^dagger d_j), integral over transverse photon polarizations.
 ∫nn^T dR/dΩ dΩ = (2 I TrM - ReM)/5.
 """
 d=cartesian_dipoles(dq)
 ops=np.einsum('aji,bjk->abik',d.conj(),d)
 M=np.einsum('...ji,abij->...ab',rho,ops)
 rate=np.trace(M,axis1=-2,axis2=-1).real
 Q=(2*rate[...,None,None]*np.eye(3)-M.real)/5
 return rate,Q

def diagnostic(ts,rho,dq,eta):
 rate,Q=emission_tensor(rho,dq)
 N=float(np.trapezoid(rate,ts));cov=eta**2*np.trapezoid(Q,ts,axis=0)
 return dict(expected_spontaneous_photons=N,free_recoil_velocity_covariance=cov.tolist(),free_recoil_rms_components=np.sqrt(np.maximum(np.diag(cov),0)).tolist(),free_recoil_rms_norm=float(np.sqrt(np.trace(cov))),max_emission_rate_error=float(np.max(np.abs(rate-np.trace(rho[:,5:,5:],axis1=1,axis2=2).real))),min_instantaneous_tensor_eigenvalue=float(np.linalg.eigvalsh(Q).min()))
