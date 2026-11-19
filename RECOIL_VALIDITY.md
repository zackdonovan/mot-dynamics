# recoil budget and computational limits

## why this step

The coherent solver now includes optical phases, internal-state memory, motion and spatial fields. The remaining omission is not cosmetic: spontaneous photons impart random momentum. Before scaling up deterministic capture calculations, this study measures the size of that missing effect and the cost of longer integrations.

Eight spatial trajectories compare 2000 and 8000 Γ⁻¹ (52.5 and 210 μs) for six/eight beams at the origin and (6,3,2) mm, using the same initial phase seed and frequency-offset architecture as v7. They remain deterministic. The new calculation evaluates a **free spontaneous-recoil noise budget along each deterministic path**, not a stochastic trajectory ensemble or final velocity distribution.

## angular recoil tensor

Let d_i be the Cartesian dipole lowering matrices, and M_ij=Tr(ρ d_i†d_j). Summing the two transverse photon polarizations gives an angular emission weight proportional to Tr M−nᵀ Re(M)n. Integrating the unit-sphere moments yields, for Γ=1:

```
R = Tr M = excited-state population
Q = (2 R I − Re M)/5
C_free = η² ∫ Q(t) dt
η = ħk²/(mΓ)
```

Q is the emission-direction second moment **times rate**, not a convention-dependent Fokker–Planck coefficient D. Thus C_free is directly the accumulated free velocity covariance in dimensionless units. It retains dipole angular anisotropy and excited-state coherences. For isotropic emission Q=R I/3; the eigenvalues of Q/R lie between 1/5 and 2/5. Trace C_free=η² times the expected emitted-photon count.

The tensor is derived here from the transverse-dipole radiation law, and checked by 500,000-direction numerical angular quadrature. Atomic normalization is independently checked by R=Trρ_exc to 1.12×10⁻¹⁶. Ground-only and isotropically excited limits and a constant-rate integration test pass.

## measured scale

| duration | expected spontaneous photons | free recoil RMS in individual Cartesian components |
|---|---:|---:|
| 52.5 μs | 159–417 | approximately 0.042–0.071 m/s |
| 210 μs | 634–1554 | approximately 0.084–0.137 m/s |

Initial deterministic speed is 0.237 m/s. Even the short-run recoil scale is substantially larger than the 0.00060–0.00427 m/s full-versus-frozen spatial-model differences in v7. This does not invalidate implementation checks; it shows that making deterministic spatial integration more precise is not enough to establish physically accurate trajectory distributions.

The longer off-axis runs have a three-dimensional free recoil RMS norm near 0.23 m/s, comparable to the initial speed. At that scale, a perturbative noise budget is a warning, not a reliable prediction of the final spread.

Reducing the output quadrature interval from 2 to 0.25 Γ⁻¹ in two short cases changes integrated photon counts by at most 0.0155%. All 29 project tests pass. A selected 210 μs case is also repeated at 100× tighter ODE tolerances; the final dimensionless velocity changes by 9.59×10⁻¹⁰ and the integrated count by 7.01×10⁻⁷ photons. Its check is retained in `long_tolerance_check.json`.

## what this does not establish

C_free is **not** the final covariance of a trapped or damped atom, and is not a lower bound on that covariance. It ignores mechanical response to kicks, recoil-induced changes of optical phase/internal state, and correlations with subsequent optical forces. Damping can reduce accumulated velocity spread. Absorption and stimulated-force fluctuations contribute additional diffusion, with internal-state correlations important beyond weak excitation. Consequently, adding independent isotropic kicks to a mean-force path would not by itself validate a coherent cooling model. See the [NIST-hosted laser-cooling review](https://www.nist.gov/system/files/documents/2017/05/09/rmp1.pdf) for the role of absorption/emission fluctuations and their correlations.

No temperature, capture fraction or laboratory performance is inferred from this budget. It is a quantitative adequacy check that identifies the next necessary physics.

## computational feasibility

The four 210 μs spatial density-matrix runs took 22–25 seconds each in this single-threaded environment, including diagnostic evaluation. A naive linear extrapolation gives about 2.2–2.5 hours per 75 ms trajectory, or 47–53 serial days for 512 trajectories at one parameter setting. This is a rough workload estimate, not a benchmark guarantee: costs change with intensity, velocity, tolerances, trajectory escape, implementation and parallelism. It excludes a stochastic ensemble multiplier and optimizer evaluations.

A brute-force replacement of the original fast optimizer with long coherent trajectories is therefore not a sensible next step here. The stronger research design is multi-fidelity: retain the inexpensive declared population model for exploration; use coherent, spatial calculations on a controlled subset; quantify approximation failures; then choose a justified stochastic treatment before any claims about temperature or capture.

A physically credible next stage would require a suitable quantum-regression force-correlation treatment or validated stochastic unravelling/semiclassical diffusion model, with reference limits and convergence checks. This study stops short of presenting an arbitrary noise model as that result.

## reproduce

```sh
OPENBLAS_NUM_THREADS=1 python benchmark/recoil_study.py
OPENBLAS_NUM_THREADS=1 python benchmark/recoil_checks.py
OPENBLAS_NUM_THREADS=1 python benchmark/recoil_long_check.py
OPENBLAS_NUM_THREADS=1 python -m unittest benchmark/test_recoil_validity.py
PYTHONPATH=. python benchmark/recoil_report.py
```

Data, covariance tensors, trajectory histories, timing and plots: `results/recoil_validity/`. Dependency and physical conventions remain in the previous benchmark reports. This is an independent simulation study, not recovered internship work or experimental validation.
