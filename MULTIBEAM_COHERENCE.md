# multibeam coherence and phase sensitivity

## scope

This extends the external optical-Bloch benchmark from a single beam to the actual six axial and eight axial-plus-oblique beam directions. The calculation retains the complete F=2→F′=3 density matrix in one fixed laboratory basis. It is a separate frozen-atom benchmark, not a replacement for the trajectory solver.

The closed manifold, g factors, transition strengths and signs match the production model. Total incident saturation is 6, distributed equally over the beams. Every beam is purely circular, with the existing (+,+,−) axial-pair helicity pattern and positive oblique-pair helicity. All share detuning Δ/Γ=−2. Magnetic fields in μ_B B/(ℏΓ) units are (0,0,0), (10⁻⁸,0,0), (0.1,0,0) and (0,0,0.1). The atom is stationary; plane-wave amplitudes are evaluated at the local origin. There is no Gaussian envelope gradient, magnetic-gradient force, gravity, motion or recoil.

Equal detunings make a time-independent rotating-frame Hamiltonian possible. The optimized production configurations have unequal detunings and mixed helicities, so these are not their coherent force predictions. An incoherent helicity mixture cannot simply be replaced by two phase-coherent circular fields.

## what the phases mean

Each geometry has 128 seeded independent samples of relative optical phases, uniform on [0,2π). One common phase is fixed to zero because it is unobservable. The same phase samples are reused across fields, allowing a paired near-zero-field continuity check. Each density matrix is solved to its own steady state, then forces are averaged. This is an ensemble of stationary, phase-coherent configurations with unknown fixed phases: **quenched phase disorder**.

This is not fast phase diffusion, frequency-offset averaging or incoherent illumination. Averaging a steady-state solution and solving an averaged generator are different operations. Independent phase sampling is also not necessarily the same as averaging over the three-dimensional spatial trajectory of one fixed interferometer. Interference changes the local intensity and polarization even though the sum of incident beam powers is fixed. Different local forces therefore do not imply a normalization failure.

## numerical checks

The trace-constrained 144-component Liouvillian solve follows the previously documented PyLCP 1.0.2 Hamiltonian normalization and density-vector conventions. Every sample records residual, trace, Hermiticity, positivity, excited-state population and force. One phase realization per geometry is also integrated from equal ground populations for 12000 Γ⁻¹ with BDF (relative tolerance 10⁻⁹, absolute tolerance 10⁻¹¹). A common phase shift of 0.731 radians must leave force and populations invariant. Unit tests independently check a different common shift and near-zero-field continuity.

Raw ensembles and numerical summaries are in `results/multibeam/`; the plot compares phase-dependent x forces at B=(0.1,0,0). Standard errors describe finite sampling of this chosen phase ensemble, not uncertainty about a real apparatus or model error. The first-64 and full-128 means are saved as a limited sampling-stability check, not proof of convergence.

## measured results

There are 1,024 coherent steady states: two geometries × four fields × 128 phases. The maximum Liouvillian residual is 6.50×10⁻¹⁵, trace error 3.33×10⁻¹⁵, and Hermiticity error 7.32×10⁻¹³. The smallest density eigenvalue is −4.97×10⁻¹⁴ (numerical roundoff). Common phase shifts change the force by less than 6×10⁻¹⁵. Direct steady states agree with time-evolved density matrices within 6.89×10⁻⁹ in Frobenius norm for the two selected checks.

Using matched phases at zero and nonzero field reduces sampling noise in the field response. In force units ℏkΓ:

| geometry / field | force component | coherent mean change ± sampling SE | population-model change |
|---|---|---:|---:|
| six beams, Bx=0.1 | x | −0.01886 ± 0.00318 | −0.00286 |
| six beams, Bz=0.1 | z | +0.01309 ± 0.00266 | +0.00286 |
| eight beams, Bx=0.1 | x | −0.01552 ± 0.00282 | −0.00286 |
| eight beams, Bz=0.1 | z | +0.01228 ± 0.00343 | +0.00143 |

These are finite field changes from zero, not fitted spring constants. Their larger coherent response is specific to the declared stationary, monochromatic phase ensemble. It is not an improved capture estimate. Individual phase realizations also have substantial force offsets: at Bx=0.1 the six-beam x-force standard deviation is 0.0681, much larger than its population-model mean force magnitude of 0.00286. The plot shows this variability, not an error bar on a laboratory measurement.

Across all matched phase realizations, changing the field from zero to (10⁻⁸,0,0) changes the full force vector by at most 2.07×10⁻⁸ (six beams) or 4.13×10⁻⁸ (eight beams). The fixed-basis coherent calculation is continuous in this tested near-zero limit.

## interpretation and remaining work

The relevant question is whether a population-only model can stand in for mutually coherent beams without specifying the phases. The phase-force distributions answer that question within this monochromatic, stationary setting. The results do not quantify capture-fraction error. Forces in a moving atom require treatment of Doppler-shifted frequencies and possibly nonadiabatic internal evolution; treating every point along a trajectory as a stationary coherent equilibrium would need separate justification.

A defensible next extension would choose an experimentally meaningful laser architecture, then compare static-phase, phase-diffusing and frequency-offset illumination on a small set of controlled trajectories. Merely inserting these steady forces into thousands of capture trajectories would be computationally expensive and scientifically under-specified.

This independent research implementation is not recovered Sygaldry internship work and has no experimental validation. PyLCP comparison does not imply NIST endorsement.

## reproduce

Use the optional benchmark environment from `requirements-benchmark.txt`:

```sh
OPENBLAS_NUM_THREADS=1 python benchmark/pylcp_multibeam.py
OPENBLAS_NUM_THREADS=1 python -m unittest benchmark/test_multibeam.py
python benchmark/multibeam_report.py
```

Reference: [PyLCP framework paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9890571/).
