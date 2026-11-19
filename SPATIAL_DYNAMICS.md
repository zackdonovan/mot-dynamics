# spatially varying light and magnetic field

## model extension

The coupled density/position/velocity solver now evaluates collimated Gaussian beam amplitudes and the quadrupole field at the atom's evolving position. It retains one laboratory basis through B=0; there is no local-axis switching or population-only projection.

Parameters use the original independent model's declared choices: beam 1/e² intensity radius w=12 mm, radial magnetic gradient b=0.12 T/m, B=b(x,y,−2z), total on-axis saturation 6, and the F=2→F′=3 closed Rb87 manifold. These are synthetic design choices, not measured apparatus parameters. The pair-offset optical frequencies and pure helicities follow v5/v6.

For beam direction k̂, the electric-field envelope is exp[−|r−k̂(k̂·r)|²/w²]. Its derivative is included alongside the optical phase gradient. The force is −Tr(ρ∇H), including both optical amplitude/phase derivatives and the magnetic-gradient force from the Zeeman Hamiltonian. The latter was omitted from earlier local benchmarks because their field was uniform. Wavefront curvature and longitudinal beam expansion are omitted under the collimated-beam approximation.

## controlled study

Eight main trajectories form four matched pairs: six/eight beam geometries, starting at the field zero or at (6,3,2) mm. Every pair compares:

- full spatial evolution of beam envelopes and quadrupole field;
- a locally frozen approximation, with amplitudes and field held at their initial values and their gradients removed, while optical phase and force-feedback motion still evolve.

The same phase realization is defined locally at the initial point for each comparison. Shifting the macroscopic start point therefore does not accidentally introduce unrelated optical phases. All runs start with v=(0.05,0,0)Γ/k and equal ground populations, lasting 2000 Γ⁻¹, about 52.5 μs. The initial speed is 0.237 m/s. Initial positions and phases are probes, not a thermal ensemble.

## measured validity of the local approximation

| geometry / start | final full-vs-frozen velocity difference, m/s | maximum relative beam-intensity change during path |
|---|---:|---:|
| six / origin | 0.00114 | 0.000223% |
| six / off-axis | 0.00427 | 0.314% |
| eight / origin | 0.00192 | 0.000172% |
| eight / off-axis | 0.00060 | 0.245% |

The maximum magnetic-field change is 0.00615 in μ_B B/(ħΓ) units. For these short paths, freezing the local envelopes and field introduces a small but measurable velocity error. This gives an actual scale-limited justification for the local approximation rather than assuming it is always adequate.

It does not justify using the **same** field and intensity everywhere. Off-axis axial-beam intensity factors at the starting point range from about 0.535 to 0.835 of their on-axis values. Local settings differ materially between locations even when they vary little over one short trajectory. The off-axis coherent trajectories can reverse longitudinal motion; this is not a capture classification or a statement about cooling performance.

## verification

- Four full-duration repeats at 100× tighter integration tolerances change final velocities by ≤2.11×10⁻¹² in Γ/k units and final displacements by ≤8.62×10⁻¹⁷ m.
- Four direct PyLCP spatial coupled-motion comparisons over 100 Γ⁻¹ agree within 2.48×10⁻⁹ in Γ/k velocity units and 1.57×10⁻¹⁴ m in displacement. Both solve the same mean-force approximation, so this validates implementation rather than physical completeness.
- At ten points per full trajectory, public PyLCP Hamiltonians agree within 7.84×10⁻¹², and forces within 5.80×10⁻⁸ in ħkΓ units. PyLCP's Gaussian-beam force uses numerical spatial derivatives; our implementation uses the analytic envelope/phase derivatives.
- A separate random physical density matrix verifies the full force against finite differences of −Tr(ρH), including magnetic and optical gradients.
- Zero field gradient and infinite beam waist recover the frozen uniform-plane-wave limit.
- Maximum density trace error is 1.58×10⁻¹⁴; recorded density matrices remain positive to numerical precision.
- All 26 project tests pass.

## interpretation

This is a spatially resolved deterministic semiclassical model, not experimentally validated capture simulation. It includes spontaneous decay of the internal state but not momentum diffusion from random photon recoil. It does not simulate a loading source, collisions, finite trap lifetime, experimental laser noise or measured optical imperfections. Two positions, one phase seed and a short time interval cannot establish trap performance.

The next physical uncertainty to address is stochastic momentum evolution and its consistency with the coherent internal dynamics. Adding arbitrary independent kicks would not automatically produce a validated temperature or loading rate. Longer trajectories and broader initial conditions should follow, with explicit computational and modeling limits.

The work remains an independent project, not recovered Sygaldry code or a reproduction of company-specific apparatus.

## reproduce

```sh
OPENBLAS_NUM_THREADS=1 python benchmark/spatial_motion.py
OPENBLAS_NUM_THREADS=1 python benchmark/spatial_checks.py
OPENBLAS_NUM_THREADS=1 python -m unittest benchmark/test_spatial.py
python benchmark/spatial_report.py
```

Raw trajectories, full/frozen comparisons and checks are in `results/spatial/`. Dependencies: `requirements-benchmark.txt`. Earlier conventions and [PyLCP reference](https://pmc.ncbi.nlm.nih.gov/articles/PMC9890571/) are documented in `PYLCP_BENCHMARK.md`.
