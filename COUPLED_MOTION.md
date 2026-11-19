# coupled internal and mechanical dynamics

## what changed

The atom's motion is no longer prescribed. The full 12-state density matrix, three position coordinates and three velocity components are integrated together. The mean optical force changes the velocity, which changes the position and hence the optical phases seen by the atom. This closes the feedback loop missing from the fixed-path v5 experiment.

In dimensionless units, with time Γt, position kr and velocity kv/Γ:

```
dρ/dt = −i[H(r,t),ρ] + spontaneous-emission Lindblad decay
dr/dt = v
dv/dt = η F(r,t,ρ)
η = ħk²/(mΓ) = 0.00124352
F = −Tr(ρ ∇H), in units of ħkΓ
```

This is deterministic semiclassical **mean-force** motion. Spontaneous decay is included in the internal equation, but its random recoil is not applied to the external motion. This is not a Monte Carlo wavefunction simulation, a temperature prediction, or a capture calculation.

## controlled experiment

Eight runs use the actual six/eight beam directions, two optical-phase realizations (seeds 721 and 722), and two laser architectures: common frequency or the previously declared distinct frequency offsets per counterpropagating pair. All start at r=0, v=(0.05,0,0)Γ/k with equal ground-state populations, at total incident saturation 6 and base detuning −2Γ. The fixed local magnetic field is B=(0.1,0,0) in μ_B B/(ħΓ) units.

Each trajectory lasts 2000 Γ⁻¹, approximately 52.5 μs for the project's Rb87 constants. The initial speed is approximately 0.237 m/s. There are no Gaussian intensity envelopes, quadrupole-field variation, collisions, gravity or recoil diffusion. Paths cover only micrometers. The two phase realizations are controlled probes, not an ensemble large enough to estimate population statistics.

The coherent Hamiltonian is evaluated at the evolving position. Its exact plane-wave gradient supplies the force. The force implementation is independently compared with PyLCP's public force method at ten points along each trajectory.

## results

Final x velocities in units Γ/k, from initial vx=0.05:

| geometry / phase seed | common frequency | offset pairs |
|---|---:|---:|
| six / 721 | −0.00718 | +0.05030 |
| six / 722 | −0.01944 | +0.04719 |
| eight / 721 | −0.01514 | +0.04633 |
| eight / 722 | +0.00316 | +0.04239 |

The common-frequency cases develop substantial transverse motion and sometimes reverse vx. In these probes, the offset-pair cases remain closer to their initial longitudinal speed. This confirms that force feedback was important and that the imposed-path result could not be relabeled as a free trajectory. It does not establish that one architecture has better trapping or cooling: a velocity reversal in a coherent interference field is not a capture event, and changes in a few atoms' kinetic energies are not a temperature measurement.

The plot `results/coupled/coupled_velocity.png` shows complete velocity histories; compressed trajectory files include time, position, velocity and force. Frequency offsets also alter individual laser detunings, as documented in v5. Architectural differences cannot be attributed solely to loss of coherence.

## validation

- Eight main trajectories complete successfully.
- Maximum force difference against the public PyLCP force calculation: 6.62×10⁻¹⁵ in ħkΓ units.
- Four direct public PyLCP `evolve_motion` comparisons over 200 Γ⁻¹ agree to 8.45×10⁻¹⁵ in density norm, 1.28×10⁻¹⁶ in velocity and 2.16×10⁻¹⁴ in position. Both implementations solve the same mean-force equations; agreement is an implementation check, not independent physical validation.
- Four full-duration runs repeated with 100× tighter ODE tolerances change final velocities by at most 3.92×10⁻¹² and final positions by 2.79×10⁻¹⁰ in dimensionless units.
- Maximum trace error is 8.22×10⁻¹⁵; density positivity holds at recorded points to numerical precision.
- An independent force-integral check compares Δv with η∫Fdt. Coarse output quadrature has error up to 1.99×10⁻⁴; reducing the sampling interval from 2 to 0.25 Γ⁻¹ in four selected runs reduces the maximum error to 2.19×10⁻⁶. These are output-quadrature errors, not the adaptive solver's velocity error.
- Unit tests recover the analytic zero-light ballistic trajectory, the infinite-mass prescribed-path limit, and a public coupled-motion reference. The complete project suite has 24 passing tests.

## remaining physical scope

These are self-consistent local mean-force trajectories, not a full MOT model. A next physical extension would add spatially varying intensity and quadrupole field and study stochastic recoil or another justified momentum-diffusion treatment. Doing so should start with a small, controlled reference set and independent limits, not immediately reinterpret the earlier population-model capture scores as coherent predictions.

There is no experimental validation, no claim of company-specific geometry or parameters, and no claim that this independent project reproduces Sygaldry internship work.

## reproduce

```sh
OPENBLAS_NUM_THREADS=1 python benchmark/coupled_motion.py
OPENBLAS_NUM_THREADS=1 python benchmark/coupled_checks.py
OPENBLAS_NUM_THREADS=1 python benchmark/coupled_quadrature.py
OPENBLAS_NUM_THREADS=1 python -m unittest benchmark/test_coupled.py
python benchmark/coupled_report.py
```

Dependencies and atomic conventions are in `requirements-benchmark.txt` and `PYLCP_BENCHMARK.md`. Raw results are in `results/coupled/`. Reference framework: [PyLCP](https://pmc.ncbi.nlm.nih.gov/articles/PMC9890571/).
