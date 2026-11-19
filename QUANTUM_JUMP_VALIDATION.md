# photon-resolved stochastic validation

## finish bar

Before adding noise to moving-atom predictions, establish that a photon-resolved stochastic internal model reproduces the same atomic master equation, emitted-photon statistics and spontaneous-recoil covariance. This milestone is complete for one explicitly controlled fixed-atom configuration. It is not a validated stochastic MOT or capture simulation.

## model and algorithm

The test uses the full F=2→F′=3 manifold, six mutually coherent stationary plane waves, total incident saturation 6, common detuning −2Γ, pure circular helicities, field (0.1,0,0) in μ_B B/(ħΓ) units and phase seed 721. Each trajectory starts in one uniformly sampled ground magnetic sublevel, reproducing the equal-ground-population mixed state. The atom is held fixed for 200 Γ⁻¹, about 5.25 μs.

Between emissions, the state propagates under the non-Hermitian generator:

```
A = −iH − ½ Σ_i d_i†d_i
|ψ̃(t)⟩ = exp(A t)|ψ(0)⟩
```

A uniform random threshold determines the next jump when the squared norm falls to that threshold. Root finding locates the event time; there is no fixed-step Bernoulli approximation or artificial restriction to one photon per time bin. Spectral propagation is checked against a direct matrix exponential, with maximum error 7.14×10⁻¹⁵ over tested times and eigenvector condition number 2.00.

The emitted direction n is sampled from the transverse-dipole radiation distribution. A transverse polarization is then sampled with its conditional transition probability, and the corresponding dipole operator collapses the internal state. Emitted momentum −ħk n is accumulated in a bookkeeping vector but **does not feed back on position or internal dynamics**. Summing the jump operators recovers the excited-state decay operator to 3.85×10⁻¹⁶.

This follows the [Monte Carlo wave-function formulation](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.68.580). The implementation and benchmarks are new parts of this independent project, not original internship code.

## independent references

The deterministic master equation supplies the final density and the integrated spontaneous-emission recoil tensor derived in v8. A separate counting-field derivative calculation evolves density operators for the first and second factorial photon-count moments:

```
ρ̇₀ = Lρ₀
ρ̇₁ = Lρ₁ + Jρ₀
ρ̇₂ = Lρ₂ + 2Jρ₁
Jρ = Σ_i d_i ρ d_i†
E[N] = Trρ₁
Var(N) = Trρ₂ + E[N] − E[N]²
```

This reference does not estimate photon-count variance by assuming a Poisson process.

## measured outcome

16,384 independent seeded trajectories give:

| quantity | sampled | master-equation reference |
|---|---:|---:|
| mean emitted photons | 10.7351 | 10.7933 |
| photon-count variance | 24.4472 | 24.2796 |
| Fano factor, variance/mean | 2.2773 | 2.2495 |

The mean-count discrepancy is 1.51 estimated sampling standard errors; variance discrepancy is 0.68. The distinctly non-Poisson Fano factor demonstrates that internal-state dynamics matter for photon statistics. A Poisson count model with the same mean would miss the variance in this case.

- Final density Frobenius error decreases from 0.106 at 128 trajectories to 0.00703 at 16,384. Four independent 4096-trajectory blocks have errors 0.0115–0.0149, supporting the observed sampling convergence. Prefix estimates are nested, not independent replicates, and errors need not decrease monotonically in general.
- The spontaneous-recoil covariance differs from its integrated master-equation reference by 1.73% in relative Frobenius norm. The three marginal variance discrepancies are within 1.84 estimated sampling standard errors; mean recoil is consistent with zero.
- Exact single-decay waiting-time tests recover the exponential-law event time. Jump states normalize and land in the ground manifold. The complete project suite has 33 passing tests.

Sampling standard errors here describe this finite Monte Carlo ensemble, not physical model uncertainty. The tests establish one controlled configuration, not every parameter regime or numerical conditioning case.

## what is now stronger

The project has moved beyond deterministic trajectories with an estimated noise budget: it now contains an explicit photon-resolved internal stochastic process, with exact event timing for a constant Hamiltonian and independently verified density, count and recoil observables. The count variance provides a useful result that a simple independent-kick approximation would not capture.

## what remains unresolved

The atom is fixed and the Hamiltonian is constant. This spectral event method cannot be reused unchanged for moving atoms or multiple optical frequencies. More importantly, tracking emitted momentum alone does not supply total motional diffusion: absorption/stimulated-force fluctuations, quantum correlations, recoil feedback and the chosen treatment of external motion must be justified. A classical trajectory driven by a conditional internal expectation value plus emission kicks is not automatically equivalent to a fully validated quantum motional model.

This is therefore a stochastic **validation foundation**, not a temperature, capture or apparatus-performance claim. The next scientifically difficult step is selecting and checking a consistent external-motion treatment, rather than simply adding more Monte Carlo samples or labeling emission bookkeeping as a full trap simulation.

## reproduce

```sh
OPENBLAS_NUM_THREADS=1 python benchmark/quantum_jumps.py
OPENBLAS_NUM_THREADS=1 python benchmark/quantum_count_moments.py
OPENBLAS_NUM_THREADS=1 python -m unittest benchmark/test_quantum_jumps.py
python benchmark/quantum_jump_report.py
```

Results include terminal states, counts, recoil vectors, exact reference arrays, numerical-propagator checks, independent counting moments and statistical diagnostics under `results/quantum_jumps/`. Seeds are `700000+i`, i=0…16383. Dependencies remain in `requirements-benchmark.txt`.
