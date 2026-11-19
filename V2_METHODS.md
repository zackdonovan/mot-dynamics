# validation follow-up

This follow-up is deliberately not another optimization run. It freezes the v1 seed-11 GP and Sobol settings and asks whether their apparent advantage survives changed assumptions. It is an independent portfolio study, unrelated to any proprietary apparatus or original internship code.

## additions

- A 36-condition matrix: six- and eight-beam configurations, standard/Sobol/GP settings, velocity-component standard deviations of 4/8/12 m/s, deterministic and approximate recoil dynamics. Each starts with the same 512 new synthetic positions and normalized velocities, from a seed not used in v1. Settings are never selected on these results. Threefold source-speed variation is a sensitivity study, not a measured distribution.
- Snapshot retention at 25 and 75 ms, including atoms lost and newly entering the threshold. Passing at both snapshots is **not** proof of continuous containment or asymptotic stability. Nine radius/speed thresholds expose dependence on the arbitrary success criterion.
- Frozen-field rate-generator spectral gaps diagnose population-relaxation times and their scale relative to source motion and turning magnetic axes. This only checks initial sampled states and frozen-rate convergence, not a coherent or dynamically coupled calculation.
- A reproducible near-zero-field direction test explicitly exposes a model failure, rather than smoothing it away.
- Exact requested integration duration (short final step), integer-state conversion, and finite input validation. Earlier runs used exactly divisible timesteps, so the duration fix does not change their specified integration intervals.

## recoil model

`simulate(..., recoil=True, seed=...)` adds a Gaussian momentum increment at each midpoint-force step. Its covariance per unit time is

`C = (hbar k)^2 [sum_b R_b khat_b khat_b^T + Gamma P_exc I/3]`,

where `R_b = sum_gq W_bgq (p_g + p_e)` counts stimulated absorption and emission. The force still uses `W (p_g - p_e)`. The two event directions have the same variance and opposite means. The second term approximates spontaneous-emission directions as isotropic. The resulting covariance is positive semidefinite, and its sampled kicks are checked against the requested covariance with 24,000 independent draws.

This is a **recoil sensitivity model**, not an exact diffusion calculation: independent-event statistics neglect internal-state correlations, Gaussian kicks approximate discrete events, emission is not isotropically exact, and no coherences/sub-Doppler mechanisms are restored. The exact force-noise covariance involves temporal correlations. The PyLCP authors discuss this distinction and their own diffusion approximation in section 2.1: https://pmc.ncbi.nlm.nih.gov/articles/PMC9890571/ . Our implementation does not reproduce PyLCP's diffusion prescription and has not been compared against it.

A passing covariance unit test establishes the sampler implements the declared approximation, not that the approximation is quantitatively correct for rubidium. Never infer reliable MOT temperature from these runs.

## uncertainty and interpretation

Wilson intervals refer to sampled atoms within a fixed condition. A single stochastic seed per condition does not establish repeated-seed stability. Common input samples make comparisons paired; bootstrap differences resample paired atom IDs. No intervals include physical-model error, source calibration error, or uncertainty in a real apparatus. Early validation is now exposed and should not be reused as an untouched test set for tuning.

The local quantization axis changes with position, while the model discards coherences. At effectively zero magnetic field this approximation becomes basis dependent. The test uses one fixed velocity and five approach directions at 1 nm, comparing force to the arbitrary +z fallback at the origin. It is a diagnostic, not a worst-case bound. Accurate near-center predictions require a coherent, fixed-basis treatment or an independently validated approximation. The current package does not implement that.

## reproduce

```sh
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s tests -v
OPENBLAS_NUM_THREADS=1 python validation/compare.py
OPENBLAS_NUM_THREADS=1 python validation/compare_eight.py
OPENBLAS_NUM_THREADS=1 python relaxation_diagnostic.py
OPENBLAS_NUM_THREADS=1 python zero_field_diagnostic.py
OPENBLAS_NUM_THREADS=1 python validation_study.py
python validation_report.py
```

All v1 results remain in `results/`; new evidence is segregated in `results/v2/`. Raw per-atom positions, velocities and masks are saved to support independent recomputation.

Three predeclared six-beam recoil cases also have a half-step (12.5 µs) distributional check in `recoil_timestep_check.py`. Random paths are not Brownian-coupled across timesteps, so changes mix integration effects and stochastic sampling; these are not strong-convergence estimates or definitive convergence evidence. Full deterministic v1 independent force checks still pass after refactoring, and a 32-atom deterministic v1/v2 replay is bit-identical.
