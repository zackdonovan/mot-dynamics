# recorded simulation results

These are synthetic finite-time retention outcomes, not experimental capture probabilities. Intervals describe atom sampling only. All optimizer finalists are reported; no test-set winner selection.

## six beams

| setting | retained / 4096 | fraction | 95% interval |
|---|---:|---:|---:|
| standard | 1487 | 36.304% | 34.844%–37.789% |
| sobol_11 | 1933 | 47.192% | 45.667%–48.723% |
| gp_ucb_11 | 2208 | 53.906% | 52.377%–55.428% |
| sobol_23 | 2011 | 49.097% | 47.567%–50.628% |
| gp_ucb_23 | 2173 | 53.052% | 51.521%–54.577% |
| sobol_37 | 2024 | 49.414% | 47.884%–50.945% |
| gp_ucb_37 | 2088 | 50.977% | 49.445%–52.506% |

Surrogate held-out MAE: 0.02317; constant predictor MAE: 0.06520. Tested on 32 new parameter settings, using the training atom cloud so this isolates parameter-space interpolation rather than atom generalization.

| perturbation | standard | GP seed 11 |
|---|---:|---:|
| nominal | 36.304% | 53.906% |
| gradient_low | 46.069% | 64.722% |
| gradient_high | 28.760% | 44.312% |
| power_low | 31.934% | 49.731% |
| power_high | 40.161% | 57.251% |

| step (µs) | retained fraction (256 new atoms) | smooth objective |
|---|---:|---:|
| 200 | 56.641% | 0.553002 |
| 100 | 56.641% | 0.552483 |
| 50 | 56.641% | 0.552381 |
| 25 | 56.641% | 0.552354 |

100µs vs 25µs final-position discrepancy: median 0.0048mm, 95th percentile 2.6559mm. Trajectory differences and threshold crossings should not be hidden by aggregate scores.

### full holdout timestep confirmation

Fixed comparison at 25µs, on the same 4096 held-out atoms. This refines the reported numerical outcome; it is not an extra optimization pass.

- standard: 100µs 36.304%, 25µs 36.450%; 8 atom labels changed. Fine-retained atoms: 95th-percentile position discrepancy 0.000002mm.
- gp_ucb_11: 100µs 53.906%, 25µs 54.053%; 6 atom labels changed. Fine-retained atoms: 95th-percentile position discrepancy 0.006700mm.

The larger all-atom position discrepancies primarily affect escaped trajectories. Aggregate retention is much more stable, but individual escaped paths should not be treated as numerically converged.

## eight beams

| setting | retained / 4096 | fraction | 95% interval |
|---|---:|---:|---:|
| standard | 1272 | 31.055% | 29.656%–32.489% |
| sobol_11 | 1820 | 44.434% | 42.918%–45.960% |
| gp_ucb_11 | 2031 | 49.585% | 48.055%–51.116% |
| sobol_23 | 1959 | 47.827% | 46.300%–49.358% |
| gp_ucb_23 | 2000 | 48.828% | 47.299%–50.359% |
| sobol_37 | 1713 | 41.821% | 40.319%–43.339% |
| gp_ucb_37 | 1818 | 44.385% | 42.869%–45.911% |

Surrogate held-out MAE: 0.08066; constant predictor MAE: 0.09976. Tested on 32 new parameter settings, using the training atom cloud so this isolates parameter-space interpolation rather than atom generalization.

| perturbation | standard | GP seed 11 |
|---|---:|---:|
| nominal | 31.055% | 49.585% |
| gradient_low | 41.553% | 59.985% |
| gradient_high | 24.146% | 39.258% |
| power_low | 27.515% | 42.188% |
| power_high | 34.082% | 52.686% |

| step (µs) | retained fraction (256 new atoms) | smooth objective |
|---|---:|---:|
| 200 | 51.953% | 0.484352 |
| 100 | 51.953% | 0.483731 |
| 50 | 51.953% | 0.486190 |
| 25 | 51.953% | 0.486318 |

100µs vs 25µs final-position discrepancy: median 0.0073mm, 95th percentile 2.1962mm. Trajectory differences and threshold crossings should not be hidden by aggregate scores.

### full holdout timestep confirmation

Fixed comparison at 25µs, on the same 4096 held-out atoms. This refines the reported numerical outcome; it is not an extra optimization pass.

- standard: 100µs 31.055%, 25µs 31.128%; 5 atom labels changed. Fine-retained atoms: 95th-percentile position discrepancy 0.000013mm.
- gp_ucb_11: 100µs 49.585%, 25µs 49.658%; 7 atom labels changed. Fine-retained atoms: 95th-percentile position discrepancy 0.014645mm.

The larger all-atom position discrepancies primarily affect escaped trajectories. Aggregate retention is much more stable, but individual escaped paths should not be treated as numerically converged.

Search plots show the mean best training score over three seeds; shading spans the seed minimum and maximum, not a confidence interval.

## interpretation

The learned model predicts simulator objective values from laser controls. Its predictions guide which expensive simulation to run next; the simulator itself is not a trained neural network. The evidence here is limited to three optimizer seeds and one declared cloud distribution. Gains over the standard setting include simple parameter tuning; gains over equal-budget Sobol search isolate the value of the surrogate more directly.

The eight-beam system is not automatically superior: fixed total power is divided across more pairs, and the extra oblique pair is not optimized in geometry. This is a useful negative engineering finding, not evidence against all eight-beam MOT designs.

Next scientifically meaningful work: recoil-diffusion trajectories, nonstationary internal dynamics, coherence-aware benchmark cases, calibrated source distributions and beam profiles, then measurements. Arbitrarily expanding the optimizer or claiming production use would not resolve these physical-model limits.
