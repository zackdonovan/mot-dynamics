# multi-fidelity candidate-selection pilot

## question

Can the inexpensive instantaneous population model screen laser settings without missing a candidate preferred by coherent dynamics? How much do the scores and rankings depend on optical phase and initial position?

This is a **pilot comparison**, not a new validated optimizer. The original surrogate-search studies remain separate. We do not reuse their capture metric here because the coherent model has not established a validated stochastic capture prediction.

## declared experiment

Four base detunings, −1, −2, −3 and −4 Γ, share six beams, total incident saturation 6, pure helicities and pair-frequency offsets (0,0.12,0.29)Γ. The initial position is (6,3,2) mm, using the declared 12 mm Gaussian waist and 0.12 T/m quadrupole gradient. Both models freeze those local intensities and magnetic field, but update velocity. The coherent model additionally tracks density dynamics and optical phase as position evolves. Profile-gradient and magnetic-gradient forces are absent in both frozen models. There is no recoil diffusion in either.

Initial velocity is (0.05,0,0)Γ/k. The horizon is 1000 Γ⁻¹, about 26.2 μs. The score is:

```
J = 1 − |v_final|² / |v_initial|²
```

Larger is better; negative values mean kinetic energy increased. Using all three velocity components avoids rewarding a longitudinal reversal or transverse acceleration as if it were cooling.

Candidate settings and eight discovery phase seeds (810–817) and eight held-out phase seeds (910–917) were fixed before their coherent results were inspected. Each phase seed is shared across candidates, enabling paired comparisons. Four rate-model evaluations are compared with 64 coherent runs. Candidate selection uses only discovery means; the held-out batch evaluates that selected candidate without retuning it.

A separate 16-run stress probe uses four new phase seeds (1010–1013) at the field origin, with all four candidates. This tests transfer to another condition, not a further optimization round. It was specified separately and is not used to change the discovery selection.

## statistical interpretation

The report includes paired bootstrap intervals from 20,000 resamples of phase seeds. These describe finite phase sampling under this deterministic model, not experimental confidence, model accuracy or recoil uncertainty. Eight seeds per main batch and four per stress condition are small samples. Bootstrap winner frequencies are sampling-stability diagnostics, not posterior probabilities that a setting is physically optimal.

## measured results

| base detuning / Γ | population-model score | coherent discovery mean | coherent held-out mean |
|---|---:|---:|---:|
| −1 | 0.881 | 0.942 | 0.874 |
| −2 | −5.028 | −1.640 | −1.685 |
| −3 | 0.155 | 0.582 | 0.601 |
| −4 | 0.515 | 0.430 | 0.447 |

Both models select −1 among these four candidates. The coherent ordering is −1, −3, −4, −2 in both phase batches; the population model swaps the middle two. Spearman rank correlation is 0.8, based on only four candidates. The held-out selected-candidate mean is 0.874 with a phase-bootstrap interval [0.815,0.930]. Its paired bootstrap advantage over the best other candidate is [0.183,0.363] score units.

The cheap model therefore identifies the best tested candidate **under this conditional objective**, but its numerical scores are not calibrated to the coherent model. A retrospective rate-model top-two shortlist would retain the winner but discard the coherent runner-up. This is a screening illustration, not evidence of general search guarantees or a realized optimizer speedup.

The origin stress probe exposes the transfer limit. The off-axis selected candidate −1 has coherent mean −0.0037 at the origin across four phases, with individual scores from −0.179 to +0.155, while the population model predicts +0.205. The stress sample is too small and phase-dependent to establish a reliable new optimum, but it does not support transferring the off-axis result into a universal setting recommendation.

The main coherent evaluations consumed 219 seconds in total; the four population evaluations took 0.066 seconds. A median single coherent phase evaluation was approximately 211× slower than one population evaluation in this implementation. These wall timings include environment/workload effects and are not portable speed guarantees; estimating a phase-averaged coherent score also requires multiple such evaluations.

## conclusion

The experiment supports using a cheap model to propose a shortlist, **then checking rankings and failure modes under richer physics and changed conditions**. It does not establish an experimentally optimal detuning, a transferable cooling efficiency, or a new capture performance result. The strongest result is the combination of partial ranking agreement, large score miscalibration, condition dependence and a quantified omitted-recoil scale.

Further deterministic parameter sweeps would add precision to a conditional objective without removing its main physical limitation. A justified stochastic model and a meaningful distribution of initial conditions are necessary before elevating this into a stronger performance claim.

## scope and physical limit

The detailed model is not automatically ground truth. It retains coherences but still lacks stochastic momentum evolution, experimental source distributions and measured apparatus parameters. The short-horizon score can favor a force that happens to bring a particular initial velocity near zero at one time, rather than a setting that maximizes trap loading or sustained cooling.

For the selected off-axis candidate, the spontaneous free-recoil noise budget is about 0.075–0.078 m/s in three-dimensional RMS norm. That is comparable to, or larger than, many of its deterministic residual speeds. As established in v8, this budget is not a final spread prediction or a lower bound; it identifies a missing physical effect that can change the objective and potentially the ranking.

Thus a good held-out phase score does not license a claim of physical cooling efficiency. Phase holdout does not test recoil, initial-condition distributions, geometry, laser noise or actual hardware.

## reproduce

```sh
OPENBLAS_NUM_THREADS=1 python benchmark/multifidelity_study.py
OPENBLAS_NUM_THREADS=1 python benchmark/multifidelity_stress.py
python benchmark/multifidelity_report.py
```

The protocols, all per-run final velocities, recoil diagnostics, timing, bootstrap summaries and plot are in `results/multifidelity/`. The project's 30 tests include a new base-detuning Hamiltonian check. No changes were made to public repositories or physical hardware, and no company-specific or experimental validation is claimed.
