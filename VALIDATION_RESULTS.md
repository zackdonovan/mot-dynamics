# validation results (v2)

Independent, synthetic follow-up. No laboratory validation; no new optimizer selection.

## main findings

- At the nominal 8 m/s source spread and 75 ms, approximate-recoil retention is 31.6% standard, 37.3% Sobol, and 41.6% GP in six beams; 26.8%, 28.3%, and 37.9% in eight beams. These are conditional simulated fractions, not experimental capture rates.
- Against the standard setting, the GP advantage in mean-force dynamics is 18.0–18.2 percentage points at this source spread; with approximate recoil it shrinks to 10.0–11.1 points. The model simplification changes the size of the claimed gain.
- At the slower 4 m/s source spread, the recoil-model GP-minus-standard difference is only 0.8 points (six beams) and 2.3 points (eight), with paired intervals spanning zero. There is no convincing advantage over standard settings in that condition.
- None of these results establish superiority of eight beams, true trap temperature, real loading efficiency, or hardware-ready settings.
- Halving the timestep in three six-beam recoil cases changes nominal standard retention from 31.6% to 31.6%, nominal GP from 41.6% to 43.9%, and slower-source GP from 69.3% to 72.5%. Paths are not Brownian-coupled and only one seed was used, so this is a distributional spot check, not a convergence proof.


## complete matrix

| geometry | setting | input velocity sigma (m/s) | dynamics | retained 25 ms | retained 75 ms | lost from 25 ms cohort |
|---|---|---:|---|---:|---:|---:|
| axial6 | standard | 4 | mean force | 358/512 (69.9%) | 358/512 (69.9%) | 0 |
| axial6 | standard | 4 | recoil approximation | 349/512 (68.2%) | 351/512 (68.6%) | 12 |
| axial6 | standard | 8 | mean force | 167/512 (32.6%) | 167/512 (32.6%) | 0 |
| axial6 | standard | 8 | recoil approximation | 160/512 (31.2%) | 162/512 (31.6%) | 4 |
| axial6 | standard | 12 | mean force | 70/512 (13.7%) | 70/512 (13.7%) | 0 |
| axial6 | standard | 12 | recoil approximation | 69/512 (13.5%) | 66/512 (12.9%) | 3 |
| axial6 | sobol_11 | 4 | mean force | 379/512 (74.0%) | 379/512 (74.0%) | 0 |
| axial6 | sobol_11 | 4 | recoil approximation | 302/512 (59.0%) | 318/512 (62.1%) | 50 |
| axial6 | sobol_11 | 8 | mean force | 227/512 (44.3%) | 227/512 (44.3%) | 0 |
| axial6 | sobol_11 | 8 | recoil approximation | 179/512 (35.0%) | 191/512 (37.3%) | 31 |
| axial6 | sobol_11 | 12 | mean force | 113/512 (22.1%) | 113/512 (22.1%) | 0 |
| axial6 | sobol_11 | 12 | recoil approximation | 86/512 (16.8%) | 94/512 (18.4%) | 16 |
| axial6 | gp_ucb_11 | 4 | mean force | 427/512 (83.4%) | 428/512 (83.6%) | 0 |
| axial6 | gp_ucb_11 | 4 | recoil approximation | 350/512 (68.4%) | 355/512 (69.3%) | 62 |
| axial6 | gp_ucb_11 | 8 | mean force | 258/512 (50.4%) | 259/512 (50.6%) | 0 |
| axial6 | gp_ucb_11 | 8 | recoil approximation | 204/512 (39.8%) | 213/512 (41.6%) | 36 |
| axial6 | gp_ucb_11 | 12 | mean force | 134/512 (26.2%) | 134/512 (26.2%) | 0 |
| axial6 | gp_ucb_11 | 12 | recoil approximation | 104/512 (20.3%) | 110/512 (21.5%) | 20 |
| oblique8 | standard | 4 | mean force | 322/512 (62.9%) | 324/512 (63.3%) | 0 |
| oblique8 | standard | 4 | recoil approximation | 310/512 (60.5%) | 311/512 (60.7%) | 11 |
| oblique8 | standard | 8 | mean force | 143/512 (27.9%) | 144/512 (28.1%) | 0 |
| oblique8 | standard | 8 | recoil approximation | 141/512 (27.5%) | 137/512 (26.8%) | 9 |
| oblique8 | standard | 12 | mean force | 60/512 (11.7%) | 60/512 (11.7%) | 0 |
| oblique8 | standard | 12 | recoil approximation | 59/512 (11.5%) | 58/512 (11.3%) | 2 |
| oblique8 | sobol_11 | 4 | mean force | 349/512 (68.2%) | 360/512 (70.3%) | 0 |
| oblique8 | sobol_11 | 4 | recoil approximation | 213/512 (41.6%) | 259/512 (50.6%) | 61 |
| oblique8 | sobol_11 | 8 | mean force | 207/512 (40.4%) | 212/512 (41.4%) | 0 |
| oblique8 | sobol_11 | 8 | recoil approximation | 123/512 (24.0%) | 145/512 (28.3%) | 40 |
| oblique8 | sobol_11 | 12 | mean force | 96/512 (18.8%) | 98/512 (19.1%) | 0 |
| oblique8 | sobol_11 | 12 | recoil approximation | 54/512 (10.5%) | 66/512 (12.9%) | 15 |
| oblique8 | gp_ucb_11 | 4 | mean force | 396/512 (77.3%) | 398/512 (77.7%) | 0 |
| oblique8 | gp_ucb_11 | 4 | recoil approximation | 274/512 (53.5%) | 323/512 (63.1%) | 51 |
| oblique8 | gp_ucb_11 | 8 | mean force | 236/512 (46.1%) | 237/512 (46.3%) | 0 |
| oblique8 | gp_ucb_11 | 8 | recoil approximation | 165/512 (32.2%) | 194/512 (37.9%) | 33 |
| oblique8 | gp_ucb_11 | 12 | mean force | 116/512 (22.7%) | 117/512 (22.9%) | 0 |
| oblique8 | gp_ucb_11 | 12 | recoil approximation | 81/512 (15.8%) | 92/512 (18.0%) | 20 |

## paired comparisons at 75 ms

95% percentile bootstrap intervals resample paired atoms (5,000 draws); conditional on these controls and approximate physics. No optimization-seed or physical-model uncertainty.

| geometry | source sigma | dynamics | GP minus standard (percentage points) | GP minus Sobol (percentage points) |
|---|---:|---|---|---|
| axial6 | 4 | mean force | 13.7 [10.7, 16.6] | 9.6 [7.0, 12.3] |
| axial6 | 4 | recoil approximation | 0.8 [-3.3, 4.5] | 7.2 [3.3, 11.1] |
| axial6 | 8 | mean force | 18.0 [14.8, 21.5] | 6.2 [4.1, 8.4] |
| axial6 | 8 | recoil approximation | 10.0 [6.1, 13.9] | 4.3 [1.2, 7.4] |
| axial6 | 12 | mean force | 12.5 [9.6, 15.4] | 4.1 [2.5, 5.9] |
| axial6 | 12 | recoil approximation | 8.6 [5.5, 11.5] | 3.1 [0.8, 5.5] |
| oblique8 | 4 | mean force | 14.5 [11.3, 17.8] | 7.4 [4.9, 10.2] |
| oblique8 | 4 | recoil approximation | 2.3 [-2.0, 6.6] | 12.5 [8.2, 16.8] |
| oblique8 | 8 | mean force | 18.2 [14.8, 21.7] | 4.9 [2.7, 7.2] |
| oblique8 | 8 | recoil approximation | 11.1 [7.0, 15.2] | 9.6 [6.1, 13.1] |
| oblique8 | 12 | mean force | 11.1 [8.4, 13.9] | 3.7 [2.0, 5.7] |
| oblique8 | 12 | recoil approximation | 6.6 [3.7, 9.8] | 5.1 [2.3, 7.8] |

## model limitation exposed

At a fixed velocity (0.1, 0.2, 0.3) m/s, moving the position only 1 nm away from the origin along different directions changes the computed force by up to 33.7% for the six-beam finalists and 15.7% for the eight-beam finalists, relative to the origin fallback. This is a local-axis/incoherent-population approximation failure at vanishing field. It prevents claiming high-fidelity near-center physics; more samples or better optimization do not resolve it.

Initial-cloud population-relaxation medians span 2.3–7.1 microseconds; the slowest sampled case is 425 microseconds. In the declared initial cloud, the 95th-percentile product of relaxation time and magnetic-axis turning rate is below 0.023. That supports timescale separation for many initial states, not all trajectories or near-zero-field coherences.

See V2_METHODS.md for the independent-event recoil approximation and uncertainty scope. Snapshot retention is not proof of continuous or long-term confinement.
