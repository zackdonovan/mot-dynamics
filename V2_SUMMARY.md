# validation results (v2)

Independent, synthetic follow-up. No laboratory validation; no new optimizer selection.

## main findings

- At the nominal 8 m/s source spread and 75 ms, approximate-recoil retention is 31.6% standard, 37.3% Sobol, and 41.6% GP in six beams; 26.8%, 28.3%, and 37.9% in eight beams. These are conditional simulated fractions, not experimental capture rates.
- Against the standard setting, the GP advantage in mean-force dynamics is 18.0–18.2 percentage points at this source spread; with approximate recoil it shrinks to 10.0–11.1 points. The model simplification changes the size of the claimed gain.
- At the slower 4 m/s source spread, the recoil-model GP-minus-standard difference is only 0.8 points (six beams) and 2.3 points (eight), with paired intervals spanning zero. There is no convincing advantage over standard settings in that condition.
- None of these results establish superiority of eight beams, true trap temperature, real loading efficiency, or hardware-ready settings.
- Halving the timestep in three six-beam recoil cases changes nominal standard retention from 31.6% to 31.6%, nominal GP from 41.6% to 43.9%, and slower-source GP from 69.3% to 72.5%. Paths are not Brownian-coupled and only one seed was used, so this is a distributional spot check, not a convergence proof.

