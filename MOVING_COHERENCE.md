# moving atoms and frequency-offset beams

## question and scope

Does retaining time-dependent optical coherence along a prescribed path reproduce the population-only force? Does assigning separate optical frequencies to different beam pairs suppress the large static-phase effects seen in v4?

This is a controlled local force study, not a trap-loading calculation. It integrates the full 12-state density matrix along **prescribed straight paths**, with no force feedback on velocity. The fixed field B=(0.1,0,0), plane-wave intensities, closed Rb87 manifold, total incident saturation 6, and pure circular helicities match the declared multibeam benchmark. Gaussian envelopes, changing quadrupole fields, recoil, gravity and stochastic laser linewidth are absent.

Two explicitly specified laser architectures are compared:

- common optical frequency, detuning Δ/Γ=−2 for every beam;
- counterpropagating pairs with deterministic frequency offsets: (0,0.12,0.29)Γ for six beams, plus 0.47Γ for the oblique pair in eight beams. Each pair remains mutually coherent. Pair-relative phases drift deterministically; this is not phase diffusion or incoherent illumination.

Frequency offsets alter both interference dynamics and individual detunings. The rate comparator therefore uses each actual detuning −2+offset and the same Doppler shifts. Differences between architectures are not attributed solely to coherence while ignoring their changed spectra. These offsets are a synthetic, reproducible laser architecture, not recovered laboratory settings.

## implementation and validation

Positions follow r(t)=v t. The optical field contains phases exp[i(δ_b−k_b·v)t−iφ_b], with the package's documented field convention. A cached Hamiltonian integrates the Lindblad equation without repeatedly rebuilding atomic structure. Its time dependence is checked against PyLCP's public Hamiltonian, and four selected integrations are compared directly with PyLCP's public density evolution over 30 Γ⁻¹. No local steady-state approximation is used during motion.

The main matrix has 24 runs: two geometries, two initial phase realizations, two frequency architectures, and prescribed vx/(Γ/k)=0,0.05,0.2. Two phases are probes, **not a statistically converged phase ensemble**. Main runs last 2000 Γ⁻¹, approximately 52.5 μs using the project's Rb87 constants. Speeds 0.05 and 0.2 correspond to 0.237 and 0.946 m/s. Forces are averaged over the second half of each run. Full force time series are saved.

Adaptive DOP853 uses rtol=2×10⁻⁸, atol=2×10⁻¹⁰. Force quadrature samples every 2 Γ⁻¹. Four selected runs are repeated with 100× tighter solver tolerances; four with 4× denser force samples; four with duration 8000 Γ⁻¹ (210 μs). Two stationary cases are compared with the direct steady-state solver. Additional unit tests verify Doppler equivalence: moving through plane waves produces the same internal dynamics as a stationary atom with the corresponding per-beam Doppler frequency offsets in this uniform-field problem.

Measured checks:

- Maximum public PyLCP density difference: 1.83×10⁻¹¹.
- Maximum Hamiltonian reconstruction difference: 4.39×10⁻¹⁶.
- Maximum density trace error: 1.20×10⁻¹⁴; no negative density eigenvalues at saved samples in the main runs.
- Tightening ODE tolerances changes the selected mean-force vectors by less than 1.3×10⁻¹¹.
- Denser force quadrature changes them by at most 1.38×10⁻⁵.
- Longer averaging changes them by 2.48×10⁻⁴ to 2.85×10⁻³, much larger than ODE error. Finite-window averaging, not integrator tolerance, is the dominant numerical limitation here.
- Stationary finite-window means differ from direct steady forces by up to 4.91×10⁻⁶.

## results

At vx/(Γ/k)=0.05, phase seed 721, the second-half x-force averages are:

| geometry | common-frequency coherent | common-frequency rates | offset-pair coherent | offset-pair rates |
|---|---:|---:|---:|---:|
| six beams | −0.02341 | −0.00453 | −0.00128 | −0.00430 |
| eight beams | −0.02532 | −0.00452 | −0.00322 | −0.00536 |

Force units are ℏkΓ. These are selected finite-window results, not universal improvement factors or population estimates. The common-frequency beams retain large phase-sensitive transverse forces. Frequency offsets reduce these effects in the tested cases, but **do not make the coherent result equal to the rate model**. Eight-beam common-frequency evolution also shows slow force beating: the two halves of its averaging interval differ in full force vector by up to 0.0187. Calling these finite-window averages a converged equilibrium force would be incorrect.

The more important outcome is identifying the modeling assumptions that must be specified before claiming predictive trap dynamics: optical frequency architecture, mutual coherence, initial phase handling, motion and averaging timescale.

## why these are not capture trajectories

The paths are imposed to measure velocity-dependent internal response. They are not solutions of Newton's equation. This distinction matters numerically: the dimensionless recoil-to-motional scale is ħk²/(mΓ)=0.001244. A force of 0.05 sustained for 2000 Γ⁻¹ would change v/(Γ/k) by about 0.124, larger than the imposed 0.05 velocity. Thus the fixed-path assumption cannot simply be relabeled as a self-consistent free trajectory, even over a short physical duration.

A next step would integrate coupled density, position and velocity for a small, explicitly chosen laser architecture, checking timestep/force feedback and recoil assumptions before any capture comparison. The full thermal cloud optimizer remains the earlier approximate population model. No experimental or company-specific validation is claimed.

## reproduce

```sh
OPENBLAS_NUM_THREADS=1 python benchmark/moving_coherence.py
OPENBLAS_NUM_THREADS=1 python benchmark/moving_checks.py
OPENBLAS_NUM_THREADS=1 python benchmark/moving_quadrature.py
OPENBLAS_NUM_THREADS=1 python -m unittest benchmark/test_moving.py
python benchmark/moving_report.py
```

Raw data and summary checks: `results/moving/`. Optional benchmark dependencies remain isolated in `requirements-benchmark.txt`. Atomic conventions follow the prior [PyLCP comparison](PYLCP_BENCHMARK.md) and [PyLCP framework paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC9890571/).
