# mot design study

A new, independently implemented numerical study of laser-cooling controls for a simplified rubidium-87 magneto-optical trap. It is not recovered internship code, a reproduction of a company system, or evidence of laboratory deployment. Geometry and input distributions here are declared synthetic choices.

## photon-resolved validation (v10)

`QUANTUM_JUMP_VALIDATION.md` adds 16,384 exact-event fixed-atom quantum trajectories, verified against density evolution, independent photon-count moments and spontaneous-recoil covariance. Non-Poisson photon statistics are resolved. This validates the internal stochastic foundation, not total moving-atom diffusion or capture. All 33 tests pass.

## multi-fidelity pilot (v9)

`MULTIFIDELITY_PILOT.md` compares four locked detuning candidates across inexpensive population dynamics and 64 coherent evaluations with separate discovery/held-out phases, plus 16 condition-shift probes. The cheap model retains the best tested candidate but miscalibrates scores and does not support transferring results across positions. This is a model-dependent screening study, not validated cooling or capture performance. All 30 tests pass.

## recoil adequacy and compute limits (v8)

`RECOIL_VALIDITY.md` extends selected spatial paths to 210 μs and derives a dipole-anisotropic spontaneous recoil budget. The omitted free-recoil scale is already material, and measured runtime makes brute-force long coherent optimization impractical. This establishes a quantitative modeling limit rather than claiming an unvalidated stochastic capture result. All 29 tests pass.

## spatial dynamics (v7)

`SPATIAL_DYNAMICS.md` adds evolving Gaussian envelopes, quadrupole Zeeman fields and their optical/magnetic force gradients. Four full-versus-locally-frozen trajectory pairs quantify the local approximation over short paths. Analytic-gradient, limiting-case and public spatial-motion checks pass; all 26 tests pass. Recoil diffusion and experimental capture validation remain out of scope.

## coupled internal/mechanical dynamics (v6)

`COUPLED_MOTION.md` closes the force-feedback loop: density matrix, position and velocity evolve together in eight controlled local trajectories. Ballistic and infinite-mass limits, direct PyLCP motion comparisons, force-gradient checks and tighter-tolerance runs validate implementation. These deterministic mean-force paths still omit recoil diffusion and spatial trap profiles, so they are not capture predictions. The full suite has 24 passing tests.

## moving-atom force study (v5)

`MOVING_COHERENCE.md` adds time-dependent density evolution on prescribed straight paths, comparing common-frequency and pair-offset lasers. Twenty-four main runs plus duration, tolerance and quadrature checks show that optical frequency architecture materially changes forces, while rate equations remain an approximation. These are controlled force probes, not self-consistent capture trajectories. The complete suite now has 21 passing tests.

## multibeam coherence (v4)

`MULTIBEAM_COHERENCE.md` extends the coherent calculation to the actual six/eight beam directions, with 1,024 stationary density-matrix solutions and explicit relative-phase ensembles. Phase-dependent forces and field responses differ materially from population-only predictions. This is a controlled approximation study, not a new capture result or an experimental claim. All 19 tests pass (14 simulator, three original external benchmarks, two multibeam tests).

## external comparison (v3)

`PYLCP_BENCHMARK.md` adds 40 production-force comparisons against PyLCP's rate equations and separate coherent optical-Bloch benchmarks. Production rates agree to numerical precision, while a controlled single-beam limit independently confirms the near-zero-field approximation failure. This does not validate capture or fix the full trajectory model. Optional benchmark dependencies are isolated from the simulator.

## validation follow-up (v2)

Start with `VALIDATION_RESULTS.md` and `V2_METHODS.md`. The follow-up freezes prior controls, tests three incoming velocity spreads, extends trajectories to 75 ms, and introduces optional approximate recoil noise. It exposes a 15–34% near-zero-field force sensitivity to the arbitrary local quantization basis. This is a substantive limitation of the approximation, not a validated apparatus simulator. Original results below are retained for provenance, not superseding the follow-up.

## question

Can a learned surrogate find better laser settings than a space-filling search under the same simulation budget? Does an apparent improvement survive unseen atoms, changed power/field gradient and a smaller integration timestep?

We study a conventional six-beam geometry and an independent eight-beam geometry: three axial counterpropagating pairs plus one pair along ±(1,1,1)/√3. This eighth-beam extension is a generic design, not an inferred proprietary layout. Total on-axis saturation parameter is six in both cases. Extra beams do not imply better trapping.

## model

The internal model has five ground populations (F=2, m=-2…2) and seven excited populations (F'=3, m'=-3…3). Only Δm=q∈{-1,0,+1} transitions are allowed. Clebsch–Gordan strengths weight each transition. At each trajectory step:

1. Evaluate Gaussian beam intensity at the atom and the quadrupole field B=b(x,y,-2z).
2. Decompose each beam's polarization into q components relative to the local field.
3. Shift detuning by Doppler and differential Zeeman shifts.
4. Build stimulated absorption/emission and spontaneous-decay rates. Solve the normalized stationary 12-state population system.
5. Sum net photon momenta, add gravity, and advance position/velocity with explicit midpoint integration.

For a transition, W=(Γ/2)s C² w_q/[1+4δ²/Γ²], with δ=Δ-k·v-(μ_B/ℏ)|B|(g_e m_e-g_g m_g). Stimulated rates transfer populations in both directions; decay transfers excited population to ground states with rate ΓC². Saturation emerges from population balance: no extra saturation denominator is inserted.

The force is Σ_beams ℏk Σ_transitions W(p_g-p_e). Ground and excited Landé factors are approximated as 1/2 and 2/3. Constants use the 780nm rubidium D2 cycling transition.

## what the polarization control means

Each pair has a helicity-mixture parameter h, not a general coherent polarization ellipse. For c=k̂·B̂, weights are w_±=(1+c²±2hc)/4 and w_0=(1-c²)/2. Axial pair signs are +,+,- for restoring forces under this field and detuning convention. Opposite circular components are treated as an incoherent mixture; phase-sensitive polarization interference is absent. At B≈0 the quantization axis defaults to +z.

## experiment protocol

- Incoming cloud: independent Gaussian coordinates with 6mm standard deviation and Gaussian velocity components with 8m/s standard deviation. These are assumptions, not a measured source distribution or thermal vapor flux.
- 25ms evolution; default step 100µs. Success proxy: final radius <2mm **and** speed <0.5m/s. This is finite-time retention, not a measured loading rate or guarantee of long-term capture.
- Search objective: average exp[-|r/3mm|²-|v/(1m/s)|²], a smooth alternative to a discontinuous count.
- Fixed training cloud: 192 atoms. Each search gets 48 simulator evaluations, with 16 identical starting points. Three optimizer seeds.
- Baselines: a declared standard setting and scrambled Sobol space-filling search. Learned method: a Gaussian-process surrogate with a Matérn kernel, choosing from 1024 candidate controls using mean+1.5×standard deviation. This is Bayesian-style black-box optimization, not gradient descent through atom dynamics.
- Six-beam controls: three pair detunings and three helicity mixtures. Eight-beam controls: four pair detunings and four mixtures. Power is fixed, preventing gains from silently increasing the resource budget.
- Final validation: 4096 new atoms, all six search finalists reported, no winner chosen using validation outcomes. Wilson intervals quantify finite atom-sampling uncertainty only.
- Separate surrogate assessment: 32 held-out parameter vectors. Compare regression MAE against a constant-prediction baseline.
- Perturbation study:±20% power and field gradient; seed 11 GP candidate specified before evaluation. Timestep refinement uses a further independent 256-atom cloud.

## run

Python3.12 was used. A single BLAS thread avoids oversubscription.

```sh
python -m pip install -r requirements.txt
python -m unittest discover -s tests -v
OPENBLAS_NUM_THREADS=1 python validation/compare.py
OPENBLAS_NUM_THREADS=1 python validation/compare_eight.py
OPENBLAS_NUM_THREADS=1 python study.py
OPENBLAS_NUM_THREADS=1 python eight_beam.py
OPENBLAS_NUM_THREADS=1 python confirm_timestep.py
python report.py
python visualize.py
```

Results are saved as JSON so claims can be inspected without rerunning. See `WALKTHROUGH.md` for a conceptual tour and `RESULTS.md`, `results/`, and `validation/AUDIT.md`. No external training dataset is required: simulator queries generate the surrogate's training targets.

## limits that matter

- This is a population rate model, not a coherent optical-Bloch model. It neglects interference and sub-Doppler mechanisms. V1 omits recoil; v2 optionally adds an independent-event recoil approximation, not a validated quantum diffusion calculation. It cannot predict a reliable final temperature.
- Internal populations relax instantaneously to a stationary distribution. This can fail during weak illumination or rapid state changes.
- The F=2→F'=3 manifold is closed. Off-resonant hyperfine leakage and a real repump system are omitted. Magnetic dipole-gradient forces are also omitted.
- If total stimulated rate falls below 10⁻⁸Γ, populations are assigned equal ground weights instead of solving an effectively singular dark-state problem. Momentum exchange there is negligible on this study's timescale; no artificial capture boundary is imposed.
- Small statistical error bars do not bound physical-model error. Threshold choice, source distribution and force approximations can dominate.
- Better simulation settings do not establish better experimental capture. Real apparatus calibration, beam measurements and independent experimental data would be needed before suggesting hardware settings.

## sources and provenance

Equations were independently implemented, not copied from a third-party codebase. The independent audit uses Wigner coefficients and a different stationary solver to cross-check the production implementation.

- Daniel Steck, [Rubidium87 D Line Data](https://steck.us/alkalidata/rubidium87numbers.pdf): atomic constants, transition strengths and angular momentum structure.
- Eckel et al., [PyLCP: A Python package for computing laser cooling physics](https://pmc.ncbi.nlm.nih.gov/articles/PMC9890571/): optical-Bloch, rate-equation and heuristic modeling hierarchy and limitations.
- [PyLCP rate-equation documentation](https://python-laser-cooling-physics.readthedocs.io/en/latest/_modules/pylcp/rateeq.html): independent methodological reference for stimulated and spontaneous rates.

Referenced papers and data remain subject to their own terms; they are linked, not bundled.
