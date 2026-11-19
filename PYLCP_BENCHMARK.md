# external force benchmarks

The independent force calculation was checked against NIST-developed PyLCP 1.0.2. This checks implementation and approximation limits; it does **not** validate the trap's capture fractions, initial cloud, recoil model, or a real apparatus.

## 1. production 12-state rate model

Forty new randomized cases compare the existing F=2→F′=3 force with PyLCP rate equations: twenty six-beam and twenty eight-beam settings, with randomized detuning, intensity, helicity mixture, position and velocity. Gaussian intensity is evaluated at each test position, then held constant in PyLCP's local plane-wave calculation. Each helicity mixture is split into two weighted circular beams, which is valid for this incoherent **rate** comparison. It is not silently treated as a coherent OBE beam.

Maximum relative force difference is **9.03 × 10⁻¹⁴**. This corroborates the production implementation of the declared population model, including signs and normalization. Both implementations discard coherences, so agreement cannot establish adequacy near zero field.

## 2. simple F=0→F′=1 model with optical coherences

A separate four-state toy uses a stationary atom and one or two counterpropagating circular plane waves. There is no trajectory integration, gravity, magnetic-gradient force, or Gaussian envelope in this test.

- Decay rate Γ=1, wave number k=1 and ℏ=1; force is reported in ℏkΓ.
- Detuning Δ/Γ=−2 is placed on the excited Hamiltonian; laser frequency offsets are zero.
- For the pair, saturation per beam is 0.05 or 1; signed field is −2, −1, −0.1, 0, 0.1, 1 or 2. The complete configuration is rotated onto x, y and z, giving 42 cases.
- Fifteen single-beam checks use saturation 0.01, 0.1 or 1 and fields −3, −2, 0, 1 or 3. The exactly solvable cycling transition has force `s/[2(1+s+4(delta-beta)^2)]`.
- PyLCP integrates the optical Bloch equations for 100 Γ⁻¹ from the ground state, with relative tolerance 10⁻⁹ and absolute tolerance 10⁻¹¹. Force is checked at 90 and 100 Γ⁻¹; three selected pair cases are repeated at 200 Γ⁻¹.

The independent scalar rate formula agrees with PyLCP to **4.2 × 10⁻¹⁷** absolute force. Single-beam OBE agrees with the exact two-level formula to **3.4 × 10⁻¹⁰**. The largest steady-equation residual is **5.8 × 10⁻⁹**; the largest tail-force change is **7.1 × 10⁻¹⁰**. Density trace and positivity checks pass within numerical tolerance.

Rate equations and OBE need not agree exactly for the pair. At field 0.1, they differ by **0.30%** at saturation 0.05 and **7.05%** at saturation 1, relative to OBE force. The coherent calculation contains excited-state coherences absent from the population-only model. Agreement in a weakly excited limit is not proof of agreement at stronger intensity.

## 3. isolating the near-zero-field failure

A second controlled toy uses the production angular momenta F=2→F′=3 but **only one** stationary, circularly polarized plane-wave beam along +z, saturation 1 and detuning−2. It is not the six/eight-beam trap. It uses a fixed-z-basis density matrix, retaining all 144 complex entries and the same decay branches as PyLCP.

At zero field, optical pumping leads to the stretched cycling transition, with force 1/36=0.02777778 in normalized units. For a field of magnitude 10⁻⁸ in five different directions:

- The local-eigenbasis rate model yields forces between **0.01358981 and 0.02777778**.
- The fixed-basis OBE yields **0.02777777753–0.02777777778**.

The rate result has a direction spread equal to **51.1% of the zero-field reference**, while the coherent result approaches a common limit. This independently demonstrates that the basis sensitivity is an approximation issue, not just a coding error in our force function. It does **not** mean the full trap's capture score has a 51% error.

The coherent steady state is solved from the trace-constrained Liouvillian. Selected cases at zero field and field 0.1 alongx are independently evolved with PyLCP's OBE integrator for 12,000 Γ⁻¹ from equal ground populations. The largest density difference between the two procedures is **5.5 × 10⁻⁹**. Across 26 frozen-field cases, steady residuals are below 6 × 10⁻¹⁵; small negative eigenvalues are numerical roundoff below 3 × 10⁻¹³ in magnitude.

## conventions and numerical details

For the production comparison, positions and velocities are converted using `Bbar=mu_B B/(hbar Gamma)` and `vbar=k v/Gamma`; output force is multiplied by `hbar k Gamma`. Ground/excited Landé factors are 1/2 and 2/3. Magnetic-gradient forces and gravity are disabled in both implementations. Intensities and optical linewidth conventions match before comparisons are made.

PyLCP's OBE uses a coupling ΓE/4 with its beam-field normalization, whereas the public full-Hamiltonian builder uses E/2. For Γ=1 the independently assembled Hamiltonian therefore receives E/2. Its density vector is indexed column-major, while its exposed solution array is reshaped row-major; the physical density matrix is transposed for conventional commutator checks. The package itself was not modified. Single-beam analytic limits and actual time evolution check these conversions rather than assuming the API conventions match.

## implications

The production rate model is numerically well corroborated. It is still **not a high-fidelity coherent simulator**. A fixed-basis coherent solver is a promising next step, but the present benchmark has not replaced the trajectory model. Full multi-frequency beams, relative phases, realistic repumping, recoil correlations and experimental calibration remain separate problems. No proprietary Sygaldry code, data or geometry is used.

## reproduce

Install the regular requirements and `requirements-benchmark.txt` in a separate environment. Then, from the project root:

```sh
OPENBLAS_NUM_THREADS=1 python benchmark/pylcp_production_rates.py
OPENBLAS_NUM_THREADS=1 python benchmark/pylcp_toy.py
OPENBLAS_NUM_THREADS=1 python benchmark/pylcp_zero_field.py
OPENBLAS_NUM_THREADS=1 python benchmark/report.py
OPENBLAS_NUM_THREADS=1 python -m unittest discover -s benchmark -p 'test_*.py' -v
```

The exact environment is recorded in `results/pylcp/environment.json`. Raw results, residuals and controls are in `results/pylcp/`; the figure is `coherent_comparison.png`. These benchmark scripts are optional; the production simulator does not depend on PyLCP.

## primary sources

- Eckel et al., *PyLCP: A Python package for computing laser cooling physics*: https://pmc.ncbi.nlm.nih.gov/articles/PMC9890571/
- Authors' F=0→F′=1 OBE MOT example, including axis rotations: https://python-laser-cooling-physics.readthedocs.io/en/latest/examples/MOTs/02_F0_to_F1_1D_MOT_OBE.html
- PyLCP package 1.0.2: https://pypi.org/project/pylcp/1.0.2/

The documentation displays version 1.0.0; these runs used the pinned 1.0.2 release. The benchmark is our comparison against published software, not endorsement by its authors or NIST.
