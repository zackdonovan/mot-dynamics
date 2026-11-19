# independent physics audit

## checked

- 15 absorption strengths independently computed with SymPy Wigner Clebsch–Gordan coefficients. Agreement within floating precision; all seven excited-state spontaneous-decay branch sums equal one.
- 30 randomized positions, velocities, detunings, intensities and helicity mixtures compared with a separately written scalar generator and SVD nullspace population solver. Largest relative force difference: 3.73 × 10^-14. Populations agree to absolute tolerance 2 × 10^-12.
- Radiation pressure has restoring signs along all three axes at 1 mm displacement and damping signs at 1 m/s.
- 100 randomized simultaneous position/velocity inversions yield maximum |F(r,v)+F(-r,-v)| = 1.58 × 10^-35 N, with gravity disabled.
- Light-off test returns gravitational force only.

These are implementation checks within the stated rate-equation approximation, not experimental validation or independent validation of all physical assumptions.

An additional 20 randomized eight-beam force comparisons passed against the independent reference; maximum relative difference 2.68×10^-14.

## conventions

Delta = laser frequency minus unshifted atomic resonance, in angular frequency. Effective detuning is Delta - k dot v - (mu_B/hbar)|B|(g_e m_e - g_g m_g). Consequently, for red light and B=b(x,y,-2z), restoring helicities are +,+,- relative to each beam's own propagation direction. Reversing propagation does not reverse the scalar helicity parameter within each pair.

For the F=2 to F'=3 manifold, the absorption CG-squared sums over q at fixed ground m are 7/5, not one. Spontaneous branches sum to one at fixed excited m'. Stimulated pumping is symmetric between ground and excited populations; the resulting population solve provides saturation. Inserting a separate 1+s denominator into the stimulated rates would double count saturation.

## important scope boundaries

- This is a closed-manifold population rate model, not an optical Bloch solver. It ignores optical coherences, laser interference, polarization-gradient/sub-Doppler cooling and coherent dark-state physics.
- The polarization parameter describes an incoherent mixture of circular components, not an arbitrary coherent polarization ellipse or phase-controlled beam geometry.
- The model assumes internal populations instantly reach the local steady state. That assumption weakens where pumping is slow, near dark regions and in rapidly changing fields. B=0 uses an arbitrary quantization-axis fallback; neglecting coherences can make nearby behavior axis dependent.
- Linear Zeeman shifts use approximate g_F values. Hyperfine mixing, off-resonant leakage and the actual repump dynamics are excluded. Real Rb systems need appropriate repumping.
- Mean-force deterministic trajectories omit photon recoil diffusion. They cannot establish equilibrium temperature, true long-time stability or experimental loading rate. Terminal finite-time retention of a declared synthetic initial cloud must not be labeled measured vapor capture probability.
- Magnetic dipole gradient force is omitted; gravity is included separately. Collisions, multiple scattering, reabsorption, light-assisted loss and atom–atom interactions are omitted.
- No company apparatus, beam geometry, original code or measurements are represented. The six-axis and added oblique-pair geometries are generic.

## primary references

- Daniel Steck, Rubidium 87 D Line Data: https://steck.us/alkalidata/rubidium87numbers.pdf (atomic constants and line strengths).
- Eckel et al., PyLCP paper: https://arxiv.org/abs/2011.07979 (rate-equation framework and limitations).
- Authors' executable rate-equation documentation: https://python-laser-cooling-physics.readthedocs.io/en/latest/_modules/pylcp/rateeq.html . The pumping formula and symmetric stimulated terms are independently corroborated at lines 268–291 of the retrieved source; our reference implementation is independently written, not copied.
