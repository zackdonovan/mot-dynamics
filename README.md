# laser cooling and trap dynamics

a numerical study of laser cooling and trapping in a simplified rubidium-87 magneto-optical trap. it compares six- and eight-beam geometries and asks whether a Gaussian-process optimizer finds better laser settings than a Sobol search at the same simulation budget.

the geometry, fields, and incoming atom distributions are synthetic. this is an independent simulation, not a model of a particular lab or company apparatus.

## run

```sh
python -m pip install -r requirements.txt
python -m unittest discover -s tests -p 'test_physics.py' -v
```

see [methods](METHODS.md) for full run commands, assumptions, and cited references. comparisons against PyLCP use the separate benchmark dependencies.

## what the study found

the rate-equation model is useful for fast search, but its behavior near zero magnetic field depends on the choice of quantization axis. coherent dynamics and photon recoil matter enough that the simpler model should not be treated as a capture predictor. the [validation results](VALIDATION_RESULTS.md) report how proposals fare on new simulated atoms and changed conditions. the [photon-resolved check](QUANTUM_JUMP_VALIDATION.md) validates a fixed-atom stochastic model, not moving-atom capture.

scripts and compact result summaries are included. large generated trajectories and figures can be regenerated using the commands in the study notes. no experimental capture efficiency or hardware performance is claimed.
