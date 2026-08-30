# Phase-Ramsey simulation

`phase_ramsey_simulation.py` simulates a phase-scanned Ramsey experiment on
Ba-137 transitions from BVA state `0` to states with different measured magnetic
field sensitivities.

`phase_ramsey_simulation.ipynb` contains the same workflow as a fully executed,
sectioned notebook. Its cells separately cover BVA configuration loading,
calibration lookup, propagation, phase-fringe fitting, Voigt fitting, user
controls, simulation, exports, and plotting.

For each state and wait time, the simulated sequence is:

1. a phase-zero pi/2 pulse;
2. free evolution;
3. a pi/2 pulse scanned over at least 10 phases;
4. measurement of the coupled-state population.

Each fringe is fitted to

```text
P(phi) = offset + amplitude * sin(phi + phase)
```

The output table records the fitted offset, amplitude, phase, their propagated
uncertainties, residual RMS, and contrast. Contrast is defined as
`amplitude / offset`, so an ideal population fringe has contrast one.

For each transition, the contrast is then fitted to the time-domain
characteristic function of a Voigt frequency distribution:

```text
C(t) = C0 * exp[-(t/T2_G)^2 - t/T2_L]
```

The reported effective `T2*` is the positive time at which the fitted envelope
falls to `C0/e`. Its uncertainty includes the fitted covariance between the
Gaussian and Lorentzian time constants, matching the Ramsey-analysis scripts.

## BVA model provenance

The script reads the active scalar noise values directly from
`BVA_analysis_scripts/BVA_qudit_analysis_script.ipynb` every time it runs. It
imports the noise distributions and sequence conventions from
`General_Simulations/sequence_simulation_helpers.py`, which implements the same
model used by the BVA simulation:

- truncated Voigt frequency-calibration noise;
- Gaussian magnetic-field noise multiplied by the transition sensitivity;
- truncated Voigt laser-frequency noise;
- a Gaussian shot-to-shot Rabi-frequency scale shared by both Ramsey pulses;
- the fitted deterministic AC-line signal and optional stochastic field offset;
- quasi-static stochastic detunings shared across the complete Ramsey sequence.

The sensitivity, transition strength, reference pi times, fitted line-signal
parameters, and retained BVA pi-time scale are taken from the same inputs as the
BVA notebook. Repository copies are used only if the notebook's `Z:` paths are
not mounted.

## Run

For interactive work, open and run:

```text
phase_ramsey_simulation.ipynb
```

The standalone command-line version can be run from this folder with:

```powershell
& C:\Users\iamga\anaconda3\python.exe phase_ramsey_simulation.py
```

The default run uses 16 scan phases, 17 wait times from 0 to 2 ms, the BVA shot
count, and five transitions spanning the measured absolute magnetic
sensitivities.

To choose states or sampling explicitly:

```powershell
& C:\Users\iamga\anaconda3\python.exe phase_ramsey_simulation.py --state-indices 8 6 2 7 1 --num-phases 16 --num-wait-times 17 --max-wait-us 2000 --shots 200
```

`--max-wait-us` is intentionally limited to 2000 us. `--num-phases` must be at
least 10.

## Outputs

The `results` folder contains:

- `phase_ramsey_scans.csv`: mean population and noise-shot standard deviation
  for every state, wait time, and scan phase;
- `phase_ramsey_fits.csv`: fitted amplitude, phase, contrast, and residual for
  every state and wait time;
- `phase_ramsey_t2_star.csv`: Voigt Gaussian and Lorentzian time constants,
  effective `T2*`, propagated uncertainties, and fit statistics for each state;
- `phase_ramsey_metadata.json`: exact noise values, input paths, transitions,
  timing, fitted AC-line parameters, and random seed used for the run;
- `phase_ramsey_contrast.png` and `.pdf`: contrast data with the fitted Voigt
  characteristic functions and extracted `T2*` values.
