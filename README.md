# Synthetic coherency demonstration for GNSS scintillation

This repository contains a MATLAB script that generates the synthetic dataset and figure used to illustrate the coherency-based measurement framework developed in:

**A Coherency-Based Measurement Framework for GNSS Scintillation in a Magnetized Ionosphere**

The goal of the experiment is not to simulate a full GNSS receiver chain at the raw-sample or tracking-loop level, but to provide a controlled **coherency-level demonstration** of the paper’s central measurement-theoretic claims.

## Purpose of the synthetic experiment

The script demonstrates three key ideas:

1. **Receiver dependence of scalar observables**  
   The same underlying propagated coherency state can produce different scalar receiver outputs when observed through different receiver polarization response matrices.

2. **Receiver de-embedding**  
   If the receiver response matrix and noise covariance are known, the receiver-dependent distortions can be removed to recover a receiver-normalized estimate of the propagated second-order field statistics.

3. **Invariant behavior under unitary evolution versus mixing**  
   Basis-invariant quantities such as degree of polarization (DoP) and polarization entropy remain constant under unitary polarization evolution, but change under effective depolarization or mixing.

## Files

- `synthetic_coherency_demo.m`  
  Main MATLAB script that generates the synthetic dataset and the multi-panel figure.

- `output_figure4/`  
  Output folder created automatically by the script. It contains:
  - editable MATLAB figure (`.fig`)
  - vector PDF
  - high-resolution PNG
  - EPS
  - SVG (if supported by MATLAB version)
  - TIFF (if supported)
  - `.mat` file containing the synthetic data used in the figure

## Conceptual structure of the experiment

The synthetic propagated coherency state is divided into three regimes:

### Regime I: Unitary rotation
The normalized Stokes vector evolves smoothly while its magnitude remains constant. This corresponds to unitary polarization evolution in which polarization state changes but spectral invariants such as DoP and polarization entropy remain fixed.

### Regime II: Same field, different receivers
The underlying propagated coherency state remains governed by the same unitary evolution, but two distinct coherent receiver response matrices are used to generate measured covariance data. This shows that scalar receiver products can differ even when the underlying field is the same.

### Regime III: Mixing / depolarization
The propagated coherency state is constructed as a convex mixture of two locally unitary states with different polarization states. This produces effective depolarization: DoP decreases and polarization entropy changes. This regime is intended to represent the effect of unresolved path superposition, multipath, scattering, or other mixing processes at the level of second-order statistics.

## Description of the figure panels

### Panel (a): Synthetic propagated coherency state
This panel shows the **normalized Stokes components**
\[
p_i = S_i/S_0, \quad i=1,2,3,
\]
derived from the synthetic propagated coherency matrix. Normalized Stokes components are used rather than raw Stokes parameters so that polarization-state evolution is shown independently of total intensity fluctuations.

- In Regimes I and II, the state evolves under unitary rotation with constant DoP.
- In Regime III, the state enters a mixing/depolarization regime.

### Panel (b): Receiver-dependent scalar observables
This panel shows a **windowed scalar \(S_4\)-proxy** computed from the first receiver-channel power for two distinct receiver response matrices.

The purpose of this panel is to show that scalar observables remain receiver-dependent projections of the same underlying vector field. The quantity `mean |Δ|` is the mean absolute difference between the two scalar proxy curves over the plotted interval.

### Panel (c): De-embedded receiver-normalized recovery
This panel compares the true propagated intensity
\[
S_0 = \mathrm{tr}(\Jmat_{\mathrm{out}})
\]
with the receiver-normalized recovered intensities
\[
\hat S_0^{(A)} = \mathrm{tr}(\hat{\Jmat}^{(A)}), \qquad
\hat S_0^{(B)} = \mathrm{tr}(\hat{\Jmat}^{(B)}).
\]

The close agreement of the recovered A and B curves shows that de-embedding removes the receiver-specific distortion. The remaining mismatch to the synthetic truth is intentional and arises from the inclusion of common scalar weighting and additive noise. The RMSE values quantify deviation from the synthetic truth.

### Panel (d): Invariant diagnostics: true and recovered
This panel shows the true and recovered **degree of polarization** and **polarization entropy**.

- In Regimes I and II, these invariants remain constant, as expected for unitary polarization evolution.
- In Regime III, they change as a result of mixing / depolarization.

The close overlap between true and recovered invariant curves demonstrates that receiver-normalized invariant quantities are recovered consistently across the two receiver models.

## Important modeling note

This experiment is carried out at the **coherency level**. It is not intended as a full end-to-end simulation of GNSS modulation, tracking loops, or receiver firmware. Its purpose is instead to isolate and demonstrate the measurement-theoretic structure emphasized in the paper:

- propagation acts on the coherency matrix,
- receivers observe that state through instrument-dependent mappings,
- and de-embedding can recover receiver-normalized second-order quantities.

## How to run

Open MATLAB, navigate to the folder containing `synthetic_coherency_demo.m`, and run:

```matlab
synthetic_coherency_demo
