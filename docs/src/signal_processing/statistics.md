# `src/signal_processing/statistics.rs`

## Purpose

Provides shared local baseline, noise, corrected-amplitude, SNR, and metric-rounding primitives for observation-only signal processing.

## Responsibilities

- Estimate one channel baseline as the median of a validated sample span.
- Estimate one channel noise sigma from first-difference MAD with the existing one-unit floor.
- Convert one raw analyzed-channel value into a non-negative baseline-corrected amplitude.
- Convert a corrected amplitude and local noise sigma into SNR.
- Provide deterministic median helpers for integer/floating evidence summaries.
- Round emitted or thresholded signal metrics to six decimal places.

## Non-responsibilities

No rolling-window geometry, per-call context selection, peak selection, channel mutation, smoothing, calling, quality scoring, or serialization.

## Formula

```text
baseline = median(samples)
noise_sigma = max(1, MAD(first_difference(samples)) / (0.67448975 × sqrt(2)))
corrected_amplitude = max(0, raw_height - baseline)
snr = corrected_amplitude / noise_sigma
```

## Invariants and errors

- Baseline estimation rejects an empty sample span.
- Noise estimation requires at least two samples.
- The noise floor guarantees finite non-negative SNR for finite corrected amplitudes.
- Rolling-window and per-call calculations share these exact primitives.

## Tests

Tests cover odd/even integer and floating medians, invalid non-finite float input, the quantization floor, first-difference MAD scaling, negative-value correction, finite SNR, and six-decimal rounding.

## Status

Implemented.
