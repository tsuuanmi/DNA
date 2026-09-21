# `src/signal_processing/features.rs`

## Purpose

Calculates bounded rolling baseline, noise, and peak-SNR observations from analyzed A/C/G/T samples.

## Responsibilities

- Build full-width, stride-one base-call windows from retained basecalling sample intervals.
- Reuse shared median-baseline and first-difference-MAD noise primitives for each
  channel.
- Apply a one-channel-unit noise floor so every result remains finite.
- Baseline-correct and rank the four selected channel peaks per call.
- Retain minimum primary SNR, maximum secondary SNR, and the configured
  candidate-noisy decision in each internal window. Pipeline logging computes its
  aggregate maximum secondary SNR from those windows, while compact v5 serializes
  neither rolling windows nor maximum secondary SNR.
- Round output and comparison metrics to six decimal places.

## Non-responsibilities

No channel mutation, smoothing, peak selection, per-call primary-event metrics,
calibrated quality, region merging, trimming, or variant filtering.

## Formula

For one channel in one rolling sample span:

```text
baseline = median(samples)
noise_sigma = max(1, MAD(first_difference(samples)) / (0.67448975 × sqrt(2)))
peak_snr = max(0, peak_height - baseline) / noise_sigma
```

A window is candidate-noisy only when its rounded minimum primary SNR is strictly below `minimum_primary_snr`. The run-length setting is applied by `regions`; it does not alter individual window flags.

## Errors and tests

A read shorter than `window_size_bases`, an invalid sample interval, or a shared
statistics failure returns `Error::DNAProcessing`. Unit tests cover full
windows, rolling primary/secondary SNR regression values, finite metrics, exact
threshold equality, and short reads. Statistical primitive tests live in
`statistics`.

## Status

Implemented.
