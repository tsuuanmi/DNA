# Windowed SNR Method

Part of the canonical [signal-processing method](README.md).

Configuration chooses `window_size_bases` in `5..=10`, a positive finite `minimum_primary_snr`, and `minimum_noisy_windows` of at least `2`. The default window is 10 bases. Windows have that complete width and stride one; short partial windows are never emitted. A read with fewer calls than one window is therefore not processed at all: read processing applies the one length rule (`signal_processing::check_length`) before basecalling, so `analyze` and `basecall` fail typed in the `signal_processing` stage, and `sample` and `consensus` read that failure (`ShortRead::of`) as a rejection with reason `too_short` ([[design/sample-evidence/read-registry|read registry]]). The window of at least 5 calls also covers the two loci that locus-window geometry needs. A noisy interval is emitted only when a consecutive run contains at least the configured number of candidate-noisy windows.

For each channel in a rolling sample span:

```text
baseline = median(samples)
noise_sigma = max(minimum_noise_sigma, MAD(first_difference(samples)) / (0.67448975 × sqrt(2)))
peak_snr = max(0, selected_peak_height - baseline) / noise_sigma
```

The floor, `signal_processing.minimum_noise_sigma` (one unit by default), reflects signed-short quantization and prevents NaN or infinity. Within each call, baseline-corrected selected peaks are ranked deterministically by value and then A/C/G/T order. A window records its minimum primary SNR and maximum secondary SNR internally. Values are rounded to six decimal places before threshold comparison; only each merged region's minimum primary SNR is serialized.

A window is `candidate_noisy` only when its minimum primary SNR is strictly below the configured threshold. Overlapping or adjacent candidate windows are unioned; clean gaps are never filled. Secondary SNR participates only in internal observation and does not make a window noisy because a strong secondary peak may be real mixed dna. Compact analysis v9 and basecalls v3 omit individual windows and secondary-SNR values.
