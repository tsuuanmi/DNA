# Observational DNA-Quality Analysis

## Scope

DNA reads the analyzed ABIF `DATA.9`–`DATA.12` arrays in canonical A/C/G/T order. These are instrument-analyzed fluorescence channels, not raw detector channels. The current ABIF boundary does not retain a spectral matrix, mobility model, or raw-channel baseline metadata.

The signal-processing stage is deliberately observational. It retains `dna.windowed_snr/v1` noisy-window behavior, derives internal basecall-independent `LocusEvidence` / `EvidenceProfile`, and derives concise whole-trace integrity evidence. Public JSON emits trace-integrity observations plus merged candidate-noisy regions; signal processing itself does not smooth channels, re-call bases, trim internal sequence, mutate an alignment, classify dye blobs, or remove a variant. Reference alignment may consume the immutable evidence profile under ADR-0029.

## Coordinate domains

- **Sample indexes** address A/C/G/T channel values and PLOC positions.
- **Call indexes** address base calls, quality records, and variant mappings.

Both are 0-based. Window and region intervals are half-open. Shared PLOC geometry defines one midpoint-derived locus window per vendor locus. Basecalling and signal evidence consume that same geometry without one stage re-deriving the other's classification.

## Method

Configuration chooses `window_size_bases` in `5..=10`, a positive finite `minimum_primary_snr`, and `minimum_noisy_windows` of at least `2`. The default window is 10 bases. Windows have that complete width and stride one; short partial windows are never emitted. A noisy interval is emitted only when a consecutive run contains at least the configured number of candidate-noisy windows.

For each channel in a rolling sample span:

```text
baseline = median(samples)
noise_sigma = max(1, MAD(first_difference(samples)) / (0.67448975 × sqrt(2)))
peak_snr = max(0, selected_peak_height - baseline) / noise_sigma
```

The one-unit floor reflects signed-short quantization and prevents NaN or infinity. Within each call, baseline-corrected selected peaks are ranked deterministically by value and then A/C/G/T order. A window records its minimum primary SNR and maximum secondary SNR internally. Values are rounded to six decimal places before threshold comparison; only each merged region's minimum primary SNR is serialized.

A window is `candidate_noisy` only when its minimum primary SNR is strictly below the configured threshold. Overlapping or adjacent candidate windows are unioned; clean gaps are never filled. Secondary SNR participates only in internal observation and does not make a window noisy because a strong secondary peak may be real mixed dna. Compact analysis v7/basecalls v2 omit individual windows and secondary-SNR values.

## Trace integrity evidence

DNA retains the validated PLOC series as the current method's event-anchor
authority, but optional vendor PBAS/PCON cardinality is treated separately.
A vendor series may be shorter or longer without creating or deleting DNA
calls. Public integrity evidence records the PLOC count and, when present, each
vendor-series count.

For adjacent PLOC values DNA records minimum, median, and maximum spacing in
trace-sample units. One-locus inputs have no spacing summary.

Analyzed DATA values are signed 16-bit samples. Values exactly equal to the
signed-16-bit extrema are counted as `clipped_channel_samples`. This is an exact
representation-boundary observation, not a general artifact classifier.

For each PLOC locus, let:

~~~text
event_signal = sum(A/C/G/T corrected amplitudes at the refined event sample)
~~~

Among loci with positive event signal, DNA records:

~~~text
maximum_to_median_event_signal_ratio =
    max(event_signal) / median(event_signal)
~~~

using the existing six-decimal metric rounding. This ratio is deliberately not
thresholded into a dye-blob/high-amplitude-artifact label. It is evidence for
validation and review until a classifier has its own specification and truth
data.

Trace-integrity observations never change calls, SNR windows, quality, trim,
alignment, or variants.

## Locus evidence and evidence profile

For each PLOC-defined locus, DNA selects a deterministic rolling context and estimates A/C/G/T local baseline and first-difference-MAD noise with the same primitives used by windowed SNR. It then forms the total non-negative baseline-corrected A+C+G+T signal series inside the locus window and identifies positive local maxima. Event refinement selects the local maximum nearest PLOC; equally distant candidates prefer greater total corrected signal and then the lower sample coordinate. If no positive total-signal local maximum exists, the validated PLOC sample is used directly. This keeps event placement nucleotide/basecall-independent while preventing a stronger neighboring base event elsewhere in the midpoint window from taking over the locus profile.

At that one event sample, `LocusEvidence` retains raw A/C/G/T channel values, local baseline/noise, corrected amplitudes, and per-channel SNR. `EvidenceProfile` normalizes only the positive corrected signal mass:

```text
weight[channel] = corrected_amplitude[channel] / sum(corrected_amplitudes)
```

If the total corrected amplitude is zero, the profile is absent. DNA does not inject a uniform profile, reference base, or caller-derived fallback.

This profile is intentionally independent of primary base, ambiguity/IUPAC code, selected basecall peaks, qualifying-channel membership, and `secondary_peak_ratio`. The primary caller remains unchanged. Reference-guided alignment consumes the retained post-trim profile sequence through the fixed-point profile-aware Gotoh method in ADR-0029; persistent mixed-signal interpretation remains a separate downstream method.

## Interpretation limits

The SNR is an uncalibrated local feature, not a Phred score or error probability. Fixed thresholds may not transfer across instruments, chemistries, or runs.
Because first differences are measured across the complete span, real peak edges,
broad peaks, and homopolymers can inflate the estimated noise; this v1 metric is
not a detector-background measurement. A merged region is the union of low-SNR
windows, not a per-call noise classification, so it can include calls that were
not individually weak. Candidate-noisy regions therefore have no
variant-eligibility authority in v1.

Phred demonstrates that trace features require empirical calibration before becoming error probabilities: [Ewing et al. 1998](https://doi.org/10.1101/gr.8.3.175) and [Ewing & Green 1998](https://doi.org/10.1101/gr.8.3.186).

## Deferred cleaning

No disabled transform or compatibility branch is included. A later behavior-changing method must preserve the decoded trace and produce a separate processed projection. Candidate methods include peak-preserving Savitzky–Golay smoothing ([Savitzky and Golay 1964](https://doi.org/10.1021/ac60214a047)), asymmetric baseline correction ([Eilers 2003](https://doi.org/10.1021/ac034173t); [airPLS](https://doi.org/10.1039/b922045c)), and wavelet soft-thresholding ([Donoho 1995](https://doi.org/10.1109/18.382009)).

Before any transform or noisy-region filter affects calls, validation must use approved truth-labeled traces and synthetic major/secondary peaks, baseline drift, impulses, compressed peaks, homopolymers, and read ends. It must measure secondary-peak retention and both false-positive and false-negative variants.
