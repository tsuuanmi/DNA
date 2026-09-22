# Observational Signal Processing

This is the canonical entry point for production signal-processing methods.

## Scope

DNA reads the analyzed ABIF `DATA.9`–`DATA.12` arrays in canonical A/C/G/T order. These are instrument-analyzed fluorescence channels, not raw detector channels. The current ABIF boundary does not retain a spectral matrix, mobility model, or raw-channel baseline metadata.

The signal-processing stage is deliberately observational. It retains `dna.windowed_snr/v1` noisy-window behavior, derives internal basecall-independent `LocusEvidence` / `EvidenceProfile`, and derives concise whole-trace integrity evidence. Public JSON emits trace-integrity observations plus merged candidate-noisy regions; signal processing itself does not smooth channels, re-call bases, trim internal sequence, mutate an alignment, classify dye blobs, or remove a variant. Reference alignment may consume the immutable evidence profile under ADR-0029.

## Coordinate domains

- **Sample indexes** address A/C/G/T channel values and PLOC positions.
- **Call indexes** address base calls, quality records, and variant mappings.

Both are 0-based. Window and region intervals are half-open. Shared PLOC geometry defines one midpoint-derived locus window per vendor locus. Basecalling and signal evidence consume that same geometry without one stage re-deriving the other's classification.

## Methods

- [Windowed SNR](windowed-snr.md): rolling local signal-quality observations and candidate-noisy regions.
- [Trace integrity](trace-integrity.md): PLOC/vendor cardinality, spacing, clipping, and whole-trace event-signal observations.
- [Locus evidence](locus-evidence.md): basecall-independent per-locus A/C/G/T evidence and normalized evidence profiles.

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

No disabled transform or alternate legacy branch is included. A later behavior-changing method must preserve the decoded trace and produce a separate processed projection. Candidate methods include peak-preserving Savitzky–Golay smoothing ([Savitzky and Golay 1964](https://doi.org/10.1021/ac60214a047)), asymmetric baseline correction ([Eilers 2003](https://doi.org/10.1021/ac034173t); [airPLS](https://doi.org/10.1039/b922045c)), and wavelet soft-thresholding ([Donoho 1995](https://doi.org/10.1109/18.382009)).

Before any transform or noisy-region filter affects calls, validation must use approved truth-labeled traces and synthetic major/secondary peaks, baseline drift, impulses, compressed peaks, homopolymers, and read ends. It must measure secondary-peak retention and both false-positive and false-negative variants.
