# Locus Evidence and Evidence Profile

Part of the canonical [signal-processing method](README.md).

For each PLOC-defined locus, DNA selects a deterministic rolling context and estimates A/C/G/T local baseline and first-difference-MAD noise with the same primitives used by windowed SNR. It then forms the total non-negative baseline-corrected A+C+G+T signal series inside the locus window and identifies positive local maxima. Event refinement selects the local maximum nearest PLOC; equally distant candidates prefer greater total corrected signal and then the lower sample coordinate. If no positive total-signal local maximum exists, the validated PLOC sample is used directly. This keeps event placement nucleotide/basecall-independent while preventing a stronger neighboring base event elsewhere in the midpoint window from taking over the locus profile.

At that one event sample, `LocusEvidence` retains raw A/C/G/T channel values, local baseline/noise, corrected amplitudes, and per-channel SNR. `EvidenceProfile` normalizes only the positive corrected signal mass:

```text
weight[channel] = corrected_amplitude[channel] / sum(corrected_amplitudes)
```

If the total corrected amplitude is zero, the profile is absent. DNA does not inject a uniform profile, reference base, or caller-derived fallback.

This profile is intentionally independent of primary base, ambiguity/IUPAC code, selected basecall peaks, qualifying-channel membership, and `secondary_peak_ratio`. The primary caller remains unchanged. Reference-guided alignment consumes the retained post-trim profile sequence through the fixed-point profile-aware Gotoh method in ADR-0029; persistent mixed-signal interpretation remains a separate downstream method.
