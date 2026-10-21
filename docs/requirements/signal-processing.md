---
implements: [SRS-CALL, SRS-SAMPLE]
decided-by: [ADR-0013, ADR-0028, ADR-0031, ADR-0046, ADR-0088]
verified-by:
  - crates/dna-sanger/src/locus.rs
  - crates/dna-sanger/src/model/locus_evidence.rs
  - crates/dna-sanger/src/signal_processing
---

# Observational Signal Processing Requirements

**Requirement namespace:** `SRS-SIG-*`

These requirements are part of the canonical [[requirements/SRS|DNA SRS]].

- **SRS-SIG-001:** DNA processing MUST read the immutable analyzed A/C/G/T channels and basecalling evidence; it MUST NOT replace decoded channels or selected peaks.
- **SRS-SIG-002:** Configuration MUST require a rolling window in `5..=10` bases, a finite positive minimum primary SNR, and a minimum noisy-window run length of at least 2.
- **SRS-SIG-003:** Every internally calculated window MUST have the complete configured width and stride one with explicit 0-based half-open call and sample intervals. Compact JSON MUST emit merged noisy regions only, not individual windows.
- **SRS-SIG-004:** Baseline, first-difference MAD noise, primary SNR, secondary SNR, noise flooring at `signal_processing.minimum_noise_sigma`, rounding, and threshold comparison MUST follow the documented `dna.windowed_snr/v1` formula.
- **SRS-SIG-005:** Overlapping or adjacent candidate-noisy windows MUST be unioned without bridging clean gaps, and a noisy interval MUST require at least the configured minimum run length.
- **SRS-SIG-006:** Candidate-noisy window/region annotations (`dna.windowed_snr/v1`) MUST NOT alter quality scores, trim bounds, alignment, warning totals, or variant eligibility. `LocusEvidence` / `EvidenceProfile` MAY be consumed by reference alignment only under the explicit [[requirements/alignment|alignment requirements]] and by read callability only under [[requirements/callability|`SRS-CALL-*`]], and the Sanger evidence MAY withhold the calls inside an electrophoretic spike (SRS-SIG-012) only under SRS-CALL-027; no other consumer may turn a signal feature into a trim, placement, or eligibility decision. No signal feature may be described as Phred, error probability, genotype, or heteroplasmy evidence.
- **SRS-SIG-007:** A read with fewer calls than the configured signal window MUST fail `analyze` and `basecall` with a typed signal-processing error before basecalling, rather than emit partial windows; in `sample` and `consensus` it MUST be recorded as rejected instead (SRS-SAMPLE-030).
- **SRS-SIG-008:** DNA processing MUST derive one internal `LocusEvidence` record per validated canonical Sanger locus directly from analyzed A/C/G/T channels and shared locus-window geometry. Event refinement MUST consider positive local maxima of total non-negative baseline-corrected A/C/G/T amplitude inside the locus window and select the candidate nearest the canonical locus position; equally distant candidates MUST prefer greater total corrected signal and then the lower sample coordinate. If no positive total-signal local maximum exists, the validated canonical locus sample MUST be used. Event selection MUST remain independent of primary/ambiguity calls, selected channel identity, qualifying-channel membership, alignment, and reference sequence.
- **SRS-SIG-009:** `EvidenceProfile` MUST normalize positive corrected A/C/G/T amplitudes at the refined event sample and MUST NOT depend on primary/ambiguity calls, selected basecall peaks, qualifying-channel membership, or `secondary_peak_ratio`. A locus with zero total corrected signal MUST have no profile; no uniform or reference-guided fallback is permitted.
- **SRS-SIG-010:** DNA processing MUST derive one immutable Sanger-integrity record containing canonical locus count; optional PBAS/PCON counts; minimum/median/maximum adjacent locus spacing when at least two loci exist; the count of analyzed channel samples exactly equal to signed-16-bit extrema; and, when positive locus event signal exists, the six-decimal ratio of maximum to median total corrected event dna.
- **SRS-SIG-011:** Sanger-integrity observations MUST NOT rewrite decoded channels, add/remove canonical loci, change base calls, alter quality/trim, steer reference placement, or change variant eligibility. Vendor-length mismatch and exact clipping MAY contribute warning counts; the event-signal ratio MUST remain unthresholded observation evidence until a separately specified artifact classifier exists (the spike detector of SRS-SIG-012 does not read it).
- **SRS-SIG-012:** DNA processing MUST find a read's electrophoretic spikes ([[ADR-0088]]) from its analysed channels and its own calls only, never a reference position, alignment, or truth. A spike MUST be a set of four local channel maxima, one per channel, all within `signal_processing.spike_span` median call spacings of the earliest, each at least `signal_processing.spike_height` times its local level: the median primary peak (the highest selected peak of a call) of the `signal_processing.spike_level_calls` calls on each side of the call whose locus window holds the maximum, that call left out, and at least one. Of the qualifying maxima within the span, each channel MUST contribute its highest, and the next spike MUST start after the latest of them. A spike's samples MUST run from the earliest of its maxima's half-height starts to the latest of their half-height ends, and its calls MUST be those whose primary peak lies among them. A read without a median spacing has no spike. Spikes MUST NOT rewrite channels, loci, calls, or selected peaks.


## Implemented in

<!-- docs:begin backlinks -->
- **Implemented in:** [design/callability](../design/callability.md), [design/signal-processing](../design/signal-processing/README.md), [design/signal-processing/spikes](../design/signal-processing/spikes.md), [reference/configuration](../reference/configuration.md), [requirements/callability](callability.md), [requirements/sample-evidence/callability](sample-evidence/callability.md)
<!-- docs:end backlinks -->
