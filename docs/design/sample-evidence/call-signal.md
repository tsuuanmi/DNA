# Reference-oriented Call Signal Projection

Part of the canonical [sample evidence aggregation method](README.md).

Sample aggregation works on modality-neutral `CalledRead` records
([ADR-0069](../../decisions/adr/0069-plugin-first-modality-core-post-calling.md)).
Every call-backed observation keeps its source call index and the call's
optional basecall-independent `EvidenceProfile` from the read's `ReadEvidence`.
The profile is projected into reference orientation:
- forward reads keep A/C/G/T order;
- reverse reads complement it into reference-oriented A/C/G/T.

A zero-signal call remains profile-less. A deletion observation has no call
index and no profile.

Sanger signal evidence stays in the read's Sanger attachment:
- baseline-corrected amplitudes and per-channel SNR from `LocusEvidence`;
- merged candidate-noisy-region membership;
- relative quality and primary-event peaks.

The sample report joins it by read identity and call index. It publishes the
reference-oriented profile, noisy-region membership, and quality at retained
differential loci, and peaks and quality for variant calls. Corrected
amplitudes and per-channel SNR feed only operational metrics.

None of this evidence alters placement, overlap admission, callability,
variant eligibility, or consensus weighting.
