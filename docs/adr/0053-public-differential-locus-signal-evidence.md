# ADR-0053: Expose concise differential-locus signal evidence in sample output

- **Status:** Accepted
- **Date:** 2026-09-20

## Context

DNA already preserves basecall-independent reference-oriented A/C/G/T
`EvidenceProfile` values, candidate-noisy-region membership, and factorized
support topology for every retained differential sample locus.

Those values currently remain internal while `dna.sample_evidence/v7`
serializes only each covering read's state, observed base, and relative quality.
That loses the most useful Tracy-derived lesson at the public review boundary:
primary sequence should remain one interpretation of the chromatogram rather
than the only representation available to downstream sample review.

The validation tooling can recover richer evidence because it uses internal
types directly, but normal `signal sample` output should not require a separate
validation path merely to distinguish reproducible mixed signal from clean
cross-read disagreement.

## Decision

Replace `dna.sample_evidence/v7` with
`dna.sample_evidence/v8`.

Every serialized `locus_differences[]` record MUST additionally expose the
existing factorized locus support topology:

- total reads;
- forward/reverse reads;
- reference/alternate/unresolved/deletion reads;
- profile-bearing reads;
- profile-bearing forward/reverse reads.

Every call-backed locus observation MUST additionally expose:

- `in_noisy_region`, copied from the already computed sample call-signal
  evidence;
- optional normalized reference-oriented A/C/G/T `profile` when a real
  basecall-independent `EvidenceProfile` exists.

A zero-signal call remains profile-less. A deletion has no source call signal
and therefore serializes neither profile nor noisy-region membership.

Reverse-read profiles MUST already be complemented into reference A/C/G/T
orientation before publication. The report layer MUST project existing evidence
only; it MUST NOT recompute signal features or infer missing profiles.

## Compactness boundary

v8 intentionally does **not** serialize:

- corrected A/C/G/T amplitudes;
- per-channel SNR;
- mean aggregate profiles;
- within/between heterogeneity geometry;
- directional profile distance;
- call index, PLOC, or event sample coordinates.

Those remain internal/validation evidence until a concrete reviewer or
production need justifies their public cost.

The public profile is the smallest useful signal-preserving representation for
one differential-locus observation.

## Consequences

- Reviewers can inspect the A/C/G/T evidence shape behind a differential call
  without reconstructing it from primary bases or raw peak thresholds.
- Reference support, alternate support, unresolved calls, deletions, selected
  orientation, and profile availability remain explicit at the same locus.
- Candidate-noisy context is visible without becoming an exclusion or weight.
- The sample contract remains sparse: routine all-reference loci are still
  omitted.
- No consensus state or winning nucleotide is introduced.
- v7 is not emitted as a compatibility alias.

## Tracy-phase implication

This promotion closes the remaining high-value Tracy lesson around preserving
chromatogram evidence through sample reconciliation and reviewer-facing output.

Further Tracy-derived work is intentionally separate research:

- persistent mixed-signal / phase-shift detection;
- poly-C phase-recovery modeling;
- calibrated contribution weighting;
- final sample-consensus policy;
- large-reference candidate search;
- VCF/BCF projection.

Those items require separate validation or product scope and are not required to
keep the Tracy research phase open.

## Non-goals

This ADR does not define a consensus base, sample-level variant verdict,
heteroplasmy interpretation, quality weighting, artifact classifier, phase-shift
detector, assay metadata contract, or additional read-admission rule.
