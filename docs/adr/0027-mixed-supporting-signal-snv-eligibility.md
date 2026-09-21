# ADR-0027: Mixed supporting signal is not a clean SNV

- **Status:** Accepted
- **Date:** 2026-09-19

## Context

Tracy research highlighted a high-ROI failure mode in primary-sequence variant
reporting: the strongest channel can differ from the reference while meaningful
co-localized secondary signal is present at the same nucleotide event.

DNA already solves the harder low-level part in `dna.peak_recall/v3`.
A secondary channel is counted as qualifying only when it satisfies the configured
ratio both at its selected peak and at the uniquely strongest primary-event sample.
Therefore `BaseCall.qualifying_channels` already distinguishes a clean canonical
call from a call with meaningful co-localized competing dna.

Before this decision, a canonical strongest primary could still produce a reported
SNV even when two or three channels qualified. That output looked like an ordinary
clean substitution even though the chromatogram evidence was explicitly mixed.

## Decision

DNA promotes the variant eligibility method to
`dna.primary_difference/v4`.

For normalized SNVs:

- the supporting call must continue to pass the configured peak and relative-quality
  gates;
- the supporting call must also have exactly one qualifying channel under the
  authoritative basecalling method;
- if more than one channel qualifies, the normalized SNV remains in
  `ObservedVariant` / sample evidence but is not clean-report eligible;
- the stable exclusion reason is `mixed_supporting_dna`.

The mixed-signal rule is intentionally limited to SNVs. It does not apply to
insertions or deletions because mixed-length products require a separate
persistent phase-shift / length-mixture method rather than a point-SNV heuristic.

The rule reuses existing co-localization and `secondary_peak_ratio` behavior.
It introduces no new threshold and does not reinterpret the evidence as genotype,
heteroplasmy, contamination, or allele fraction.

Because sample evidence exposes exclusion reasons through a closed enum,
the public sample contract advances from `dna.sample_evidence/v2` to
`dna.sample_evidence/v3`. The v2 schema/example are removed rather than
retained as compatibility output. `dna.analysis/v6` keeps the same serialized
shape; only scientific eligibility behavior changes.

## Consequences

- A mixed canonical call can still be aligned and normalized, preserving where the
  strongest primary differs from the reference.
- Single-read `analysis/v6` no longer presents that observation as an ordinary
  clean SNV.
- Sample evidence can retain and reconcile the observation with explicit
  `eligible=false` and `mixed_supporting_dna`.
- Existing reviewer-facing A/C/G/T peak evidence remains available to inspect the
  competing signal directly.
- Future richer mixed-base or length-mixture analysis can build on the retained
  observation without undoing this conservative gate.

## Non-goals

This decision does not estimate mixture fraction, call heteroplasmy, add genotype
semantics, change basecalling thresholds, alter alignment, or detect persistent
post-indel phase shifts.

## Follow-up

ADR-0030 advances the current sample contract from `dna.sample_evidence/v3`
to `dna.sample_evidence/v4` to expose Tracy-derived pairwise overlap/admission
evidence. The `mixed_supporting_dna` eligibility semantics defined here remain
unchanged.
