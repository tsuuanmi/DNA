# ADR-0026: Prefer reviewer-facing signal evidence over implementation call coordinates

- **Status:** Accepted
- **Date:** 2026-09-19

## Context

Earlier compact analysis contracts retained original call index, mapped call position,
and ABIF PLOC for every variant-associated call. Those fields were useful while
validating coordinate plumbing, but they duplicate implementation coordinates that
a routine reviewer rarely needs once the normalized variant already provides its
biological position.

Review of real multi-read output showed a more important question for manual
assessment: whether the called allele is supported by the chromatogram at that
call, and how strong the call quality is.

The existing basecalling model already retains two different peak concepts:

1. one strongest selected peak per A/C/G/T channel inside the call window; these
   peaks can occur at different sample coordinates;
2. `PrimaryPeakEvidence.channel_heights`, which samples all four analyzed
   channels together at the uniquely strongest primary-event coordinate.

The second representation is the appropriate compact reviewer view because the
four values describe one shared signal event.

Reverse reads add another usability problem: raw trace channel labels are on the
trace strand, while variants are reported on the reference strand.

## Decision

DNA replaces `dna.analysis/v5` with `dna.analysis/v6`.

Each public variant-associated call contains only:

- `role`: `supporting` or `flanking`;
- `base`: called base projected onto the reference strand;
- `peaks`: co-located A/C/G/T raw analyzed channel heights at the primary-event
  coordinate, with channel labels projected onto the reference strand;
- `quality`: the existing uncalibrated relative quality score under a concise
  reviewer-facing field name.

Original call index, mapped call position, ABIF PLOC, trace-strand
primary/ambiguity, selected-peak positions/sources, and maximum-peak-only summaries
remain internal.

The public field name `quality` does not change the quality algorithm or claim
Phred calibration. Configuration and exclusion reasons may continue to use
explicit `relative_quality` terminology because they describe the underlying
method.

For deletions, DNA emits real flanking call signal/quality and never fabricates
evidence at a deleted reference base. Insertions retain supporting inserted calls
plus available flanks.

At adoption, sample evidence v2 used the same call-evidence shape. Its public `read` references
use unique filename stems rather than numeric registry indexes; SHA-256 remains the
scientific content identity and filenames remain non-authoritative for placement or
reconciliation.

## Consequences

- Reviewers can compare the normalized alternate allele directly with four-channel
  evidence without resolving implementation coordinates.
- Reverse-read evidence is immediately comparable with reference-oriented alleles.
- The contract is smaller and more semantically focused.
- Exact call/PLOC mappings remain available internally for tests and scientific
  implementation invariants but are no longer public API.
- Analysis v5 schema/example are removed instead of retained as compatibility
  output.

## Non-goals

This decision does not calibrate quality to Phred, infer allele fraction or
heteroplasmy, change basecalling, alter trimming/alignment, or use filenames as
scientific keys.

## Follow-up

ADR-0027 advanced sample evidence to v3 for mixed-SNV eligibility while retaining the same reviewer-facing `role/base/peaks/quality` call-evidence shape. ADR-0030 advances the current sample contract to v4 with pairwise overlap/admission evidence and leaves this call-evidence shape unchanged.
