# `src/sample/overlap.rs`

## Purpose

Builds deterministic pairwise overlap/admission evidence between independently
placed reads before any future sample consensus.

## Responsibilities

- Consume the SHA-sorted `ReadObservation[]` produced by sample aggregation.
- Reconstruct each read's covered reference coordinates from the selected
  alignment columns without using filename, HV/amplicon, primer, or declared
  direction metadata.
- Materialize an unordered pair only when the two reads share at least one
  reference coordinate.
- Count `shared_positions` for every shared reference coordinate.
- Count `comparable_bases` only when both aligned query observations are
  canonical A/C/G/T.
- Count canonical base/base agreements and conflicts and derive the optional
  agreement fraction.
- Apply `SampleReconciliationConfig.minimum_comparable_bases` and
  `minimum_overlap_agreement` in deterministic rule order.
- Retain exact overlap exclusion reasons without mutating either read.

## Non-responsibilities

No read placement, F/R pairing, raw-trace processing, profile-to-profile
alignment, consensus voting, variant eligibility changes, genotype/heteroplasmy
inference, or gap-quality synthesis.

## Invariants

Non-overlapping reads have no pair edge and remain valid sample evidence.
Unresolved query symbols and deletions do not enter the nucleotide agreement
denominator. Insertions have no reference-coordinate column and therefore do not
create synthetic overlap positions. Gap/indel evidence remains owned by the
existing locus/variant layers.

## Tracy mapping

The method adopts Tracy's explicit pre-consensus minimum-overlap and
minimum-match gates (25 bases and 0.5 by default) while adapting them to DNA's
N-read reference-coordinate architecture. It deliberately does not copy Tracy's
pair-first domain boundary or quality-blind multi-trace majority consensus.

## Tests

Unit tests cover non-overlapping pairs, admitted canonical overlap, exclusion of
unresolved/gap observations from the nucleotide denominator, and deterministic
multi-rule exclusion ordering.

## Status

Implemented.
