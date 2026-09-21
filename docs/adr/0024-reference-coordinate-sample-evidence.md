# ADR-0024: Reference-coordinate sample evidence contract

- **Status:** Superseded in part by ADR-0025
- **Date:** 2026-09-18

## Context

ADR-0023 established `ReadObservation` as the one-read boundary and required future
sample reconciliation to operate in mapped reference-coordinate/variant space rather
than canonical F/R pairs. Real mtDNA sample output confirms the need for that
boundary: reverse reads retain raw trace-strand bases that differ from their
reference-strand normalized variants, non-canonical assay pairs can overlap, and
absence from a read's reported variant list does not prove reference support.

## Decision

DNA adds a separate sample operation:

`signal sample <sample-id> <trace.ab1>... --reference <reference.fasta>`

Every trace is independently processed through the same `pipeline::observation`
path used by single-read `analyze`. Only completed `ReadObservation` values enter
sample aggregation.

Sample evidence has two complementary layers:

1. **Locus observations** come directly from selected alignment columns on the
   supplied reference strand. Covered loci are classified as `reference`,
   `alternate`, `unresolved`, or `deletion`, with original call index and
   relative quality retained when a query call exists.
2. **Normalized variants** aggregate every canonical normalized read-level observation by biological `(position, reference, alternate, kind)` identity and retain every contributing read basename, SHA-256, derived orientation, configured eligibility, exclusion reasons, and concise original-call mappings. Read-level filtering changes eligibility, not whether the normalized observation existed.

A read absent from a locus contributes no observation and never counts as
reference support. Inserted query columns have no reference locus and are handled
through normalized variant evidence. Duplicate trace content is rejected so copied
files cannot be counted twice.

Read basenames are retained strictly as reviewer-facing provenance. Stable identity, duplicate rejection, placement, overlap, and variant reconciliation remain content/alignment driven. Each read also carries a concise selected-alignment summary so placement quality can be reviewed without serializing gapped rows or traceback columns.

The public contract is `dna.sample_evidence/v1`, published atomically as
`results/<sample-id>.sample.json`. The sample identifier affects naming and
provenance only; it never affects placement or scientific reconciliation.

## Consequences

- F/R, HV, primer, replicate, and filename labels are unnecessary for placement
  and variant merging.
- Forward and reverse observations are directly comparable because
  `AlignmentColumn.query_base` is already oriented to the reference strand.
- A normalized variant observed by a read remains sample evidence even when read-level region/signal thresholds make that read ineligible for single-read reporting. Eligibility and exact exclusion reasons stay attached to that read's variant support.
- The evidence contract is intentionally not a consensus or sample-variant
  interpretation contract.

## Non-goals

No consensus sequence, majority voting, conflict adjudication, haplogroup
interpretation, heteroplasmy/genotype inference, metadata-driven placement,
canonical F/R pair object, or compatibility alias is introduced.

## Follow-up

ADR-0025 replaces the dense public loci table with a SHA-sorted read registry and sparse differential loci in dna.sample_evidence/v2. ADR-0024 remains the historical basis for independent read placement, normalized variant aggregation, filtered-observation retention, and the no-consensus boundary.
