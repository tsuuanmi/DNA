# `src/variant_calling/extract.rs`

## Purpose

Walks a completed alignment and extracts normalized primary-sequence SNVs and
small indels, mapping each to its original trace calls.

## Responsibilities

- Walk the alignment columns, identifying substitutions, insertions, and
  deletions.
- Exclude unresolved alleles and indels longer than the configured maximum while
  retaining kind, position when available, and every structural rejection reason.
- Build supporting and flanking call mappings for each reported difference.
- Sort, deduplicate, and merge reported variants with their call mappings.

## Non-responsibilities

No heteroplasmy inference, breakpoint detection, two-allele decomposition,
allelic fraction fitting, serialization, or addition of fixed report-only
classification/normalization labels.

## Key types and functions

- `call(alignment, reference, config) -> Result<VariantCallingResult>`: the
  extraction entry point consumed by the parent module before filtering.
- `previous_reference`/`next_reference`, `previous_flank`/`next_flank`,
  `optional_pair`, `is_canonical`: helpers.

## Invariants and errors

- Variants use validated alleles and explicit coordinates.
- Indels longer than `max_indel_length` or with unresolved alleles are excluded
  with one allele-free diagnostic; candidates that fail both rules retain both
  reasons.
- Reported variants are sorted by `(contig, position, ref, alt)` and deduplicated;
  same-locus variants merge their call mappings.
- A missing reference coordinate on an SNV column returns `Error::Variant`.
- Insertions gather only supporting inserted calls plus optional flanking calls;
  deletions use flanking calls only.

## Dependencies

- `config` for `VariantCallingConfig`.
- `normalize` and `mapping`.
- `model::alignment`, `model::reference`, `model::variant`.
- `error` for `Error`/`Result`.

## Biological semantics

The alignment differences are the primary-sequence variants: single-nucleotide
substitutions and small insertions/deletions relative to the reference. Unresolved
or oversized differences are excluded rather than reported. Each reported variant
keeps a direct mapping to its original trace calls so the observed evidence is not
lost. A missing aligned left flank is passed through so normalization can derive
the actual reference predecessor.

## Tests

- `leading_alignment_deletion_uses_reference_predecessor`: verifies a leading
  deletion without an aligned left flank still derives the reference predecessor.
- `unresolved_primary_difference_is_excluded`: verifies an `N` difference keeps
  its SNV kind, position, and noncanonical-allele reason.
- `leading_alignment_insertion_keeps_inserted_and_flanking_mappings`: verifies
  inserted and flanking call mappings are retained.

## Status

Implemented.
