# `src/variant_calling/mapping.rs`

## Purpose

Derives the original trace-call mappings for a variant from the selected alignment
columns and validates their shape.

## Responsibilities

- Build a `VariantCallMapping` for a supporting or flanking alignment column.
- Sort and deduplicate the collected call mappings.
- Validate that SNV, insertion, and deletion evidence has the correct roles and
  reference positions.

## Non-responsibilities

No extraction, normalization, or serialization.

## Key types and functions

- `supporting(column) -> Result<VariantCallMapping>`: a difference-bearing query
  column mapping to an original call; carries the column's reference position
  when present.
- `flanking(column) -> Result<VariantCallMapping>`: a reference-aligned query
  column used to bound an indel; requires a reference position.
- `sort_dedup(calls) -> Vec<VariantCallMapping>`: sorts and deduplicates stable
  call mappings.
- `validate_snv(calls, position) -> Result<()>`: every SNV call must support the
  substituted reference position.
- `validate_insertion(calls) -> Result<()>`: at least one supporting inserted call
  with no reference position, plus optional aligned flanks.
- `validate_deletion(calls) -> Result<()>`: deletion evidence consists only of
  reference-aligned flanking calls.

## Invariants and errors

- A column missing an original call index returns `Error::Variant`.
- An inserted call must not carry a reference position; a flanking call must carry
  one; otherwise `Error::Variant`.
- An insertion with no supporting call or a deletion with no flanks returns
  `Error::Variant`.

## Dependencies

- `model::alignment` for `AlignmentColumn`.
- `model::variant` for `VariantCallMapping` and `VariantCallRole`.
- `error` for `Error`/`Result`.

## Biological semantics

A mapped call is either `Supporting` (its observed base contributes the alternate
allele) or `Flanking` (a reference-aligned call that bounds an indel). Inserted
query calls have no reference base, so their biological reference position is
absent; flanking calls always carry one.

## Tests

No dedicated unit tests; behavior is exercised through `extract`, `normalize`,
and the integration tests.

## Status

Implemented.
