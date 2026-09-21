# `src/alignment/orient.rs`

## Purpose

Builds post-trim forward/reverse evidence-profile queries, selects one uniquely
best orientation, and projects circular reference coordinates.

## Responsibilities

- Validate signal-locus and QC call cardinality before alignment.
- Slice profiles by the QC trim interval.
- Build the reverse candidate by reversing profile order and complementing A/T and C/G weights.
- Map oriented query indexes back to original call indexes.
- Compare forward/reverse candidates by fixed-point profile score only.
- Reject exact orientation-score ties and modulo-distinct placement ties.
- Apply the existing primary-sequence callable-base and identity admission gates.
- Derive linear/circular reference segments and origin-wrap state.

## Important semantic separation

Profile evidence determines placement score. The retained primary sequence remains
attached to traceback columns and still defines `callable_columns`,
`callable_identity`, exact/mismatch counts, and unresolved-query counts.

Primary-sequence metrics do not break an exact forward/reverse profile-score tie.

## Circular origin

Circular references still align to a doubled working sequence. Selected
coordinates are projected modulo the original reference and split into two
segments when crossing the origin.

## Traceability

ADR-0029; `SRS-ALN-005`, `SRS-ALN-006`, `SRS-ALN-010`, and `SRS-ALN-011`.

## Status

Implemented.
