# Reads and Coverage

Part of the canonical [sample evidence contract](README.md).

## Reads and post-trim coverage

`reads[]` is the one registry of contributing reads. Records are sorted by
SHA-256, so CLI argument order does not change the scientific document. Reads
that could not be analyzed are listed separately in `rejected_reads[]` (see
[Rejected reads](#rejected-reads)).

Each record contains:

- `name`: the AB1 filename stem, for example
  `D11_20260404_LN_26_AB0442_HV1F_11`;
- `sha256`: stable content identity;
- `integrity`: the same concise PLOC/vendor cardinality, PLOC-spacing,
  exact-clipping, and event-signal-scale observations retained by the one-read
  pipeline; this evidence remains read-local and does not by itself admit/reject
  a read;
- `callability`: the read's `dna.read_callability/v1` view — callable span,
  ordered phase-state segments with repeat attribution and the shadow offsets of
  dephased segments, and masked-call count —
  exactly as the one-read contracts publish it ([method](../../design/callability.md));
- `alignment`: the evidence-derived orientation, callable-base count and
  identity, unresolved-base count, masked-base count, gap-open count, mapped
  reference segments, callable reference segments, and origin-wrap state.

The scientific pipeline trims each read to its callable span, plus a few
dephased context calls, before alignment. `reference_segments` therefore
describe where the **retained post-trim sequence** aligned on the reference,
not the untrimmed raw call span. `masked_bases` counts aligned calls that the
read's callability masks: dephased ones align with their call and profile,
every other masked call as unresolved, and none of them counts in
`unresolved_bases`. `callable_reference_segments` are the parts of the mapped
segments observed by unmasked calls and by deletions between them; a covered
position outside them is covered by a masked call.

`reference_segments` and `callable_reference_segments` are 0-based half-open.
For a segment
`{"start": S, "end": E}`, the covered 1-based biological positions are
`S + 1` through `E`, inclusive.

For a circular reference, a read can cross the reference origin. In that case
`wraps_origin` is `true` and `reference_segments` contains two segments:
one from the mapped start to the end of the reference and one from reference
position 0 to the mapped end. A normal non-crossing read has
`wraps_origin: false`.

## Rejected reads

`rejected_reads[]` lists reads whose callability left fewer than
`callability.minimum_callable_calls` unmasked calls (ADR-0067). Records are
sorted by SHA-256 and contain `name`, `sha256`, `integrity`, and `callability`
as in `reads[]`, plus `rejection`:

- `reason`: `callable_calls_below_minimum`;
- `callable_calls`: the read's unmasked call count;
- `minimum_callable_calls`: the configured minimum.

A rejected read has no `alignment` and contributes nothing to `coverage[]`,
`overlaps[]`, `locus_differences[]`, `variants[]`, or `notation`. Names are
unique across both registries. The array is empty when every read was
analyzed; a sample whose every read is rejected fails without a document.

Read names are unique within one emitted sample document because they are used as
human-readable references from overlap, locus, and variant evidence. SHA-256 remains the
scientific content identity. Filename semantics are never used as placement or
merge keys.

## Coverage topology

`coverage[]` is a run-length encoded summary of selected post-trim reference
coverage. Each item contains a 0-based half-open `reference` interval plus:

- `read_depth`: number of independently placed reads covering every coordinate
  in the interval;
- `forward_depth`: covering reads whose selected orientation is forward;
- `reverse_depth`: covering reads whose selected orientation is reverse.

For every item, `read_depth = forward_depth + reverse_depth`. Adjacent intervals
with identical depth tuples are merged even when the identity of the covering
read changes at the boundary; exact read identities and segments remain in
`reads[]`.

Coverage counts all independently placed reads. It does **not** remove a read
because a pairwise overlap is ineligible, and it does not mean canonical-base
agreement, nucleotide comparability, consensus confidence, or biological strand
independence. Deletion columns remain reference-coordinate coverage; insertions
do not create extra reference coordinates.

Circular origin-spanning reads contribute through their two explicit
`reference_segments`, so the linearized JSON coverage map can contain runs near
both reference ends without an implicit wrapped interval.
