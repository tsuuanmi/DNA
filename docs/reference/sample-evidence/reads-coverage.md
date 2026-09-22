# Reads and Coverage

Part of the canonical [sample evidence contract](README.md).

## Reads and post-trim coverage

`reads[]` is the one registry of contributing reads. Records are sorted by
SHA-256, so CLI argument order does not change the scientific document.

Each record contains:

- `name`: the AB1 filename stem, for example
  `D11_20260404_LN_26_AB0442_HV1F_11`;
- `sha256`: stable content identity;
- `integrity`: the same concise PLOC/vendor cardinality, PLOC-spacing,
  exact-clipping, and event-signal-scale observations retained by the one-read
  pipeline; this evidence remains read-local and does not by itself admit/reject
  a read;
- `alignment`: the evidence-derived orientation, callable-base count and
  identity, unresolved-base count, gap-open count, mapped reference segments, and
  origin-wrap state.

The scientific pipeline trims each read before alignment. `reference_segments`
therefore describe where the **retained post-trim sequence** aligned on the
reference, not the untrimmed raw call span.

`reference_segments` are 0-based half-open. For a segment
`{"start": S, "end": E}`, the covered 1-based biological positions are
`S + 1` through `E`, inclusive.

For a circular reference, a read can cross the reference origin. In that case
`wraps_origin` is `true` and `reference_segments` contains two segments:
one from the mapped start to the end of the reference and one from reference
position 0 to the mapped end. A normal non-crossing read has
`wraps_origin: false`.

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
