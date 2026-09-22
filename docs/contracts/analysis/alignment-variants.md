# Analysis Alignment and Variants

Part of the canonical [analysis contract](../analysis.md).

## Alignment summary

The selected alignment reports:

- `orientation` (`forward` or `reverse`);
- `callable_bases` and callable `identity`;
- `unresolved_bases` and `gap_opens`;
- one or two 0-based half-open `reference_segments`;
- `wraps_origin`.

Gapped query/reference rows, operation runs, fixed-point profile score, exact-match/mismatch redundancy, and traceback columns remain internal. Placement consumes the retained post-trim `EvidenceProfile` sequence, while the retained primary sequence remains attached to traceback columns for callable/identity metrics and downstream primary-sequence variant extraction. `reference_segments` therefore describe post-trim mapped coverage.

## Variant calls

Each normalized variant contains only `position`, `reference`, `alternate`,
`kind`, and direct `calls`. The variant position is 1-based on the supplied
reference strand.

Each call is optimized for scientific review rather than implementation
traceability:

- `role`: `supporting` or `flanking`;
- `base`: called base projected onto the reference strand;
- `peaks`: raw analyzed A/C/G/T channel heights sampled together at the unique
  primary-event coordinate and projected to reference orientation;
- `quality`: the existing uncalibrated relative score under a concise public
  field name.

Original call indexes, ABIF PLOC coordinates, mapped call positions, selected-peak
positions/sources, penalties, calibration flags, and vendor scores remain internal.

For reverse reads, both `base` and `peaks` are reference-oriented, so the
reviewer can compare the normalized alternate allele directly with the signal
channels.

### SNVs

An SNV contains one or more supporting calls. The top-level variant `position`
already identifies the biological coordinate, so that coordinate is not repeated
inside each call.

### Insertions

An insertion contains one supporting call per inserted base plus available
flanking calls. Supporting-call bases and peaks show the observed inserted
sequence directly.

### Deletions

A deletion has no trace signal at the deleted reference base. It therefore
contains available flanking calls with their real peak/quality evidence and never
fabricates deleted-base dna.

For repeat-associated indels, normalization can move the reported allele
representation away from the observed alignment gap. The normalized
`position/reference/alternate` remain authoritative for the biological variant.
