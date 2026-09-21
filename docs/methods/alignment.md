# Alignment Method

Part of the canonical [DNA pipeline](pipeline.md).

Aligns the retained basecall-independent evidence-profile sequence to the reference with affine-gap Gotoh dynamic programming. The retained primary sequence stays attached to traceback coordinates for admission metrics and downstream primary-sequence variant extraction. Alignment is semi-global: the retained query is fully consumed while unaligned reference flanks are allowed.

### Substep 5.1 — Orientation candidates

The retained query evidence is aligned in both orientations:

- **forward:** retained profile order and retained primary sequence as-is;
- **reverse:** reverse profile order with A↔T/C↔G profile complementation, plus the reverse-complemented retained primary sequence for traceback character/provenance mapping.

For a circular reference, the reference is duplicated (concatenated with itself)
so the query may wrap across the origin; the working reference length is the
modulo length. A traceback may consume at most one reference length, so a query
whose required reference span is longer than the circle is unsupported.

### Substep 5.2 — Gotoh scoring

Three dynamic-programming matrices track match, insertion, and deletion states. All score deltas use fixed-point scale 1024. For canonical reference base `r`, profile support is quantized as `u = round(weight[r] × 1024)` and substitution score is `u × match_score + (1024-u) × mismatch_score`. A missing profile or non-canonical reference base receives `1024 × ambiguous_score`. Gap open/extension deltas use the same scale, preserving `open + k × extension` semantics and preserving the clean one-hot method ordering. Endpoint candidates are ranked from the last query row, allowing free reference flanks. For a circular reference, DNA selects the highest-scoring traceback whose consumed reference span is at most one circle rather than letting an invalid unbounded candidate mask a valid placement. Allocation is bounded by a compiled cell cap.

### Substep 5.3 — Traceback

The traceback internally decodes the selected path into equal-length gapped query and
gapped reference strings, an operation-run string (e.g. `5M`, `3M1I1M`), and
alignment metrics. Compact analysis v7 emits only the selected alignment summary and
reference segments, not the rows, operation runs, or score. When multiple paths tie, a documented state order
(match > deletion > insertion) makes the result deterministic. Metrics are:

- `exact_matches`, `mismatches`, `gap_opens`;
- `callable_columns` (columns where both bases are canonical);
- `callable_identity` = `exact_matches / callable_columns` (0 when no callable
  columns);
- `unresolved_query_bases` (query `N` columns).

### Substep 5.4 — Orientation selection

The forward and reverse candidates are compared by fixed-point profile score only. The strictly better orientation is selected; an exact score tie is an error and primary-sequence exact/mismatch metrics do not break it. After placement, `callable_columns`, `callable_identity`, exact/mismatch counts, and unresolved-query count are still computed from the retained primary sequence on the selected traceback. The selected orientation must meet the existing `minimum_callable_bases` and `minimum_identity` primary-sequence admission gates, otherwise analysis fails.

### Substep 5.5 — Reference segments

For a linear reference, the alignment maps to one half-open reference segment.
For a circular reference, the aligned span is projected back onto the reference;
if it crosses the origin it is split into two segments and `wraps_origin` is
`true`.
