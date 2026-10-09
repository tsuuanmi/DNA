# Quality Control Method

Part of the canonical [DNA pipeline](pipeline.md).

**Method identifiers:** `dna.apollo_relative_quality/v1`, `dna.callable_span_trim/v1`

Computes one bounded, uncalibrated quality value per call and derives one
retained interval from the read's [callability](callability.md). Masked
segments inside the interval stay in the retained sequence; alignment keeps
dephased calls with their evidence profile and presents every other masked call
as unresolved ([alignment](alignment.md)).

### Substep 5.1 — Per-call penalty

For each call `i`, the window of `penalty_window_size` calls that starts
`penalty_window_size / 2` calls before `i` (shifted right at the read start and
truncated at the read end) is examined. The penalty is the sum of two components:

- **Ambiguity penalty:** the count of calls in the window whose ambiguity symbol
  is not a canonical A/C/G/T.
- **Spacing penalty:** with `mean_spacing` the average distance between adjacent
  basecall positions across the whole read, and `min`/`max` the minimum and
  maximum adjacent spacing inside the window, the spacing penalty is
  `floor((|max - mean| + |min - mean|) / 2)`.

The penalty is `ambiguity + spacing_penalty`.

### Substep 5.2 — Relative quality score

Scores are uncalibrated and bounded. Let `max_penalty` be the largest penalty in
the read. If `max_penalty <= 0`, every call receives
`max_relative_quality_score`. Otherwise each call receives

```text
floor(max_relative_quality_score * (1 - penalty / max_penalty))
```

clamped to `[0, max_relative_quality_score]`. These scores are **not** Phred
calibrated; `phred_calibrated` is always `false`.

### Substep 5.3 — Trim interval

The trim interval starts from the callable span `[s, e)`, the hull of the
read's unmasked calls. A masked segment that touches the span on either side
extends the interval by up to `variant_calling.read_end_margin` of its calls
when it is `dephased`, because a dephased call still reads the main ladder and
anchors the alignment next to the span; a `mixed`, `weak`, or `irregular`
neighbour carries no alignment information and is not kept. The read-end
margin of [variant eligibility](variant-calling.md) is measured from this
interval, so a call close to an uninformative end cannot support a variant.
The retained sequence is `primary_sequence[trim_start..trim_end]`.

A read with fewer than `callability.minimum_callable_calls` unmasked calls fails
in the callability stage before this substep; the span is therefore never
empty here, and a span outside the read is a typed failure.

### Substep 5.4 — Per-call record

Each call records its penalty, relative quality score, and optional vendor
quality. `vendor_quality_applies` is true only when a vendor quality exists and
the vendor primary agrees with the signal primary. Retention is represented once
by the global trim interval rather than duplicated per call.
