# Quality Control Method

Part of the canonical [DNA pipeline](pipeline.md).

**Method identifiers:** `dna.apollo_relative_quality/v1`, `dna.apollo_end_trim/v1`

Computes one bounded, uncalibrated quality value per call and selects one
retained interval. It never removes internal sequence regions.

### Substep 4.1 — Per-call penalty

For each call `i`, a window of `trim_window_size` calls centered on `i` is
examined. The penalty is the sum of two components:

- **Ambiguity penalty:** the count of calls in the window whose ambiguity symbol
  is not a canonical A/C/G/T.
- **Spacing penalty:** with `mean_spacing` the average distance between adjacent
  basecall positions across the whole read, and `min`/`max` the minimum and
  maximum adjacent spacing inside the window, the spacing penalty is
  `floor((|max - mean| + |min - mean|) / 2)`.

The penalty is `ambiguity + spacing_penalty`.

### Substep 4.2 — Best section

The best contiguous section is the window of length
`max(1, floor(call_count * best_section_fraction))` with the minimum summed
penalty. Its average penalty is recorded.

### Substep 4.3 — Relative quality score

Scores are uncalibrated and bounded. Let `max_penalty` be the largest penalty in
the read. If `max_penalty <= 0`, every call receives
`max_relative_quality_score`. Otherwise each call receives

```text
floor(max_relative_quality_score * (1 - penalty / max_penalty))
```

clamped to `[0, max_relative_quality_score]`. These scores are **not** Phred
calibrated; `phred_calibrated` is always `false`.

### Substep 4.4 — End trimming

The trim threshold is `trim_stringency * best_average * trim_window_size`.
Starting from the best section, the algorithm walks outward and stops when a
window's summed penalty exceeds the threshold, producing `trim_start` and
`trim_end`. The retained interval must contain at least
`minimum_retained_bases` calls, otherwise analysis fails. The retained sequence
is `primary_sequence[trim_start..trim_end]`.

### Substep 4.5 — Per-call record

Each call records its penalty, relative quality score, and optional vendor
quality. `vendor_quality_applies` is true only when a vendor quality exists and
the vendor primary agrees with the signal primary. Retention is represented once
by the global trim interval rather than duplicated per call.
