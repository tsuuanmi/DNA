# Review

Owns the review list of one sample, `dna.review/v1` (PROP-0004): the `review`
post-calling plugin, whose version is the labelling rules' version, and the
rules that label the evidence the sample's analysis and call already hold.
It changes no call and reads no truth.

Output contract: [`dna.review/v1`](../../../docs/reference/schemas/review-v1.schema.json)
([review list](../../../docs/reference/review.md)).

- `mod.rs` — the plugin descriptor, building the document (`build`) and the
  in-process review list (`list`), and the rules over a sample's sites,
  reads, second peaks, run lengths, and items.
- `labels.rs` — label order within an item, an item's priority, and the
  summary's counts over the items.
- `ranges.rs` — reference ranges: 1-based inclusive `[first, last]` pairs
  from alignment segments, merged, clipped to the reportable regions, and
  tested for a position.
- `test_support.rs` — helpers the unit tests share.
- `zones.rs` — the zones a reviewer opens whole: maximal stretches of one
  kind (repeat window, not callable, dephased, masked, clean) and the check
  items each touches.

The label table (names, ranks, and checks) is `LabelKind` in
[`model::review_result`](../../model/README.md). The run-length rules of the
profile's anchored repeat windows are in [`report::length`](../README.md).
