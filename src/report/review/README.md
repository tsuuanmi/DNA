# Review

Owns the review list of one sample, `dna.review/v1` (PROP-0004): the `review`
post-calling plugin, whose version is the labelling rules' version, and the
rules that label the evidence the sample's analysis and call already hold.
It changes no call and reads no truth.

Output contract: [`dna.review/v1`](../../../docs/reference/schemas/review-v1.schema.json)
([review list](../../../docs/reference/review.md)).

- `mod.rs` — the plugin descriptor, building the document (`build`) and the
  in-process review list (`list`), the sample under review (`Sample`), and
  its items.
- `labels.rs` — label order within an item, an item's priority, and the
  summary's counts over the items.
- `lengths.rs` — the runs of the profile's anchored repeat windows and
  length tracts, judged from the reads' tract decodes and in-phase readings,
  for the `estimated_length`, `length_unconfirmed`, `length_not_recorded`,
  and `tract_length` labels; a length tract's run is noted when its
  dominant length is determined (SRS-API-020).
- `ranges.rs` — reference ranges: 1-based inclusive `[first, last]` pairs
  from alignment segments, merged, clipped to the reportable regions, and
  tested for a position.
- `reads.rs` — what each read shows over an item: its read labels (one
  read, one strand, noisy signal, phase shift, read edge, poor signal) and
  evidence, the sequence it shows between two clean flanks, whether it
  disagrees with the decision, and the stretches where it is masked or
  unreliable.
- `second.rs` — second peaks, the Sanger part of the review: a call's
  second peak judged against its run's signal (the sample's own calls taken
  out) and the configured floor, its cause, `reference_second`, and a
  possible heteroplasmy both strands confirm, `explained` when the modality
  explains every qualifying read's second peak.
- `sites.rs` — the labels of a decided site: reads that disagree across or
  within a strand, a site resting on one strand (run deletion, run indel,
  edge deletion, one-strand insertion), a long run, an undecided site, and a
  difference beside a region end not analysed.
- `test_support.rs` — helpers the unit tests share.
- `zones.rs` — the zones a reviewer opens whole: maximal stretches of one
  kind (repeat window, not callable, dephased, masked, clean) and the check
  items each touches.

The label table (names, ranks, and checks) is `LabelKind` in
[`model::review_result`](../../model/README.md).
