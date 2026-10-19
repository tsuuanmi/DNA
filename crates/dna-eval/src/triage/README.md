# triage

Draft stage-1 triage for phase 0 of
[PROP-0004](../../../../docs/proposals/0004-review-triage.md): research-grade
detectors over the evidence DNA already produces, scored against the truth to
decide whether production triage is worth building. It runs only under
`dna-eval run --triage`. It is not a DNA contract, nothing in DNA reads it,
and its rules are drafts, never approved or active.

It reads the compact triage evidence that the runner takes per sample
(`../evidence.rs`, stage E) and the compared rows of `../report.rs`.

- `detect.rs` — stage R: the trace and region detectors, the primer
  footprints, and the region, trace, and sample outcome under one policy.
  Pure and deterministic.
- `repeats.rs` — the reference geometry of the one-strand exception: the
  repeat windows (the measures' table) with their margin, and the
  reference's homopolymer runs. Pure.
- `minor.rs` — the `region.minor_signal` rules: a kind (how a read's minor
  channel counts and how reads must agree) and a threshold; a policy joins
  one or more of them, any of which marks. Pure and deterministic.
- `score.rs` — stage C: A/B/C, leaks with Clopper–Pearson bounds, catches
  and burden per detector, region leak cases, the policy sweep with every sample's outcome under
  each policy, and the `triage*` output files; `score/minor.rs` measures each policy's minor
  rules against the policy without them (`triage_minor.tsv`, [minor score](score/README.md)).

The detectors, outcomes, measures, and outputs are specified in the
[evaluation harness](../../../../docs/validation/evaluation-triage.md#draft-triage)
document. Tests use synthetic evidence only.
