# evidence

Children of `../evidence.rs`, the compact triage evidence of
`dna-eval run --triage`.

- `minor.rs` — the minor-channel evidence of the draft `region.minor_signal`
  detectors: per region position, each read that calls a canonical base
  without a mask, its strongest other channel (kept from 0.05), whether that
  base is the base the read calls just before or after (a shadow), how many
  of the read's nearby calls show each of those shadows, and its local
  noise. It reads every covered position through
  `SampleAnalysis::observations` (once over the span of the regions) and the
  read's own calls from its Sanger evidence. When they cannot be read, the
  stage-1 evidence is kept without minor channels and `Evidence::minor_failure`
  says why; a policy with a minor rule then reviews the would-pass regions.

The evidence is specified in the
[evaluation harness](../../../../docs/validation/evaluation-harness.md#draft-triage)
document. Tests use synthetic observations only.
