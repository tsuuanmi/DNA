# evidence

Children of `../evidence.rs`, the compact triage evidence of
`dna-eval run --triage`.

- `minor.rs` — the minor-channel evidence of the draft `region.minor_signal`
  detectors: per region position, each read that calls a canonical base
  without a mask, its strongest other channel (kept from 0.05), whether that
  base is the base the read calls just before or after (a shadow), and its
  local noise. It reads every covered position through
  `SampleAnalysis::loci` and the read's own calls from its Sanger evidence.

The evidence is specified in the
[evaluation harness](../../../../docs/validation/evaluation-harness.md#draft-triage)
document. Tests use synthetic observations only.
