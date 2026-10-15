# score

Children of `../score.rs`, the scoring of the draft triage.

- `minor.rs` — for the reported policy and each policy with another
  `region.minor_signal` rule, against the same policy without the rule
  (decided once per such policy, like any other): corpus-weighted A
  without and with the rule, the samples it alone moves from pass to review
  and the burden among them, leaks under every measure, and its recall of
  truth IUPAC positions; `triage_minor.tsv` and the minor-signal section of
  `triage.md`.

The measures are specified in the
[evaluation harness](../../../../../docs/validation/evaluation-harness.md#draft-triage)
document. Tests use synthetic evidence only.
