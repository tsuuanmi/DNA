# Second-peak cause

Owns why a call shows a second peak (SRS-CALL-021): `second_base`, the one
computation of a call's second base and its share (SRS-CALL-015), which the
`ReadEvidence` adapter also uses, and `classify`. `classify` labels each call
whose second share exceeds `signal_processing.second_floor`, or whose
callability state is mixed, from the peak's shape against the call's
neighbours. The labels are `noise`, `neighbour_tail`, `off_centre`,
`dye_blob`, `phase_shift`, `second_template`, and `mixed_stretch`.
`run_signal::refine` turns a `second_template` into `context` given a run
signal. That function sits in [run_signal](../run_signal/README.md), so this
module does not depend on it.

The labels decide nothing: read processing keeps them on the read's
`SangerAttachment` (`second`), beside the core's records, and no verdict
reads them.

The module reads the basecalls, the locus evidence's peak widths, the
callability mask, and the `[signal_processing]` configuration. It depends on
[model](../model/README.md) and [signal_processing](../signal_processing/README.md),
and [read_processing](../read_processing/README.md) calls it. It reads no
reference and no other read.

See [design](../../../../docs/design/second-peak-cause.md) and
[dna-sanger](../README.md).
