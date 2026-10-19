# Run signal

Owns the second-channel signal a run's reads share (PROP-0009, SRS-CALL-019):
`RunSignal`, its `Key` (the read-strand context of the read's own primary
calls around a call, `run_signal_flank` calls on each side, the call's base,
and its second base), the per-key histograms of the second-channel share in
`run_signal_resolution` steps, and the one call selection,
`contributions`, that building a run signal and taking a read out of it both
apply (`SELECTION` versions it). `read` processes one trace as a sample
command does up to its core evidence and returns its contributions, or none
for a read the sample commands reject. A run signal records the flank and
resolution it was made with, and `RunSignal::check` refuses one whose values
differ from the configuration's.

It also owns the one second-channel threshold (SRS-CALL-020):
`threshold(key, signal, config)` is `signal_processing.second_floor`, raised
to the key's run level (`RunSignal::level`) when a run signal holds the key
with enough reads and calls. `Key::of` gives a call's key from its read's
evidence. Every rule that judges a second peak share uses it: the review's
`possible_heteroplasmy`. The Sanger mixed verdict keeps its own ratio: judging
it against the run signal is deferred (PROP-0009 increment iii, roadmap).

The module sits above [read_processing](../read_processing/README.md): its
selection reads the callability mask and the signal edges of the read's
`ReadEvidence`. It reads no reference and no other read.

See [dna-sanger](../README.md).
