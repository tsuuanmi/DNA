# Run signal

Owns the second-channel signal a run's reads share (PROP-0009, SRS-CALL-019):
`RunSignal`, its `Key` (the read-strand 5-mer of the read's own primary calls
around a call, the call's base, and its second base), the per-key histograms
of the second-channel share in thousandths, and the one call selection,
`contributions`, that building a run signal and taking a read out of it both
apply (`SELECTION` versions it). `read` processes one trace as a sample
command does up to its core evidence and returns its contributions, or none
for a read the sample commands reject.

The module sits above [read_processing](../read_processing/README.md): its
selection reads the callability mask and the signal edges of the read's
`ReadEvidence`. It reads no reference and no other read.

See [dna-sanger](../README.md).
