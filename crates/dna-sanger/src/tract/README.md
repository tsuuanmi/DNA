# Tract

Owns a read's homopolymer tract decodes (SRS-CALL-025,
[PROP-0012](../../../../docs/proposals/0012-trace-crate.md) increment 2):
the Sanger side of the work that
[dna-trace](../../../dna-trace/src/README.md) computes
([ADR-0084](../../../../docs/decisions/adr/0084-trace-crate.md)).

`tracts` takes each call where callability sees the read leave phase (the
start of a segment not `in_phase`). It hands `dna-trace` the read's analysed
channels at its calls' locus samples, each row divided by its sum, and turns
each decoded tract into the kernel's `TractLengths`. On each marginal it
sets the margin state: `near_equal` when the deflated margin is below
`trace.near_equal_z`, `undefined` (never clear) when two lengths show
without a finite margin, and `clear` otherwise.
A tract whose fit holds no signal is left out. The `tracts_decoded` log
event counts the onsets, the tracts, and the tracts where the iteration or
joint-length cap binds.

`config.rs` owns the `[trace]` configuration section, which the Sanger plugin
declares: its raw record, validated record, and rules.

The decodes decide nothing: the `ReadEvidence` adapter attaches them, and no
verdict, mask, call, or review item reads them. The module reads the
chromatogram, the basecalls, and the callability segments. It depends on
[model](../model/README.md), and [read_processing](../read_processing/README.md)
calls it. It reads no reference and no other read.

See [design](../../../../docs/design/length-mixture.md) and
[dna-sanger](../README.md).
