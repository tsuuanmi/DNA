# Decompose

Owns the phase mixtures of a read (SRS-CALL-023,
[ADR-0082](../../../../docs/decisions/adr/0082-phase-mixtures.md)): `mixtures`
fits every dephased segment that a repeat run opens, from the trim
interval's start on, as the read's primary
calls superposed with copies of themselves offset by whole calls, by
non-negative least squares (`nnls`), and asserts each segment's components
(offset, fraction), its residual, its unexplained calls, and whether the
components explain it. The read's own in-phase calls set the noise that adds
an offset and the level that explains a call or the mixture; the only key is
`signal_processing.decompose_unexplained_quantile`.

A mixture decides nothing: read processing passes it to the `ReadEvidence`
adapter, which puts it on the kernel contract, and no verdict, mask, or call
reads it.

The module reads the basecalls, the locus evidence's normalised profiles,
the callability segments and mask, and the trim interval's start. It depends on
[model](../model/README.md), and [read_processing](../read_processing/README.md)
calls it. It reads no reference and no other read.

See [design](../../../../docs/design/length-mixture.md) and
[dna-sanger](../README.md).
