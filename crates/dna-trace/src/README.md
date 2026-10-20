# dna-trace

Owns trace signal engineering
([ADR-0084](../../../docs/decisions/adr/0084-trace-crate.md),
[PROP-0012](../../../docs/proposals/0012-trace-crate.md)): numerical work on
the trace windows a modality flags, from plain data. It depends on
`dna-kernel` only, and `dna-sanger` alone depends on it. It is a library, not
a plugin: it has no plugin descriptor and no configuration section, and the
modality's plugin version covers its output.

- `decompose.rs` — the phase mixtures of a read (SRS-CALL-023,
  [ADR-0082](../../../docs/decisions/adr/0082-phase-mixtures.md)):
  `mixtures` fits every dephased call range that a repeat run opens as the
  read's primary calls superposed with copies of themselves offset by whole
  calls, by non-negative least squares, and gives each stretch's components
  (offset, fraction), its residual, its unexplained calls, and whether the
  components explain it. Its input, `TraceCall`, is one call's primary base,
  whether one channel qualifies, its normalised A/C/G/T profile, and whether
  the modality keeps it in phase.

The crate reads no instrument format, reference, alignment, or truth.

Output contract: none; the modality turns its results into
`ReadEvidence::mixtures`, published in
[`dna.read_evidence/v1`](../../../docs/reference/schemas/read-evidence-v1.schema.json).

See [design](../../../docs/design/length-mixture.md) and
[crates](../../README.md).
