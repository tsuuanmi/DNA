# dna-trace

Owns trace signal engineering
([ADR-0084](../../../docs/decisions/adr/0084-trace-crate.md),
[PROP-0012](../../../docs/proposals/0012-trace-crate.md)): numerical work on
the trace windows a modality flags, from plain data. It depends on
`dna-kernel` only, and `dna-sanger` alone depends on it. It is a library, not
a plugin: it has no plugin descriptor and no configuration section of its
own (the `sanger` plugin's `[trace]` section holds the tract decode's
constants), and the modality's plugin version covers its output.

- `decompose.rs` — the phase mixtures of a read (SRS-CALL-023,
  [ADR-0082](../../../docs/decisions/adr/0082-phase-mixtures.md)):
  `mixtures` fits every dephased call range that a repeat run opens as the
  read's primary calls superposed with copies of themselves offset by whole
  calls, by non-negative least squares, and gives each stretch's components
  (offset, fraction), its residual, its unexplained calls, and whether the
  components explain it. Its input, `TraceCall`, is one call's primary base,
  whether one channel qualifies, its normalised A/C/G/T profile, and whether
  the modality keeps it in phase.
- `locate.rs` — where a tract opens (SRS-CALL-025, PROP-0012 increment 2):
  the latest run of equal primaries before a call where the read leaves
  phase, extended back over a single-call mark between two runs.
- `tract.rs` — a tract's length decode (SRS-CALL-025): `analyse` locates each
  tract, `choose` picks its structure (a single run, or two runs around a
  mark) by the lowest residual over one shared window, and `decode` fits it
  as a mixture of molecules sharing one downstream sequence, by NNLS (`dna_kernel::nnls`) and an
  exact dynamic programme from a start per candidate length. Only the first
  run's and the total's marginals are identifiable, so only they are
  reported, with standard errors and autocorrelation-deflated margins. The shared downstream sequence is reported as far as
  the decode determines it (`TractDecode::decoded`).
  `signature`, `leakage`, and `refit` clean the minor shares with the read's
  own base signatures and neighbour leakage. Its constants come from the
  modality (`TractParams`, `LeakageParams`, `LocateParams`); the iteration
  and joint-length caps are reported when they bind.

The crate reads no instrument format, reference, alignment, or truth.

Output contract: none; the modality turns its results into
`ReadEvidence::mixtures` and `ReadEvidence::tracts`, published in
[`dna.read_evidence/v1`](../../../docs/reference/schemas/read-evidence-v1.schema.json).

See [design](../../../docs/design/length-mixture.md) and
[crates](../../README.md).
