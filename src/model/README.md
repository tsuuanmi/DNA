# Model

Owns the facade's records:
- `read_observation`: a Sanger read after the core, pairing the `CalledRead`
  with its `SangerAttachment`, and a rejected read with its cause
  (`RejectedRead`, `RejectionCause`, public through `sample_analysis`);
- the serializable result contracts: `result` and `basecall_result` for
  analysis and basecalls, `evidence_result` (`dna.read_evidence/v1`), `sample_result`, `run_signal_result` (`dna.signal/v1`), `variants_result`,
  `notation_result`, `consensus_result`, and `review_result`. `result` also
  holds the records the other contracts share (intervals, identities,
  callability, a call's placement, reference column, and trust scope).

Shared, core, and Sanger vocabulary lives in the plugin crates' `model`
modules ([crates](../../crates/README.md)).

See [system architecture](../../docs/architecture/overview.md) and
[system invariants](../../docs/architecture/invariants/README.md).
