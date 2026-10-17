# Model

Owns the facade's records:
- `read_observation`: a Sanger read after the core, pairing the `CalledRead`
  with its `SangerAttachment`, and a rejected read with its cause
  (`RejectedRead`, `RejectionCause`, public through `sample_analysis`);
- the serializable result contracts: `result` and `basecall_result` for
  analysis and basecalls, `sample_result`, `variants_result`,
  `notation_result`, `consensus_result`, and `review_result`.

Shared, core, and Sanger vocabulary lives in the plugin crates' `model`
modules ([crates](../../crates/README.md)).

See [system architecture](../../docs/architecture/overview.md) and
[system invariants](../../docs/architecture/invariants/README.md).
