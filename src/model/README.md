# Model

Owns the facade's records:
- `read_observation`: a Sanger read after the core, pairing the `CalledRead`
  with its `SangerAttachment`, and a rejected read with its cause
  (`RejectedRead`, `RejectionCause`, public through `sample_analysis`);
- `sanger_read`: one read's Sanger evidence as the documents report it
  (`SangerRead`: trim, callability, signal quality, warning counts, and each
  call's scan, window peaks, primary-event peaks, relative quality, and
  second peak with its cause), and a rejected read as they report it
  (`Rejection`: its cause and, when it reached callability, its integrity and
  callability), projected from a trace's Sanger products (the one projection
  of Sanger products into the documents' records: callability, signal
  quality, integrity, warning counts, call signals, rejections) or read back
  from a read evidence document;
- the serializable result contracts: `result` and `basecall_result` for
  analysis and basecalls, `evidence_result` (`dna.read_evidence/v1`), `sample_result`, `run_signal_result` (`dna.signal/v1`), `variants_result`,
  `notation_result`, `consensus_result`, and `review_result`. `result` also
  holds the records the other contracts share (intervals, identities,
  callability, a call's placement, reference column, and trust scope).

Shared, core, and Sanger vocabulary lives in the plugin crates' `model`
modules ([crates](../../crates/README.md)).

See [system architecture](../../docs/architecture/overview.md) and
[system invariants](../../docs/architecture/invariants/README.md).
