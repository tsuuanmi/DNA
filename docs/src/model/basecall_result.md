# `src/model/basecall_result.rs`

## Purpose

Defines the serializable `dna.basecalls/v2` reference-free result contract.

## Responsibilities

- Represent input and configuration identities; software/build provenance is
  intentionally deferred.
- Represent complete primary, ambiguity, and retained sequences with call count
  and trim bounds.
- Reuse shared interval and signal-quality result types, including trace-integrity evidence.
- Expose unresolved-primary and multi-channel-unresolved counts, vendor
  disagreement counts, PLOC/vendor cardinality mismatch counts, and exact
  clipped-channel-sample counts.

## Non-responsibilities

No scientific computation, filesystem access, serialization logic, reference,
alignment, variant, per-call peak, vendor, or calibrated-quality output.

## Dependencies

- `model::result` for shared input, interval, and signal-quality result records.
- `serde` for deterministic serialization.

## Tests

`tests/basecall.rs` and the basecall schema validator exercise the contract.
