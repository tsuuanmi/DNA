# Report

Owns projection of completed scientific state into versioned public contracts,
deterministic JSON serialization, and atomic no-overwrite publication.

Key children separate analysis, basecall, sample, signal, variant,
`per_base_decimal` notation rendering (`notation.rs`), serialization, and atomic
publication concerns. Analysis and sample provenance record the target-profile
identity.

Scientific decisions remain upstream.

See [contracts](../../docs/reference/README.md),
[output requirements](../../docs/requirements/output.md), and
[output/operations invariants](../../docs/architecture/invariants/output-operations.md).
