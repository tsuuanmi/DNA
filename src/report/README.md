# Report

Owns projection of completed scientific state into versioned public contracts,
deterministic JSON serialization, and atomic no-overwrite publication.

Key children separate analysis, basecall, sample, signal, variant,
`per_base_decimal` notation rendering (`notation.rs`), serialization, and atomic
publication concerns. Analysis and sample provenance record the target-profile
identity and the workflow's plugins.

Core records carry no Sanger evidence
([ADR-0069](../../docs/decisions/adr/0069-plugin-first-modality-core-post-calling.md)).
`call.rs` projects core-only calls into `dna.variants/v1`,
and `notation.rs` renders and projects the per-read notation shared by sample,
call, and notation documents and builds `dna.notation/v1`. `sanger_call.rs` joins a call's reference-strand peaks and
relative quality, and
`sample.rs` joins each read's Sanger attachment (integrity, callability,
quality, noisy regions) by read identity and call index.

Scientific decisions remain upstream.

See [contracts](../../docs/reference/README.md),
[output requirements](../../docs/requirements/output.md), and
[output/operations invariants](../../docs/architecture/invariants/output-operations.md).
