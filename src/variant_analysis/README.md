# Variant Analysis

Owns the public high-level capability that converts supported sequencing input
into canonical variant evidence.

The initial adapter is Sanger AB1 via `analyze_sanger`. It reuses the existing
validated read/alignment/variant scientific path but does not own CLI logging,
JSON projection, or result-file publication.

Public result types belong to this capability boundary and must not expose
private pipeline/report DTOs. New input modalities such as NGS should adapt into
the same canonical result semantics rather than teach downstream consumers about
source-specific implementation types.

See [Rust API contract](../../docs/reference/rust-api.md),
[interface architecture](../../docs/architecture/interfaces.md), and
[ADR-0058](../../docs/decisions/adr/0058-canonical-contracts-and-modular-analysis-composition.md).
