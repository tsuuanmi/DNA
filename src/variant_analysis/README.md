# Variant Analysis

Owns the public high-level capability that converts supported sequencing input
into typed evidence-backed variant results.

The initial public entry point is Sanger AB1 via `analyze_sanger`. Sanger
filesystem loading is owned by `input::sanger`; Variant Analysis consumes the
validated trace/reference/configuration models and owns the reference-guided
one-read observation path: shared read processing, alignment, variant calling,
warning accounting, and assembly of the internal `ReadObservation`. Shared
reference-free read processing is provided by the crate-internal
`read_processing` module.

CLI and sample pipelines consume this capability through a crate-private
observation seam. Variant Analysis does not own file-backed CLI logging, JSON
projection, or result-file publication.

Public result types belong to this capability boundary and must not expose
private pipeline/report DTOs. The canonical called-variant contracts
(`Variant`, `VariantKind`, `ReferenceIdentity`, `CalledVariantSet`) live in the
core `dna::variant` module (ADR-0064); `VariantAnalysisResult::called_variants()`
projects the current Sanger result into that cross-modality boundary. New input modalities such as NGS should expose
compatible called-variant semantics without teaching downstream consumers about
source-specific implementation types.

Optional post-calling representation normalization is owned by
`variant_normalization`, and target nomenclature by `variant_nomenclature`,
as separate capabilities under ADR-0060.

See [Rust API contract](../../docs/reference/rust-api.md),
[interface architecture](../../docs/architecture/interfaces.md), and
[ADR-0058](../../docs/decisions/adr/0058-canonical-contracts-and-modular-analysis-composition.md).
