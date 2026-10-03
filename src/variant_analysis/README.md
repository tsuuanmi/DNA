# Variant Analysis

Owns the public high-level capability that converts supported sequencing input
into typed evidence-backed variant results.

The public Sanger boundary supports both one-shot `analyze_sanger` and a
reusable immutable `SangerAnalyzer`. `SangerAnalyzer::load` validates and
materializes one shared reference/configuration context once; each later
`analyze(trace)` validates/decodes only that ABIF trace and reuses the same
scientific context.

Sanger filesystem loading is owned by `input::sanger`; Variant Analysis
consumes the validated trace/reference/configuration models and owns the
reference-guided one-read observation path: shared read processing, alignment,
variant calling, warning accounting, and assembly of the internal
`ReadObservation`. Shared reference-free read processing is provided by the
crate-internal `read_processing` module.

`SangerAnalyzer` owns no thread pool, operational logger, or result publisher.
It is `Send + Sync` so outer validation/batch orchestration may share one
prepared context across independent trace workers without introducing
parallelism into the scientific kernel.

CLI and sample pipelines consume this capability through a crate-private
observation seam. Variant Analysis does not own file-backed CLI logging, JSON
projection, or result-file publication.

Public result types belong to this capability boundary and must not expose
private pipeline/report DTOs. `VariantAnalysisResult::called_variants()`
projects the current Sanger result into the implemented `CalledVariantSet`
cross-modality boundary. New input modalities such as NGS should expose
compatible called-variant semantics without teaching downstream consumers about
source-specific implementation types.

Optional post-calling representation normalization is owned by
`variant_normalization`. Target nomenclature remains a separate future
capability under ADR-0060.

See [Rust API contract](../../docs/reference/rust-api.md),
[interface architecture](../../docs/architecture/interfaces.md), and
[ADR-0058](../../docs/decisions/adr/0058-canonical-contracts-and-modular-analysis-composition.md).
