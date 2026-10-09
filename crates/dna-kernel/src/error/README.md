# Error

Owns the typed `Error` boundary, the shared `Result<T>` alias, and one failure
vocabulary per stage (`AbifError`, `FastaError`, `ConfigError`, `ProfileError`,
`BasecallingError`, `SignalError`, `CallabilityError`, `QualityControlError`, `AlignmentError`,
`VariantError`, `SampleError`, `ReportError`, and the representation,
normalization, and nomenclature errors), plus shared `LocusWindowError` and
`CallEvidenceError`.

`Error` names the failing stage; the wrapped stage failure carries observed
values as structured data, names fixed rules as distinct variants (some with
static rule text), and gives arithmetic overflow and internal-consistency guards
a static description. Stage failures render inline as
`"<stage prefix>: <failure>"`, which is the CLI and operational-log text.

This module is a dependency leaf: stage modules depend on it, never the
reverse, so limits travel as variant data rather than imported constants.
Errors preserve structured failure context without printing, exiting, or
deciding recovery policy.

See [Rust public API](../../../../docs/reference/rust-api.md#errors),
[system architecture](../../../../docs/architecture/overview.md), and
[quality requirements](../../../../docs/requirements/quality-attributes.md).
