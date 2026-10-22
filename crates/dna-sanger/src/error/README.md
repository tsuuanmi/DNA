# Error

Owns the Sanger modality's typed failures: `SangerError`, the plugin's
`Result<T>` alias, and one failure vocabulary per Sanger stage (`AbifError`
with its `Tag`, `LocusWindowError`, `BasecallingError`, `SignalError`,
`CallabilityError`, `QualityControlError`), plus `ProductsError` for one
read's Sanger products disagreeing while its `ReadEvidence` is assembled or
a reported call's Sanger signal is looked up, `RunSignalError` for a
`dna.signal/v1` document, and `ReadEvidenceError` for a `dna.read_evidence/v1`
document read back (PROP-0011). Configuration and evidence
validation failures stay the kernel's and travel as `SangerError::Kernel`.

Every variant renders with the stage prefix it had before the modality owned
it (`"invalid ABIF input: …"`, `"alignment failed: …"` for a products
mismatch), so CLI messages, operational logs, and evaluation failure text do
not depend on which crate owns a failure. The `dna` facade wraps
`SangerError` beside the kernel's error (PROP-0008 step 4, part 2c).

Like the kernel's error module, this module is a dependency leaf of the
plugin: stage modules depend on it, never the reverse.

See [Rust public API](../../../../docs/reference/rust-api.md#errors) and the
[kernel's errors](../../../dna-kernel/src/error/README.md).
