# ADR-0008: Emit a Single JSON-Only, Auditable Analysis Document

- **Status:** Superseded in part by ADR-0011
- **Date:** 2026-08-22

## Context

The MVP must produce a result that is deterministic, machine-readable, and
auditable. Earlier planning (ADR-0005) anticipated both JSON and VCF 4.2 text
outputs. The implemented pipeline, however, produces exactly one output file,
`<prefix>.analysis.json`, and no VCF. The JSON document carries the full decoded
trace, every stage record, the effective configuration, and all warnings, so a
consumer can reconstruct and audit the analysis without any other artifact.

## Options

1. Emit JSON plus VCF 4.2 text as planned in ADR-0005.
2. Emit JSON only, with the complete auditable record in one document.

## Decision

Choose option 2. `signal analyze` writes a single `dna.analysis/v1` JSON
document validated against `docs/schemas/analysis-v1.schema.json`. The document
includes the full A/C/G/T channel arrays, per-call basecalling evidence, quality
and trim records, the selected alignment and both orientation summaries,
normalized primary-sequence variants, the effective configuration and checksums,
and all warnings. No VCF or other output file is written.

## Consequences

Consumers get one self-contained, schema-validated artifact that is sufficient
to audit every stage. There is no VCF surface to keep in sync, and no risk of
mislabeling text as BCF. Downstream tools that require VCF must convert from the
JSON document. If VCF is later required, it should be derived from the same
typed result and added as a new output contract, not as a parallel re-derivation.

## Supersession

This ADR supersedes the VCF portion of ADR-0005 for the MVP. A future ADR may
add a VCF output derived from the same completed analysis result.
