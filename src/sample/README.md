# Sample

Owns deterministic aggregation of independently processed reads into sample-level
coverage, overlap, locus, nucleotide-profile, and normalized-variant evidence.

The entry point is [`aggregate/`](aggregate/README.md). Its sibling modules
separate coverage, overlap, locus aggregation, call evidence, contribution
policy, profile geometry, nucleotide support, and variant aggregation.

The module is modality-neutral
([ADR-0069](../../docs/decisions/adr/0069-plugin-first-modality-core-post-calling.md)).
It aggregates `CalledRead` records (read identity, `ReadEvidence`, alignment,
variants) and imports no Sanger type. Locus observations and variant calls keep
their source call index, so the report can join Sanger quality, peaks, and
noisy-region context. A call whose evidence carries a mask observes a locus as
`masked`: it never retains the locus by itself and adds no nucleotide mass.

This module does not discover input files, infer samples from filenames, or emit
final reports.

`config.rs` owns the `[sample_reconciliation]` configuration section: its raw record, validated
record, and rules (ADR-0069).

See [sample requirements](../../docs/requirements/sample-evidence/README.md),
[sample methods](../../docs/design/sample-evidence/README.md), and
[sample invariants](../../docs/architecture/invariants/sample-boundaries.md).
