# Sample

Owns deterministic aggregation of independently processed reads into sample-level
coverage, overlap, locus, nucleotide-profile, and normalized-variant evidence.

Key children separate coverage, overlap, locus aggregation, call evidence,
contribution policy, profile geometry, nucleotide support, and variant aggregation.

This module does not discover input files, infer samples from filenames, or emit
final reports.

See [sample requirements](../../docs/srs/sample-evidence/README.md),
[sample methods](../../docs/methods/sample-evidence/README.md), and
[sample invariants](../../docs/architecture/invariants/sample-boundaries.md).
