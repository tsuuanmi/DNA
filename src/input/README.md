# Input Adapters

Owns source-specific loading boundaries for external sequencing inputs.

The current adapter is `sanger`, which validates and loads AB1 traces, FASTA
references, and explicit DNA configuration into the existing validated internal
models. It does not own CLI output paths, overwrite protection, log naming,
serialization, or result publication.

No generic input-adapter trait exists yet. A shared trait should be introduced
only when another modality such as NGS creates a real substitutability need.

See [input requirements](../../docs/requirements/input.md),
[interface architecture](../../docs/architecture/interfaces.md), and
[ADR-0058](../../docs/decisions/adr/0058-canonical-contracts-and-modular-analysis-composition.md).
