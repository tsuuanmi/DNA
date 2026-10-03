# Input Adapters

Owns source-specific loading boundaries for external sequencing inputs.

The current adapter is [`sanger`](sanger/README.md), which validates and loads
Sanger sequencing evidence, FASTA references, and explicit DNA configuration
into validated internal models. The current Sanger source format is ABIF;
`.ab1` is a common filename/extension rather than the scientific contract. It does not own CLI output paths, overwrite protection, log naming,
serialization, or result publication.

No generic input-adapter trait exists yet. A shared trait should be introduced
only when another modality such as NGS creates a real substitutability need.

See [input requirements](../../docs/requirements/input.md),
[interface architecture](../../docs/architecture/interfaces.md), and
[ADR-0058](../../docs/decisions/adr/0058-canonical-contracts-and-modular-analysis-composition.md).
