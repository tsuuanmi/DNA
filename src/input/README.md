# Input Adapters

Owns source-specific loading boundaries for external sequencing inputs.

The current adapter is [`sanger`](sanger/README.md), which validates and loads
Sanger sequencing evidence, FASTA references, explicit DNA configuration, and
the target profile it names into validated internal models; a reference that is
not the profile's pinned sequence fails here. The current Sanger source format is ABIF;
`.ab1` is a common filename/extension rather than the scientific contract. It does not own CLI output paths, overwrite protection, log naming,
serialization, or result publication.

[`sequence.rs`](sequence.rs) is the sequence modality
([ADR-0069](../../docs/decisions/adr/0069-plugin-first-modality-core-post-calling.md)).
It loads reviewed consensus FASTA records under the same configuration, profile,
and reference rules, and builds each record's `ReadEvidence`: one-hot profiles
for bases, shared-weight profiles for IUPAC codes, no masks or vetoes, and
vouched read ends. The configuration, profile, and reference loaders in
`mod.rs` are shared by both adapters.

[`variants.rs`](variants.rs) reads a `dna.variants/v1` document for the
post-calling plugins: the sample, the reference identity, and each read's
eligible variants, under the same configuration, profile, and reference rules.

No generic input-adapter trait exists yet. A shared trait should be introduced
only when another modality such as NGS creates a real substitutability need.

See [input requirements](../../docs/requirements/input.md),
[interface architecture](../../docs/architecture/interfaces.md), and
[ADR-0058](../../docs/decisions/adr/0058-canonical-contracts-and-modular-analysis-composition.md).
