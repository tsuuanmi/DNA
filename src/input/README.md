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
and reference rules (from files, or with the same record rules from an
in-memory consensus FASTA for the in-process sample call), and builds each
record's `ReadEvidence`: one-hot profiles
for bases, shared-weight profiles for IUPAC codes, no masks or vetoes, and
vouched read ends. The configuration, profile, and reference loaders in
`mod.rs` are shared by both adapters.

[`consensus.rs`](consensus.rs) reads a `dna.consensus/v1` document for `call`
(PROP-0006): the sample, the reference identity, the segments, and the sites
with their runs, back into the core's consensus model, under the same
configuration, profile, and reference rules; the segments are named as the
consensus FASTA names them.

[`evidence.rs`](evidence.rs) reads a `dna.read_evidence/v1` document for
`analyze --evidence` (PROP-0011, SRS-API-018) under the same configuration,
profile, and reference rules: closed shape, at most 16 MiB, written under the
loaded configuration by this build's evidence plugins, calls in call order.
It maps every reason label onto the Sanger plugin's own reasons
(`dna_sanger::read_processing::reason`; an unknown label fails typed) and
rebuilds the kernel's `ReadEvidence` through its constructors, with the
`SangerRead` the analysis document reports (trim, callability, signal quality,
warnings, each call's primary-event peaks and relative quality), checked
against the calls.

[`variants.rs`](variants.rs) reads a `dna.variants/v1` document for the
post-calling plugins: the sample, the reference identity, and each read's
eligible variants, under the same configuration, profile, and reference rules.

No generic input-adapter trait exists yet. A shared trait should be introduced
only when another modality such as NGS creates a real substitutability need.

See [input requirements](../../docs/requirements/input.md),
[interface architecture](../../docs/architecture/interfaces.md), and
[ADR-0058](../../docs/decisions/adr/0058-canonical-contracts-and-modular-analysis-composition.md).
