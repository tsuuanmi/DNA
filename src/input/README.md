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

[`run_context.rs`](run_context.rs) reads and validates a run context
(`dna.run_context/v1`, PROP-0014): size, schema, the length method, noise,
edges, windows, runs 1 and 2 with histograms that sum to their pairs, and
read identities as SHA-256, recording the document's own SHA-256.

[`signal.rs`](signal.rs) reads and validates a run signal (`dna.signal/v1`,
PROP-0009) for `review --signal`: at most 64 MiB, the schema, the Sanger
plugin's call selection, every key's context of `2 × flank + 1` canonical
bases centred on its base with another second base, shares within the
resolution, and counts that add up to three times the calls, recording the
document's own SHA-256. The Sanger plugin refuses a signal whose flank or
resolution differs from the configuration's
(`dna_sanger::run_signal::RunSignal::check`).

[`consensus.rs`](consensus.rs) reads a `dna.consensus/v1` document for `call`
(PROP-0006): the sample, the reference identity, the segments, and the sites
with their runs, back into the core's consensus model, under the same
configuration, profile, and reference rules; the segments are named as the
consensus FASTA names them.

[`evidence.rs`](evidence.rs) reads `dna.read_evidence/v1` documents for
`analyze --evidence` (one read, SRS-API-018) and for `sample`, `consensus`,
and `review` with `--evidence` (a sample's reads, SRS-API-019), under the
same configuration, profile, and reference rules: closed shape, at most
16 MiB, written by this build's evidence plugins under a configuration
whose sections of those plugins equal the loaded configuration's (each
plugin's configuration digest, not the file's SHA-256), calls in call order. A document records a read or a read the Sanger
plugin set aside (`Recorded`). For a read it maps every reason label onto the
Sanger plugin's own reasons (`dna_sanger::read_processing::reason`; an
unknown label fails typed) and rebuilds the kernel's `ReadEvidence` through
its constructors, with the `SangerRead` the documents report (trim,
callability, signal quality, warnings, and each call's scan, window peaks,
primary-event peaks, relative quality, and second peak with its cause),
checked against the calls; the read is named by the trace file name the
document records. For a set-aside read it reads the `Rejection` the sample
documents report, checked against the loaded configuration's minimums; the
analysis of such a document fails with the Sanger error the analysis of its
trace fails with.

[`variants.rs`](variants.rs) reads a `dna.variants/v1` document for the
post-calling plugins: the sample, the reference identity, and each read's
eligible variants, under the same configuration, profile, and reference rules.

No generic input-adapter trait exists yet. A shared trait should be introduced
only when another modality such as NGS creates a real substitutability need.

See [input requirements](../../docs/requirements/input.md),
[interface architecture](../../docs/architecture/interfaces.md), and
[ADR-0058](../../docs/decisions/adr/0058-canonical-contracts-and-modular-analysis-composition.md).
