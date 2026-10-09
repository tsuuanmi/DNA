# ADR-0063: Target knowledge in versioned profiles

- **Status:** Accepted
- **Date:** 2026-10-08
- **Related decisions:** [ADR-0060](0060-separate-variant-canonicalization-nomenclature.md),
  [ADR-0062](0062-read-callability.md)

## Context

Human-mtDNA knowledge was spread through code and configuration: the rCRS
checksum, the HVS-II/HVS-III/HVS-I windows and their `C`/`T` bases lived in
`variant_nomenclature::mtdna`, the right-alignment policy was named
`MtDnaRightAligned`, sample notation was switched on by recognising the rCRS
checksum, and the reportable regions and reference topology sat in the method
configuration. Supporting another target (a nuclear amplicon, another
organism's mitochondrion, an NGS panel) would have meant editing Rust in
several modules.

Two kinds of setting were mixed:

- **method parameters**: how a modality is processed (peak ratios, scoring,
  trimming, read callability); and
- **target knowledge**: what is being sequenced (reference identity and
  topology, reportable regions, representation conventions).

## Decision

1. Target knowledge lives in a **target profile**: a strict, versioned TOML file
   (`schema_version = 1`) that is validated completely at load. It declares an
   `id`; the reference `topology` and optionally the `sequence_sha256` it is
   validated against; reportable `variant_calling.regions`; and optionally
   `nomenclature.windows` and the sample notation chain
   (`normalization.indel_placement` with `notation.style`, declared together).
2. The scientific configuration names exactly one profile (`profile`, resolved
   against the configuration file's directory). It keeps only method
   parameters. The shipped profile is
   `config/profiles/human-mtdna-rcrs.toml`.
3. Nomenclature windows are data. Each declares a name, a 1-based start, its
   reference sequence, a structure (`anchored_homopolymer` with repeat and
   anchor bases, or `tandem_repeat` with its motif), ordered rules, and the
   data a rule needs (`canonical_haplotype`). Validation proves each rule is
   compatible with its structure, so the engine cannot meet an invalid pairing.
   `variant_nomenclature` owns only the generic engine and its
   haplotype-preservation proof.
4. When a profile declares `sequence_sha256`, every operation that uses the
   profile with a reference fails closed on another sequence; a reference that
   does not carry every window sequence fails closed too. Sample notation
   is produced exactly when the profile declares `notation`.
5. Results record the profile `id` and the SHA-256 of the profile file
   (`provenance.profile` in `dna.analysis/v8` and `dna.sample_evidence/v9`),
   beside the configuration checksum. Basecalling does not read the profile.

## Alternatives

- **Keep mtDNA in code behind a trait per target.** Every new target still needs
  Rust, a release, and review of scientific code paths.
- **Embed the profile in `dna.toml`.** One checksum would cover both, but the
  same target knowledge could not be shared by configurations with different
  method parameters, and a profile could not be reviewed or versioned alone.
- **A general rule language.** Arbitrary user rules cannot be proven
  haplotype-preserving at load. The rule vocabulary stays a closed, validated
  set; new kinds of rule are code changes with tests.

## Consequences

- A new target with the existing structures and rules needs only a profile file.
  A new kind of structure or rule is still a reviewed code change.
- Scientific output for the human-mtDNA target is unchanged: on the 44-sample
  Sequencher set all 430 real outputs match the previous build except for the
  added profile provenance, the analysis schema version, and the notation
  `style` field.
- The configuration checksum no longer covers target knowledge; provenance must
  carry the profile checksum too.
- Release packages ship `config/profiles/` with checksums.
