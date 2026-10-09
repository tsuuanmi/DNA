# Changelog

All notable user-visible and production-readiness changes are recorded here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
DNA uses semantic versions for tagged releases, while scientific JSON contracts
remain independently versioned and are never silently changed in place.

## [Unreleased]

### Added

- Post-calling notation of a variants document (ADR-0069 phase 4, PROP-0002,
  SRS-IN-014, SRS-OUT-011, SRS-NOM-016, SRS-PRF-008). The command is
  `dna notation <sample-id> <variants.json> --reference <reference.fasta>`.
  - It reads a `dna.variants/v1` document and runs normalization and
    nomenclature over each read's eligible variants.
  - It writes `results/<sample-id>.notation.json` (`dna.notation/v1`): the
    per-read notation (equal to the in-process notation), the source
    document's identity, and the findings of the new `conformance` plugin.
  - That plugin reports represented calls that break notation conventions
    declared in a new profile `[conformance]` section
    (`insertion_at_run_end`, `insertion_matches_run`); it never changes a
    call.
  - The shipped human-mtDNA profile declares both rules, so its file SHA-256
    recorded in result provenance changes.
  - **Rust API:** `Error` gains `Variants(VariantsError)` and
    `VariantsParse`.
- Core-only calls from reviewed consensus sequences (ADR-0069 phase 3,
  PROP-0002, SRS-IN-013, SRS-OUT-010, SRS-VAR-015). The command is
  `dna call <sample-id> <sequences.fasta>... --reference <reference.fasta>`.
  - Each FASTA record becomes modality-neutral `ReadEvidence` through the new
    `sequence` modality plugin: one-hot profiles, shared-weight profiles for
    IUPAC codes, no masks or vetoes, and vouched read ends.
  - The core caller alone places each record and calls its variants into
    `results/<sample-id>.variants.json` (`dna.variants/v1`), with per-read
    alignment summaries, observed variants with eligibility, and optional
    notation.
  - Sanger results are unchanged.
  - **Rust API:** `Error` gains `Sequence(SequenceError)`, and
    `analyze_sanger` opens a `read_evidence` stage span, for the Sanger
    evidence adapter, between `quality_control` and `alignment`.
- Signal-derived read callability (ADR-0067, SRS-CALL-001 to SRS-CALL-008): a
  new `callability` stage derives, from each read's own signal in trace order,
  phase-state segments (`in_phase`, `dephased`, `mixed`, `weak`, `irregular`),
  a per-position mask, and the callable span, and publishes them as
  observation. Trim bounds, alignment, warnings, and variant eligibility are
  unchanged in this increment. **Breaking:** `dna.basecalls/v2` → `v3`,
  `dna.analysis/v8` → `v9`, and `dna.sample_evidence/v9` → `v10` gain the
  required `read.callability` / `reads[].callability` view; configuration
  schema 6 → 7 adds `[callability]` (`window_calls`, `onset_defect_fraction`,
  `exit_defect_fraction`, `minimum_main_share`, `maximum_far_share`,
  `minimum_shadow_share`, `weak_amplitude_fraction`) and
  renames `quality_control.trim_window_size` to `penalty_window_size`.
  **Breaking (Rust API):** `Error` gains `Callability(CallabilityError)`, and
  `analyze_sanger` opens a `callability` stage span between
  `signal_processing` and `quality_control`.
- Target profiles (ADR-0063, SRS-PRF-001 to SRS-PRF-007): target knowledge —
  reference identity and topology, reportable regions, indel placement,
  nomenclature windows and their rules, notation style — moves out of code and
  configuration into strict versioned TOML profiles. The shipped profile is
  `config/profiles/human-mtdna-rcrs.toml`; a profile that pins a reference fails
  closed on another sequence. Scientific output for human mtDNA is unchanged.
  **Breaking:** configuration schema 6 gains a root `profile` key and loses
  `[reference]` and `variant_calling.regions`; `dna.analysis/v8` and
  `dna.sample_evidence/v9` record `provenance.profile` (`id`, file `sha256`);
  sample `notation.policy` becomes `notation.style = "per_base_decimal"`;
  the variant exclusion reason `outside_configured_region` becomes
  `outside_target_region`.
  **Breaking (Rust API):** new `dna::profile::{Profile, ProfileIdentity}`;
  `variant_nomenclature::mtdna::apply_control_region` is replaced by
  `variant_nomenclature::apply(reference, &profile, input)`;
  `NormalizationPolicy::MtDnaRightAligned` is renamed `RightAligned`;
  `VariantAnalysisResult` gains `profile`; `Error` gains `Profile` and
  `ProfileParse`; `ConfigError::RegionOutOfBounds` moves to `ProfileError`;
  `NomenclatureError` window names are `String`s.
- Read callability (ADR-0062, SRS-VAR-013): variants with a call near either
  read end or right after a long homopolymer are ineligible with reasons
  `read_end` / `post_homopolymer`. **Breaking:** configuration schema 6 adds
  `variant_calling.read_end_margin`, `homopolymer_min_length`, and
  `post_homopolymer_window`.
- `dna.sample_evidence/v9`: when the target profile declares notation, an
  optional `notation` view publishes each read's eligible calls after the
  profile's right alignment and nomenclature windows (HVS-II, HVS-III, HVS-I for
  human mtDNA), rendered per base (`73G`, `249DEL`, `309.1C`) with supporting
  reads (SRS-NOM-010 to SRS-NOM-012). All other v8 fields are unchanged.
- Production-oriented CI/CD and supply-chain verification, including strict cargo-shear dependency/source hygiene and an explicit Ubuntu 24.04 runner baseline.
- Declare the crate proprietary (`license = "LicenseRef-Proprietary"`) so SBOMs
  and audits record its ownership.
- Crate-wide Rust lint policy in `Cargo.toml` (`missing_docs`, `unreachable_pub`,
  Clippy `pedantic`) and a tuned release profile (thin LTO, one codegen unit).
- Rust/Python source-boundary enforcement.
- Dependency policy, dependency review, CodeQL, OpenSSF Scorecard, fuzzing,
  self-contained explicit-target auditable release bundles with post-strip metadata verification, release SBOMs, and artifact attestations.

### Changed

- Plugin-first architecture, phase 5a (PROP-0002): the sources are untangled
  along the future crate boundaries. Each plugin owns its configuration
  sections and its plugin descriptor; the facade composes the configuration
  envelope, the plugin registry, and the workflow compositions. The module
  validator enforces the crate map (kernel, core, sanger, post, dna) instead of
  layers. Configuration files, result documents, and the public Rust API are
  unchanged.

- Plugin-first architecture, phase 2b (ADR-0069, PROP-0002): sample
  aggregation is modality-neutral. It consumes each read's `CalledRead`
  (identity, evidence, alignment, variants), and the sample report joins the
  read's Sanger attachment by read identity and call index. Every result
  document is unchanged.
- Plugin-first architecture, phase 2 (ADR-0069, PROP-0002): a static plugin
  registry (`sanger` modality, `core`, and the `normalization` and
  `nomenclature` post-calling plugins) and the plugin composition of each
  workflow are validated at compile time. Each configuration section has one
  owning plugin. Scientific results are unchanged.
  **Breaking (unreleased, revised in place):**
  - configuration schema 7 moves `minimum_peak_height` and
    `relative_quality_threshold` from `[variant_calling]` to a new
    Sanger-owned `[sanger_evidence]` section;
  - `dna.analysis/v9`, `dna.basecalls/v3`, and `dna.sample_evidence/v10`
    require `provenance.plugins`: the workflow's plugins in execution order,
    each with `id`, `family`, and method `version`.
- Plugin-first architecture, phase 1 (ADR-0069, PROP-0002, SRS-VAR-014):
  alignment and variant calling consume a modality-neutral per-read evidence
  contract (`ReadEvidence`) instead of Sanger types. The Sanger adapter in
  `read_processing` owns the peak, relative-quality, and mixed-signal vetoes
  and the phase-state mask reasons. The module validator rejects any dependency
  from a neutral module on a Sanger module. Every result document is unchanged.
  **Breaking (Rust API):** `Error` gains `Evidence(EvidenceError)`.

- Sample variants publish opposition evidence (ADR-0068, SRS-SAMPLE-028): the
  admitted reads whose callable reference segments cover the variant's evidence
  span without supporting it, with orientation counts. It is evidence only;
  variant calling and notation are unchanged. **Breaking (unreleased, revised
  in place):** `variants[].opposition` is required in
  `dna.sample_evidence/v10`.
- The read-callability mask acts (ADR-0067 increment 2, SRS-CALL-009 to
  SRS-CALL-012). The trim interval is the callable span plus up to
  `read_end_margin` adjacent dephased calls (`dna.callable_span_trim/v1`
  replaces the best-section walk-out). Alignment keeps dephased calls with
  their profile and aligns other masked calls as unresolved. Masked evidence
  calls exclude a variant with the reason of their segment, replacing the
  sequence-only `post_homopolymer` window, and a masked call is no SNV
  candidate. Masked sample observations are kept as `masked` and never retain a
  locus. A mask that starts in the window after a repeat run is pulled back to
  the run's end. `read_end` is measured from every uninformative call,
  including internal unresolved masked segments, and the shipped
  `weak_amplitude_fraction` rises to 0.2. On the local corpus precision rises
  from 0.949 to 0.979, recall from 0.970 to 0.977, and the reads that failed
  the identity gate now
  place. **Breaking (unreleased, revised in place):** configuration schema 7
  removes `quality_control.best_section_fraction`, `trim_stringency`,
  `minimum_retained_bases`, and `variant_calling.post_homopolymer_window`,
  moves `variant_calling.homopolymer_min_length` to
  `callability.repeat_min_length`, adds `callability.minimum_callable_calls`,
  and ships `read_end_margin = 12`; `dna.analysis/v9` and
  `dna.sample_evidence/v10` alignment summaries gain `masked_bases` and
  `callable_reference_segments`; sample evidence gains the `masked` locus
  state, `support_topology.masked_reads`, the exclusion reasons
  `dephased_signal`, `mixed_signal`, `weak_signal`, and `irregular_spacing`, and
  the required `rejected_reads[]` array.
  A `sample` read with too few callable calls is recorded in the new
  `rejected_reads[]` registry, contributes nothing else, and no longer fails
  the operation. **Breaking (Rust API):** `QualityControlError` replaces
  `TooFewCalls` and `RetainedTooShort` with `InvalidCallableSpan`,
  `CallabilityError` gains `TooFewCallableCalls`, and `SampleError` gains
  `NoAdmittedReads`.
- Read callability separates dephased from mixed double-peak segments with a
  shadow model (ADR-0067 calibration amendment): a non-negative least-squares
  fit of each segment's normalized amplitudes to the read's own primary calls
  shifted by −3..+3 calls replaces the single modal shift offset, so two-sided
  slippage is dephased rather than mixed. Segment boundaries, the mask, and the
  callable span are unchanged. **Breaking (unreleased, revised in place):**
  `callabilitySegment` in `dna.basecalls/v3`, `dna.analysis/v9`, and
  `dna.sample_evidence/v10` gains `shadow_offsets`, present exactly on dephased
  segments; configuration schema 7 replaces `callability.shift_coherence` with
  `minimum_main_share`, `maximum_far_share`, and `minimum_shadow_share`.
- **Removed:** the Python batch runner `tools/python/scripts/analyze_samples.py`,
  its `SRS-BAT-*` requirements, and the batch runbook. Python in DNA is limited
  to repository checks, tests, measurement, and research (ADR-0066); running
  DNA over a corpus is the job of downstream pipelines.
- Result comparison with external call sets, reviewer truth, or other tools is
  out of DNA's scope (ADR-0065, SRS-COMPAT-003): downstream pipelines convert
  DNA's published results into their canonical sample representation and
  compare there.
- **Breaking (Rust API):** the canonical called-variant contracts move from
  `dna::variant_analysis` to the core module `dna::variant` (`Variant`,
  `VariantKind`, `ReferenceIdentity`, `CalledVariantSet`). Modules now form
  five acyclic layers (ADR-0064), enforced in CI by
  `validate_module_layers.py`, so a later crate split needs no contract redesign.
- **Breaking (Rust API):** `variant_nomenclature::mtdna::apply_hv2_polyc` is
  replaced by `apply_control_region`, which adds HVS-II `311T 315.1C`,
  `310C 315DEL` and EMPOP `315.1C`, HVS-III `513A 523DEL 524DEL`, and HVS-I
  `16183C 16184A 16189C` / `16189C 16193DEL` (SRS-NOM-013 to SRS-NOM-015).
  `NomenclatureError` window failures name their window.
- **Breaking (Rust API):** the empty public `dna::config` and `dna::model`
  modules are now private; library-produced results (`VariantAnalysisResult`,
  `VariantNormalizationResult`, `VariantNomenclatureResult`) and the enums
  `VariantKind`, `NormalizationPolicy`, `cli::Command`, and `error::Error` are
  `#[non_exhaustive]`.
- Operational logging uses `tracing`: scientific stages emit structured events
  and per-stage spans, and a DNA layer renders them in the unchanged
  per-operation record format (ADR-0061). Library callers can observe stage
  progress through their own `tracing` subscriber.
- **Breaking (Rust API):** `dna::error::Error` stage variants now wrap typed,
  `#[non_exhaustive]` per-stage failure enums (`AbifError`, `FastaError`,
  `ConfigError`, `AlignmentError`, and so on) instead of `String` messages;
  `DNAProcessing` is renamed `Signal`; `Path::reason` is `&'static str`; and
  `ConfigParse`/`Serialize` no longer expose `toml`/`serde_json` types. CLI and
  log error text is unchanged, except two internal-invariant messages that
  valid input cannot reach: `VariantError::NoSupportingCalls` uses the
  uppercase variant label, and a non-ASCII rendered allele reports
  `RepresentationError::InvalidAllele`.
- CI documentation now matches the workflows: CodeQL, fuzzing, MSRV, release
  build, and RustSec audit run on `main` and on schedules, not on pull requests.
- Corrected the minimum supported Rust version to match language features used
  by the codebase.

### Fixed

- An indel flanking call whose primary base is unresolved (`N`) is omitted from
  public calls (SRS-VAR-006). Previously a tied flank aborted `analyze` or a
  whole `sample` operation, and a mixed-signal flank was published with the
  schema-invalid base `N`.
- PLOC loci one sample apart, which SRS-IN-003 accepts, no longer fail
  basecalling with an empty locus window.

## [0.1.0] - Unreleased

Initial pre-production development baseline.
