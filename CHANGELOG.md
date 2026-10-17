# Changelog

All notable user-visible and production-readiness changes are recorded here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
DNA uses semantic versions for tagged releases, while scientific JSON contracts
remain independently versioned and are never silently changed in place.

## [Unreleased]

Plugin versions rise with every change to a plugin's output, also before a
release (`docs/governance/versioning.md`). Two changes are recorded without
a bump of their own (2026-10-17, user: no retroactive bump): the
`nomenclature` window form resting only on a wholly read window (#77) shares
the `nomenclature` 2 of #76, and the bounded deletion trade (#79) shares the
`normalization` 2 of #78.

### Added

- The EMPOP query (SRS-NOM-018, PROP-0007 increment 2): when the profile
  declares `[empop]`, as the human-mtDNA profile does, `dna notation` adds an
  offline `empop` query to `dna.notation/v1`, and `SampleAnalysis::call`
  gives it as `ConsensusCall::empop`. It is sent nowhere. The query has
  interpretation ranges over the positions exactly one record calls (no `N`,
  no excluded variant), the represented calls inside them with an insertion
  written by its anchor, as EMPOP numbers it, and the record positions cut,
  with their reason. New post-calling plugin `empop` 1 in `dna-empop`.
  **Breaking:** `dna.notation/v1` adds the optional `empop`, revised in place
  (unreleased); the profile gains `[empop]`. The evaluation harness writes
  each query to `scores.tsv` (`empop_ranges`, `empop_mitotype`), with the
  positions inside its ranges that the consensus left undecided or decided by
  a tie (`empop_undecided`, `empop_ties`), and summarises them.

- Every `technology` finding of the evaluation harness's discrepancy causes
  is deferred, to signal processing in `dna-sanger`: the core is judged on
  good input, or input the lab made poor. A defect heads a sample before a
  deferred item, so a sample with both a `lab` and a `technology` finding is
  the lab's. The summary's cause table lists defects first, then deferred
  items, each in cause order (`technology (deferred)` replaces
  `technology`). Statuses and classes are unchanged.

- Report conventions (SRS-NOM-017, PROP-0007 increment 1): a profile's
  `[[nomenclature.conventions]]` labels the rendered calls the reviewers
  write differently, and keeps them. The human-mtDNA profile declares
  `hvs1_c_length` (`not_recorded`): the length calls of the uninterrupted
  HVS-I C run that `16189C` creates (C10 and longer). Every notation view
  carries the label as `convention` on a labelled call, and
  `SampleAnalysis::call` carries it on `NotationCall`. The `nomenclature`
  plugin is version 3. **Breaking:** `dna.variants/v1`,
  `dna.sample_evidence/v10`, and `dna.notation/v1` add the optional
  `notation.calls[].convention`, revised in place (unreleased). The
  evaluation harness gives a `sequence` difference made only of such calls,
  with no truth-only call, the cause `lab`, deferred
  (`convention_not_recorded`): 31 samples on development panel v5, all
  previously `undetermined`. The summary's cause table adds a
  `lab (deferred)` column.

- The evaluation harness's haplotype measures compare only inside the
  truth's analysed intervals: a DNA call where the reviewers did not analyse
  no longer counts against DNA there (the primary, strict, and
  IUPAC-compatible measures, the lab's rules, are unchanged). A sample whose
  only primary differences are such calls takes the new class `unanalysed`,
  cause `undetermined` (`outside_truth_intervals`)
  ([measures](docs/validation/evaluation-harness.md#measures)).
- `local_shadow` minor kind and minor rule sets in the draft triage
  (evaluation harness only, #70): a read's neighbour-compatible minor channel
  counts as a shadow only where at least 0.05 of the read's calls within 10
  show the same shadow, and a policy's `minor` is a set of rules any of which
  marks (`upstream@100+local_shadow@80`); the pre-registration parses
  `+`-joined rules.
- The `haplotype_lod` measure (evaluation harness only, #71): the haplotype
  measure, except that a truth two-base IUPAC substitution the reads do not
  observe at the 10 % limit of detection compares IUPAC-compatibly.
  `scores.tsv` adds its status, code, and blocks and `below_lod`; it is the
  approval safety endpoint
  ([measures](docs/validation/evaluation-harness.md#measures)).
- Readiness panels (evaluation harness only, #72): `dna-eval run` reads a
  `readiness` panel as well as `dev`, and `power` reads a development or
  readiness run; any other (hold-out) panel is still refused.

- Discrepancy causes in the evaluation harness: `dna-eval run` gives each
  sample with a discrepancy a cause, `dna`, `lab`, `technology`, or
  `undetermined`, a scope (`defect`, or `deferred` for IUPAC codes DNA does
  not write by design), and its findings (`cause`, `cause_scope`,
  `cause_evidence` in `scores.tsv`), read from the reads' evidence at each position, the
  limit-of-detection signals, and how the reviewers write each window
  sequence over the whole truth. The summary counts causes per class and
  samples per cause and evidence
  ([discrepancy causes](docs/validation/evaluation-harness.md#discrepancy-causes)).
  A sample whose every trace DNA rejects now fails as `rejected`. A
  coverage gap read cleanly yet untrusted is placed by where the call
  stands, with the margins and thresholds read from DNA's configuration:
  beside phase loss or a read end DNA attributes to a repeat run
  (technology), at the read's own end or with a weak support veto (lab), a
  weak veto on a call that reads the reference (DNA, `veto_out_of_scope`),
  mixed channels (undetermined), or otherwise DNA's callability policy; a
  truth IUPAC code a read shows is judged as under heteroplasmy. Each run
  writes the written-form index as `forms.tsv`, keyed by the truth
  database's SHA-256, and `--forms` reads an earlier run's instead of the
  whole database.
- `dna-eval convention` (evaluation harness only): runs every IUPAC-free
  truth sample's calls through DNA's normalization and nomenclature windows
  and counts the samples that come back unchanged, with the change patterns
  of the rest ([notation convention](docs/validation/evaluation-harness.md#notation-convention)).
- Discrepancy classes in the evaluation harness: `dna-eval run` gives every
  sample one class from its statuses (`failed`, `sequence`, `heteroplasmy`,
  `representation`, `coverage`, `not_comparable`, or `none`), naming the part of DNA that owns
  its difference from the truth, with the windows involved (`class`,
  `class_regions` in `scores.tsv`). The summary counts the classes per window
  and, against a baseline, the samples whose class changed
  ([discrepancy classes](docs/validation/evaluation-harness.md#discrepancy-classes)).
- One-strand exception policies in the draft triage (evaluation harness
  only): `region.one_strand` can apply everywhere, nowhere, or only inside
  the repeat windows and near reference homopolymers (R6m, R6m-h4), and a
  single-strand minor-signal rule applies where one strand may pass. The
  approval pre-registration accepts them, and the marginal endpoint reports
  apart the leaks inside the exception.
- `dna-eval approve` and `dna-eval power` (evaluation harness only). `approve`
  runs the pre-registered approval of triage pass rules on a hold-out panel
  once: it checks every frozen input before DNA runs, decides the fixed
  sequence of hypotheses, and bounds each gate's leak rate three ways
  (count-based Clopper–Pearson, Korn–Graubard, batch-cluster bootstrap).
  `power` projects the passes and attainable bounds from a development run
  without the hold-out truth. `dna-eval run` now refuses a hold-out panel
  ([approval](docs/validation/evaluation-harness.md#approval)).
- Haplotype measures in the evaluation harness (PROP-0005).
  - **Harness:** `dna-eval run` also scores every sample under `haplotype`
    and `haplotype_iupac_compatible`: the primary and IUPAC-compatible rules
    applied only to alleles that do not describe the same sequence, so a
    haplotype written or placed differently (`16190.1C` against `16193.1C`,
    `248DEL` against `249DEL`) is not a difference. Scores list the blocks
    where the sequences differ; estimates, anchors, the baseline diff, the
    summary, and the draft triage report both measures
    ([measures](docs/validation/evaluation-harness.md#measures)). The primary,
    strict, and IUPAC-compatible measures are unchanged. The harness now
    loads the reference it passes to DNA, and a sample DNA analyses on another
    reference sequence fails (`reference`).
- Draft triage measurement for PROP-0004 phase 0 (evaluation harness only).
  - **Harness:** `dna-eval run --triage` takes a compact triage evidence
    record per sample from the same analysis, runs draft stage-1 detectors
    (`trace.too_short`, `trace.weak`, `trace.unplaced`, `trace.short`,
    `region.unreliable`, `region.one_strand`, `region.gap`), decides pass,
    review, or rerun per region, trace, and sample, and scores A/B/C, leaks
    with Clopper–Pearson bounds, catches per detector, and a policy sweep
    against the truth (`triage*` outputs,
    [draft triage](docs/validation/evaluation-harness.md#draft-triage)).
    Without `--triage` a run is unchanged. The detectors are research drafts,
    not DNA rules; their results are in
    [PROP-0004](docs/proposals/0004-review-triage.md#phase-0-results).
  - **Rust API:** `SampleAnalysis::regions()` returns the reportable regions
    of the target profile the analysis ran with (SRS-API-010).
  - **Harness:** draft `region.minor_signal` detectors for point
    heteroplasmy over each clean read's minor channel at every covered
    position: cross-strand rules with and without a shadow test, a noise
    ratio, single-read and one-strand variants, swept at several thresholds
    and measured against the same policy without them (`triage_minor.tsv`).
    The draft applies `upstream@100`; results are in
    [PROP-0004](docs/proposals/0004-review-triage.md#minor-signal-results).
  - **Rust API:** `SampleAnalysis::observations(first_1based, last_1based)`
    gives every admitted read's observation of each covered position of an
    interval, as the sample evidence records it at a differential locus,
    whether or not a read differs there, without support aggregates
    (SRS-API-013). A zero or reversed interval is a typed
    `SampleError::InvalidInterval`; a span across the circular origin is two
    calls. `Orientation::reference_base` (`dna-core`) is public.
  - **Rust API:** `sample_analysis::analyze_sanger_outcome` runs the same pass
    as `analyze_sanger` but returns `SampleOutcome::Rejected` with every
    rejected read and its cause when none is admitted, instead of the
    `NoAdmittedReads` error, so a caller can name the traces to rerun
    (SRS-API-009).
- Sample analysis in one placement pass (PROP-0005, SRS-API-009/010/011).
  - **Rust API:** `dna::sample_analysis::analyze_sanger(traces, reference,
    config, workers)` places each trace once on up to `workers` threads and returns a typed `SampleAnalysis`: the
    core sample evidence, the consensus, each admitted read's Sanger evidence,
    and each rejected read with its cause. A caller that needs both no longer
    places every trace twice.
  - `dna sample` and `dna consensus` run the same pass; their documents and
    FASTA are unchanged.

- Evaluation harness `dna-eval` (PROP-0005 harness v0, ADR-0075,
  SRS-COMPAT-003).
  - **Crate:** `crates/dna-eval`, a workspace binary outside the `dna` binary
    and the release. It is the only place for truth formats and comparison;
    no crate depends on it.
  - **Run:** it verifies a panel manifest, runs `consensus` then `call` for
    every sample in process on a thread pool, reads the reviewer truth from
    the lab's statistics database (read-only, `rusqlite`), and classifies each
    sample `Y`, `N`, or `N/A` with the lab's comparison rules, plus a strict
    measure that counts DNA's coverage gaps and an IUPAC-compatible measure.
  - **Panel versions:** `dna-eval panel` derives a new panel version against
    the database, keeping every row, retiring samples without truth, and
    recording per-sample and panel truth digests; a run reports the samples
    whose truth changed since the panel was built.
  - **Outputs:** per-sample comparisons, stratified weighted estimates with
    95 % intervals, anchors, the diff against a baseline run, a Markdown
    summary, and a run header
    ([evaluation harness](docs/validation/evaluation-harness.md)).
  - **Rust API:** `SampleAnalysis::call(sample_id)` calls a sample
    analysis's consensus as `dna call` reads the consensus FASTA, completing
    the recommended `consensus` then `call` workflow in process. It returns a
    `ConsensusCall`: the callable segments and structured `NotationCall`s
    (SRS-API-011, SRS-API-012).
- Sample consensus (PROP-0003 increment 1, ADR-0073, SRS-CONS).
  - **Command:** `dna consensus <sample-id> <trace>... --reference
    <reference.fasta>` adjudicates between a sample's placed reads.
  - **Outputs:** `results/<sample-id>.consensus.json` (`dna.consensus/v1`:
    segments, notable sites with their supporting, opposing and
    uninformative reads, and decision counts) and
    `results/<sample-id>.consensus.fasta`. The FASTA is ready for `dna call`.
  - **Rule:** only clean observations decide. Differences are decided as
    stretches of whole reference runs (as a whole, run by run, or position by
    position). One read against one read takes the reference and is marked
    `contested`.
  - **Measured on two disjoint 320-sample held-out subsets** (fair scale):
    exact samples rise from 288 to 290 and from 290 to 292 against the union
    of reads, and precision from 0.9924 to 0.9967 and from 0.9934 to 0.9964.
  - `dna sample` and `dna call` are unchanged.
  - **Rust API:** `dna_core::consensus` and `dna_core::model::consensus` are
    new; `dna_kernel::plugin::Contract` gains `Consensus`.

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

- The evaluation harness's own measures (strict, IUPAC-compatible, and the
  haplotype measures, the approval endpoint `haplotype_lod` included) judge
  insertions and deletions at 309, 455, 463, and 573 like any other: DNA is
  target-agnostic, and only the primary measure, the lab's rules, keeps the
  lab's exception for these positions. The haplotype measures no longer
  build each interval with the other side's special-position indels.
  Development panel v5: the primary measure is unchanged; 738 samples become
  discordant on `haplotype_lod` (695 at 309, 44 at 573, 1 at 455/463), its
  corpus not-concordant rate goes from 0.0650 to 0.1112, and the `sequence`
  class from 346 to 1,125 samples, mostly the HVS-II 303-309 run DNA writes
  at the reference length without a label.

- The evaluation harness's discrepancy causes read the core's own trust
  verdict on each call (`SampleAnalysis::trust`, `CallTrust`: mask, vetoes,
  `read_end` trust, ADR-0076 exemption, and what withholds trust) instead of
  re-deriving it from four configuration keys. A read's end counts as phase
  loss only where DNA attributes it so (a phase-loss mask label, or a segment
  after a repeat run); development panel v5 moves `beside_phase_loss` from
  141 to 43 findings, mostly to `read_edge` (45 to 112), and statuses and
  classes are unchanged.
- Right-aligned normalization takes an insertion across a base that a
  substitution of the same haplotype makes equal to its first inserted base,
  until the insertion follows a reference run of that base: with `16172C`, an
  inserted C is written `16174.1C` rather than `16171.1C`, as the lab's
  reviewers write it, while `16262.1C 16263C` stays. A deletion that cannot
  shift further moves up to 20 bases 3' across substitutions when the same
  haplotype needs no more edits, every base it crosses is read, and it
  crosses no nomenclature-window boundary (`16293DEL 16298C` becomes
  `16293C 16298DEL`). A candidate shift onto a substituted base is rejected
  rather than failing the read. The `normalization` plugin version rises to 2.
- Window representation now applies per window rather than per read
  (SRS-NOM-007, SRS-NOM-010). An edit that crosses one nomenclature window
  leaves that window in its normalized form, but no longer prevents the
  read's other windows from taking their forensic form: for example a
  deletion of 512-513 beside HVS-III no longer leaves HVS-II's
  `309.1C 315.1C` as `309.1C 309.2T 310C`. **Breaking for library
  consumers:** `dna::variant_nomenclature::apply` no longer fails with
  `NomenclatureError::WindowCrossing`, which is removed. The `nomenclature`
  plugin version rises to 2.
- The nomenclature windows follow the lab's reviewer convention (ADR-0077,
  SRS-NOM-005, SRS-NOM-014, SRS-NOM-015). HVS-II writes any C/T haplotype as
  the left run's length at 309 and the rest left-anchored with C insertions
  after 315 (`312T 315.1C 315.2C 315.3C`); HVS-III now starts at 512 and
  writes whole AC-copy losses and gains at its 3' end (`513A 523DEL 524DEL`,
  `524.1A 524.2C`); HVS-I now starts at 16180 and gains the `anchored_tract`
  rule for the A-tract and C-tract (`16183DEL`, `16188.1C`,
  `16182C 16183C 16189C`, `16189C 16190T 16193.1C 16193.2C`). In `sample`,
  `call`, and `notation`, a window keeps its normalized form when its
  representation would place a call where the read reads nothing, or rests
  on a base the read does not read (for the phylogenetic HVS-I form, any base
  of the window). Every window sequence of at least 5
  reviewed corpus samples now takes the form most reviewers wrote; on
  development panel v5 the primary measure gains 60 samples, strict 56, and
  IUPAC-compatible 80, with no regression on any measure.
- **Breaking for library and log consumers:** the post-calling crate
  `dna-post` is split into `dna-normalize` (representation, normalization,
  nomenclature) and `dna-empop` (conformance) (ADR-0077). Rust paths under
  `dna_post::` move to `dna_normalize::` or `dna_empop::`, and the operation
  log target `dna_post` becomes `dna_normalize` and `dna_empop`. CLI outputs,
  document contracts, the `dna::variant_normalization` and
  `dna::variant_nomenclature` re-exports, and plugin identifiers and versions
  are unchanged.
- Bounded runs before a post-repeat phase loss are exempt from `read_end`
  (ADR-0076, SRS-VAR-013, SRS-CONS-003). An insertion or deletion of a run the
  read bounds with unmasked calls, and the consensus's observations of that
  run's matching and inserted calls and of its reference-matching bounds, are
  no longer distrusted when the margin's edge is the read losing phase right
  after a repeat run. Substitutions, and bounds that read another base, keep
  the margin. `ReadEvidence` masks carry a `MaskOnset`; the `core` plugin
  version rises to 4, the `sanger` plugin version to 3, and the `consensus`
  plugin version to 3.
- Faster per-sample analysis with byte-identical results (SRS-ALN-014).
  - **Seeded pruning** now covers reads with masked or mixed calls (wildcard
    loci) and windows across the circular origin, and an orientation
    evaluated exactly sets the threshold that rejects the other. On 1,197
    real reads, 96 % avoid full-reference Gotoh, up from none. A read whose
    seed probes cover most of the reference goes straight to Gotoh.
  - **Gotoh** fills its matrix in one pass per row with one traceback write
    per cell; **signal statistics** take medians by selection.
  - `dna consensus` drops from about 4 s to 0.13 s per four-trace sample.
- Faster sample commands with byte-identical results (SRS-SAMPLE-001,
  SRS-OUT-008).
  - **Signal statistics** of each sample span are estimated once per read and
    shared by the windowed SNR features and the locus contexts: signal
    processing takes about half the time.
  - **Reads run concurrently** in `dna sample` and `dna consensus`, up to the
    available parallelism. The first failing read in trace order still decides
    the failure, and the operation log still lists each read's records in
    trace order, through that read.
  - On 2,618 panel samples every `consensus` JSON and FASTA and every `sample`
    JSON is byte-identical to the previous release. On a shared server, mean
    `dna consensus` wall time per sample drops from about 134 ms to 72 ms, and
    from 123 ms to 102 ms when pinned to one CPU.
- `consensus` method version 2 decides a stretch by its run structure and
  labels how each run length is known (ADR-0074, SRS-CONS-003/004/008).
  - **Stretches.** A stretch also takes in a neighbouring run that a
    substituted base or a deletion joins, as `16189C` does with the HVS-I
    C-tract.
  - **Composition, then lengths.** The ordered run bases are decided first.
    Run lengths follow, voted as whole sequences when reads show every run's
    length. A length comes only from reads that show both ends of the run
    (ADR-0071); an insertion or deletion of a run the read does not bound is
    never clean.
  - **A run no read shows both ends of** takes one of two lengths:
    - `reference_frame`: the reference positions between the ends that reads
      locate, when every base between is read cleanly;
    - otherwise `phase_loss`: the in-phase length before a read loses phase.
  - **Measured** on the fair scale against method version 1:
    - exact samples rise from 290 to 297 and from 292 to 299 on the two
      320-sample subsets;
    - recall rises from 0.9866 to 0.9899 and from 0.9901 to 0.9928;
    - precision moves from 0.9967 to 0.9970 and from 0.9964 to 0.9967;
    - reviewer-consensus calls identical on 145 of 160 samples (was 142);
    - the local Sequencher set has 28 exact samples (was 26).
    - The false `16182.1C` that a phase-lost read added in HVS-I is gone, and
      HVS-II dominant lengths (`309.1C`, `309.2C`) are kept.
  - **Breaking:** `dna.consensus/v1` is revised in place. Each site gains
    `runs` (`base`, `length`, `length_evidence`), and the summary gains
    `phase_loss_runs` and `reference_frame_runs`.
  - `dna sample` and `dna call` are unchanged apart from the `core` version.

- The `core` plugin method version rises to 2. It records the output changes
  of `run_boundary` (ADR-0071) and run-structure re-expression (ADR-0072),
  which earlier entries in this cycle made without raising it
  ([versioning](docs/governance/versioning.md#plugins)).

- A read the core cannot place no longer fails `dna sample` (SRS-SAMPLE-029,
  PROP-0003 increment 0).
  - The covered cases are too few callable columns, identity below
    `alignment.minimum_identity`, more than one best placement, and tied
    orientations.
  - Such a read is recorded in `rejected_reads[]` with its own reason and
    measured values, and the other reads are still aggregated.
  - On a fresh 320-sample held-out subset, this lets 5 samples that failed
    before produce evidence.
  - **Breaking:** the `rejection` object of `dna.sample_evidence/v10` gains
    four reasons (revised in place), and the error for a sample whose every
    read is rejected now reads "all N reads were rejected and none could be
    analyzed".

- The core re-expresses masked repeat stretches of the selected alignment
  (ADR-0072, SRS-ALN-016).
  - A stretch qualifies when its read calls spell the reference's runs in
    order, with every changed run longer or every changed run shorter, and
    it holds a masked call.
  - It is rewritten as one length edit at the 3' end of each changed run.
  - This recovers HVS-II length insertions that the optimal path split into
    a substitution next to a masked insertion.
  - On the 320-sample held-out set, 6 true `309.1C`/`309.2C` calls are
    gained and one false `310C` is removed; precision rises from 0.9894 to
    0.9897 and recall from 0.9800 to 0.9816.
  - Consensus calls from `dna call` are unchanged.
  - Alignment summaries (`identity`, `gap_opens`) and observed evidence
    change for reads with a rewritten stretch.

- A variant that changes the length of a run of one base is ineligible with
  the new core reason `run_boundary` when the read does not show where that
  run ends on the side the edit touches (ADR-0071, SRS-VAR-016).
  - This targets Sanger reads that dephase right after a poly-C run and
    miscount it by one base. A dephased call that still reads the next base
    resolves the run end.
  - On the 320-sample held-out set, 9 false length calls are removed and
    none is gained; precision rises from 0.9870 to 0.9894, and recall stays
    at 0.9800. Consensus calls from `dna call` are unchanged.
  - **Breaking:** the exclusion-reason vocabulary of
    `dna.sample_evidence/v10` and `dna.variants/v1` gains `run_boundary`
    (revised in place).

- Plugin-first architecture, phase 5b (ADR-0070, PROP-0002): DNA is a Cargo
  workspace.
  - `dna-kernel` (shared contracts), `dna-core` (the core caller),
    `dna-sanger` (the Sanger modality), and `dna-post` (post-calling plugins)
    live under `crates/`; the `dna` facade at the root composes them.
  - Each plugin crate builds and tests alone.
  - The public Rust API (`dna::error`, `dna::variant`, `dna::profile`,
    `dna::variant_normalization`, `dna::variant_nomenclature`,
    `dna::variant_analysis`, `dna::cli`, `dna::run`) and every result
    document are unchanged.
  - The operation log's module path for events of moved code now names its
    crate.
  - Cargo commands take `--workspace`, and builds take `-p dna`.

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

- Alignment places every gap at its 3'-most sequence-equivalent placement
  whatever the read's profile scores (ADR-0047 §1 amended, SRS-ALN-012): the
  optimal score chooses the explanation, and a move is kept when the
  primary-base score is unchanged. With profile-weighted scores two reads of
  one sequence placed an AC-repeat insertion at `513.1C 513.2A` and at
  `524.1A 524.2C`, and the consensus lost it (118 contested AC sites on
  development panel v5, now 4). The one-base move is shared with variant
  normalization (`dna_kernel::model::placement`). An insertion into a repeat
  the read bounds is judged for support vetoes on its placements along the
  repeat, in variant calling and the consensus alike (SRS-VAR-020), so a
  weak call that happens to be taken as the inserted one no longer decides.
  The `core` plugin version rises to 8 and `consensus` to 4. Development
  panel v5 against main: 232 more samples concordant on haplotype LoD and
  227 on the primary measure; the one LoD regression is a truth C-shadow
  listed for re-check.
- `dna call` and `SampleAnalysis::call` no longer report a position DNA
  excluded a call at as callable (SRS-VAR-019). A record's callable
  reference segments leave out every position of an excluded variant (a
  substitution's position, a deletion's deleted positions, both positions
  beside an insertion's junction), so an excluded call is undecided rather
  than read as the reference, and a record excluded for `placement_conflict`
  contributes no callable segments at all; per-read sample evidence is
  unchanged. The `core` plugin version rises to 7. Development panel v5:
  102 more samples concordant on haplotype LoD and 97 on the primary
  measure, no regression; 287 samples lose 441 callable positions inside the
  HV regions (at most 7).
- `call` writes an AC-repeat deletion that a consensus decided at `523–524`
  there, rather than as `513DEL 514DEL` beside the unresolved bases the
  consensus left at `513–516`. On a reviewed consensus record, whose ends are
  vouched, a deletion beside an `N` without evidence now moves 3' across it
  and the repeat beyond (SRS-ALN-017): the `N` scores the same against any
  base, so the traceback's placement there was arbitrary and could delete a
  base the consensus kept. A variant with a mapped call read as `N` is no
  longer reported, with the new core reason `unresolved_flank`
  (SRS-VAR-018); the exclusion-reason vocabulary of
  `dna.sample_evidence/v10` and `dna.variants/v1` gains it (revised in
  place), and the `core` plugin version rises to 6. The same reason also
  removes three deletions beside an `N` that no truth has (`16197DEL`,
  `533–535DEL`). Development panel v5: 29 more samples concordant on the
  primary measure and 31 on haplotype IUPAC-compatible, no regression.
- `call` no longer reports phantom variants from a consensus record placed
  away from the interval its name declares (SRS-VAR-017). A mostly
  undecided consensus segment named for HVS-I's 5' flank could align onto the
  real HVS-I segment and call its few resolved bases there; on the
  development panel (v5) 24 samples carried 115 such calls, which triage
  could not see. Every variant of such a record is now excluded with the core
  reason `placement_conflict`; the `core` plugin version rises to 5.
- **Breaking:** callable reference segments (`callable_reference_segments`)
  no longer include positions read by an unresolved call (`N`, or an IUPAC
  code a modality presents as `N`), masked or not (SRS-ALN-015). A reviewed
  or sample consensus that writes `N` where no read observes a position
  cleanly therefore no longer reports those positions as callable, so a
  comparison counts truth alleles there as uncovered instead of missed: a
  reclassification of DNA's abstentions, not a change of any call. A read
  with an unresolved call in a variant's evidence span no longer counts as
  opposing it.
  `dna.analysis/v9`, `dna.sample_evidence/v10`, and `dna.variants/v1` are
  edited in place (unreleased cycle); the `core` plugin version rises to 3
  ([versioning](docs/governance/versioning.md)). Calls are unchanged.
- A trace with fewer calls than one signal window no longer fails a whole
  `dna sample` or `dna consensus` (SRS-SIG-007, SRS-SAMPLE-030).
  - Read processing checks the call count before any stage runs. Such a read
    is recorded in `rejected_reads[]` with the new reason `too_short`
    (`calls`, `minimum_calls`), and the other reads continue.
  - `analyze` and `basecall` still fail typed with "N calls are fewer than
    window_size_bases W". A one-locus trace now fails the same way, instead
    of in locus-window geometry during basecalling.
  - On the 2620-sample dev panel, the 2 samples that failed now produce
    results (3 reads rejected as `too_short`); every other consensus
    document and FASTA differs only in the `sanger` plugin version.
  - The `sanger` plugin method version rises to 2, because a short read
    that produced an error now produces a rejection
    ([versioning](docs/governance/versioning.md#plugins)). Every document
    that lists plugins in its provenance records it.
  - **Breaking:** `dna.sample_evidence/v10` and `dna.consensus/v1` are
    revised in place. The rejection reasons gain `too_short`, and a
    `too_short` record in sample evidence has no `integrity` or
    `callability`.
  - **Rust API:** `dna_sanger::model::signal::ShortRead` (with
    `ShortRead::of`, which reads the short-read error) is new;
    `SignalError::TooFewLoci` is removed, because the one length check runs
    before signal processing;
    `dna_sanger::model::attachment::SangerRejection` no longer carries the
    trace name and SHA-256. In `dna::sample_analysis` (new in this cycle),
    `RejectedRead` carries `input_name` and `input_sha256` beside its
    `cause`, and `RejectionCause` gains `TooShort(ShortRead)`; its
    `Callability` and `Placement` variants carry the read's `sanger`
    evidence and their `rejection`, which `RejectionCause::sanger` returns.
- A trace whose `PLOC.2` repeats one position (seen at the first or last
  locus of production traces) no longer fails the whole `analyze`,
  `basecall`, or `sample` operation (SRS-IN-003, ADR-0031). Decode merges the
  repeated position into its predecessor and drops the vendor entries at its
  index; a decreasing position is still rejected, now with "PLOC.2 positions
  must not decrease". **Breaking:** trace integrity in `dna.analysis/v9`,
  `dna.basecalls/v3`, and `dna.sample_evidence/v10` gains the required
  `duplicate_ploc_positions` count (revised in place).
- An indel flanking call whose primary base is unresolved (`N`) is omitted from
  public calls (SRS-VAR-006). Previously a tied flank aborted `analyze` or a
  whole `sample` operation, and a mixed-signal flank was published with the
  schema-invalid base `N`.
- PLOC loci one sample apart, which SRS-IN-003 accepts, no longer fail
  basecalling with an empty locus window.

## [0.1.0] - Unreleased

Initial pre-production development baseline.
