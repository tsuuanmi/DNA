# dna-eval

The evaluation harness
([ADR-0075](../../../docs/decisions/adr/0075-evaluation-harness-in-workspace.md)):
the `dna-eval` binary runs DNA over a test panel and compares the calls with
reviewer truth. It is the only crate that reads truth formats or compares
with them (SRS-COMPAT-003). It depends on the `dna` facade, whose
`sample_analysis` capability it runs, on `dna-align`, `dna-sample`, and
`dna-consensus` only to name the evidence types the discrepancy causes read,
on `dna-sanger` for the plugin version that keys the evidence cache, on
the facade's `dna::plugins()` for the version key of the regression ledger,
on
`dna-kernel` for the reference the haplotype measures apply alleles to and
the profile, variant, and nucleotide types it reads, and on `dna-trace` with
its `simulate` feature to render the simulated bench's reads; its tests use
`dna-normalize`'s normalization and nomenclature as the oracle of sequence
equivalence. No crate depends on it, and it is not part of the release.

Output contract: `dna-eval run` writes
[`dna.eval_run/v1`](../../../docs/reference/schemas/eval-run-v1.schema.json):
which files a run directory holds, `run.json`, and every table's columns
([evaluation harness](../../../docs/validation/evaluation-harness.md#outputs),
[PROP-0011 increment 3](../../../docs/proposals/0011-step-contracts.md#increment-3-evaluation-contracts)).
`dna-eval regress check` writes
[`dna.eval_regress/v1`](../../../docs/reference/schemas/eval-regress-v1.schema.json)
([regression set](../../../docs/validation/regression-set.md#outputs)).
The outputs of `diff`, `bench`, and `sim-bench`, and the manifests `panel`
and `subset` write, have no schema; the evaluation harness document
describes them.

- `main.rs` — command line (`run`, `diff`, `panel`, `subset`, `cache`,
  `bench`, `sim-bench`, `regress`) and run orchestration: the options that
  run DNA (`DnaOptions`, shared by `run`, `bench`, and `regress check`) and
  count the review workload (`ReviewOptions`, shared by `run` and `regress
  check`), a panel's inputs loaded and checked before DNA runs (`Loaded`),
  and DNA's outcomes over them with the caches and run contexts it read
  (`Ran`), which `run` and `regress check` share.
- `rundir.rs` — the run directory: its file names, the run header
  `run.json` as one typed record (`RunHeader`) that `run` writes and `diff`
  reads back, the tables `run` writes beside the report's (`checks.tsv`,
  `notes.tsv`, `signals.tsv`, `run_context.tsv`, `components.tsv`,
  `tracts.tsv`), the output directory, and the strict table reader, which
  refuses a table without a column it needs, with a column twice, or with
  a row of another width.
- `bench.rs` — `dna-eval bench`: the artefact regression bench. It parses a
  case list (adjudicated artefacts and truth errors, each with a span and
  the outcome DNA must keep there), judges each case from DNA's alleles,
  check items, and callable intervals for its sample, and reports `BENCH
  FAIL` when a guarded case regresses; known failures are reported, not
  gated. It runs DNA with `run`'s options, the run contexts included. The
  case list names samples, so it stays under `results/`.
- `regress/` — `dna-eval regress`: the versioned regression set, its golden
  output, intents, verdict, and ledger ([regress](regress/README.md)).
- `simbench.rs` — `dna-eval sim-bench`: the simulated tract-decode bench
  ([PROP-0015](../../../docs/proposals/0015-forward-model-tract-decoder.md)).
  It builds each case's templates from the reference, renders them
  with `dna-trace`'s renderer, writes the read as ABIF, runs DNA's read
  evidence stage on it, and scores the tract's decode against the simulated
  truth. It reports the measures per tract and noise level against the
  committed baseline of the set's size: `baselines/sim-bench.tsv` for the
  gating set (60 reads per cell), `baselines/sim-bench-smoke.tsv` for the
  smoke set (6).
- `manifest.rs` — panel manifests: header fields, body hash, sample rows, and
  the optional `retired`, `truth_digest`, and `batch_date` columns; the
  validated `Month` (`YYYY-MM`) a batch date and `--current-since` parse
  into once.
- `panel.rs` — a new panel version from a parent and the truth database.
- `truth.rs` — reviewer truth from the lab's statistics database, opened
  read-only with `rusqlite`; the sample and panel truth digests and the check
  of a panel's recorded truth.
- `runner.rs` — the scoped thread pool that removes repeated traces and runs
  DNA per sample, with typed failures, in sample order, and takes the
  positions where a read shows a second base at the limit of detection. With
  `run --run-signal` it first builds a run signal per batch and per batch
  month from every sample's admitted reads, then runs again with each
  sample's review judged against its group's signal. With `run
  --run-context` it reads each batch's `dna.run_context/v1` from the
  directory (`contexts`) and rebuilds each sample's consensus under its
  batch's context; `main.rs` writes the coverage (`run_context.tsv`).
- `geometry.rs` — the target geometry, read once from the configuration's
  target profile: the reportable regions (named `HV1`, `HV2`, … in the
  profile's order) that bound the comparison, the callable yield, the
  review windows, and the zones; the nomenclature windows the zones are
  placed around; and the windows a discordant position is classified by
  (the nomenclature windows, then the length tracts outside them). Nothing
  else in the harness holds a reference coordinate, but the simulated
  bench's frozen cases.
- `conventions.rs` — the conventions the harness reproduces, with their
  sources (the lab's special positions and region names, the limit of
  detection, the faint second base), and the defaults of a run's options.
- `compare.rs` — concordance of DNA's alleles with the truth under the primary
  (the lab's rules), strict, IUPAC-compatible, haplotype, and haplotype LoD
  measures, and the comparable intervals and edits of the haplotype measures.
- `lod.rs` — the limit of detection of point heteroplasmy: from an analysis,
  where a read shows a second base at 10 % or more; from the truth, the
  two-base IUPAC positions no read shows, which haplotype LoD compares
  compatibly.
- `haplotype.rs` — the haplotype equivalence of two edit sets: both applied
  to the reference over each comparable interval and cut into blocks where
  their lengths agree;
  built once per sample and read under exact or compatible symbols. An open
  block knows whether its interval ends where the truth's does, so an
  insertion past the truth's interval is truth N (`Block::past_interval`).
- `classes.rs` — the class of a status across scoring conventions, and the
  discrepancy class of a sample from its statuses: the part of DNA that owns
  its difference from the truth.
- `cause/` — the cause of a sample's discrepancy (DNA, the lab, or the
  technology) from the reads' evidence and the truth's written forms
  ([cause](cause/README.md)).
- `estimate.rs` — stratified weighted rates with their variance and 95 % intervals.
- `workload.rs` — the review workload: windows over the HV regions, each
  sample's positions to review and unreviewed difference blocks from `dna
  review`'s check items (`SampleAnalysis::review_list`), forced-window
  lists, the Wilson bound of `windows.tsv`, and what passes without review
  (samples, regions, and zones no check item touches, and the passed ones
  with an unreviewed block by the class of their blocks: non-het, truth N,
  het, re-check, or length note) for the summary and `run.json`. The
  counting is `PassSample`'s, which `diff` counts a run's files with too.
- `recheck.rs` — the re-check list of `run --recheck`: the difference blocks
  whose truth is under re-check, per sample and span, that the pass
  measures count in their own class.
- `unrecorded.rs` — the run-length alleles a profile convention leaves
  unrecorded, left out of every comparison on both sides.
- `cache.rs` — the read evidence cache of `run --evidence-cache`: each
  trace's `dna.read_evidence/v1` document, gzip-compressed, in a generation
  keyed by the configuration's SHA-256, the `sanger` version, and the read
  stage's code identity (git object ids recorded by `build.rs`), written
  atomically, a corrupt entry being a miss; `cache prune` keeps the newest
  generations.
- `components.rs` — `components.tsv`: every phase-mixture component
  (ADR-0082) with the reference run clusters it covers, measured as
  haplotypes against the truth, and an unexplained mixture's offsets and
  fractions.
- `tracts.rs` — `tracts.tsv`: every read's tract decodes (SRS-CALL-025,
  PROP-0012 increment 2), with the read's reason for masking the onset call,
  whether the read stays out of phase to its end, the reference position of
  each tract's first call, its chosen structure, and its first-run and
  total marginals; not compared with the truth.
- `zones.rs` — the review zones of the target profile's regions (repeat
  window, after a window, amplicon start and end, interior), the zones a
  check item touches, and the zone instance that holds a position.
- `report.rs` — scores, estimates, anchors, the unreviewed blocks
  (`unreviewed.tsv`), and the Markdown summary; a row's unreviewed blocks
  and cause cells as these files write them, which the regression set's
  golden output reads too.
- `subset.rs` — `dna-eval subset`: a deterministic stratified quick subset
  of a panel as a manifest of its own, every stratum kept, for a first-pass
  `diff`; a pre-filter, never the gate.
- `rundiff.rs` — `dna-eval diff`: a run against a base run of the same
  panel, from their directories alone, read strictly. It refuses runs that
  measure different things, a sample scored twice, or files that count
  other passes than their `run.json` records over all samples or the
  current data; it counts the passes per group (all, period, source) and
  zone with the review workload's counting, and gives the verdict
  that gates a probe on non-het blocks in the current period alone (all
  samples when the runs record no current-data month); het blocks, the
  earlier period, and the sources are reported only.
- `error.rs` — failures that stop a run.
- `contract.rs` — test support: the `dna.eval_run/v1` schema, read so the
  tests check each table's header and cells, and `run.json`'s fields, as the
  harness writes them.

`build.rs` records the git revision the harness is built from, and whether
tracked files differed from it, for `run.json`; the read and placement
stages' code identities, for the caches; and the result contracts the
published schemas declare, for the regression ledger's version key.

Formats, comparison rules, and outputs are specified in the
[evaluation harness](../../../docs/validation/evaluation-harness.md) document.
Tests use synthetic manifests, traces, and SQLite databases built from SQL in
the test.
