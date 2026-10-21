# dna-eval

The evaluation harness
([ADR-0075](../../../docs/decisions/adr/0075-evaluation-harness-in-workspace.md)):
the `dna-eval` binary runs DNA over a test panel and compares the calls with
reviewer truth. It is the only crate that reads truth formats or compares
with them (SRS-COMPAT-003). It depends on the `dna` facade, whose
`sample_analysis` capability it runs, on `dna-align`, `dna-sample`,
`dna-consensus`, and `dna-sanger` only to name the evidence types the draft triage and the discrepancy causes read, on `dna-kernel` for the
reference the haplotype measures apply alleles to and the profile, variant,
and nucleotide types it reads, and on `dna-normalize`
to run the truth's calls through DNA's notation (`convention`); its tests use
`dna-normalize`'s normalization and nomenclature as the oracle of sequence
equivalence. No crate depends on it, and it is not part of the release.

Output contract: `dna-eval run` writes
[`dna.eval_run/v1`](../../../docs/reference/schemas/eval-run-v1.schema.json):
which files a run directory holds, `run.json`, and every table's columns
([evaluation harness](../../../docs/validation/evaluation-harness.md#outputs),
[PROP-0011 increment 3](../../../docs/proposals/0011-step-contracts.md#increment-3-evaluation-contracts)).
The outputs of `approve` and `convention`, and the manifest `panel` writes,
have no schema; the evaluation harness document describes them.

- `main.rs` — command line (`run`, `panel`, `approve`, `power`, `convention`), run orchestration, output
  directory, and the `run.json` header.
- `manifest.rs` — panel manifests: header fields, body hash, sample rows, and
  the optional `retired` and `truth_digest` columns.
- `panel.rs` — a new panel version from a parent and the truth database.
- `truth.rs` — reviewer truth from the lab's statistics database, opened
  read-only with `rusqlite`; the sample and panel truth digests and the check
  of a panel's recorded truth.
- `runner.rs` — the scoped thread pool that removes repeated traces and runs
  DNA per sample, with typed failures, in sample order, and takes the
  positions where a read shows a second base at the limit of detection; with `run --triage` it
  also takes each sample's triage evidence, from the analysis even when the
  call then fails, or from the rejected reads when none is admitted. With
  `run --run-signal` it first builds a run signal per batch and per batch
  month from every sample's admitted reads, then runs again with each
  sample's review judged against its group's signal.
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
- `classes.rs` — the discrepancy class of a sample from its statuses: the
  part of DNA that owns its difference from the truth.
- `cause/` — the cause of a sample's discrepancy (DNA, the lab, or the
  technology) from the reads' evidence and the truth's written forms
  ([cause](cause/README.md)).
- `estimate.rs` — stratified weighted rates with their variance and 95 % intervals.
- `workload.rs` — the review workload: windows over the HV regions, each
  sample's positions to review, the zones its check items touch and those
  only fast check items touch, and unreviewed difference blocks from `dna
  review`'s check items (`SampleAnalysis::review_checks`), forced-window
  lists, the Wilson bound of `windows.tsv`, and what passes without review
  (samples, regions, and zones no check item touches, and the passed ones
  with an unreviewed block by the class of their blocks: non-het, truth N,
  het, or re-check) for the summary and `run.json`.
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
- `baseline.rs` — the statuses of an earlier run and their classes.
- `report.rs` — scores, estimates, anchors, the baseline diff, the
  unreviewed blocks (`unreviewed.tsv`), and the Markdown summary.
- `subset.rs` — `dna-eval subset`: a deterministic stratified quick subset
  of a panel as a manifest of its own, every stratum kept, for a first-pass
  `diff`; a pre-filter, never the gate.
- `rundiff.rs` — `dna-eval diff`: a run against a base run of the same
  panel, from their directories alone. It refuses runs that measure
  different things or whose files disagree with their `run.json`, counts the
  passes per group (all, period, source) and zone, and gives the verdict
  that gates a probe on non-het blocks; het blocks are reported only.
- `evidence.rs` — with `run --triage`, the compact triage evidence the
  runner takes from each sample's analysis: per region position, the clean
  strands, masked reads, and consensus marks; per trace, its primer, state,
  and callable span. It names `dna_align`, `dna_sample`, `dna_consensus`, and `dna_sanger` evidence types
  (`Orientation`, `SiteState`, `RunEvidence`, `ReadCallability`) that the
  facade exposes only as fields. `evidence/minor.rs` adds each clean read's
  minor channel per region position, read through `SampleAnalysis::observations`
  ([minor evidence](evidence/README.md)).
- `triage/` — the draft stage-1 triage of PROP-0004 phase 0 over that
  evidence ([triage](triage/README.md)).
- `approve/` — the pre-registered approval run on a hold-out panel and its
  power projection ([approve](approve/README.md)).
- `convention.rs` — `convention`: every IUPAC-free truth sample's calls
  through DNA's normalization and nomenclature windows, and the change
  patterns where DNA writes them differently.
- `error.rs` — failures that stop a run.
- `contract.rs` — test support: the `dna.eval_run/v1` schema, read so the
  tests check each table's header and cells, and `run.json`'s fields, as the
  harness writes them.

`build.rs` records the git revision the harness is built from, and whether
tracked files differed from it, for `run.json`.

Formats, comparison rules, and outputs are specified in the
[evaluation harness](../../../docs/validation/evaluation-harness.md) document.
Tests use synthetic manifests, traces, and SQLite databases built from SQL in
the test.
