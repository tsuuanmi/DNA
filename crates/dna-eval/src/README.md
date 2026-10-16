# dna-eval

The evaluation harness
([ADR-0075](../../../docs/decisions/adr/0075-evaluation-harness-in-workspace.md)):
the `dna-eval` binary runs DNA over a test panel and compares the calls with
reviewer truth. It is the only crate that reads truth formats or compares
with them (SRS-COMPAT-003). It depends on the `dna` facade, whose
`sample_analysis` capability it runs, on `dna-core` and `dna-sanger` only
to name the evidence types the draft triage reads, and on `dna-kernel` to load
the reference the haplotype measures apply alleles to; its tests use
`dna-post`'s normalization and nomenclature as the oracle of sequence
equivalence. No crate depends on it, and it is not part of the release.

- `main.rs` — command line (`run`, `panel`, `approve`, `power`), run orchestration, output
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
  call then fails, or from the rejected reads when none is admitted.
- `compare.rs` — concordance of DNA's alleles with the truth under the primary
  (the lab's rules), strict, IUPAC-compatible, haplotype, and haplotype LoD
  measures, and the comparable intervals and edits of the haplotype measures.
- `lod.rs` — the limit of detection of point heteroplasmy: from an analysis,
  where a read shows a second base at 10 % or more; from the truth, the
  two-base IUPAC positions no read shows, which haplotype LoD compares
  compatibly.
- `haplotype.rs` — the haplotype equivalence of two edit sets: both applied
  to the reference over each comparable interval, once per way of taking
  the special-position indels, and cut into blocks where their lengths agree;
  built once per sample and read under exact or compatible symbols.
- `estimate.rs` — stratified weighted rates with their variance and 95 % intervals.
- `baseline.rs` — the statuses of an earlier run and their classes.
- `report.rs` — scores, estimates, anchors, the baseline diff, and the
  Markdown summary.
- `evidence.rs` — with `run --triage`, the compact triage evidence the
  runner takes from each sample's analysis: per region position, the clean
  strands, masked reads, and consensus marks; per trace, its primer, state,
  and callable span. It names `dna_core` and `dna_sanger` evidence types
  (`Orientation`, `SiteState`, `RunEvidence`, `ReadCallability`) that the
  facade exposes only as fields. `evidence/minor.rs` adds each clean read's
  minor channel per region position, read through `SampleAnalysis::observations`
  ([minor evidence](evidence/README.md)).
- `triage/` — the draft stage-1 triage of PROP-0004 phase 0 over that
  evidence ([triage](triage/README.md)).
- `approve/` — the pre-registered approval run on a hold-out panel and its
  power projection ([approve](approve/README.md)).
- `error.rs` — failures that stop a run.

`build.rs` records the git revision the harness is built from, and whether
tracked files differed from it, for `run.json`.

Formats, comparison rules, and outputs are specified in the
[evaluation harness](../../../docs/validation/evaluation-harness.md) document.
Tests use synthetic manifests, traces, and SQLite databases built from SQL in
the test.
