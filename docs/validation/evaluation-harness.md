---
implements: [SRS-ALN, SRS-API, SRS-CALL, SRS-COMPAT, SRS-CONS, SRS-NOM, SRS-SAMPLE]
decided-by: [ADR-0065, ADR-0069, ADR-0075, ADR-0076, ADR-0077, ADR-0082, ADR-0087]
---

# Evaluation harness

`dna-eval` runs DNA over a fixed test panel and compares the calls with
reviewer truth ([[PROP-0005]],
[[ADR-0075]]). It is
the only place in the repository for truth formats and comparison
(SRS-COMPAT-003). Its outputs are validation evidence about a DNA revision,
not DNA result contracts; their shape is the harness's own versioned contract
([outputs](#outputs)). Source: [crates/dna-eval-lib](../../crates/dna-eval-lib/src/README.md),
under the thin [crates/dna-eval](../../crates/dna-eval/src/README.md) binary.

Measures, discrepancy classes, and causes are in [[validation/evaluation-measures]]; the review workload, its passes, estimates, and anchors are in [[validation/evaluation-workload]].

## Commands

```bash
cargo run --release -p dna-eval -- run \
  --manifest <panel.tsv> --truth <statistics.sqlite> \
  --config config/dna.toml --reference references/rCRS.fasta \
  --output <new-run-directory> \
  --threads <1-256> \
  [--review-window <bp>] [--forced-windows <file>] [--recheck <file.tsv>] \
  [--current-since <YYYY-MM>] \
  [--run-signal batch|month [--run-signal-batch-reads <n>]] \
  [--zone-after <bp>] [--zone-edge <bp>] [--zone-overrun <bp>] \
  [--evidence-cache <dir>] [--run-context <dir>]

cargo run --release -p dna-eval -- cache list --dir <dir>
cargo run --release -p dna-eval -- cache prune --dir <dir> --keep <n> [--dry-run]

cargo run --release -p dna-eval -- panel \
  --manifest <parent-panel.tsv> --truth <statistics.sqlite> \
  --version <v2> --created <YYYY-MM-DD> --output <new-panel.tsv>

cargo run --release -p dna-eval -- subset \
  --manifest <panel.tsv> --rows <n> --minimum <n> --salt <text> \
  --output <new-subset.tsv>
```

Every path is an argument. `run` needs a new or empty output directory, and
`panel` a new output file; nothing is overwritten. `panel` prints the new
body's checksum line (`<sha256>  <file> (body)`). The harness is built from
the workspace but not shipped in a release.

`run` and `bench` refuse a manifest whose header names no panel, or a panel
other than `dev` or `readiness`. A hold-out panel stays reserved: the approval path
that will read it is designed later, with endpoints per zone kind among
auto-passed zones on recent batches (evaluation validity review,
recommendation 5; [[engineering/improvement-method#6. Panels and their uses]]).

## Comparing two runs

```bash
cargo run --release -p dna-eval -- diff --base <run-directory> --run <run-directory> \
  [--output <new-directory>] [--force]
```

`diff` compares a candidate run with a base run of the same panel from
their directories alone, and gates it.

- **Comparability.** Both runs must record the same manifest body, truth
  digest, re-check list, current-data month, review window, zone widths,
  reference and profile (their SHA-256 in `identities`), run signal grouping
  and batch reads, and run contexts (their digest), each from a clean build;
  the configuration may differ, since a candidate changes it. The manifest
  the base run names must still hash to the body it recorded, since the
  periods and sources are read from it. Otherwise `diff` refuses, or with
  `--force` it
  compares and warns at the top of the report.
- **Reading.** It reads each run's `run.json`, `scores.tsv`, and
  `unreviewed.tsv` strictly: a missing or unparsable header field, a table
  without a column `diff` reads or with a column twice, a row of another
  width, a sample scored twice, an unknown block class, or a block whose
  `zone_checked` disagrees with its sample's `review_zones` stops the
  comparison.
- **Counting.** From each run's `scores.tsv` and `unreviewed.tsv`, it counts
  the passes with the review workload's own counting, as the run's summary
  does: a passed unit counts under its most severe unreviewed block. It
  first checks that these counts equal the ones the run's `run.json`
  records, the counted and auto-passed samples by class and every zone's
  units, over all samples and over the current data, and refuses a run whose
  files disagree.
- **Groups.** Every count is given for all samples, by period (the
  manifest's `batch_date` against the current-data month), and by the
  manifest's `source` column when it has one.

The report (`diff.md`) shows, before and after:

- **Samples:** samples, clear samples, samples with a non-het unreviewed
  block and those whose only unreviewed blocks are het, passed samples by
  their most severe block (length notes reported, never gated), and
  positions to review.
- **Zones:** for every zone that changed, the units needing a check and the
  passed units with a non-het block.
- **Coverage:** the truth alleles DNA's callable intervals newly leave out
  or now cover, by allele.
- **Status changes:** per measure, by status.
- **Changed samples:** the samples whose unreviewed blocks or clear state
  changed, with the blocks gone and new. The first 40 are shown;
  `samples.tsv` lists all, and `statuses.tsv` lists every status change.

The verdict reads one group, the current period: the samples whose
`batch_date` is in or after the runs' current-data month. The LAB keeps
improving, so recent data decides (user decision 2026-10-21). The earlier
period, all samples, and the sources are reported in full and never gate; a
source mixes periods, and its current samples already count in the current
period. A sample without a `batch_date` is in no period and does not gate.
When the runs record no current-data month, the verdict reads all samples.
The report names the group under **Verdict**.

The verdict is `GATE FAIL`, with exit status 2, when in that group more
passed samples carry a non-het block, more samples carry a non-het
unreviewed block, or more passed units of any zone carry a non-het block. It
also fails when the samples that failed to score in that group differ
between the runs, whatever the cause (the change, or a transient read
error): the group then holds different samples. Failures that differ only
outside it are a warning. Heteroplasmy is optional (user decision): het
blocks are reported, beside re-check and truth-`N` blocks, and never gate.
Otherwise it is `GATE PASS`.

### Same outputs

```bash
cargo run --release -p dna-eval -- same <run-directory> <run-directory>
```

`same` says whether two runs' outputs are byte-identical: every file but
`run.json`, compared by SHA-256. It is the check of a refactor, a cache, or
a build change that must not change what a run writes. It lists the files
that differ, with the first line at which each does, and the files only one
run holds. Each run's `run.json` records the digests of its other files
(`outputs`); they are checked against the files first, so a file edited
after its run is refused rather than compared. A run of an earlier revision
records none, and its files are hashed as found. The exit status is 0 when
every output is identical, 2 with `DIFFERENT` when one is not, and 1 when a
run cannot be read.

### Quick subset (pre-filter)

`subset` writes a deterministic stratified sample of a panel's staying rows
as a manifest of its own: each stratum keeps its share of `--rows`, at least
`--minimum` rows (all of them when it has fewer, so the rare strata are
always present), chosen by the smallest SHA-256 of `--salt`, a tab, and the
sample identifier, and written verbatim in the parent's order. Its header
keeps the parent's panel fields, records `subset_rows`, `subset_minimum`,
`subset_salt`, and the parent's body hash, and recomputes the row count, the
body hash, and the truth digest, so a run records them in `run.json`. Run
the base and the candidate on the same subset and compare them with `diff`.
**A subset verdict is a pre-filter, never the gate: a candidate that passes
it must still pass `diff` on the full development panel before a pull
request.**

`bench` checks adjudicated artefact cases before any panel run:
[[validation/artefact-bench]]. `sim-bench` measures the tract decoder on
simulated length mixtures against the committed baseline:
[[validation/simulated-bench]].

## Inputs

**Panel manifest.** `#` comment lines with `key=value` fields come first, then
one tab-separated column header, then one row per sample.

- `body_sha256` is required: the SHA-256 of every byte after the column header
  line. A mismatch stops the run. `rows`, when present, must equal the row
  count. `panel`, `version`, `salt`, and from version 3 `truth_digest` (the
  panel truth digest) are recorded.
- The columns `sample_id`, `batch`, `stratum`, `weight`, `weight_corpus`,
  `anchor_candidate` (`0` or `1`), and `trace_paths` (`;`-separated) are read.
  The optional `retired` column (`0` or `1`, from panel version 2) marks
  samples kept for continuity but left out of the estimates. The optional
  `truth_digest` column (from version 3, 64 lowercase hexadecimal digits)
  records each sample's truth digest when the panel was built. The optional
  `batch_date` column gives the batch's run date as `YYYY-MM-DD` (a day of
  01 to 31), or empty when unknown; another value stops the run. Its month
  is the sample's batch month everywhere the harness reads one: the
  current data, the run signal's month groups, `run_context.tsv`, and the
  summary's months. Other columns are ignored. Sample identifiers are
  unique.
- Bounds: 256 MiB, 1,000,000 rows, and 64 traces per sample.

**Truth.** The lab's statistics database (`statistics.sqlite`, at most
4 GiB), opened read-only. It is the indexed copy of the lab's merged reviewer
results (`statistics.json`) and the single truth source of the harness. Its
`info` table must declare `schema_version` 1. Per panel sample the harness
reads the batch (`samples`), the reviewed alleles (`variants`: `pos`, `ref`,
`seq`, with `-` for an inserted or deleted base; at most 10,000 per sample),
the analysed intervals (`intervals`), and the sample's status in its batch
(`sample_status`). Sample flags and counts are not used.

## Panel versions

`panel` derives a new version from a parent manifest and the database:

- membership is frozen: every parent row is kept, in order;
- a sample the database has no record of is retired (`retired = 1`, both
  weights `0`), and the weights of the other samples of its stratum are scaled
  in proportion, `w' = w · total / Σ staying w` (written with four decimals),
  so each stratum still sums to its population and keeps its differences
  within; a zero weight stays zero. When every staying weight of a stratum is
  zero, its weights stay zero and `panel` reports the stratum as having lost
  weight. Strata without a newly retired sample keep their weights as written;
- each row records its sample's truth digest (`truth_digest`), and the header
  the panel truth digest;
- the header records the panel, the new version, the salt, the creation date,
  the row and retired counts, and the parent version and body hash; the
  database path, file SHA-256 (`truth_database_sha256`), `created_at`, and
  sample count are recorded as provenance only. Header values escape `%`, `=`,
  and white space, so any path reads back.

Earlier versions stay runnable; report results on the rows both versions
share.

## Panel truth digest

The lab rebuilds the database in place, so its file hash changes with every
refresh even when no panel sample changes. A panel therefore records the
truth it was built against as digests of exactly what the harness uses:

- a sample's digest is the SHA-256 of a canonical text: `batch\t<batch>`,
  `status\t<status>` (empty when unrecorded), one `allele\t<pos>\t<ref>\t<seq>`
  line per allele with the position normalised, and one
  `interval\t<start>\t<end>` line per analysed interval, the alleles and
  intervals sorted, each line ending in a newline; a sample without a record
  digests `missing\n`;
- the panel's digest is the SHA-256 of one `<sample_id>\t<digest>` line per
  sample, in panel order.

`run` computes the digests of the current database. When the panel digest
equals the manifest's, the truth is `unchanged` and the database's file hash
and `created_at` are recorded in `run.json` for provenance only. Otherwise the
run is marked `truth changed`: the summary starts with a notice, the samples
whose digest differs, or that record none under a header that does, are
listed in `truth_changes.tsv` (written whenever the truth changed, with only
its header when no row differs; the notice then says so), and a new panel
version is due. A manifest without a digest (versions 1 and 2) reports
`not_recorded`.

## Running DNA

For each sample the harness removes repeated trace content (keeping the first
path of each content), then runs `dna::sample_analysis::analyze_sanger` with
one worker and `SampleAnalysis::call`: the `consensus` then `call` workflow in
process, with the configuration and reference given. Its notation calls become alleles (`309.1C` →
`309.1`, `-`, `C`; `249DEL` → `249`, `A`, `-`), and the consensus reads'
callable reference segments become DNA's intervals. Samples run on a pool of
scoped threads (`--threads`). A sample whose traces cannot be read, that DNA
rejects, or on which DNA panics gets a failure row; the others continue. The
truth of the whole panel and the reference load before DNA runs, so a
database that cannot be opened, has another schema, or fails a query, or a
reference that cannot be read, stops the run first. The haplotype measures
compare on that reference sequence; a sample that DNA analyses on another
sequence gets a failure row (kind `reference`).

`--evidence-cache <dir>` analyses every read from its `dna.read_evidence/v1`
document instead of its trace, with
`dna::recorded::analyze_evidence_outcome`, as `dna sample --evidence`
does ([[PROP-0011]]). Each document is read from the cache or made with
`read_evidence_document` and stored. A hit and a miss therefore analyse
alike, and the outputs equal a run from the traces: on 410 samples of panel
v5, every table was byte-identical from the traces, a cold cache, and a warm
one.

- **Key.** The read stage reads no profile and no reference ([[ADR-0069]]),
  and nothing of the core, so a document depends on the trace, the
  configuration, and the read-stage code alone. An entry is the trace's
  content SHA-256 in a generation named for the configuration file's SHA-256
  (which the document's provenance records and the reader checks), the
  `sanger` version, and the read-stage code identity.
- **Code identity.** The identity is the git object ids of every source from
  a trace to its document, which `crates/dna-eval/stages.txt` lists
  (`[read]`): the build's inputs (`Cargo.lock`, `Cargo.toml`,
  `rust-toolchain.toml`), `dna-kernel`, `dna-sanger`, `dna-trace`, and the
  facade's configuration, trace loading, evidence models, plugin
  composition, and writer, with the module files above them, recorded by
  the `dna-eval` binary's build script. A harness built with any of those
  sources differing from its commit has no identity, and the option fails.
  A plugin version does not move on a work-in-progress commit, so it alone
  would serve stale documents. `validate_module_layers.py` checks the list
  is closed: every module a listed source uses, every module file above
  it, and every workspace crate a listed crate depends on, is listed, or
  marked as used but off the stage's path (the configuration envelope
  parses the core's sections, which the read stage never reads). A file on
  the path for some items only is listed with them, and what its other
  items use is not followed.
- **Storage.** Entries are gzip-compressed and renamed into place. An entry
  that does not decompress, is not JSON, or names another trace is a miss,
  remade and replaced.
- **Record.** `run.json` records the directory, the key, the hits, misses,
  and replaced entries, and the generation's size.
- **Hygiene.** Every generation labels itself when a run opens it:
  `key.json` records its kind (`evidence`, `placement`, or `forms`) and key,
  written once, and `last_used` the time of the latest run that opened it,
  rewritten by each run, so recency is use, not the last write. A
  generation written before the labels lists as `unlabelled`, by its newest
  entry, until a run opens it: in a cache from before the labels, a
  generation still in use can look old, so run once on what you keep, and
  read `--dry-run`, before pruning it. `cache list` prints one line per
  generation (TSV: generation, kind, last use in UTC, entries, bytes, key)
  by kind and most recently used first, and counts the temporary files of
  runs no longer alive. `cache prune --keep <n>` keeps the `n` most recently used
  generations of each kind, unlabelled ones being a kind of their own, and
  removes those temporary files; a temporary file is named for its writer's
  process, and one whose process no longer runs on this host is orphaned.
  `--dry-run` reports what it would remove and removes nothing.
- **Forms cache.** The same directory caches the index of the truth's
  written forms ([[validation/evaluation-measures#Discrepancy causes|forms]]),
  which otherwise indexes the whole truth database on every run. A
  generation is named for the truth database file's SHA-256, the reference
  sequence's SHA-256, the target profile's repeat windows (each name and
  bounds, in order), and the index's code identity (`[forms]`: the build's
  inputs, the kernel, and the harness's truth reader, comparison alleles,
  geometry, and index). A run on the same
  inputs reads the index; a renamed or moved window, another reference or
  database, or a change to that code is a new generation, so a stale index
  is never read. An entry that does not decompress or parse is rebuilt and
  replaced. A harness built with those sources differing from its commit
  builds the index and stores nothing; the evidence cache still works.
  `run.json` records the key and whether it was a hit (`forms_cache`). On
  development panel v5 the tables are byte-identical with a cold forms
  cache, a warm one, and none.
- **Placement cache** ([[ADR-0087]], proposed). The same directory holds
  each evidence document's `dna.placed_read/v1` document, written by
  `dna-align`. Entries are keyed by the evidence document's SHA-256, in a
  generation named for the configuration, the placement stage's code
  identity, and the reference. The code identity (`[placement]`) is the git
  object ids of the build's inputs, the kernel, `dna-align`, and the
  facade's configuration, evidence reading, and writer. Every read is then placed from its document,
  which `dna-align` re-checks against the read's evidence and the reference.
  A change to the placement stage is a new generation. `run.json` records it
  as `placement_cache`. On development panel v5 the tables are
  byte-identical with a cold placement cache, a warm one, and none. A
  consensus-only probe takes 29-36 % less CPU.

`--run-signal` judges each sample's review against a run signal
([[PROP-0009]], SRS-CALL-020) built in
process, the way `dna review --signal` would with a per-run file. A first
pass analyses every sample and adds its admitted reads' contributions
(`SampleAnalysis::signal_contributions`) to its groups' run signals: its
manifest `batch` and its batch month (`batch_date`, `unknown` without one).
A second pass runs DNA again and reviews each sample against its group's
signal, the sample's own reads taken out. With `batch`, a sample uses its
batch's signal when that batch has at least `--run-signal-batch-reads` reads
(default 40), otherwise its month's; with `month`, always its month's. The
manifest's batch stands for the laboratory's run, which the panel does not
record. Without the flag, no run signal is built and the floor alone sets
every threshold.

`--run-context <dir>` calibrates each sample's consensus by its batch's
[[reference/run-context|run context]] ([[PROP-0014]], SRS-CONS-014):
`<dir>/<batch>.run_context.json`, as `dna run-context` writes it over the
batch's full folder, its `run_id` the manifest batch. A panel holds too few
of a batch's samples to build one in process (about 16, or 2.5 run-1 pairs),
so the directory is built beforehand and serves as the cache. A sample whose
batch has no file is analysed with the configured τ; one whose reads are
only partly in its batch's context fails. `run.json` records the directory,
how many of the panel's batches had a context, and the contexts' digest
(the SHA-256 of one `<batch>\t<file SHA-256>` line per context, in batch
order), and `run_context.tsv` the coverage per batch month. With
`--run-signal batch`, `run.json` also records `run_signal_batch_reads`.

## EMPOP queries

When the profile declares `[empop]`, each sample's EMPOP query
(`SampleAnalysis::call`, SRS-NOM-018) is written to `scores.tsv`. The harness
checks it against the consensus. A range must not cover a position the
consensus left undecided (`empop_undecided`):
- a site without a call counts its positions inside the ranges;
- a call of the reference's length counts its `N`s inside them;
- a call that changes the length counts its `N`s less the site positions the
  ranges leave out.

The export follows the variants document, so a nonzero count means that `call`
placed a consensus `N` elsewhere: a consensus-to-`call` interface loss, not an
export defect. A tie the consensus resolved is exported as DNA's decision;
`empop_ties` counts its positions. The summary gives the samples with a query
and the positions and samples of each count.

## Outputs

`dna-eval run` writes one contract, `dna.eval_run/v1`
([PROP-0011](../proposals/0011-step-contracts.md#increment-3-evaluation-contracts)).
Its schema,
[`eval-run-v1.schema.json`](../reference/schemas/eval-run-v1.schema.json), is
authoritative for which files a run writes, every table's columns in order,
and each column's format and meaning; it also says how a directory reads as
one instance. The synthetic
[example](../reference/examples/eval-run-v1.example.json) is such an
instance. `validate_result_schemas.py --eval-run <directory>` checks a run
directory against it, and the harness's tests check the headers and
enumerated cells it writes. A change to a file's columns, formats, or meaning,
or to which files a run writes, is a new contract version
([[governance/versioning#evaluation-contracts|versioning]]).

| File | Written | Content |
|---|---|---|
| `scores.tsv` | always | per sample: the comparison under each [measure](evaluation-measures.md#measures), the [discrepancy class](evaluation-measures.md#discrepancy-classes) and [cause](evaluation-measures.md#discrepancy-causes), the [EMPOP query](#empop-queries), the callable yield per HV region, the traces, the failure, and the [review workload](evaluation-workload.md#review-workload); read back by `dna-eval diff` |
| `components.tsv` | always | per phase-mixture component (ADR-0082) and reference run cluster it covers: the run it enters in phase against the truth's |
| `tracts.tsv` | always | per tract decode (SRS-CALL-025, [[PROP-0012]] increment 2): the read's reason for masking the onset call and whether it stays out of phase to its end, where the tract's first call lies on the reference, its chosen structure, and its first-run and total dominant lengths, margins, and margin states; evidence only, not compared with the truth |
| `checks.tsv` | always | per sample and `dna review` check item: its positions and labels |
| `notes.tsv` | always | per sample, `dna review` item, and run its `tract_length` label notes: the tract, the run, its status, DNA's length, the reading, and the convention's reading |
| `unreviewed.tsv` | always | per sample and unreviewed `haplotype_lod` difference block: its anchors and text, its class as the pass measures split it (non-het, het, truth `N`, re-check, length note), the zone that holds its first anchor, and whether a check item touches that zone; `dna-eval diff` counts the passes from it |
| `windows.tsv` | always | per review window: the samples it reaches, checks, and unreviewed differences |
| `estimates.tsv` | always | per weighting, measure, and metric: the stratified rates and interval |
| `anchors.tsv` | always | the anchor candidates under each measure |
| `forms.tsv` | always | how the reviewers write each sequence over the whole truth database, keyed by the database file's SHA-256, without sample identifiers; cached with `--evidence-cache` |
| `run.json` | always | the run header: harness and source revision, inputs and identities, options, measurement constants, counts per measure, what passes without review (the [passes](evaluation-workload.md#review-workload) of the review workload), timings, and the SHA-256 of every other file (`outputs`, written last, which [`same`](#same-outputs) compares) |
| `summary.md` | always | the run in aggregate: the panel checked against the truth (truth digest state, changed, retired samples, samples without a record, batch differences, status not `OK`), outcomes and estimates under each measure, primary discordance, discrepancy classes (per class and window, retired samples left out, the causes per class and the samples per cause and evidence), the EMPOP queries, the callable yield per HV region (samples, analysed and callable positions, their share, and the samples short of full yield), coverage gaps, IUPAC alleles, representation (discordant samples a haplotype measure finds concordant, by window, and the blocks of the rest), anchors, and the review workload |
| `truth_changes.tsv` | when the panel's truth changed | the samples whose truth digest differs from the manifest's record |
| `signals.tsv` | with `--run-signal` | per run-signal group: reads, calls, and keys |
| `run_context.tsv` | with `--run-context` | per batch month, then `all`: samples, analysed, with a run context, and in how many the context rescaled run 1's and run 2's τ (the rest fail closed); and the month's batches with a context and their run-1 and run-2 pairs, which show why (a batch below `consensus.length_context_minimum_pairs` calibrates nothing) |

The Markdown reports are for people: the contract requires them, not their text.

Every file except `run.json` depends only on the panel, the truth, DNA's
calls, and the configuration, reference, profile, and re-check list the summary
names, so it is byte-identical whatever the thread count and whatever revision
or database file it ran on; `forms.tsv` alone also records the database
file's SHA-256 in its first line. Run outputs and
panel manifests hold sample identifiers and stay outside the repository under
the [[governance/data|data policy]]; DNA documents quote aggregates only.
