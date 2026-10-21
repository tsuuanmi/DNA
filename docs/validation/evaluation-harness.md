---
implements: [SRS-ALN, SRS-API, SRS-CALL, SRS-COMPAT, SRS-NOM, SRS-SAMPLE]
decided-by: [ADR-0065, ADR-0069, ADR-0075, ADR-0076, ADR-0077, ADR-0082, ADR-0087]
---

# Evaluation harness

`dna-eval` runs DNA over a fixed test panel and compares the calls with
reviewer truth ([[PROP-0005]],
[[ADR-0075]]). It is
the only place in the repository for truth formats and comparison
(SRS-COMPAT-003). Its outputs are validation evidence about a DNA revision,
not DNA result contracts; their shape is the harness's own versioned contract
([outputs](#outputs)). Source: [crates/dna-eval](../../crates/dna-eval/src/README.md).

Measures, discrepancy classes, and causes are in [[validation/evaluation-measures]]; review workload, estimates, and the draft triage are in [[validation/evaluation-triage]].

## Commands

```bash
cargo run --release -p dna-eval -- run \
  --manifest <panel.tsv> --truth <statistics.sqlite> \
  --config config/dna.toml --reference references/rCRS.fasta \
  --output <new-run-directory> [--baseline <earlier-run-directory>] \
  --threads <1-256> [--triage] [--forms <earlier-run-directory>/forms.tsv] \
  [--review-window <bp>] [--forced-windows <file>] [--recheck <file.tsv>] \
  [--current-since <YYYY-MM>] \
  [--run-signal batch|month [--run-signal-batch-reads <n>]] \
  [--zone-after <bp>] [--zone-edge <bp>] [--zone-overrun <bp>] \
  [--evidence-cache <dir>]

cargo run --release -p dna-eval -- cache prune --dir <dir> --keep <n>

cargo run --release -p dna-eval -- panel \
  --manifest <parent-panel.tsv> --truth <statistics.sqlite> \
  --version <v2> --created <YYYY-MM-DD> --output <new-panel.tsv>

cargo run --release -p dna-eval -- subset \
  --manifest <panel.tsv> --rows <n> --minimum <n> --salt <text> \
  --output <new-subset.tsv>

cargo run --release -p dna-eval -- power \
  --dev-run <development-or-readiness-run-directory> --manifest <hold-out-panel.tsv> \
  [--base <sweep-policy> --candidate <sweep-policy> \
   [--true-rate <rate> --bound <rate>...]]

cargo run --release -p dna-eval -- approve \
  --prereg <committed-prereg.toml> --output results/<new-directory> \
  --threads <1-256>

cargo run --release -p dna-eval -- convention \
  --truth <statistics.sqlite> --profile config/profiles/human-mtdna-rcrs.toml \
  --reference references/rCRS.fasta --output <new-directory>
```

Every path is an argument. `run` needs a new or empty output directory, and
`panel` a new output file; nothing is overwritten. `panel` prints the new
body's checksum line (`<sha256>  <file> (body)`). The harness is built from
the workspace but not shipped in a release.

`run` refuses a manifest whose header names a panel other than `dev` or
`readiness`: a hold-out panel is read only by `approve`, once, under a
committed pre-registration ([approval](#approval)). `power` reads a
development or readiness run.

## Comparing two runs

```bash
cargo run --release -p dna-eval -- diff --base <run-directory> --run <run-directory> \
  [--output <new-directory>] [--force]
```

`diff` compares a candidate run with a base run of the same panel from
their directories alone, and gates it.

- **Comparability.** Both runs must record the same manifest body, truth
  digest, re-check list, current-data month, review window, and zone widths,
  each from a clean build. Otherwise `diff` refuses, or with `--force` it
  compares and warns at the top of the report.
- **Counting.** From each run's `scores.tsv` and `unreviewed.tsv`, it counts
  the passes as the run's own summary does: a passed unit counts under its
  most severe unreviewed block. It first checks that these counts equal the
  ones the run's `run.json` records, and refuses a run whose files disagree.
- **Groups.** Every count is given for all samples, by period (the
  manifest's `batch_date` against the current-data month), and by the
  manifest's `source` column when it has one.

The report (`diff.md`) shows, before and after:

- **Samples:** samples, clear samples, samples with a non-het unreviewed
  block and those whose only unreviewed blocks are het, passed samples by
  their most severe block, and positions to review.
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

## Approval

`dna-eval approve` is the single approval run of the [approval
protocol](../proposals/0005-evaluation-feedback-loop.md). Its input is a
pre-registration (TOML) committed to the repository:

- **`[inputs]`** — the manifest and its body SHA-256, the panel truth digest,
  a pinned read-only copy of the truth database and its SHA-256, the
  configuration and its SHA-256, the reference, and the DNA revision.
- **`[analysis]`** — `bootstrap_draws` and `bootstrap_seed`.
- **`[utility]`** — `form` (`point` or `lower_bound`) and
  `minimum_pass_rate`.
- **`[[gates]]`** — a measure and an `upper_bound` each; every tested
  hypothesis must meet every gate.
- **`[[hypotheses]]`** — the fixed sequence in order: a name, `requires`
  (`first`, then `previous_pass`: the sequence stops at the first failure, so
  the family-wise alpha holds), the full policy
  (`short_calls`, `causes`, `one_strand` as `everywhere`, `off`,
  `outside_repeats`, or `outside_repeats_h<k>`, `proxy`, `minor` as one or
  more `<kind>@<t>` joined by `+` (any of which marks) or `off`, and
  `single_minor` as `<kind>@<t>` or `off`), and an optional `marginal`
  (`against` an earlier hypothesis, a measure, an `upper_bound`, the
  `projected_n` that `dna-eval power` gave at its lower bound, and
  `projected_from`: the development `run`, the SHA-256 of its
  `triage_policies.tsv`, and the `base` and `candidate` sweep policies
  projected).

The run stops before DNA when any value is still `"pending"`, a policy lacks
the `n` cause, a marginal `upper_bound` is below the count-based bound at 0
leaks among its `projected_n` samples (it cannot be met), the file is not
committed unmodified, the harness build is not the frozen revision or has
local changes, or the manifest, configuration,
truth copy, or panel truth digest differs from the frozen value. The output
directory must be new and lie under `results/`.

Per tested hypothesis, a sample passes when DNA called it and the policy
decides `pass`; failed samples never pass, and retired ones are left out.
Each gate bounds the leak rate (passed samples not concordant under the
gate's measure) three ways, at one-sided alpha 0.05:

- the count-based Clopper–Pearson upper bound;
- the Korn–Graubard upper bound of the corpus-weighted rate: the variance is
  the stratified variance of `estimate.rs` over the passed samples, each
  stratum's passed samples its sample, which gives the effective sample size
  (a zero rate uses the passed count); the degrees-of-freedom correction is
  omitted, and the bootstrap is what covers clustering by batch;
- the 95th percentile of a bootstrap that resamples the panel's batches with
  replacement, each draw the combined ratio `Σ w·leak / Σ w·pass` with fixed
  corpus weights, seeded by `SplitMix64(seed ^ draw · 0x9E3779B97F4A7C15)`; a
  draw without a passed sample counts as rate 1.

A gate passes only when all three are within its bound; at zero leaks the
bootstrap bound is 0 and proves nothing alone. The marginal endpoint bounds
the leaks under its measure among the samples a hypothesis passes and
`against` does not (count-based bound), and reports apart the leaks with a discordant position that the one-strand exception lets one strand read. The utility floor applies to the
corpus-weighted pass rate or its one-sided lower bound. A hypothesis passes
when its gates, marginal, and utility pass; the next one is tested only when
the previous one passed.

The outputs are `approval.json` (identities, every endpoint and verdict,
with each marginal's `projected_n` beside the realised newly passed count),
`approval.md` (aggregates only), `leaks.tsv` (every leak of every gate with
an empty `verdict` column for the adjudicators; the only file with sample
identifiers), and `prereg.toml`. Adjudication never changes the verdict.

`dna-eval power` projects an approval run before the freeze. It reads the
development run's `triage.tsv` and the hold-out manifest's stratum sizes,
never the hold-out truth. Per stratum it reports the development pass rate
and the passes it projects. For the whole panel it reports:

- the stratified rate: the development rates weighted by the panel's
  stratum sizes, with the variance `Σ W_h² p_h(1 − p_h) / n_h` over the
  development samples `n_h`;
- its Korn–Graubard one-sided 95 % lower bound, the gates' method: the
  Clopper–Pearson lower bound at the effective sample size;
- the expected passes at each, with the count-based bound attainable at
  0–10 leaks.

A panel stratum without development samples counts at rate zero. The
development rates come from the panel the rules were tuned on, so even the
lower bound may be optimistic.

With `--base` and `--candidate`, two policy names of the development run's
sweep, it projects the samples the candidate passes and the base does not,
the samples a marginal endpoint against the base counts. It reads them from
`triage_policies.tsv` and reports the same projection, with the marginal
bound attainable at 0–3 leaks. With `--true-rate` and one or more `--bound`,
it also reports, per bound, the most leaks among the lower-bound samples
that still meet it and the binomial chance of at most that many at the true
rate. The pre-registration records the lower-bound projection as
`projected_n`.

## Notation convention

`dna-eval convention` measures how far DNA's notation agrees with the way the
lab's reviewers write their own calls, without traces. Every truth sample
whose calls are plain A/C/G/T substitutions, deletions, and insertions is
rebuilt as anchored variants (runs of deleted positions as one deletion,
insertion ordinals `1, 2, …` without a gap as one insertion) and goes through
DNA's right-aligned [[design/variant-normalization|normalization]] and
the profile's nomenclature windows, as DNA's own calls do, with a deletion
trade bounded by the windows and every base taken as read. A sample whose
written calls come back unchanged is a fixed point of DNA's notation; a
changed one shows a convention DNA writes differently. A sample with an
IUPAC call, or calls DNA cannot apply, is skipped.

The output directory, new or empty, holds `convention.tsv`: per change
pattern, most samples first and at most 200, its sample count, the
truth-only calls, and the DNA-only calls. It names no sample. The command
prints the samples checked, unchanged (with their share), changed, the
patterns, and the samples skipped. The truth is only read; it is never
rewritten or scored as DNA output.

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
  records each sample's truth digest when the panel was built. Other columns
  are ignored. Sample identifiers are unique.
- Bounds: 256 MiB, 1,000,000 rows, and 64 traces per sample.

**Truth.** The lab's statistics database (`statistics.sqlite`, at most
4 GiB), opened read-only. It is the indexed copy of the lab's merged reviewer
results (`statistics.json`) and the single truth source of the harness. Its
`info` table must declare `schema_version` 1. Per panel sample the harness
reads the batch (`samples`), the reviewed alleles (`variants`: `pos`, `ref`,
`seq`, with `-` for an inserted or deleted base; at most 10,000 per sample),
the analysed intervals (`intervals`), and the sample's status in its batch
(`sample_status`). Sample flags and counts are not used.

**Baseline.** The `scores.tsv` of an earlier run. `sample_id` and `status`
are required; `strict_status` and `iupac_compatible_status` are compared when
present.

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

`--triage` also takes each sample's triage evidence from the same analysis
and runs the [draft triage](#draft-triage); without it the run and its
outputs are unchanged.

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
  a trace to its document: `Cargo.lock`, `dna-kernel`, `dna-sanger`,
  `dna-trace`, and the facade's configuration, trace loading, evidence
  models, and writer, recorded by `dna-eval`'s build script. A harness built
  with any of those sources differing from its commit has no identity, and
  the option fails. A plugin version does not move on a work-in-progress
  commit, so it alone would serve stale documents.
- **Storage.** Entries are gzip-compressed and renamed into place. An entry
  that does not decompress, is not JSON, or names another trace is a miss,
  remade and replaced.
- **Record.** `run.json` records the directory, the key, the hits, misses,
  and replaced entries, and the generation's size. `cache prune --keep <n>`
  removes every generation but the `n` most recently used.
- **Placement cache** ([[ADR-0087]], proposed). The same directory holds
  each evidence document's `dna.placed_read/v1` document, written by
  `dna-align`. Entries are keyed by the evidence document's SHA-256, in a
  generation named for the configuration, the placement stage's code
  identity, and the reference. The code identity is the git object ids of
  `Cargo.lock`, the kernel, `dna-align`, and the facade's configuration,
  evidence reading, and writer. Every read is then placed from its document,
  which `dna-align` re-checks against the read's evidence and the reference.
  A change to the placement stage is a new generation. `run.json` records it
  as `placement_cache`. On development panel v5 the tables are
  byte-identical with a cold placement cache, a warm one, and none. A
  consensus-only probe takes 29-36 % less CPU.
- **Limits.** Not with `--triage`, which reads the traces' own Sanger
  products. `approve` always runs from the traces.

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
| `scores.tsv` | always | per sample: the comparison under each [measure](evaluation-measures.md#measures), the [discrepancy class](evaluation-measures.md#discrepancy-classes) and [cause](evaluation-measures.md#discrepancy-causes), the [EMPOP query](#empop-queries), the callable yield per HV region, the traces, the failure, and the [review workload](evaluation-triage.md#review-workload); read back as a later run's `--baseline` |
| `components.tsv` | always | per phase-mixture component (ADR-0082) and reference run cluster it covers: the run it enters in phase against the truth's |
| `tracts.tsv` | always | per tract decode (SRS-CALL-025, [[PROP-0012]] increment 2): the read's reason for masking the onset call and whether it stays out of phase to its end, where the tract's first call lies on the reference, its chosen structure, and its first-run and total dominant lengths, margins, and margin states; evidence only, not compared with the truth |
| `checks.tsv` | always | per sample and `dna review` check item: its positions and labels |
| `unreviewed.tsv` | always | per sample and unreviewed `haplotype_lod` difference block: its anchors and text, its class as the pass measures split it (non-het, het, truth `N`, re-check), the zone that holds its first anchor, and whether a check item touches that zone; `dna-eval diff` counts the passes from it |
| `windows.tsv` | always | per review window: the samples it reaches, checks, and unreviewed differences |
| `estimates.tsv` | always | per weighting, measure, and metric: the stratified rates and interval |
| `anchors.tsv` | always | the anchor candidates under each measure, against the baseline |
| `forms.tsv` | always | how the reviewers write each sequence over the whole truth database, keyed by the database file's SHA-256, without sample identifiers; read back with `--forms` |
| `run.json` | always | the run header: harness and source revision, inputs and identities, options, measurement constants, counts per measure, what passes without review (the [passes](evaluation-triage.md#review-workload) of the review workload), and timings |
| `summary.md` | always | the run in aggregate: the panel checked against the truth (truth digest state, changed, retired samples, samples without a record, batch differences, status not `OK`), outcomes and estimates under each measure, primary discordance, discrepancy classes (per class and window, retired samples left out, the causes per class and the samples per cause and evidence, and with a baseline the samples whose class changed), the EMPOP queries, the callable yield per HV region (samples, analysed and callable positions, their share, and the samples short of full yield), coverage gaps, IUPAC alleles, representation (discordant samples a haplotype measure finds concordant, by window, and the blocks of the rest), anchors, and the baseline |
| `diff.tsv` | with `--baseline` | the samples whose status class changed, per measure the baseline reports |
| `truth_changes.tsv` | when the panel's truth changed | the samples whose truth digest differs from the manifest's record |
| `signals.tsv` | with `--run-signal` | per run-signal group: reads, calls, and keys |
| `triage.tsv`, `triage_regions.tsv`, `triage_traces.tsv`, `triage_leaks.tsv`, `triage_detectors.tsv`, `triage_positions.tsv`, `triage_sweep.tsv`, `triage_policies.tsv`, `triage_minor.tsv` | with `--triage` | the [draft triage](evaluation-triage.md#draft-triage) per sample, region, trace, leak, detector, and region position; per policy of the sweep and of the minor-signal measurement; and each sample's class under every policy, which `dna-eval power` reads |
| `triage.md` | with `--triage` | the triage in aggregate: A/B/C overall and by stratum, leaks and what the leaked samples differ in, region outcomes (counts and weighted shares) and region leaks, leak classes, detectors, traces and primer footprints, the sweep, and the minor signal |

The Markdown reports are for people: the contract requires them, not their text.

Every file except `run.json` depends only on the panel, the truth, DNA's
calls (and with `--triage` its evidence), and the configuration, reference, profile, baseline, and re-check list the summary
names, so it is byte-identical whatever the thread count and whatever revision
or database file it ran on; `forms.tsv` alone also records the database
file's SHA-256 in its first line. Run outputs and
panel manifests hold sample identifiers and stay outside the repository under
the [[governance/data|data policy]]; DNA documents quote aggregates only.
