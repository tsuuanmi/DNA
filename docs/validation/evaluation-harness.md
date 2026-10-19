---
implements: [SRS-ALN, SRS-API, SRS-CALL, SRS-COMPAT, SRS-NOM, SRS-SAMPLE]
decided-by: [ADR-0065, ADR-0075, ADR-0076, ADR-0077]
---

# Evaluation harness

`dna-eval` runs DNA over a fixed test panel and compares the calls with
reviewer truth ([[PROP-0005]],
[[ADR-0075]]). It is
the only place in the repository for truth formats and comparison
(SRS-COMPAT-003). Its outputs are validation evidence about a DNA revision,
not DNA contracts. Source: [crates/dna-eval](../../crates/dna-eval/src/README.md).

Measures, discrepancy classes, and causes are in [[validation/evaluation-measures]]; review workload, estimates, and the draft triage are in [[validation/evaluation-triage]].

## Commands

```bash
cargo run --release -p dna-eval -- run \
  --manifest <panel.tsv> --truth <statistics.sqlite> \
  --config config/dna.toml --reference references/rCRS.fasta \
  --output <new-run-directory> [--baseline <earlier-run-directory>] \
  --threads <1-256> [--triage] [--forms <earlier-run-directory>/forms.tsv] \
  [--review-window <bp>] [--forced-windows <file>] [--current-since <YYYY-MM>] \
  [--run-signal batch|month [--run-signal-batch-reads <n>]] \
  [--zone-after <bp>] [--zone-edge <bp>] [--zone-overrun <bp>]

cargo run --release -p dna-eval -- panel \
  --manifest <parent-panel.tsv> --truth <statistics.sqlite> \
  --version <v2> --created <YYYY-MM-DD> --output <new-panel.tsv>

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

| File | Content |
|---|---|
| `scores.tsv` | per sample, sorted by `sample_id`: panel and truth batch, stratum, weights, anchor and retired flags; primary status (`concordant`, `discordant`, `not_applicable`, `failed`) and code (`Y`, `N`, `N/A`), allele counts, unique alleles of each side, rules, positions, `iupac_only`, windows; strict status, code, and `uncovered` alleles; IUPAC-compatible status and code; haplotype, haplotype IUPAC-compatible, and haplotype LoD status, code, and differing blocks of the way that decides each interval that hold a discordant position (`<start>-<end>[+]:<DNA>/<truth>`, `+` when open, `-` for an empty sequence, a sequence over 24 symbols as its first and last 10 and its length, `CCCTCCCCCC..TTCTGGCCAC(40)`); `below_lod`, the truth IUPAC positions below the limit of detection; the [discrepancy class](#discrepancy-classes) and its windows (`class`, `class_regions`); its [cause](#discrepancy-causes), scope, and findings (`cause`, `cause_scope`, `cause_evidence`); the [EMPOP query](#empop-queries) (`empop_ranges`, `empop_mitotype`) and the positions inside its ranges that the consensus left undecided or decided by a tie (`empop_undecided`, `empop_ties`); DNA's callable yield per HV region (`callable_yield`, `HV1=330/342 HV2=268/268 HV3=139/139`: the truth's analysed positions there that DNA's callable intervals cover); distinct and rejected traces; failure (`<kind>: <message>`); the [review workload](#review-workload): check items, positions to review, whether the sample is clear (`Y` with no check item), unreviewed blocks, and with `--forced-windows` the blocks still unreviewed, then the zones to check and their number (`review_checks`, `review_positions`, `review_clear`, `unreviewed`, `unreviewed_forced`, `review_zones`, `review_zone_count`) |
| `checks.tsv` | per sample and check item, in sample then item order: `first` and `last` (1-based) and the item's label names, space-separated |
| `signals.tsv` | with `--run-signal`, per group (`batch:<batch>`, `month:<YYYY-MM>`): reads and selected calls, keys, and keys with at least 20, 40, and 100 calls |
| `windows.tsv` | per review window, by region and start: the counted samples whose analysed intervals reach it, those with a check item in it, those without one that have an unreviewed difference in it, that share, and its Wilson 95 % upper bound |
| `estimates.tsv` | per weighting, measure, and metric: each stratum, then `all` with the interval |
| `anchors.tsv` | anchor candidates with their status under each measure, their baseline status, and the primary change |
| `diff.tsv` | with a baseline: per measure the baseline reports, the samples whose status class changed, including added and removed ones |
| `truth_changes.tsv` | whenever the panel's truth changed: each sample whose truth digest differs, with its stratum, retired flag, and recorded (`-` when none) and current digests |
| `summary.md` | the run in aggregate: the panel checked against the truth (truth digest state, changed, retired samples, samples without a record, batch differences, status not `OK`), outcomes and estimates under each measure, primary discordance, discrepancy classes (per class and window, retired samples left out, the causes per class and the samples per cause and evidence, and with a baseline the samples whose class changed), the EMPOP queries, the callable yield per HV region (samples, analysed and callable positions, their share, and the samples short of full yield), coverage gaps, IUPAC alleles, representation (discordant samples a haplotype measure finds concordant, by window, and the blocks of the rest), anchors, and the baseline |
| `triage.tsv` | with `--triage`, per sample under the draft policy: outcome, region outcomes, findings with their position counts (and `minor_unavailable=…` when the minor channels could not be read), poor traces, the traces a rerun asks for, and the status under each measure |
| `triage_regions.tsv` | with `--triage`, per sample and region: outcome, positions per detector, rerun traces, and the discordant positions under each measure |
| `triage_traces.tsv` | with `--triage`, per input trace: primer, state, callable span, findings, and the regions whose rerun asks for it (`no_admitted_read` for a sample without an admitted read) |
| `triage_leaks.tsv` | with `--triage`, every discordant position inside a passed region under each measure: the sample's outcome, region, position, window, kind, class, the alleles of each side, and the detectors marking it |
| `triage_detectors.tsv` | with `--triage`, per region detector: whether the draft applies it, samples, regions, and positions flagged, catches per measure, burden, and hits; unique catches and sole burden for an applied detector, `if_applied_*` catches and burden for one not applied, `-` in the other set |
| `triage_positions.tsv` | with `--triage`, per region position: the analysed samples each detector marks it in |
| `triage_policies.tsv` | with `--triage`, per sample (retired ones included, with `retired`): its stratum and its class under each policy of the sweep, one column per policy named as in `triage_sweep.tsv` (`pass`, `review`, `rerun`, or `failed`) |
| `triage_minor.tsv` | with `--triage`, per policy of the minor-signal measurement: its rule, A without and with it and the cost, moved samples and burden, IUPAC sites and those caught (in regions that pass without the rule, and in any region), and per measure the passed samples and leaks without and with the rule, the upper bound, and the corpus rate |
| `triage_sweep.tsv` | with `--triage`, per policy: its parameters, A/B/C counts (with the `failed` samples) and corpus-weighted shares, and per measure the passed samples, leaks, rate, upper bound, and corpus rate |
| `triage.md` | with `--triage`, the triage in aggregate: A/B/C overall and by stratum, leaks and what the leaked samples differ in, region outcomes (counts and weighted shares) and region leaks, leak classes, detectors, traces and primer footprints, the sweep, and the minor signal |
| `forms.tsv` | how the reviewers write each sequence over the whole truth database, keyed by its SHA-256: per repeat window and window sequence the samples per written form, and the samples per set of calls outside the windows (`-` for an empty one); no sample identifiers |
| `run.json` | the run header: harness version, the source revision and whether tracked files differed from it, recorded when the harness was built, threads, manifest fields and hash, the panel truth digest, the recorded one, the truth state and the number of changed samples, the database's file hash and `info` (provenance only), paths, DNA's configuration, reference, and profile identities, the review zone widths (`zones`: `after`, `edge`, `overrun`), the harness's measurement constants (`constants`: limit of detection, faint and minor shares, noise calls and ratio, local shadow share, upper-bound level, interval quantile, triage window margin and footprint), counts per measure, timings, and with `--triage` the draft policy, the number of policies, and the triage time |

Every file except `run.json` depends only on the panel, the truth, DNA's
calls (and with `--triage` its evidence), and the configuration, reference, profile, and baseline the summary
names, so it is byte-identical whatever the thread count and whatever revision
or database file it ran on. Run outputs and
panel manifests hold sample identifiers and stay outside the repository under
the [[governance/data|data policy]]; DNA documents quote aggregates only.
