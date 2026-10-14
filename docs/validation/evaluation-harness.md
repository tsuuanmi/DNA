# Evaluation harness

`dna-eval` runs DNA over a fixed test panel and compares the calls with
reviewer truth ([PROP-0005](../proposals/0005-evaluation-feedback-loop.md),
[ADR-0075](../decisions/adr/0075-evaluation-harness-in-workspace.md)). It is
the only place in the repository for truth formats and comparison
(SRS-COMPAT-003). Its outputs are validation evidence about a DNA revision,
not DNA contracts. Source: [crates/dna-eval](../../crates/dna-eval/src/README.md).

## Commands

```bash
cargo run --release -p dna-eval -- run \
  --manifest <panel.tsv> --truth <statistics.sqlite> \
  --config config/dna.toml --reference references/rCRS.fasta \
  --output <new-run-directory> [--baseline <earlier-run-directory>] \
  --threads <1-256>

cargo run --release -p dna-eval -- panel \
  --manifest <parent-panel.tsv> --truth <statistics.sqlite> \
  --version <v2> --created <YYYY-MM-DD> --output <new-panel.tsv>
```

Every path is an argument. `run` needs a new or empty output directory, and
`panel` a new output file; nothing is overwritten. `panel` prints the new
body's checksum line (`<sha256>  <file> (body)`). The harness is built from
the workspace but not shipped in a release.

## Inputs

**Panel manifest.** `#` comment lines with `key=value` fields come first, then
one tab-separated column header, then one row per sample.

- `body_sha256` is required: the SHA-256 of every byte after the column header
  line. A mismatch stops the run. `rows`, when present, must equal the row
  count. `panel`, `version`, `salt`, and `truth_sha256` are recorded.
- The columns `sample_id`, `batch`, `stratum`, `weight`, `weight_corpus`,
  `anchor_candidate` (`0` or `1`), and `trace_paths` (`;`-separated) are read.
  The optional `retired` column (`0` or `1`, from panel version 2) marks
  samples kept for continuity but left out of the estimates. Other columns are
  ignored. Sample identifiers are unique.
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
  weights `0`), and the other samples of its stratum share the stratum's
  weight total (written with four decimals), so each stratum still sums to its
  population; strata without a newly retired sample keep their weights as
  written;
- the header records the panel, the new version, the salt, the creation date,
  the row and retired counts, the parent version and body hash, and the
  database path, SHA-256, `created_at`, and sample count as the truth.

Earlier versions stay runnable; report results on the rows both versions
share.

## Running DNA

For each sample the harness removes repeated trace content (keeping the first
path of each content), then runs `dna::sample_analysis::analyze_sanger` with
one worker and `SampleAnalysis::call`: the `consensus` then `call` workflow in
process, with the configuration and reference given. Its notation calls become alleles (`309.1C` →
`309.1`, `-`, `C`; `249DEL` → `249`, `A`, `-`), and the consensus reads'
callable reference segments become DNA's intervals. Samples run on a pool of
scoped threads (`--threads`). A sample whose traces cannot be read, that DNA
rejects, or on which DNA panics gets a failure row; the others continue. The
truth loads on its own thread while DNA runs.

## Measures

The truth side is compared as stored in the database; DNA's side uses its
notation. Every sample is compared under three measures.

**Primary.** The lab's concordance (`mtdna_raw`, `src/core/comparison.py`,
`SampleComparator.is_concordant`), with DNA as tool A (the automated caller)
and the truth as tool B:

- an allele is `(position, reference, alternate)`. Positions are normalised
  (`16193.0` → `16193`, thousands separators removed), and alleles whose
  reference equals the alternate are dropped;
- `N/A` when either side has no allele;
- insertions and deletions at 309, 455, 463, and 573, and at their `.n`
  insertion positions, are left out ("special" positions); substitutions
  there count;
- `N` when a position both sides call holds different allele sets, except at
  309.1, 309.2, and 573.1;
- `N` when DNA alone calls an allele inside HV1 (16024–16365), HV2 (73–340),
  or HV3 (438–576);
- `N` when the truth alone calls an allele inside an HV region and inside
  DNA's callable intervals, or DNA reports none;
- `Y` otherwise. Interval and region tests use an allele's anchor base, so
  `16193.1` lies in `[16024, 16193]`.

A truth allele where DNA has no callable base is therefore not counted, and
a truth IUPAC allele differs from any definite base DNA calls.

**Strict.** The primary rules, and a truth-only allele in an HV region also
counts inside the truth's analysed intervals: where DNA has no callable base
it is missed (rule `uncovered`). This makes DNA's coverage gaps visible.

**IUPAC-compatible.** The primary rules, and DNA's base matches a truth IUPAC
code that admits it, like `iupac_bases_compatible`: two alleles with the same
position and reference match when their called bases share a possible base
(`16189C` and `16189Y`), and where DNA calls nothing inside its callable
intervals its base is the reference, which matches a truth IUPAC allele whose
code admits the reference (`16189Y` over a reference `T`).

For each sample the harness also reports, on the primary measure, the alleles
unique to each side (special positions left out), every rule that applies,
the positions those rules name, whether every such position holds a truth
IUPAC allele (`iupac_only`), and their windows: `hv2_polyc` 302–316,
`hv1_polyc` 16180–16195, `hv3_polyc` 560–576, `ac_repeat` 513–526, `other`.
It lists the `uncovered` truth alleles as well.

## Estimates

For each weighting (`corpus` = `weight_corpus`, `pool` = `weight`), each
measure with the metrics `not_concordant` (anything but `Y`, failures
included) and `discordant`, and `failed` (measure `all`), over samples of
positive weight that are not retired:

- rate `p = Σ w·y / Σ w`, per stratum and over all;
- stratum size `N_h = Σ_h w`, `W_h = N_h / N`;
- variance `Σ_h W_h² · p_h(1 − p_h) / n_h · max(0, 1 − n_h / N_h)`;
- 95 % interval `p ± 1.96·√variance`, and the unweighted rate.

## Anchors and baseline

Statuses fall into classes so that runs of other conventions stay comparable:
`exact` and `concordant` are good; `error` and `discordant` are bad;
`not_applicable` and `failed` are their own classes.

- Anchor candidates are the manifest's `anchor_candidate = 1` samples, listed
  under each measure. With a baseline, the anchor set is the candidates whose
  baseline `status` was good.
- Changes against the baseline, per measure the baseline reports:
  `regression` (good before, not now), `fix` (the reverse), `changed` (another
  change of class), `added`, `removed`. The summary also tabulates the
  baseline `status` against the primary status.

## Outputs

| File | Content |
|---|---|
| `scores.tsv` | per sample, sorted by `sample_id`: panel and truth batch, stratum, weights, anchor and retired flags; primary status (`concordant`, `discordant`, `not_applicable`, `failed`) and code (`Y`, `N`, `N/A`), allele counts, unique alleles of each side, rules, positions, `iupac_only`, windows; strict status, code, and `uncovered` alleles; IUPAC-compatible status and code; distinct and rejected traces; failure (`<kind>: <message>`) |
| `estimates.tsv` | per weighting, measure, and metric: each stratum, then `all` with the interval |
| `anchors.tsv` | anchor candidates with their status under each measure, their baseline status, and the primary change |
| `diff.tsv` | with a baseline: per measure, the samples whose status class changed |
| `summary.md` | the run in aggregate: the panel checked against the truth (retired samples, samples without a record, batch differences, status not `OK`), outcomes and estimates under each measure, primary discordance, coverage gaps, IUPAC alleles, anchors, and the baseline |
| `run.json` | the run header: harness version, source revision, threads, manifest fields and hash, the database's hash and `info` and whether the manifest names it, paths, DNA's configuration, reference, and profile identities, counts per measure, timings |

Every file except `run.json` depends only on the panel, the truth, and DNA's
calls, so it is byte-identical whatever the thread count. Run outputs and
panel manifests hold sample identifiers and stay outside the repository under
the [data policy](../governance/data.md); DNA documents quote aggregates only.
