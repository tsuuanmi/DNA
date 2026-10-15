---
id: PROP-0005
type: proposal
status: draft
owners: []
created: 2026-10-14
related-requirements: [SRS-NFR-003, SRS-NFR-006]
related-decisions: [ADR-0018, ADR-0065, ADR-0066, ADR-0069, ADR-0075]
implementation: []
---

# Proposal: Evaluation method and fast feedback loop

## Problem

[PROP-0004](0004-review-triage.md) lets a rule act only after ground truth
has validated it. The lab has about 100,000 samples with raw traces and
reviewer results. Such a corpus can support strong claims, but only if two
things exist:

- **A method.** Fixed labels, splits, metrics, and reports, so that a number
  means the same thing in every run.
- **A fast loop.** Changing a threshold and seeing its effect must take
  seconds to minutes, not a day.

DNA is Rust and runs on a 256-core server. One run is not the bottleneck;
recomputing everything on every change is.

**Measured on 2026-10-14** (release build of main, one thread per process):

| Command | Wall time | Peak memory |
|---|---|---|
| `basecall` on one trace | 0.07 s | 5 MB |
| `analyze` on one trace | 1.9 s | 34 MB |
| `consensus` on one sample (4 traces) | 4.7 s | 52 MB |
| `sample` on the same sample | 4.6 s | 52 MB |

- **Placement dominates.** Reference placement takes about 95 % of the time.
- **Placement runs twice.** `sample` and `consensus` each place every trace
  again.
- **Scale.** 100,000 samples take about 130 CPU-hours per command. That is
  about 30 minutes on 256 cores at best, and today every rule change would
  repeat it.

## Goals

- **One method.** An evaluation method that defines labels, splits, metrics,
  and reports for detection, diagnosis, rerun, and per-rule approval.
- **Reproducible reports.** The same DNA revision, policy version, and dataset
  manifest always give the same report.
- **Loop budgets** on the server:

| Change | Scope | Budget |
|---|---|---|
| rule parameter or rule code | development panel | under 10 s |
| rule parameter or rule code | validation panel | under 1 min |
| evidence algorithm (DNA itself) | validation panel | under 10 min |
| any, before a release | full corpus | a scheduled run, not part of the loop |

## Non-goals

- Truth formats or comparison anywhere outside the harness crate. The `dna`
  binary and the plugin crates stay free of them.
- A review user interface.
- Biological rules. They are the subject of the loop, not of this proposal.

## Proposed design

### Three stages with a cache between them

```mermaid
flowchart LR
  RAW[AB1 traces] -->|E: evidence, slow| EC[(evidence cache)]
  EC -->|R: rules, fast| OUT[findings and outcomes]
  TRUTH[(reviewer truth, rerun history)] --> C
  OUT --> C[C: compare and score]
  C --> REP[report and diff against the baseline]
```

- **Stage E (evidence).** Decoding, base calling, callability, placement, and
  consensus, per trace and per sample. This stage is expensive. Its output is
  cached under a key made of:
  - the SHA-256 of the trace bytes;
  - the evidence version (plugin versions, the hash of the evidence-relevant
    configuration, the reference checksum).

  Rule parameters are not part of the key, so a rule change never invalidates
  the cache.
- **Stage R (rules).** The pure rules of PROP-0004 run over cached evidence.
  They take microseconds to milliseconds per sample and run in parallel
  across samples.
- **Stage C (compare).** Outcomes are compared with the truth, then scored,
  reported, and diffed against a stored baseline run.

The loop then splits by what changed:

| What changed | Stages re-run |
|---|---|
| rule parameter or rule code | R and C only |
| DNA evidence code | E, R, and C |
| truth adjudication | C only |

### What DNA provides

The `dna` facade and the plugin crates stay free of truth and comparison. To
make the loop possible they provide:

- **Library calls.** A public call that computes a sample's triage evidence
  in process, and a separate call that runs the rule layer on such evidence
  (SRS-API).
- **An evidence record.** A versioned, serialisable record that holds exactly
  what the rules read, so a cache can store it. Whether this is a new compact
  record or the existing `dna.consensus` and `dna.sample_evidence` documents is
  open.
- **Placement once per trace.** One in-process call yields both the sample
  evidence and the consensus.
- **Faster placement.** Seed the placement with reference k-mers, then align
  in a band. Maintained implementations are evaluated first (SRS-NFR-006).
  This is a measured change: placements must stay identical on the corpus, or
  the changes must be explained.

### Where the harness lives

| Option | For | Against |
|---|---|---|
| **A. Separate Rust repository**, depending on DNA's crates by path or git revision | Keeps ADR-0065 and ADR-0066. Calls DNA in process with no process start or JSON parse per stage. Owns truth adapters, manifests, cache, the parallel runner, and reports. | DNA and the harness are versioned separately. Local work needs a path dependency on the DNA checkout. |
| **B. Crate in the DNA workspace** (chosen), outside the `dna` binary | One change spans both, and the loop is the tightest. | Brings lab truth formats and comparison into DNA, so it needs an ADR that supersedes ADR-0065. |
| **C. `mtdna_raw` (Python) drives the CLI** | Nothing new to build. | A process per sample, JSON on every stage, and every rule change re-runs the CLI, so the loop budgets cannot be met. |

### Evaluation method

**Truth and normalisation.**

- **What counts.** The reviewer variants of each sample, compared by the
  lab's concordance rules (primary measure), with a stricter measure on the
  truth's analysed intervals and an IUPAC-compatible measure beside it
  ([evaluation harness](../validation/evaluation-harness.md#measures)).
- **Notation.** The truth side is compared as stored in the database,
  following `comparison.py`'s rules; DNA's side uses its notation.

**Labels.**

| Label | Definition |
|---|---|
| site error | false call, missed call, wrong allele, or wrong length at a truth site or DNA call |
| sample exact | no site error in the sample |
| error caught | the site lies inside an unreliable region or a region not passed |
| leak | a site error inside a passed region |
| rerun positive | a trace that the lab's rerun history replaced |
| rerun negative | a trace of a sample whose reviewed result needed no rerun |

The rerun history is the set of rerun lists that record which trace each
rerun replaced. The lab's own rerun decisions are noisy, so these labels are
a proxy.

**Metrics.**

| Area | Metrics |
|---|---|
| detection (stage 1) | leaks; passed-sample error rate with its Clopper–Pearson 95 % upper bound; burden as unreliable bases and flagged sites per sample; A/B/C at sample and region level |
| trade-off | risk–coverage curve: passed-sample error rate against A over a threshold sweep |
| diagnosis (stage 2) | confusion matrix of assigned against inspected cause, on an inspected sample of caught errors |
| rerun | precision and recall against the rerun history |
| per rule | unique catches, false flags, and the change in every metric when the rule is removed (ablation) |
| stability | samples whose outcome changed against the baseline, listed |

**Test panels.** The loop does not run on all 100,000 samples. It runs on
small, fixed panels, each with a manifest and its hash. A panel holds two
kinds of samples:

- **Anchors.** Samples whose result is exact and stays exact from run to run.
  They catch regressions: an anchor that changes is a failure to explain.
- **Problem samples.** Enriched by stratum, so each known difficulty has
  enough cases. Examples:
  - disagreement between the older pipeline and the truth;
  - IUPAC truth sites;
  - length variants in poly-C windows;
  - rare variants;
  - rerun history;
  - known DNA errors.

The panels:

| Panel | Size | Use |
|---|---|---|
| development | about 1,000–3,000 | every change, seconds |
| validation | about 10,000 | tuning and confirmation, minutes |
| hold-out | whole batches, including the most recent | rule approval only; never used to tune |
| full corpus | all samples with raw traces | scheduled release and drift runs |

Panel rules:

- **Deterministic selection.** Within each stratum, samples are taken in
  order of the hash of the sample id.
- **Weights.** Each stratum carries the inverse of its sampling share. The
  enriched panel can then still estimate corpus rates (A, error rate of
  passed results) without bias.
- **Refreshing.** A refresh adds a new manifest version and keeps the
  previous version runnable, so numbers stay comparable across versions.

**Panel v1 (built 2026-10-14, outside the repository).**

Inventory:

- the truth holds 103,418 samples;
- 97,306 of them have raw traces, and 96,841 have all of HV1, HV2, and HV3.

Each sample goes to the first stratum it matches. The strata use only
information that needs no DNA run, so their weights hold for the whole
corpus.

| Stratum (first match) | Corpus | Dev | Baseline not exact |
|---|---:|---:|---:|
| partial traces | 465 | 120 | 0.93 (mostly missing traces, not call errors) |
| reviewer flag or short range | 1,720 | 120 | 0.29 |
| HVS-I length variant | 462 | 150 | 0.60 |
| unusual HVS-II length | 1,196 | 150 | 0.17 |
| IUPAC truth site | 8,371 | 200 | 0.17 |
| rare variant | 351 | 150 | 0.22 |
| older pipeline disagrees, point | 6,582 | 200 | 0.30 |
| older pipeline disagrees, repeat window | 9,629 | 150 | 0.18 |
| resequenced | 4,018 | 150 | 0.25 |
| 309.2C | 12,263 | 150 | 0.23 |
| 16189C | 14,943 | 150 | 0.15 |
| clean | 37,306 | 300 | 0.01 |

What the strata show:

- **They separate difficulty.** Clean samples fail 1 % of the time; the
  enriched strata fail 15–60 %.
- **Corpus estimate.** Weighted, the baseline not-exact rate is 0.135 ± 0.016
  (95 %). Unweighted, the panel gives 0.258.

The panels:

- **Development panel.** 1,990 stratified samples, plus 630 earlier
  evaluation samples kept at zero weight.
- **Anchors.** Clean samples on which the older pipeline agrees with the
  truth. 294 of 297 were exact in the baseline. Two repeated runs gave
  identical results.
- **Validation and hold-out panel.** 115 whole batches with 10,109 samples:
  - all batches from the most recent reviewed weeks that were not used for
    tuning;
  - batches drawn by hash.

  It is not run until approval.

Measured run time for `consensus` plus `call`, about 6.3 CPU-seconds per
sample:

| Set | CPU-hours | 64 workers | 200 workers |
|---|---:|---:|---:|
| development (2,620) | 4.7 | 4.5 min (measured) | about 1.5 min |
| validation (10,109) | 18 | about 17 min | about 6 min |
| full traced corpus (97,306) | 175 | about 2.8 h | about 55 min |

**Faster placement (2026-10-14, PR #51, in review).** Four exact changes
keep every result byte-identical to the current release on 299 samples
(both `consensus` and `sample`) and on a separate check of 40 panel samples:

- a single-pass Gotoh kernel (about 4× per cell);
- seeded q-gram pruning that tolerates masked and mixed loci and crosses the
  circular origin;
- an orientation bound taken from the first orientation;
- selection instead of sorting for signal medians.

Results:

- **Fast tier usage.** Today 100 % of real reads fall back to full Gotoh in
  both orientations; with the prototype 96 % need none.
- **One sample** takes about 0.13 s instead of about 3.8 s on one thread.
- **Many samples.** One process with a thread pool reaches about 280–430
  samples per second on the shared server: about 25–40 s for 10,000 samples
  and 4–6 min for the full corpus.

The pruning change alters a proof (SRS-ALN-014) and needs scientific review
before it lands. With it, the validation panel fits the rule-change budget
even without the evidence cache. The cache still matters for the
seconds-scale development loop and for the full corpus.

Caveats:

- Strata built on the older pipeline inherit its biases.
- The truth is not certain.
- Some traces are copied into two batch folders. These copies are removed,
  but they weaken the batch hold-out.

**Statistical power.** With zero leaks in *n* passed samples, the 95 % upper
bound on the error rate is about 3/*n*. A panel can show which rule is
better. An absolute error-rate claim needs the hold-out or the full corpus:
10,000 hold-out samples with 60 % passed give about 0.05 %.

**Adjudication.** Every leak is inspected.

- A truth-origin case leaves the error count only with a recorded reason.
- A method-origin case goes on DNA's work list.
- The adjudication list is versioned with the harness.

**Report.** Every run reports:

- the DNA revision, policy version, and manifest hash;
- the stage timings;
- the summary metrics and the per-rule table;
- the list of leaks;
- the diff against the baseline.

The report is written both as machine-readable TSV/JSON and as Markdown.

### Validation of the harness itself

- **Byte-identical.** Running twice gives a byte-identical report.
- **Synthetic cases.** Synthetic cases with known labels cover every label
  and metric.
- **Known numbers.** On the 320-sample evaluation subset, the harness
  reproduces the numbers measured so far.

### Approval protocol (draft)

This protocol approves pass rules on the hold-out (validation panel v3: 115
whole batches, 10,102 active samples, the same corpus weights as the
development panel). It is a draft. **The numeric bounds below await the
user's decision.**

**Freeze and pre-register.**
- **Pre-registration note.** Before any approval run, a pre-registration note
  is committed under `docs/validation/preregistration/` that fixes:
  - the DNA revision, `config/dna.toml`, the harness version, and the truth
    digest;
  - the exact policy: its detectors, the minor-signal kind and threshold, and
    the one-strand rule;
  - the measures, endpoints, bounds, and decision rule;
  - the expected number of passed samples.
- **One run.** The run happens once, and nothing is tuned after it. A rerun
  needs a new revision and a written reason.
- **Truth changes.** If the truth digest differs from the frozen one, the
  run aborts.
- **A spent hold-out.** Once a rule is approved, the hold-out is spent for it.
  A later rule change needs fresh batches, such as those received after the
  freeze.
- **Measurement fixes first.** Fixes that change what "concordant" means
  land before the freeze. These include DNA's callable intervals excluding
  unresolved consensus bases, and the haplotype measures deciding per
  interval.

**Rules under test.** The rules are tested as a fixed sequence, which keeps
the overall α without a Bonferroni split:

1. **H-base.** The stage-1 triage with the minor-signal gate (its kind and
   threshold are fixed at pre-registration), requiring both strands
   everywhere.
2. **H-os1**, tested only if H-base is approved. H-base plus the one-strand
   exception R6m-h4: one strand is enough outside the repeat windows ± 3 and
   away from reference runs of 4 or more bases. It is paired with a
   single-strand minor rule on the positions it allows, because the
   cross-strand minor rules cannot fire where only one strand reads. R6m may
   replace R6m-h4 only before the freeze.
3. **H-os2**, tested only if H-os1 fails. H-base plus R1: an edge stretch of
   fewer than 5 positions.

**Endpoints.**

- **Primary safety: the leak rate among passed samples.** This is
  P(not concordant | pass) under the haplotype IUPAC-compatible measure,
  reported three ways, each with its upper bound:
  - as a one-sided 95 % Clopper–Pearson bound on the unweighted passed
    samples;
  - as the corpus-weighted rate with the stratified interval;
  - as the 95th percentile of a cluster bootstrap over batches (10,000 draws,
    fixed seed). Leaks cluster by batch, through chemistry, run, and
    reviewer.
- **Secondary safety, reported but not gating.**
  - The same rate under the primary, strict, and IUPAC-compatible measures.
  - Leaks by window and by stratum.
  - A sensitivity analysis without the enriched strata.
- **Marginal one-strand endpoint** (H-os rules). Among the samples a rule
  newly passes, the share that differ inside the positions it allows, at
  sample level with its upper bound.
- **Utility.** The corpus-weighted A with its interval; B and C are
  reported.

**Acceptance, pending the user's decision.**

| Rule | Approved when |
|---|---|
| H-base | all three upper bounds of the primary safety endpoint are at most 1 % *(pending)* |
| H-os1, H-os2 | the rule still meets H-base's bound on all passed samples, its marginal upper bound is at most 2 % *(pending)*, and adjudication finds no leak class unseen on the development panel |
| any | optionally, A is at least 0.15 corpus-weighted *(pending)*, so that a rule passing almost nothing is not approved by default |

**Power.** With a one-sided 95 % Clopper–Pearson bound:
- **H-base.** A bound of 1 % needs 303 passed samples at 0 leaks, 474 at 1,
  630 at 2, 1,057 at 5, and 1,718 at 10. With an A near the development
  panel's, the hold-out should pass about 1,500–2,500 samples. That gives,
  for example, 0.52 % for 5 leaks in 2,000 passes.
- **H-os rules.** These pass far fewer new samples. A bound of 2 % needs 150
  newly passed samples at 0 leaks, and 237 at 1.
- **Check at freeze time.** A power check from the development panel's pass
  counts and the hold-out's size is recorded in the pre-registration.

**Adjudication in an approval run.**
- Every leak on the hold-out is inspected and classified as a DNA error, a
  truth error, a representation the measure missed, or ambiguous.
- Unlike in development runs (see *Adjudication* above), the gating endpoint
  counts every leak, and no case leaves it after the run.
- The classification is reported next to the decision and feeds the next
  development iteration.
- Sample identifiers stay in the adjudication list, never in repository
  documents.

**Exclusions.** A sample DNA fails counts as not passed and is reported.

**After approval.**
- **Monitoring.** A random 2–5 % *(pending)* of auto-passed production
  samples is reviewed each month, with a running Clopper–Pearson bound.
- **Stop rule.** Auto-pass is suspended when that bound exceeds 2 %
  *(pending)*.

## Alternatives considered

- **No cache; make DNA fast enough.** Rejected as the only path. Even a ten
  times faster placement re-runs about 13 CPU-hours for each threshold
  change. Faster placement is still worth having for evidence changes.
- **Score only exact samples.** Rejected. It hides where errors sit and
  cannot score a region or a rule.

## Rollout

1. **Inventory and panel v1.** Done, outside the repository (see *Panel
   v1*). The panel scripts move into the harness crate.
2. **Harness v0.** Comparison and the run are done (see *Implementation*):
   stage E calls DNA in process instead of reading CLI result JSON, and a
   draft stage R gives the phase-0 numbers of PROP-0004 (`run --triage`).
3. **DNA library calls and evidence record.** Placement once per trace.
4. **Placement speed-up.** PR #51, measured (see *Panel v1*); review and
   land it, then run a combined `sample` and `consensus` pass, share signal
   statistics between stages, and add a pipeline-level thread pool.

## Observability and operations

The harness logs stage timings, cache hit rates, and failures per sample.

## Security and privacy

The corpus holds lab and personal sample data. The harness and its reports
stay on the lab server under the [data policy](../governance/data.md). DNA
documents quote aggregates only.

## Compatibility and migration

DNA gains library calls and possibly an evidence record. Existing contracts
are unchanged unless the record reuses them.

## Open questions

1. Should the evidence record be a new compact contract, or the existing
   consensus and sample documents?
2. Which batches form the hold-out, and who may touch it?
3. Who adjudicates leaks, and where is the list kept?
4. How much of the server may the loop use: all 256 cores, or a share?

## Decision

Decided 2026-10-14: the harness is a crate in the DNA workspace
(`crates/`), outside the `dna` binary (option B). This reverses the
earlier choice of a separate repository.

The crate landed with
[ADR-0075](../decisions/adr/0075-evaluation-harness-in-workspace.md), which
supersedes ADR-0065 in part, the matching revision of SRS-COMPAT-003, and the
crate-graph rule in `validate_module_layers.py`: the harness may depend on
`dna` and the plugin crates, and nothing may depend on the harness.

ADR-0066 is unaffected; the harness is Rust. The rest of the proposal is not
yet reviewed.

## Implementation

- **Placement once per trace** (rollout step 3, SRS-API-009..011):
  `dna::sample_analysis::analyze_sanger` returns the typed sample evidence and
  consensus from one placement pass; `sample` and `consensus` run the same
  pass. The evidence record (open question 1) is not decided.
- **Rollout step 4:** signal statistics are shared between the windowed
  features and the locus contexts, and a sample's reads run concurrently
  with the operation log kept in trace order. On 2,618 development-panel
  samples every `consensus` and `sample` output is byte-identical to the
  previous release.

**Harness v0** (`crates/dna-eval`,
[evaluation harness](../validation/evaluation-harness.md)):

- **Scope.** It runs DNA, compares the calls with reviewer truth, and
  reports weighted estimates, anchors, and the diff against a baseline run.
  The evidence cache and rerun labels are not built yet.
- **Draft stage R** (`run --triage`, 2026-10-15). Draft stage-1 detectors of
  PROP-0004 over a compact per-sample evidence record taken in stage E and
  held in memory; every policy of a sweep runs on that record without DNA.
  Results and timings are in
  [PROP-0004](0004-review-triage.md#phase-0-results).
- **Stage E.** Instead of reading CLI result JSON (rollout step 2), it calls
  the library: `dna::sample_analysis::analyze_sanger` with one worker per
  sample, then `SampleAnalysis::call`, which completes `consensus` then `call`
  in process; the harness runs samples on its own scoped thread pool.
- **Truth.** The lab's statistics database (`statistics.sqlite`, read-only
  through `rusqlite`), the indexed copy of the merged reviewer results; it is
  the harness's single truth source.
- **Comparison.** It adapts the lab's comparison rules (`mtdna_raw`
  `comparison.py`), DNA as the automated tool A and the truth as tool B, as
  the primary measure. A strict measure also counts truth alleles inside the
  truth's analysed intervals where DNA has no callable base, and an
  IUPAC-compatible measure lets DNA's base match a truth IUPAC code that
  admits it. This replaces the panel's research scoring ("fair scale"); the
  differences of the primary measure:
  - the research scoring kept only tokens inside the truth's analysed
    intervals on both sides; the lab's rules count every DNA-only allele in
    HV1–HV3, and a truth-only allele only inside DNA's callable intervals;
  - the research scoring left out every token at a truth IUPAC position; the
    lab's rules compare IUPAC alleles like any other, so DNA's major-allele
    call is discordant;
  - the lab's rules leave out insertions and deletions at 309, 455, 463, and
    573 and never compare 309.1, 309.2, and 573.1;
  - the lab's rules key an allele by position, reference, and alternate, so a
    truth reference that differs from rCRS (25 truth rows) is a difference;
  - an empty side is `N/A` rather than exact; insertions belong to their
    anchor base.
- **Panel v2** (`dna-eval panel`, 2026-10-15). It keeps every v1 row and
  names the statistics database and its SHA-256 as the truth. Samples the
  database lacks are retired and left out of the estimates; their strata share
  the stratum weight among the remaining samples. Dev: 2,620 rows, 2 retired.
  Validation: 10,109 rows, 7 retired. No batch differs from the database. The
  v1 files are unchanged.
- **Panel v3** (`dna-eval panel`, 2026-10-15). It adds each sample's truth
  digest and the panel truth digest, so a database refresh that changes no
  panel sample needs no new version; a run lists the samples whose truth
  changed and marks itself `truth changed` otherwise. Dev and validation keep
  the v2 rows and weights; v1 and v2 files are unchanged. On dev v3 every
  per-sample result equals the v2 run.
- **Dev panel v2** (2026-10-15, revision 56424be plus this change; the
  same per-sample results as on revision 348cf13 before the one-pass sample
  analysis, and as on the database refresh of 2026-10-14 17:31, which changed
  no panel allele):
  - **Panel against the truth.** The 2,618 active samples have the same
    alleles as in the truth file v1 was built from; 2 are not `OK` in their
    batch.
  - **Primary.** 1,958 `Y`, 656 `N`, 2 `N/A`, 4 failed (2 traces too short
    for signal processing, 2 retired without truth). Corpus-weighted
    not-concordant rate 0.1783 ± 0.0138 (95 %), unweighted 0.2852. 295 of 297
    anchor candidates are concordant.
  - **Strict.** 1,818 `Y`, 796 `N`; rate 0.2061 ± 0.0152. 201 samples have
    666 truth alleles where DNA has no callable base, and 140 samples
    concordant on the primary measure are not on the strict one; partial
    traces go from 0.375 to 0.933 not concordant.
  - **IUPAC-compatible.** 2,219 `Y`, 395 `N`; rate 0.1013 ± 0.0143. It
    turns exactly the 261 samples whose primary discordance lies only at
    truth IUPAC alleles concordant; the IUPAC stratum goes from 1.000 to
    0.130.
  - **Against the research baseline** (0.1350 ± 0.0159 not exact): 278
    exact samples are discordant, 247 of them only at truth IUPAC alleles (the
    whole IUPAC stratum is discordant) and 29 at DNA alleles outside the
    truth's analysed intervals. 172 not-exact samples are concordant: 107
    differ only at truth alleles outside DNA's callable intervals, 32 only at
    special-position indels, and 33 at both.
  - **Same rules as the lab's code.** Run through `SampleComparator` on the
    same DNA calls and truth, every one of the 2,616 compared samples gets the
    same `Y`, `N`, or `N/A` on the primary measure.
- **Dev panel v3 after the too-short rejection** (revision d171ce3 plus this
  change, truth unchanged): the 2 samples that failed on a trace too short
  for signal processing now set that trace aside and are concordant on all
  three measures; no other sample changes. Primary 1,960 `Y`, 656 `N`, 2
  `N/A`, 2 failed (the retired samples), rate 0.1779 ± 0.0138; strict 1,820
  `Y`, rate 0.2057 ± 0.0152; IUPAC-compatible 2,221 `Y`, rate 0.1010 ±
  0.0143.
- **Speed.** The dev panel takes 5.5–7 s wall on 64 or 128 threads of the
  shared server (about 320–800 CPU-seconds, with the one-pass sample analysis
  at the low end); a first run under a load spike took 20 s.
- **Determinism.** Every output except the timed `run.json` is
  byte-identical between 64 and 128 threads.
