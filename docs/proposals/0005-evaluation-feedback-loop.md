---
id: PROP-0005
type: proposal
status: draft
owners: []
created: 2026-10-14
related-requirements: [SRS-NFR-003, SRS-NFR-006]
related-decisions: [ADR-0018, ADR-0065, ADR-0066, ADR-0069]
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
| rule parameter or rule code | development split (about 2,000 samples) | under 10 s |
| rule parameter or rule code | full corpus | under 2 min |
| evidence algorithm (DNA itself) | full corpus | under 30 min |

## Non-goals

- Moving truth formats or comparison into DNA without a decision. ADR-0065
  keeps them out today. The options below make that choice explicit.
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

DNA stays free of truth and comparison. To make the loop possible it
provides:

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
| **A. Separate Rust repository** (recommended), depending on DNA's crates by path or git revision | Keeps ADR-0065 and ADR-0066. Calls DNA in process with no process start or JSON parse per stage. Owns truth adapters, manifests, cache, the parallel runner, and reports. | DNA and the harness are versioned separately. Local work needs a path dependency on the DNA checkout. |
| **B. Crate in the DNA workspace**, outside the `dna` binary | One change spans both, and the loop is the tightest. | Brings lab truth formats and comparison into DNA, so it needs an ADR that supersedes ADR-0065. |
| **C. `mtdna_raw` (Python) drives the CLI** | Nothing new to build. | A process per sample, JSON on every stage, and every rule change re-runs the CLI, so the loop budgets cannot be met. |

### Evaluation method

**Truth and normalisation.**

- **What counts.** The reviewer variants of each sample, inside the analysed
  intervals. IUPAC sites are scored on their own.
- **Same notation on both sides.** Both sides pass through DNA's notation
  plugin, called as a library. A difference in writing alone is then never
  an error.

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

**Splits.** Splits are deterministic and recorded in a manifest with its
hash.

- **Hold-out.** Whole batches, including the most recent ones, so batch
  effects and drift are tested. Hold-out is only for approval and is never
  used to tune.
- **Development.** About 2,000 samples, chosen by sample-id hash and
  stratified by batch, for the seconds-scale loop.
- **Tuning.** All remaining samples.

**Statistical power.** With zero leaks in *n* passed samples, the 95 % upper
bound on the error rate is about 3/*n*. A hold-out of 20,000 samples with
60 % passed gives about 0.025 %.

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

## Alternatives considered

- **No cache; make DNA fast enough.** Rejected as the only path. Even a ten
  times faster placement re-runs about 13 CPU-hours for each threshold
  change. Faster placement is still worth having for evidence changes.
- **Score only exact samples.** Rejected. It hides where errors sit and
  cannot score a region or a rule.

## Rollout

1. **Inventory.** Which of the 100,000 truth samples have raw traces, in
   which batches, and with which rerun history.
2. **Harness v0.** Stage E reads existing CLI result JSON as its evidence, so
   DNA needs no change. Stage R reproduces the stage-1 detectors from that
   JSON for the phase-0 numbers of PROP-0004.
3. **DNA library calls and evidence record.** Placement once per trace.
4. **Placement speed-up.** Measured against the corpus.

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

1. Where does the harness live: option A, B, or C?
2. Should the evidence record be a new compact contract, or the existing
   consensus and sample documents?
3. Which batches form the hold-out, and who may touch it?
4. Who adjudicates leaks, and where is the list kept?
5. How much of the server may the loop use: all 256 cores, or a share?

## Decision

Not yet reviewed.

## Implementation

None.
