---
id: PROP-0009
type: proposal
status: draft
owners: []
created: 2026-10-18
related-requirements: [SRS-VAR-012, SRS-API-014]
related-decisions: [ADR-0027, ADR-0065, ADR-0066, ADR-0067, ADR-0069]
implementation: []
---

# Proposal: Run signal for second-channel judgement

## Problem

The Sanger plugin judges a call "mixed" when a second channel qualifies under a fixed ratio
(ADR-0027). The review list labels `possible_heteroplasmy` when a second channel's share
reaches a fixed `0.15`. Both thresholds ignore that secondary signal recurs. The same base
shows the same shadow in the same context in most samples of a run, and the level changes
from run to run.

Recurring shadows read as mixtures, and real minors at quiet places read as noise. As a result:

- Calls that read the reference cleanly apart from a recurring shadow have their identity
  withheld, so the consensus leaves `N` and the review lists them as `undecided`.
  - `16031`: second C at a median share of about 0.27.
  - `527`: second A at about 0.2.
- `possible_heteroplasmy` fires on 18.5 % of samples where the truth records no
  heteroplasmy.
- These recurring checks are the largest single blockers of reviewer workload on the
  development panel.

## Goals

- Judge a second channel against the level its own run shows, keyed so that the modality
  stays reference-free (ADR-0069).
- One scale and one threshold for the Sanger mixed verdict and for `possible_heteroplasmy`.
- Fail closed:
  - without a usable run signal, the configured floor applies;
  - a run signal lifts a mixed verdict only where the run itself shows that shadow as
    common, and never adds one.
- Reproducible: the run signal is an explicit, recorded input, not hidden state.

## Non-goals

- Quantitative heteroplasmy, or IUPAC output.
- Static tables learned from history that may suppress a signal.
- Position-keyed (reference) levels.
- Learning at run time across runs.
- The read-edge trust rule, which is a separate change to SRS-CALL-013/014.

## Evidence and research

These are research spikes, out of repo (ADR-0065/0066). The data are the development panel v5
(16,413 samples) at main `5c0310f` and every scored, trusted call observation in the profile
regions (21.6 million). Each sample's half is `sha256(sample_id) % 2`; rules are judged on the
test half (8,257 samples) against the truth's two-base IUPAC sites (739).

**Key.** A level keyed by the read-strand 5-mer around the call, the read-strand base and
the second base needs no reference position. It recovers most of what a position key gives.

| Rule (`possible_heteroplasmy` shape) | Recall of truth sites | Samples with a false flag |
|---|---:|---:|
| Fixed share 0.15 (today) | 23.7 % | 18.5 % |
| Position key, static, q99 + 0.03 | 23.4 % | 8.9 % |
| Read-strand 5-mer, static, q99 + 0.03 | 23.4 % | 11.0 % |
| Read-strand 3-mer, static | 19.1 % | 16.2 % |

A 7-mer behaves like the position key, because it is nearly unique in the profile regions,
and overfits the same way.

**Drift.** Levels move strongly between runs. Median clean second share by month:

- `16031` (forward, A then second C): 0.26-0.28 until 2025-11, then 0.18-0.23, then about 0
  from 2026-06.
- `527` (reverse, C then second A): about 0, then 0.20-0.23 in 2025-08..11, then about 0
  again.

A static table learned on early runs would place `16031`'s level near 0.27. On later runs
that hides a real 0.25 minor, which is fail-open.

**Time-local level.** Panel batches are subsamples, so a month stands in for a run here.

| Rule | Period | Recall | False-flag samples |
|---|---|---:|---:|
| Fixed 0.15 | all | 23.7 % | 18.5 % |
| | 2025-08..11 | 19.9 % | 40.2 % |
| Static 5-mer table | all | 25.3 % | 11.3 % |
| | 2025-08..11 | 22.0 % | 22.1 % |
| Per-month 5-mer q99 + 0.03 | all | 22.6 % | 9.0 % |
| | 2025-08..11 | 17.7 % | 13.9 % |

Small groups are a trap. Per-batch tables over 24-96 panel samples absorb real minors into a
high quantile: q99 recall falls to 16.8 %. The estimator and the minimum group size must be
chosen together.

**Ceiling.** Only 47 % of the truth's two-base IUPAC sites show a second share of at least 20 %
in any clean read, and 27 % show less than 10 %. About half of the truth heteroplasmy cannot be
seen in the traces at all. That is a truth or convention question, not one this proposal
answers.

## Proposed design

1. **A run signal command**, part of the Sanger signal processing:
   `dna signal <run-id> <trace>...` writes `results/<run-id>.signal.json`
   (`dna.signal/v1`). It reads the traces only; no reference, no alignment and no truth.
   - **Selection.** It keeps the calls that are in phase by the callability mask, not at a
     signal edge (SRS-CALL-013), canonical, and whose four neighbours are canonical. Peak
     floor and relative quality filter as today.
   - **Calls judged mixed today are kept**: their shadow is what is measured. The selection
     never depends on the verdict it informs.
   - **The rule is versioned in the file.** The first and last two calls of a read have no key.
   - **Key.** The read-strand 5-mer of the read's own primary calls around the call, the call's
     base, and its strongest other channel: at most 4^5 x 4 x 3 = 12,288 keys.
   - **Content.** For each key, a sparse histogram of the second-channel share in thousandths.
     The file also holds the run's trace hashes and its counts of reads and calls.
     Histograms rather than quantiles let the estimator change without rebuilding, and make
     leave-one-out exact.
   - A library entry, `dna_sanger::run_signal` (`run_signal::read`, `RunSignal::build`), which sits
     beside read processing rather than in signal processing, since the selection needs the
     callability mask and the signal edges, lets the evaluation harness build run signals per run (or per month) in process. There is also a pipeline
     operation with its own log.
2. **Sample commands accept `--signal <file>`** (`sample`, `consensus`, `review`). The
   public entry `analyze_sanger` takes an optional `&RunSignal`; this is a public API change to
   be specified with SRS-API.
   - **Identity.** The file's SHA-256 is recorded in provenance, beside the configuration and
     the profile.
   - **Leave-one-out.** When the sample's trace hashes are in the run signal, the command
     recomputes their contributions with the file's selection rule and subtracts them before
     estimating, so a sample is never judged against itself.
   - **Run mismatch.** A file from another run (no shared trace hash and another run id) is
     refused with a typed error. If the run id matches but the hashes are absent, the
     run signal is used and the mismatch is flagged.
   - **Two minimums.** Below a minimum of reads for the whole file, or for a key with fewer
     than `minimum_calls` observations, the run signal gives no level and the floor applies.
3. **One scale, one threshold, two uses** (amended in increment (ii)). Every rule that
   judges a call's second channel reads the peak share of SRS-CALL-015,
   `h_X / (h_primary + h_X)`, and one function owned by `dna_sanger::run_signal`,
   `threshold(key, signal, config)` (SRS-CALL-020): `t = max(floor, level)` when the run
   signal, minus the sample, holds the call's key with enough reads and calls, where
   `level = median + max(margin, mad_factor × MAD)` of the key's shares; otherwise
   `t = floor`. The floor is one configuration value, `signal_processing.second_floor`.
   - **The review's `possible_heteroplasmy`** fires when a clean call's second peak share
     exceeds `t`, with the neighbour-shadow and strand rules unchanged. It replaces today's
     `0.15` profile weight, with or without a run signal. The neighbour-shadow rule stays
     (2026-10-18, user): dropping it raised the review workload on the development panel
     (the context level alone does not absorb read-start and read-end shadows), and
     heteroplasmy is optional today. It goes when DNA handles noisy regions itself. The
     published evidence `second` moves to the same peak scale. Review output therefore
     changes without a run signal too, and the `review` plugin version says so.
   - **The Sanger mixed verdict (ADR-0027)**, in increment (iii), uses the same `threshold()`
     in place of its fixed `0.33` peak ratio, and ADR-0027 is amended then. A run signal can
     only raise the threshold above the floor, so it lifts a mixed verdict only where the run
     itself shows the shadow as common. Deferred to the [roadmap](roadmap.md) on 2026-10-19:
     lifting every mixed verdict, its ceiling, did more harm than good.
   - **Such a minor keeps its identity**: the consensus decides the primary base, and the
     review makes the position a check. Reporting it as IUPAC remains a separate, deferred
     decision.
   - **Parameters.** `second_floor`, `margin`, `mad_factor`, `minimum_calls` and the minimum
     reads sit in the existing `[signal_processing]` configuration section. One floor serves
     both uses; it is fixed on development runs and judged on a hold-out panel and on later
     runs.
4. **A sample's result now depends on its run-mates** through one recorded file. Rerunning the
   sample with the same run signal file reproduces it exactly.

The modality stays reference-free, and the core is unchanged: it reads verdicts as it does
today.

## Alternatives considered

- **Fixed thresholds (today).** These cause recurring false checks and miss quiet minors.
- **A static table learned on history.** It is fail-open under drift, and it cannot follow
  reagent or instrument changes.
- **A position-keyed level after placement.** It needs a new modality-after-alignment
  boundary and overfits like a 7-mer, for about two points of false flags over the 5-mer.
- **Calendar months as groups.** A month is a measurement stand-in, not a production unit; the
  run is the unit that shares chemistry and instrument.

## Validation plan

- **Measure on the development panel with the evaluation harness's review workload layer.**
  Report positions to review, samples with an unreviewed difference, and clear samples, plus
  strict and `haplotype_lod` better/worse, with samples grouped by run.
- **Measure the heteroplasmy work explicitly.** Report truth heteroplasmy recall and samples
  with a false `possible_heteroplasmy`, overall and by period. Report the flags the run signal
  removes and the flags it adds, separately, against today's rule, as a curve over the floor
  with and without a run signal. A run signal adds no flag below the floor; the floor sets
  the recall gain at quiet keys. The false-flag rate at quiet keys on the hold-out is
  the number to watch.
- **Show the effect on `undecided`.** Report items at recurring sites (`16031`, `527`, `320`,
  `109`) before and after.
- **Fix every parameter before judging it.** Fix the estimator, margin, floor and minimum on
  development runs, then judge on a hold-out panel and on runs later than the development
  panel's (readiness).
- **Test fail-closed behaviour** with:
  - a missing file;
  - a run below the minimum;
  - a key below the minimum;
  - a run signal from another run (refused);
  - a matching run id without the sample's hashes (flagged).
- **Prove leave-one-out exact:** a sample judged with the run's file equals the same sample
  judged with a file built without it.
- **Report the per-run spread** of the learned thresholds at the most recurrent keys.
- **Prove determinism:** byte-identical outputs for the same inputs and run signal, whatever
  the number of workers.

## Rollout and rollback

The run signal is an optional input. The `review` plugin version rises in increment (ii) and
the `sanger` version in increment (iii). Moving to one scale and one threshold changes review
labels and evidence with or without `--signal`; calls and verdicts change only in increment
(iii). With the file, its identity is in provenance. Rollback of the run signal is to stop
passing the file, which leaves the floor alone.

## Observability and operations

`dna signal` logs the run's read and call counts, the keys with and without enough calls,
and its output identity. Sample commands log:

- which run signal they used;
- whether leave-one-out applied;
- how many `possible_heteroplasmy` reads a run level judged and how many the floor alone judged;
- any refusal reason.

**Production flow.** The laboratory's run (plate) grouping lives in the review tool
(`mtdna_raw`). Its pipeline gains one step per run: build the run signal from the run's traces,
then call each sample with it. This is a handoff note for that repository, which DNA does not
change.

## Security and privacy

The run signal holds aggregate signal statistics and trace hashes, with no sample identity
beyond the run's. It is validated like any input file.

## Compatibility and migration

This proposal adds a contract (`dna.signal/v1`), an optional input, and an optional
`&RunSignal` argument to `analyze_sanger` (a public API change). The one-threshold amendment
changes review output and configuration (`second_floor`, schema version 8) whether or not
the input is given. The review tool adds a per-run signal step to use it.

## Decision

Pending.

## Implementation

- Increment (i): the run signal (`dna_sanger::run_signal`, `dna signal`, `dna.signal/v1`,
  SRS-CALL-019, SRS-API-016), with no verdict change; no command reads a run signal yet.
- Increment (ii): the canonical threshold (`run_signal::threshold`, SRS-CALL-020,
  `second_floor`), `dna review --signal` (SRS-API-014), and `possible_heteroplasmy` and the
  review evidence on the peak scale (`review` version 5). The evaluation harness builds run
  signals per batch, falling back to the month, in process (`dna-eval run --run-signal`).
