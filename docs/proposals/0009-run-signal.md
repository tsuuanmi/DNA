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
- One rule for the Sanger mixed verdict and for `possible_heteroplasmy`.
- Fail closed:
  - without a usable run signal, today's rules apply;
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
   - **Two minimums.** Below a minimum of reads for the whole file, the run signal is refused
     with a logged reason and today's rules apply. A key with fewer than `minimum_calls`
     observations falls back to today's rule for that key.
3. **One estimate, two exact uses.** For a call's key, let `t = max(floor, quantile + margin)`,
   computed from the run histogram minus the sample.
   - **The Sanger mixed verdict (ADR-0027)** holds when today's qualifying-channel rule holds
     **and** the second channel's share exceeds `t`. A run signal can therefore only lift a
     mixed verdict whose shadow the run shows as common; it never adds one. Identity is withheld
     no more often than today, so `undecided` can only fall: at recurring sites such as `16031`
     and `527`.
   - **The review's `possible_heteroplasmy`** fires when a clean call's share exceeds `t` (in
     place of today's fixed `0.15`), with the neighbour-shadow and strand rules unchanged. It
     can flag a real minor that sits above a quiet level but below today's ratio, and it
     drops recurring shadows.
   - **Such a minor keeps its identity**: the consensus decides the primary base, and the
     review makes the position a check. Reporting it as IUPAC remains a separate, deferred
     decision.
   - **Parameters.** `quantile`, `margin`, `floor`, `minimum_calls` and the minimum reads sit
     in the existing `[signal_processing]` configuration section. They are fixed by research on
     development runs and judged on a hold-out panel and on later runs.
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
  removes and the flags it adds, separately. It adds a flag where `t` is below `0.15` at a quiet
  key, which is the intended recall gain. The false-flag rate at quiet keys on the hold-out is
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

The run signal is an optional input. The `sanger` and `review` plugin versions rise with the
change, so every output's provenance changes. Without `--signal`, nothing else changes:
calls, verdicts and labels stay byte-identical apart from those versions. With the file, its
identity is in provenance. Rollback is to stop passing the file.

## Observability and operations

`dna signal` logs the run's read and call counts, the keys with and without enough calls,
and its output identity. Sample commands log:

- which run signal they used;
- whether leave-one-out applied;
- how many verdicts fell back to today's rule;
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
`&RunSignal` argument to `analyze_sanger` (a public API change). Without the input, outputs
change only in plugin versions. The review tool adds a per-run signal step to use it.

## Decision

Pending.

## Implementation

- Increment (i): the run signal (`dna_sanger::run_signal`, `dna signal`, `dna.signal/v1`,
  SRS-CALL-019, SRS-API-016), with no verdict change; no command reads a run signal yet.
