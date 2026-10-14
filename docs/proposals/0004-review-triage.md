---
id: PROP-0004
type: proposal
status: draft
owners: []
created: 2026-10-14
related-requirements: [SRS-CALL-005, SRS-CALL-009, SRS-SAMPLE-029, SRS-CONS-001]
related-decisions: [ADR-0018, ADR-0019, ADR-0065, ADR-0067, ADR-0068, ADR-0073, ADR-0074]
implementation: []
---

# Proposal: Review triage

## Problem

DNA's [purpose](../requirements/SRS.md#purpose) is to sort sequencing work
into three outcomes: pass, review, and rerun. It sorts at four levels:
region, trace, sample, and batch, worked on from the bottom up. Today DNA
reports the evidence for each decision, but never the decision.

- **Trace evidence.** Callability segments, the callable span, rejected reads
  with a typed cause, and integrity and SNR annotations.
- **Region evidence.** Masked calls, eligibility reasons, opposition, and the
  consensus marks `N`, `contested`, `phase_loss`, and `reference_frame`.
- **Sample evidence.** Coverage, overlaps, and the variants of the consensus.

So a reviewer still opens every sample. Nothing says which samples need no
look, which regions need one, and which traces should go back to the lab.

## Goals

- Build a catalogue of the errors that occur in the data. Give each error
  class one outcome:
  - errors DNA can resolve from the evidence → `pass`;
  - errors DNA cannot resolve, but a reviewer can settle from the traces →
    `review`;
  - errors nobody can settle because the trace is too poor or unreliable →
    `rerun`.
- Give every sample, trace, and region exactly one outcome: `pass`, `review`,
  or `rerun`. Attach typed reasons and locations to every outcome.
- Pass a result only when its evidence meets an explicit, versioned policy.
- Make rerun decisions from the trace and its coverage. Name the cause, so the
  lab can choose the action.
- Measure A % (pass), B % (review), C % (rerun), and the error rate of passed
  results against ground truth. Use these numbers to rank further work.

## Non-goals

- Setting the lab's rerun rate. The chemistry sets it; DNA only detects it.
- Comparing results with reviewers or other tools. Downstream pipelines own
  that (ADR-0065).
- Batch, plate, and control decisions across samples. Batch is the top
  level and comes last; downstream pipelines aggregate it from DNA's
  per-trace outcomes.
- Sample swaps. The data cannot show that a sample was swapped.
- Checking that a sample has all its expected traces. Samples normally do,
  so DNA needs no primer design for this.
- Clinical or forensic interpretation of a passed result.

## Evidence and research

Current DNA measurements on held-out reviewer truth (see
[phase-recovery](../research/phase-recovery/README.md) and
[point-heteroplasmy](../research/point-heteroplasmy/README.md)):

- **Most remaining errors are in known places.** Calls that are not masked are
  about 99.98 % right. Remaining errors cluster in dephased poly-C windows
  (HVS-I 16184–16193, HVS-II 303–315) and at sites with a single callable
  read. These places can already be told apart, so a cautious pass rule can
  route them to review.
- **Point heteroplasmy is common enough to matter.** Reviewers wrote IUPAC
  codes at about 70 sites in 69 of 640 samples. DNA's consensus never writes
  them. A pass rule must therefore check for minor signal, or a heteroplasmic
  sample passes as homoplasmic.
- **Consensus exactness bounds A.** On the 320-sample evaluation subset, the
  consensus matches reviewers exactly on 297 samples. A passed sample must be
  exact, so A is at most that share until the consensus improves.

External practice:

- **Both strands.** Forensic mtDNA guidelines that follow SWGDAM ask for
  concordant data from both strands across HVS-I and HVS-II. They accept
  single-strand data only when it is confirmed.
- **Fit the rerun to the failure.** Core facilities match the rerun to the
  cause:
  - off-scale signal → re-inject;
  - dye blobs → clean up again;
  - failed or weak reaction → fix the template and sequence again;
  - mixed template from the first base → PCR again, or extract again.

  Repeating the same reaction rarely fixes a failure.

## Brainstorm and review outcome

The draft was reviewed on 2026-10-14. The points below record what was
agreed, what was dropped, and what is still open.

### The error catalogue comes first

Not every error can be fixed, and a rerun does not fix every error. So the
first step is a catalogue of error classes, each with an outcome:

- **Resolved by DNA → pass.** For example, a call that one read gets wrong but
  clean reads from both strands oppose.
- **Settled by review.** Biological states belong here, because running them
  again gives the same picture: length heteroplasmy in the C-tracts and point
  heteroplasmy. So do technical limits that a reviewer can settle from the
  traces, such as a run length that no read shows.
- **Rerun.** A trace too poor or unreliable to settle anything: a failed or
  weak reaction, a read with too few callable calls, a read that cannot be
  placed.

A `review` reason therefore carries a class: technical, biological, or
convention. A `rerun` names its cause, and the lab chooses the action
(re-inject, clean up again, sequence again, PCR again). Core facilities match
the action to the cause, and repeating the same reaction rarely helps.

### A bad trace does not always need a rerun

A rejected read whose regions other reads already cover cleanly needs no
rerun. DNA keeps two things apart:

- `trace quality` — a diagnostic for lab QC;
- `rerun required` — the sample cannot pass or be reviewed without a better
  trace.

### Mixed samples

DNA cannot see a sample swap, so swaps are out of scope. It can see a mixed
sample when the traces of one sample disagree too much: reads that disagree
outside the repeat windows, or a read mixed from its first base. A mixed
sample goes to review or rerun with its own reason.

### Levels, bottom up

The levels form a matrix, not a tree: traces cut across regions. The work goes
from the lowest level up:

1. **Region.** Decided from all clean observations of the region: pass,
   review with reasons and positions, or rerun naming the traces that would
   close the gap. Regions are the profile's reportable regions (HVS-I, -II,
   -III), with flagged intervals and sites inside them.
2. **Trace.** Decided from one trace's signal and placement, and from whether
   its regions need it.
3. **Sample.** Pass only if every region passes and no sample-level reason
   such as mixing fires. Otherwise review or rerun, carrying the reasons up.
4. **Batch.** Many failing traces in one run, plate, or capillary point at the
   run. Downstream pipelines aggregate this. DNA can help by carrying the
   ABIF run and well identity in its trace result.

Reviewers work at site level, so B is counted as flagged sites per sample as
well as flagged samples.

### "Not called" is not "reference"

A region passes only when every target position is read cleanly. Rejected
reads and masked calls can leave a position with no clean reading. A silent
gap must never pass as "no variant".

### What "pass" must mean

A region normally needs clean, concordant data from both strands to pass, as
in forensic mtDNA practice. Exceptions are still to be discussed. A starting
policy, to be validated and never weakened silently:

- every target position in the region is in the consensus, with no `N` and
  no `contested` site;
- clean, concordant observations from both strands;
- no masked, `phase_loss`, or `reference_frame` run in the region;
- no minor-channel signal above the heteroplasmy threshold after shadow
  correction;
- no opposition from a clean observation, and no conformance finding.

Rules that cost passes but catch no errors can be relaxed later, with
evidence, one class at a time. That is how B goes down.

### Error rate of passed results

This is the main safety metric; A, B, and C come second. It is measured
against ground truth.

- **Sample size.** If no errors appear in *n* passed samples, the 95 % upper
  bound on the error rate is about 3/*n*. Showing below 1 % takes at least
  300 passed, error-free samples; below 0.1 % takes about 3,000. The current
  evaluation sets can support the first claim, not the second.
- **Imperfect truth.** Reviewer truth is not perfect. Where DNA and reviewers
  disagree, a third read is needed before the case counts against a pass.
- **Reason codes.** Every flag has a reason code. B can then be split by
  reason and worked on by count ([ROI ranking](roadmap.md)). Reviewer
  decisions on flagged regions, collected downstream, become the next
  evaluation truth.

## Proposed design

Sketch only; an ADR fixes it.

- **Outcome type.** An outcome is `pass`, `review`, or `rerun`, with a list of
  typed reasons. Each reason carries:
  - a level and a location;
  - a class (technical, biological, convention);
  - for a rerun, the cause.
- **Ownership.**
  - **Sanger plugin:** trace-intrinsic reasons that are about the chemistry.
  - **Core:** region and sample outcomes over the consensus and `ReadEvidence`,
    with no modality types.
  - **Policy:** the target profile and configuration, versioned with the
    result.
- **Contract.** Two options are open:
  - a `triage` section in `dna.consensus`;
  - a separate `dna.triage/v1` document from a `triage` command over the
    sample's results.
- **Determinism.** The same inputs and policy always give the same outcome.

## Alternatives considered

- **Leave triage to downstream pipelines.** Rejected as the main path. The
  pass decision needs DNA's evidence types (masks, run structure, profiles),
  and downstream rules would copy them. Downstream still owns batch rules,
  controls, and comparison.
- **One score with a threshold.** Rejected. A score hides the reason, cannot
  name a rerun cause, and is not calibrated (SRS-QC-002).
- **Sample-level outcomes only.** Rejected. Reviewers and the lab act on
  regions and traces.

## Validation plan

- **Phase 0, before any code.**
  - Build the error catalogue from the ground-truth disagreements and the
    rejected reads of the evaluation and fresh subsets. Give each class an
    outcome and a count.
  - A downstream research script then applies the draft pass policy to the
    existing `dna consensus` and `dna sample` results of the same subsets. It
    reports A, B, and C with flagged sites per sample, the errors that leak
    into pass, and the truth sites that a pass would hide.

  This decides whether the policy is worth building.
- **Acceptance.** Zero known errors in passed results on the ground-truth
  held-out sets.
  Report the upper bound with the set size, and keep the result stable on a
  disjoint set.
- **Rerun decisions.** Compare with lab rerun records, where available, for
  the precision and recall of `rerun`.

## Rollout and rollback

To be decided after phase 0. Outcomes are new output and do not change any
call. Rollback means not emitting them.

## Observability and operations

Outcome counts and reason counts per run go to the operational log. Batch
summaries belong to downstream pipelines.

## Security and privacy

No new trust boundary. ABIF run and well identity are lab metadata and follow
the [data policy](../governance/data.md).

## Compatibility and migration

New contract fields or a new document version; no migration is needed.

## Open questions

1. Which exceptions, if any, allow a region with single-strand data to pass?
   Both strands are the default.
2. Should a pass ever cover length heteroplasmy recorded by convention (the
   HVS-I C-tract), or must that always go to review?
3. How far must the traces of one sample disagree before the sample counts as
   mixed?
4. What error rate of passed results is acceptable for this lab, and how many
   ground-truth samples are needed to show it?

## Decision

Not yet reviewed.

## Implementation

None.
