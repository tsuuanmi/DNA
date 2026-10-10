# ADR-0074: Run-structure consensus and how a run length is known

- **Status:** Accepted
- **Date:** 2026-10-10
- **Related decisions:** [ADR-0009](0009-biological-semantics.md),
  [ADR-0021](0021-scientific-core-confidence-floor.md),
  [ADR-0065](0065-result-comparison-downstream.md),
  [ADR-0067](0067-signal-derived-read-callability.md),
  [ADR-0071](0071-run-length-edits-need-resolved-run-ends.md),
  [ADR-0073](0073-sample-consensus.md)

## Context

The human-mtDNA poly-C tracts carry length heteroplasmy. Rare in the reference
form, it is common once a substitution joins two runs. Examples are 16189T>C
in HVS-I (`CCCCCTCCCC` at 16184–16193) and the 303–315 window in HVS-II.
A sample then holds molecules with different run lengths.

A Sanger read reads every molecule in phase up to and through the run, in its
own reading direction. It loses phase only after the run's exit. The data show
this directly:
- **Forward HVS-I reads** are callable through 16193 and dephased after it.
- **Reverse reads** are callable from their 3′ end down to 16184.

So each read shows the run's entry end and its bases, but not where the run
ends. At most, a masked call beyond the run reads the next base of the
dominant molecule. Usually the read loses phase at the exit: the call there is
masked, unresolved, vetoed as mixed, or still reads the run's base.

Consensus method version 1 (ADR-0073) decided a stretch from edits and
alignment columns. Three faults followed:
- **An unbounded indel decided a length.** A read's indel decided a run
  length even when the read did not bound the run. This added a false
  `16182.1C` in 6 of 640 held-out samples.
- **An extra base placed beside its run.** The aligner placed an extra tract
  base as an insertion in a neighbouring run, as an `AC` insertion at the A/C
  junction. The stretch was also allowed to end inside a read's run. Both cases
  were decided locally and named a run that the read did not show.
- **True lengths rested on these loopholes.** The true HVS-II dominant lengths
  (`309.1C`, `309.2C`) came only through the same loopholes, from reads that
  lose phase at the run they lengthen.

Reviewers follow a convention (forensic SOP, SWGDAM style):
- **HVS-I:** read `16189C` when the T is absent, and note A-tract changes as
  transversions (`16182C`, `16183C`). Length heteroplasmy is noted, but the C
  count is not recorded; the reference `C5TC4` frame is kept.
- **HVS-II:** record the dominant length.

The reviewed corpus agrees:
- of 41,754 samples with `16189C` or `16189Y`, 98.9 % carry no HVS-I length
  call;
- the HVS-II length calls (`309.1C`, `309.2C`, `315.1C`) are routine.

## Decision

The consensus decides a stretch by its **run structure**. It reports, for
every run, how its length is known. Target conventions stay out of the core.

### Stretches

A seed also takes in a neighbouring reference run when the reads' runs merge
with it:
- a clean substituted base that equals a neighbouring reference base, such as
  `16189C` joining both C runs;
- a clean deletion that brings two runs of one base together, such as `AC` in
  an AC repeat.

### Run ends

Each read's bases over the stretch are grouped into runs. Each run end is one
of these:

| End | The call just beyond the run |
|---|---|
| in phase | is unmasked and reads another base |
| anchored | reads another base but is masked yet anchoring; or reads another base, is masked, and comes after one unresolved call between it and the in-phase run |
| phase loss | is masked, unresolved, or unclean and reads no other base: the read loses phase at the run |
| open | is missing, or is a clean call of the run's base, so the run continues |

Two rules refine these ends:
- **In-phase evidence wins.** When an anchored end sits elsewhere than where
  reads place that end in phase, it counts as phase loss.
- **A partly clean run shows its base only.** A run at the inner edge of a
  read's clean start or end may hold further runs. Neither its length nor its
  inner end counts.

### Composition, then lengths

1. **Composition (the ordered run bases).** The first of these that applies
   decides it:
   - reads that read every run cleanly and lose phase at neither edge vote on
     it;
   - otherwise, the single composition that joins one read's clean start to
     another read's clean end, where both reach a common reference position
     and the start does not begin, nor the end finish, in phase loss;
   - otherwise, the reference's composition, when a read shows its run
     structure;
   - otherwise, the stretch is decided position by position.
2. **Whole sequences.** Reads that show every run's length (unmasked run
   calls, no end open or phase-lost, ADR-0071) first vote on the whole
   sequence. Run lengths read through different run boundaries never mix this
   way. A tie between them leaves the stretch to the position-by-position
   decision, which uses only reads without a clean indel there.
3. **Run by run.** Otherwise, each run's length is voted from reads that show
   both of its ends.
4. **A run whose length no read shows** takes the first of these:
   - **`reference_frame`:** the number of reference positions between its
     ends. Each end must be placed alike by every read that shows it, and
     every position between must be read cleanly by some read. This assumes
     no length change.
   - **`phase_loss`:** the longest run a read reads in phase before losing
     phase. Within a length mixture this lies near the dominant length, so it
     is an estimate.
   - Otherwise the run is undecided.

### Reporting

Every run of a stretch decided this way is published with its
`length_evidence`: `in_phase`, `anchored_end`, `phase_loss`, or
`reference_frame`. The summary counts `phase_loss` and `reference_frame` runs.
The FASTA carries only the sequence.

No position, target, or region enters the rule (ADR-0021). A profile
convention, such as HVS-I lengths not being recorded, belongs to post-calling
policy and is roadmap work. That policy must classify a call, never drop it.

### Versions

The `consensus` plugin method version is 2. It covers this decision and the
rule that an unbounded indel is not clean (SRS-CONS-003).

## Alternatives

- **Strict `N` for every run that no read shows both ends of.** Rejected. It
  lost about 30 true HVS-II dominant lengths. Exact samples fell to 289 and
  277, against 293 and 298 for the edit-level check alone.
- **The edit-level check alone** (an unbounded indel is not clean). Not
  enough on its own. It misses an extra base placed as a multi-base insertion
  in a neighbouring run, and a stretch that ends inside a read's run. Both
  still decide a run the read does not show.
- **Fill the stretch's remaining reference length into the unread run.**
  Rejected. In HVS-II the other runs often carry real length changes
  (`315.1C`), and the remainder turned true `309.1C 309.2C` into `309DEL`.
- **Prefer the in-phase count over the reference frame when it is longer.**
  Rejected. In HVS-I it records the dominant tract length, which the
  convention does not record, and it added false calls.
- **Drop anchored ends.** Rejected. ADR-0071 keeps them, and they carry true
  `309.1C` and `309.2C`.
- **Decompose the length mixture from the signal.** Deferred (ADR-0009). The
  trace holds a reference-free survival curve at the run exit, but its
  calibration and its use need their own decision
  ([phase-recovery research](../../research/phase-recovery/README.md)).

## Consequences

**Measured downstream** (ADR-0065), on the fair scale with IUPAC truth scored
separately. The figures are exact samples / precision / recall.

| Set | Union | Method version 1 | Method version 2 |
|---|---|---|---|
| 320 evaluation | 288 / 0.9924 / 0.9924 | 290 / 0.9967 / 0.9866 | **297 / 0.9970 / 0.9899** |
| disjoint 320 | 290 / 0.9934 / 0.9926 | 292 / 0.9964 / 0.9901 | **299 / 0.9967 / 0.9928** |

- **Reviewer consensus calls** (160 samples): 145 identical, against 142 for
  method version 1 and 140 for the union.
- **Local Sequencher set:** 28 exact samples, against 26 for method version 1.
- **One sample** got one more error. It is a rare variant that a single read
  misreads.
- **`N`:** the fraction of decided positions written as `N` is unchanged at
  1.1–1.3 %.

**Where the labels fall.**
- **`reference_frame`:** 145 of 151 runs are in the HVS-I tract and 1 is in
  HVS-II. This is where the convention keeps the reference frame.
- **`phase_loss`:** 28 of 37 runs are in the HVS-II window, where reviewers
  record the dominant length.

The split comes from the reads' signal, not from a regional rule.

**Remaining limits.**
- A dominant HVS-I length that a read does show (an anchored end) is still
  published as a length change, although reviewers do not record it. That
  awaits the profile policy.
- A run that no read reads stays `N`.
- A mixed call at the HVS-I A/C boundary stays `N`; point heteroplasmy is not
  reported (SRS-VAR-009).

**Contract.** `dna.consensus/v1` sites gain `runs`, and the summary gains
`phase_loss_runs` and `reference_frame_runs`. The contract is revised in place
before its first release (Breaking). `dna sample` and `dna call` are
unchanged.
