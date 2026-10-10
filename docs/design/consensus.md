# Sample Consensus Method

Part of the canonical [DNA pipeline](pipeline.md). Requirements:
[SRS-CONS](../requirements/consensus.md). Decisions:
[ADR-0073](../decisions/adr/0073-sample-consensus.md),
[ADR-0074](../decisions/adr/0074-run-structure-consensus.md).

`dna consensus` runs the `sample` read path (read processing, placement,
rejected reads), aggregates the admitted reads to fix their order, and then
calls `dna_core::consensus::build`. The module reads only `CalledRead`
records, so any modality that produces `ReadEvidence` can feed it.

## Step 1 — Observations

`consensus/observe.rs` walks each read's alignment columns in reference order.

| Column | Observation | Clean when |
|---|---|---|
| call column | the base on its position | the call is unmasked, trusted (`read_end`), canonical, and raises no support veto |
| deletion column | a deleted position | the nearest calls on both sides are clean |
| inserted calls before a position | the inserted bases at the junction after the previous position | every inserted call and both neighbouring calls are clean |

Each observed base also records its trust level, `observe::Trust`, and its
reference position (none for an inserted base). Each level implies the ones
below it:

| Level | Meaning |
|---|---|
| `Clean` | the base itself counts |
| `Unmasked` | an in-phase call that counts towards a run's length |
| `Anchoring` | masked but anchoring; it still shows where a neighbouring run ends (ADR-0071) |
| `Nothing` | none of these |

Indels longer than `variant_calling.max_indel_length` are never clean.

An insertion or deletion that changes the length of a run is not clean either
when the read does not resolve both ends of that run. This is the
`run_boundary` test of variant calling (ADR-0071), shared through
`ReadEligibility::unbounded_run`. Such a read shows the run's base but not its
length. It seeds no stretch through that indel, and position by position it
still contributes the bases its alignment places on the run.

## Step 2 — Stretches

**Seeds.**
- Every position where a read cleanly differs from the reference.
- Every junction where a read has a clean inserted base; such a junction
  seeds both of its positions.
- A clean substituted base that equals a neighbouring reference base also
  seeds that neighbour. A clean deletion that brings two runs of one base
  together seeds both of them. In both cases the reads' runs merge with the
  neighbouring run.

Each seed grows to the whole reference run containing it. Overlapping or
touching intervals merge into stretches, and a stretch never crosses the
origin.

## Step 3 — Deciding a stretch by its run structure

This step follows [ADR-0074](../decisions/adr/0074-run-structure-consensus.md).

**A read's runs.** Each read that covers the stretch and both neighbouring
positions contributes its bases over the stretch, insertions included, and
the base on either side (the frame).
- **Unresolved calls (`N`) are dropped.** The bases next to a dropped call
  lose all trust, because the `N` may extend either run.
- **Exception: a lone `N` at a phase-loss exit.** This is a lone `N` between
  two different bases, exactly one of them unmasked. There the unmasked base
  keeps counting and the masked one still shows the run end, both marked
  `beside_unresolved`.
- **A frame that is a lone `N`.** It is replaced by the masked anchoring base
  beyond it, marked the same way.

The remaining bases are grouped into runs. Each run end is classified from
the call just beyond it:

| End | The call beyond the run |
|---|---|
| `InPhase` | unmasked, reads another base |
| `Anchored` | reads another base and is masked but anchoring, or is `beside_unresolved` |
| `PhaseLoss` | masked, unresolved, or unclean, and reads no other base |
| `Open` | missing, or a clean call of the run's base |

Each read also has a clean start (`prefix`) and a clean end (`suffix`). These
are the whole clean runs from each edge, plus the next run if its first call
is clean. A clean start counts only when the read does not lose phase at the
stretch's start, and a clean end likewise.

**1. Composition** (the ordered run bases), from the first that applies:
- A vote among reads that read every run cleanly and lose phase at neither
  edge.
- The single composition that joins one read's clean start to another read's
  clean end. The two must reach a common reference position, and the start's
  last runs must match the end's first runs in exactly one way.
- The reference composition, when a read without an `N` shows it.
- Otherwise, Step 3c.

**Mapping reads onto the composition.**
- A read with the composition's bases maps run by run.
- Otherwise, only its clean start and end map.
- A partly clean run at their inner edge shows its base only: it neither
  counts nor shows its inner end.
- An `Anchored` end placed elsewhere than where reads show that end
  `InPhase` becomes `PhaseLoss`.

**2. Lengths.**
- **Whole sequences.** Reads that show every run's length vote on the whole
  sequence: unmasked run calls, and no end `Open` or `PhaseLoss`. A tie goes
  to Step 3c.
- **Run by run.** Otherwise each run is voted from reads that show both of its
  ends. Its label is `in_phase` when a supporting read shows both ends
  `InPhase`, and `anchored_end` otherwise.
- **A run whose length no read shows** takes one of these:
  - `reference_frame`: the reference positions between its ends. Every read
    that shows an end must place it alike, and every position between must be
    read cleanly by some read.
  - `phase_loss`: otherwise, the longest run read in phase by a read that
    shows one end and loses phase at the other.
- **Undecided runs.**
  - At the reference composition, an undecided run is written as `N` over its
    reference length.
  - Otherwise, the stretch goes to Step 3c.

The composition decision is merged into the result. Reads that disagree on it
oppose the stretch, and a contested composition contests the stretch.

**3c. Position by position.** This uses reads without a clean indel inside
the stretch.

**Combining run or position votes.**
- An undecided run or position is written as `N` over its reference length.
- The stretch is `contested` if any part is.
- It is `unresolved` only if every part is.

## Step 4 — Decision rule, version 1 (`decide.rs`)

| Clean observations | Result | State |
|---|---|---|
| all equal | that sequence | `called` |
| a strict majority | the majority | `called` |
| one against one, one is the reference | the reference | `contested` |
| any other tie | none (`N`) | `contested` |
| none | none (`N`) | `unresolved` |

Positions outside stretches are decided with the same rule, one by one.

## Step 5 — Segments

1. Decided positions are walked in order.
2. A run of up to `alignment.minimum_callable_bases` positions with no clean
   observation stays inside a segment as `N`. A longer run of such positions,
   or an uncovered position, ends the segment, and undecided positions at a
   segment's end are dropped.
3. A segment with fewer resolved (non-`N`) bases than that minimum is
   dropped, because the core could not place it.

Each segment becomes one FASTA record, named
`<sample-id>_<first>-<last>` in 1-based reference positions.

## Step 6 — Variants

The FASTA is the input of `dna call`. It reads each segment as a reviewed
sequence (vouched ends; `N` is unresolved) and calls it with the core, so
`dna call` remains the only variant path.

## Interpretation boundary

- The consensus sequence is DNA's adjudication between reads. The evidence
  behind every non-trivial decision is listed in `sites[]`.
- It reports no heteroplasmy. A site that reads disagree on is `contested` or
  `N`, never an IUPAC code (SRS-VAR-009).
