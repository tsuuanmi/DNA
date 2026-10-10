# Sample Consensus Method

Part of the canonical [DNA pipeline](pipeline.md). Requirements:
[SRS-CONS](../requirements/consensus.md). Decision:
[ADR-0073](../decisions/adr/0073-sample-consensus.md).

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

Each observed base also records whether it *resolves* a neighbouring run's
end, the ADR-0071 notion. That means it is canonical, informative, and
unmasked or anchoring.

Indels longer than `variant_calling.max_indel_length` are never clean.

## Step 2 — Stretches

Every position where a read cleanly differs from the reference is a seed. So
is every junction where a read has a clean inserted base; such a junction
seeds both of its positions.

Each seed grows to the whole reference run containing it. Overlapping or
touching intervals merge, and an interval never crosses the origin. The
result is a set of stretches. Inside a stretch, reads that describe one
haplotype through different alignments still produce the same sequence: for
example `310C` plus an inserted `T`, against an insertion at the run end.

## Step 3 — Deciding a stretch

**Read's view.** Each read that covers the stretch and both neighbouring
positions contributes:
- its bases over the stretch, insertions included;
- the bases on either side.

Unresolved calls (`N`) are dropped. The bases next to a dropped call are
neither clean nor resolving, because the `N` may extend either neighbour.

**Order of preference.**
1. **Whole.** At least one read observes the stretch cleanly (every base
   clean, no `N`, both outer bases resolving) with a run structure other than
   the reference's. Whole-stretch sequences are then voted.
2. **Run by run.** Otherwise each reference run is voted separately, using
   reads whose run structure equals the reference's. A read's vote for a run
   is that run's length, if every base in the run is clean and the bases just
   before and after it resolve.
3. **Position by position.** This applies when no read has the reference
   structure, or when run-by-run decides nothing. It uses reads without a
   clean indel inside the stretch.

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
