---
implements: [SRS-CALL, SRS-CONS, SRS-VAR]
decided-by: [ADR-0071, ADR-0073, ADR-0074, ADR-0080, ADR-0083, ADR-0085]
covers:
  - crates/dna-consensus/src/
  - crates/dna-consensus/src/model/consensus.rs
  - src/model/consensus_result.rs
  - src/pipeline/consensus.rs
  - src/report/consensus.rs
---

# Sample Consensus Method

Part of the canonical [[design/pipeline|DNA pipeline]]. Requirements:
[[requirements/consensus|SRS-CONS]]. Decisions:
[[ADR-0073]],
[[ADR-0074]].

`dna consensus` runs the `sample` read path (read processing, placement,
rejected reads), aggregates the admitted reads to fix their order, and then
calls `dna_consensus::build`. The module reads only `CalledRead`
records, so any modality that produces `ReadEvidence` can feed it.

## Step 1 — Observations

`consensus/observe.rs` walks each read's alignment columns in reference order.

| Column | Observation | Clean when |
|---|---|---|
| call column | the base on its position | the call is unmasked, trusted (`read_end`), canonical, and raises no support veto |
| deletion column | a deleted position | the nearest calls on both sides are clean |
| inserted calls before a position | the inserted bases at the junction after the previous position | every inserted call and both neighbouring calls are clean; an insertion into a repeat the read bounds is judged for support vetoes on the repeat's placements, as in variant calling (SRS-VAR-020) |

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
`Eligibility::unbounded_run`. Such a read shows the run's base but not its
length. It seeds no stretch through that indel, and position by position it
still contributes the bases its alignment places on the run.

**Reference confirmations** (ADR-0080, SRS-CONS-010). A call its modality
marks `reference_only` is withheld for its base, but its base observation is
clean when the base is the reference's: it confirms the reference at its
position. The modality's own edge verdict does not stand against it, but
the read's unvouched ends do (`Eligibility::clear_of_ends`). Within that
margin of an unvouched end, a call the modality vouches for its base
(scored, off its edges, withheld for nothing) confirms the reference too,
when it shows no second base standing out and no unexplained second peak:
the margin keeps indels and run ends off a read end's calls, and a
confirmation shows neither. Its trust level
stays `Unmasked`, so it counts and bounds runs as before. Only the
position-by-position vote outside stretches changes; a stretch is decided as
if the confirmation were not there. The confirmation is withdrawn:
- where a read of the confirming call's strand has `second_unexplained` at
  the position, or a read of the other strand has it with no second base
  standing out, or with one that a read of the confirming strand shows
  standing out too: an unexplained second on one strand alone is not the
  two-strand signature of a second sequence, and the confirming strand reads
  its base cleanly apart from an explained artefact;
- where a read of the other strand shows the confirming call's own second
  base standing out (`SecondBase.stands_out`);
- at or beside a position where any read observes something other than the
  reference, cleanly or not (a read shifted in a repeat shows the reference
  beside its own difference);
- where no read observes the position cleanly otherwise and it forms a
  reference run of at least `consensus.seeded_run` bases with a neighbour
  across a junction no read observes cleanly (a length may change there);
- within or beside a decided site that differs from the reference (a site
  written partly as `N` included), after which the positions are decided
  again until no confirmation is withdrawn.

So a confirmation never decides a position a second sequence may explain, and
never changes what a variant's neighbours were. The positions only
confirmations observe cleanly are published as `reference_only`. Calling the
decisions treats the call there as showing its base only, so it bounds no
neighbouring run: a change that needed a run end there stays `run_boundary`.

**Poor read starts** ([[ADR-0085]], SRS-CONS-012). A position that some
read's call covers inside its unreliable span (the modality's poor start,
SRS-CALL-024), and that no read shows soundly outside its own span, is
published as `unreliable_only`; a deleted position counts when the calls on
both sides of it do. It decides nothing here. Calling the decisions walks each
reportable region inwards from its first and last position over these
positions and leaves out those up to the last one not settled as the
reference: they are not analysed, not callable, and a change there is
excluded.

**Dissent set aside** (ADR-0083, SRS-CONS-011). A site decided as the
reference keeps the reads that oppose it, and the review lists it as a check
where they are on the other strand. Two kinds of dissent are artefacts that a
read of the other strand contradicts, and are set aside. In both, a
supporting read on a strand no opposing read is on must observe every
position of the site cleanly, without a second base standing out, and far
from its own ends:
- every opposing read shows a second base standing out there, and its
  modality explains every second peak there (no `second_unexplained`): the
  signal's own shape, not a second sequence. The supporting read lies at
  least `variant_calling.read_end_margin` calls from its ends.
- every opposing read differs there only by insertions or deletions within
  `consensus.dissent_edge` calls (40) of its own ends, where a read's indels
  are mostly its own artefacts. The supporting read lies at least that far
  from its ends. A substitution is never set aside this way: near a read end
  it can be a run mix-up between reads.

The first kind also applies at a site or stretch called by a strict majority
as a value other than the reference, to an opposing read that differs only
in a run's length: every call it shows there is resolved, and it shows the
decided runs in the same order. Such a read is typically a late call that
counts one base too many in a run and shows the run's tail beside it (the
G run 66-71 in HV1F, whose own `73G` agrees). The second kind stays at
reference decisions, and a substitution against a called variant is never
set aside.

A set-aside read shows nothing at the site, not even a run's length inside a
stretch, and the positions are decided again with the confirmations' loop.
The decided value never changes: a reference decision loses only reads
against the reference, a contested site becomes called and a called one
loses its dissent, and a called variant keeps the majority that called it.

A read's decoded downstream sequence after a homopolymer tract adds
reference confirmations where its own calls are masked (SRS-CONS-010): see
[[design/length-mixture]] "Downstream confirmations".

## Step 2 — Stretches

**Seeds.**
- Every position where a read cleanly differs from the reference.
- Every junction where a read has a clean inserted base; such a junction
  seeds both of its positions.
- A clean substituted base that equals a neighbouring reference base also
  seeds that neighbour. A clean deletion that brings two runs of one base
  together seeds both of them. In both cases the reads' runs merge with the
  neighbouring run.
- A reference homopolymer of at least four bases that no read shows cleanly
  together with both neighbouring positions, while some read differs on it or
  its neighbours (a substitution, a deletion, or an insertion, clean or not).
  Decided position by position, such a run would take the reference length
  unlabelled; as a stretch, its length gets a label (Step 3). Shorter runs
  are left out: their disputes are mostly mixed signal at one base.

Each seed grows to the whole reference run containing it. Overlapping or
touching intervals merge into stretches, and a stretch never crosses the
origin.

## Step 3 — Deciding a stretch by its run structure

This step follows [[ADR-0074]].

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
  edge. Any other read, covering the stretch or only one of its edges, that
  reads that edge cleanly opposes the composition when the runs of its clean
  bases there are not the composition's first or last runs, with an
  observation no other read shares. A forward read losing phase in the
  HVS-I C-tract can read the tract through 16194 as C, anchored only by its
  dephased tail, while the reverse read shows the exit `…CA` in phase from
  the other side; that read now contests the tract instead of being left
  out.
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
- A clean start or end that holds every run of the composition while the read
  reads further runs over the stretch maps nothing: its outer run would stand
  for a run that starts or ends elsewhere. A read reading `C7 T C T C` over
  `CCCCCTCCCC` has the clean end `C T C`, whose inserted `C` would otherwise
  give the first C run a length of one.
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
  - `confirmed_end`: the count two reads agree on, each entering the run
    `InPhase` from an opposite end, every call of the run unmasked, and
    losing phase (`PhaseLoss`) at the other end on a call that does not read
    the run's base. A same-base call there means the run goes on, so its
    count is only a lower bound. Every read entering this way must agree.
  - `reference_frame`: the reference positions between its ends. Every read
    that shows an end must place it alike, and every position between must be
    read cleanly by some read.
  - `phase_loss`: otherwise, the longest run read in phase by a read that
    shows one end and loses phase at the other.
- **Two strands** (SRS-CONS-013). After every stretch is decided, where
  its call keeps the reference's run composition, each pair of runs of one
  base around a single-base anchor is a window taken from the stretch
  itself. A run of it still `anchored_end`, `phase_loss`, or
  `reference_frame` takes the length the reads' tract decodes put on top on
  both strands, as `two_strand`, when `length::decide` decides it
  ([[design/length-mixture]]).
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
2. A run of up to `consensus.maximum_undecided_gap` positions with no clean
   observation stays inside a segment as `N`. A longer run of such positions,
   or an uncovered position, ends the segment, and undecided positions at a
   segment's end are dropped.
3. A segment with fewer than `consensus.minimum_resolved_bases` resolved
   (non-`N`) bases is dropped, because the core could not place it (both 20,
   the core's `alignment.minimum_callable_bases`).

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
