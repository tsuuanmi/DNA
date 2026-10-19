---
id: PROP-0006
status: implemented
created: 2026-10-17
summary: "`call` reports the edits the consensus decided instead of aligning its FASTA again, which loses where `N`s sit (placement, run bounds, AC deletions), and every read places an equivalent indel at one canonical position before the consensus votes, with the laboratory's conventions applied after it."
implements: [SRS-API-012, SRS-VAR-015, SRS-VAR-017, SRS-ALN-015]
decided-by: [ADR-0047, ADR-0073, ADR-0074, ADR-0077]
---

# Proposal: Call the consensus's decisions without re-aligning them

## Problem

The recommended workflow is `consensus` then `call`. `SampleAnalysis::call`
runs it in process (SRS-API-012): it renders the consensus FASTA, and `call`
aligns each record to the reference again and extracts variants from that
new alignment. The consensus already placed every base on the reference and
decided every edit (`segments[].reference`, `sites[]`). Re-deriving them from
the FASTA loses what the FASTA cannot carry: where each base sits and why a
base is `N`. An `N` scores alike against any base and observes nothing, so
around undecided bases the re-alignment is free to choose differently from
the consensus.

The [[validation/discrepancy-catalogue|discrepancy catalogue]] counts
three classes with this one cause on development panel v5:

| Class | Samples | Current handling |
|---|---:|---|
| A record placed away from the interval its consensus segment covers | 24 | guarded: every variant excluded (`placement_conflict`, SRS-VAR-017) |
| A run length decided `in_phase` whose bound the FASTA writes as `N` (`315.1C`, a `309` insertion) | 8 discordant (30 affected) | open: `call` excludes the length change (`run_boundary`) |
| An AC-unit deletion decided at `523–524` placed beside `N`s at `513–516` | 31 | fixed by a tie rule for vouched reads (SRS-ALN-017, SRS-VAR-018) |

Each fix so far patches one symptom of the same gap. The guard and the tie
rule make `call` safer; they cannot make it say what the consensus decided.

The same gap exists before the consensus. Reads vote on their alignment
columns, and an equivalent edit placed differently by two reads splits the
vote. On development panel v5 a gained `AC` unit is placed at `513.1C 513.2A`
by one read and `524.1A 524.2C` by another: the consensus sees two
one-against-one ties, keeps the reference at both, and loses an insertion that
every read shows (116 contested AC sites, 112 with the truth `524.1A 524.2C`).

## Goals

- `SampleAnalysis::call` and the recommended workflow report exactly the
  edits the consensus decided, at the positions it decided them, without
  placing a record on the reference again.
- One authoritative path from consensus decisions to called variants and
  notation.
- `dna call` on a FASTA sequence from elsewhere (a reviewer-edited consensus,
  another tool) keeps working, with its current guards.

## Non-goals

- Changing how the consensus decides (ADR-0073, ADR-0074).
- Heteroplasmy, IUPAC output, or run-length estimates; a `phase_loss` or
  `reference_frame` length stays labelled as it is.
- Notation rules; the decided variants go through `dna-normalize` as now
  (ADR-0077 §5).
- Trusting weak calls: making a weak call that reads the reference count
  towards a decision was measured and set aside, because such a call also
  decides run lengths (development panel v5: 21 regressions against 11 fixes
  in the narrowest variant).

## Evidence and research

- Counts above: development panel v5 at `7931e9e`, discrepancy classes from
  `dna-eval run` (`class`, `class_regions`).
- Worked case (AC repeat): the consensus decides site `522–525` `CACC→CC`
  `in_phase`, leaves `513`, `514`, and `516` unresolved, and writes
  `…CCCA NN A N ACACACC G…`. With the `N`s neutral, every deletion placement
  from `513` to `523` scores the same.
- Run bounds: [[validation/known-limitations|known limitations]],
  "Run lengths lost between consensus and call".
- AC placement: profile-aware alignment scoring sets an insertion in a tandem
  repeat where the read's own call purities are lowest, and the canonical
  3'-most shift of ADR-0047 is refused because it changes that score. The
  post-calling `normalization` plugin would place it at `524.1A 524.2C`, but it
  runs after the vote. Reads' placements on panel v5: reverse reads mostly
  after 513, forward reads after 524, 513, 523, or 522; deletions agree.

## Proposed design

0. **Order: one canonical form before the vote, the lab's notation after.**
   - Each read's alignment places every indel at its single sequence-equivalent
     3'-most position before any read is observed: a shift is taken when every
     moved column stays a match and the read's bases are unchanged, whatever
     the profile score. The score still chooses which edits explain the read;
     it no longer chooses among equivalent placements (an amendment of
     ADR-0047 §1). One implementation of the equivalent shift serves alignment
     and post-calling normalization.
   - The consensus votes on those canonical observations. It also grows a
     stretch over whole tandem-repeat units, not only homopolymer runs, so a
     placement that still differs cannot split a vote.
   - Conventions of the laboratory stay after the consensus, in
     `dna-normalize`: the deletion trade, the nomenclature windows, and
     convention labels. They decide how a decided sequence is written, never
     what it is (ADR-0077 §5). The deletion trade crosses substitutions, so it
     is not an equivalent shift and never runs before the vote.

1. **Decided edits come from the consensus.** For each consensus segment, the
   called variants come from its `called` sites whose call differs from the
   reference, in reference coordinates, so no record is placed again:
   - a site decided by its run structure (`runs`, ADR-0074) gives one edit
     per run, as ADR-0072 writes it: a length change as an insertion or
     deletion at the run's 3' end, a changed base as a substitution;
   - any other site gives the edits of an alignment of its call against its
     reference interval alone, placed by the same canonical shift as in item
     0. That alignment is bounded to the site and holds no `N`, since an
     undecided run is undecided whole.

   A `contested` site between one read and one other contributes the
   reference, as now. A contested site without a sequence, and an
   `unresolved` one, write `N`: they contribute nothing and are not
   callable. The callable reference segments are the segment's interval less
   those positions (SRS-ALN-015).
   The consensus applied the run-bound rule (ADR-0071) when it decided a
   length, so the document path does not apply `run_boundary` again: an
   unresolved bound only removes that position from the callable bases.
   Expected outcome for the run-bound class: the decided `315.1C` or `309`
   insertion is reported, and `310` is not callable.

   A length labelled `phase_loss` or `reference_frame` is an estimate or an
   assumption, not an observation (ADR-0074). It is reported as today, where
   it already reaches `call` through the FASTA, and its label stays visible in
   `sites[].runs`; triage keeps routing it to review. The document path does
   not make it a firmer call than the FASTA path does.
2. **Notation as now.** The variants go through normalization and the
   nomenclature windows with those callable segments as the readable bases.
3. **Entry points.**
   - `SampleAnalysis::call` builds them from the analysis's consensus in
     memory.
   - `dna call` accepts the consensus document (`dna.consensus/v1`) as input,
     beside a FASTA sequence; given the document it takes this path.
   - A FASTA sequence keeps today's path (alignment, SRS-VAR-017, SRS-ALN-017,
     SRS-VAR-018) because nothing else knows its placement.
4. **Contract.** SRS-API-012 changes from "exactly as `call` reads the
   consensus FASTA" to "exactly the consensus's decided edits". The
   `dna.variants/v1` reads of a document input are its segments, each with
   its interval and callable segments and the decided variants, all
   eligible. Whether that is a revision of `dna.variants/v1` in place or a
   new version is decided with the contract change.

An ADR records the decision once this proposal is accepted.

## Alternatives considered

- **Keep re-aligning and guard each symptom.** The current state. Every new
  case needs its own rule, and the run-bound class has no guard that recovers
  the decided length.
- **Carry placement in the FASTA** (an edit string or coordinates in the
  record header). The FASTA would stop being a plain sequence other tools
  read, and `call` would trust a header over the sequence.
- **Write a base instead of `N` at a decided bound.** The consensus would
  assert a base no clean read shows, against the readability rule
  (ADR-0077 §5).
- **Constrain the re-alignment to the segment's interval.** It removes the
  placement class only; inside the interval the `N`s still leave the edits
  undecided.

## Validation plan

- Item 0, first, on development panel v5 against main: no new
  `haplotype_lod` discordance and at least as many fixes as regressions; the
  contested AC sites, the samples with the truth `524.1A 524.2C`, the
  insertions whose eligibility changes because other calls become the
  inserted ones, and the HVS-II insertions at `309`–`315`.

- Equivalence: for every development-panel consensus segment without an `N`,
  the document path and the FASTA path give the same variants and notation.
- Development panel v5: the three classes above resolve to the consensus's
  decisions; no sample regresses on any measure; the class counts move from
  `sequence` as expected.
- Unit and integration tests for each site state, run labels, insertions at
  segment edges, and a segment crossing the circular origin.

## Rollout and rollback

One change set: the document path, the `dna call` input, the contract and
requirement changes, and the tests. Rollback reverts it; the FASTA path is
unchanged throughout.

## Observability and operations

The operation log records which input path `call` took and the counts it
already reports.

## Security and privacy

No new data leaves the machine; the input is a document DNA wrote.

## Compatibility and migration

The current specification requires no compatibility: results of the
recommended workflow change where the re-alignment disagreed with the
consensus, which is the fix. Callers that pass a FASTA see no change.

## Decision

- 2026-10-17, user: item 0. The canonical placement in alignment (option A)
  is approved for implementation and measurement, and this order (canonical
  placement before the consensus, the laboratory's conventions after it) is
  the design direction.
- 2026-10-17, user: the contract change is approved. SRS-API-012 becomes "exactly the consensus's decided edits", and `dna call` accepts the consensus document (`dna.consensus/v1`) as an input beside a FASTA sequence. Whether the variants document of a consensus input is a revision in place or a new version is settled with the implementation, under the versioning policy.

## Implementation

ADR-0078 records the decision; SRS-VAR-021 and SRS-API-012 hold the rules.
Item 0's canonical placement landed earlier (ADR-0047 §1 amendment). Its
tandem-repeat stretch growth was measured and parked: it gains fewer strict
concordances than it costs.

Items 1 to 4: `dna_core::consensus::decided` (plugin `decided` 1),
`dna call <id>.consensus.json`, and `SampleAnalysis::call` on the same path.
`dna.variants/v1` is revised in place: the reads of a consensus input are its
segments.

Two choices differ from the text above, both from measurement:
- **The record rules stay.** Item 1 proposed not applying `run_boundary`
  again. Waiving it, even only for stretches whose every run length a read
  showed, reported lengths the consensus could not bound and measured worse.
  The decided variants are judged by the record rules. A length change of the
  run at a segment's first or last position is also excluded, since no read
  shows that run's outer end.
- **Ties in a site's own alignment** prefer a pair, then a deletion, then an
  insertion. The deletion trade of normalization settles a deletion beside a
  substitution (`16297C 16298DEL`); its gap was fixed separately (#103).

Measured on development panel v5 against main (`haplotype_lod`):
- 14 samples become concordant and none become discordant;
- primary: 13 better, 1 worse (a stretch decided position by position next
  to an `N`, `16027C`, which PROP-0008 step 3 covers);
- 8 samples that failed to call now call;
- no exported EMPOP range covers a position the consensus left undecided (16
  positions in 13 samples before, 0 now);
- callable yield is unchanged to the fourth decimal.

Equivalence: of the 7,144 consensus segments without an `N`, 7,135 give the
same notation as the FASTA path. In the other 9:
- the FASTA path misplaced a record that starts inside the HVS-II run
  (`301C 302C`);
- or the FASTA path's alignment chose a substitution where the consensus
  decided a deletion at the HVS-I read start (`16037G` against `16038DEL`,
  `16031C` against `16031DEL`);
- or a length at a segment-edge run that the document path now excludes.
