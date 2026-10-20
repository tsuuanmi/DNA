---
implements: [SRS-SAMPLE, SRS-VAR]
decided-by: [ADR-0071, ADR-0073, ADR-0074, ADR-0076, ADR-0080]
verified-by: [tests/consensus.rs]
---

# Sample consensus requirements

**Requirement namespace:** `SRS-CONS-*`

Design: [[design/consensus|consensus method]]. Contract:
[[reference/consensus|consensus result]]. Decision:
[[ADR-0073]], [[ADR-0074]].

- **SRS-CONS-001:** `dna consensus <sample-id> <trace>... --reference <fasta>` MUST process the traces as `dna sample` does, with rejected reads (SRS-SAMPLE-027, SRS-SAMPLE-029, SRS-SAMPLE-030). It MUST write `results/<sample-id>.consensus.json` (`dna.consensus/v1`) and `results/<sample-id>.consensus.fasta`. It MUST fail without output when every read is rejected or when either target exists.
- **SRS-CONS-002:** The consensus MUST be modality-neutral. It MUST read admitted reads only as `CalledRead` records (`ReadEvidence`, alignment) and MUST NOT use filenames, primer labels, or profile windows.
- **SRS-CONS-003:** Only clean observations MUST decide. A clean observation's calls are unmasked, trusted under `variant_calling.read_end_margin` or exempt from it as the matching or inserted calls, or the reference-matching bounds, of a bounded run before a post-repeat phase loss (ADR-0076, SRS-VAR-013), and carry no modality reason against that use, read by scope as variant calling reads it: an observed base against its base, an inserted base against an indel, and a call beside an insertion or deletion against neither (the inserted calls of an insertion into a repeat the read bounds are judged on the repeat's placements, SRS-VAR-020); a substituted call, or a bound that reads another base, is never exempt. Indels longer than `variant_calling.max_indel_length` are never clean, and neither is a deletion whose neighbouring calls are not trusted. An insertion or deletion that changes the length of a run whose ends the read does not resolve (ADR-0071) MUST NOT be clean. A call its modality vouches for by its base only (scored, its bound withheld) MUST show only its base: it MUST NOT count towards a run's length, show a run's end, make an insertion clean, or bound a clean deletion, and its base MUST NOT be clean where SRS-VAR-018 finds that it may sit one position off, or where it lies on a reference homopolymer of at least four bases whose two bounds the read does not show with calls trusted for more than their base: such a position would carry a run length no read shows.
- **SRS-CONS-004:** Every place where a read cleanly differs from the reference MUST be grown to whole reference runs, together with each neighbouring reference run that a clean substituted base or a clean deletion joins. A reference homopolymer of at least four bases MUST also be a stretch when no read shows it cleanly together with both neighbouring positions while some read differs from the reference on it or its neighbours, cleanly or not; a shorter run MUST NOT, since its disputes are mostly mixed signal. The length is `consensus.seeded_run` (4); it matches the conformance minimum run length of the human-mtDNA profile, which the consensus does not read. Overlapping or touching places MUST merge into one stretch, and a stretch MUST NOT cross the reference origin. A stretch MUST be decided by its run structure (ADR-0074):
  - its composition (the ordered run bases) from reads that read every run cleanly and lose phase at neither edge, where any other read that reads an edge of the stretch cleanly, covering the whole stretch or not, MUST oppose a composition whose first or last runs its clean bases there contradict (a dephased read can swallow a run's exit that an in-phase read of only that end shows); otherwise from the single composition joining one read's clean start to another read's clean end where both reach a common reference position; otherwise the reference's, when a read shows its run structure; otherwise position by position;
  - its run lengths first as whole sequences, from reads that show every run's length, a tie leaving the stretch to the position-by-position decision; otherwise run by run, from reads that show both ends of a run with unmasked run calls;
  - position by position from reads without a clean indel inside it.

  A run end is shown only by a call beyond the run that reads another base: unmasked (`in_phase`), or masked but anchoring, possibly beyond one unresolved call (`anchored_end`). An anchored end placed elsewhere than where reads show that end in phase MUST NOT count. A partly clean run at the inner edge of a read's clean start or end MUST contribute only its base, and a clean start or end that holds every run of the composition while the read reads further runs over the stretch MUST contribute nothing, since it places a run end of the composition inside the stretch.
- **SRS-CONS-005:** Every decision MUST follow rule version 1:
  - all clean observations agree: that sequence;
  - a strict majority wins;
  - one read against one read: the reference sequence, `contested`;
  - any other tie: no sequence, `contested`;
  - no clean observation: `unresolved`.

  A change to the rule MUST raise the `consensus` plugin's method version.
- **SRS-CONS-006:** Segments MUST be maximal runs of decided positions. Undecided positions MUST be written as `N`. A gap of more than `consensus.maximum_undecided_gap` undecided positions MUST split a segment, and a segment with fewer than `consensus.minimum_resolved_bases` resolved bases MUST be dropped.
- **SRS-CONS-007:** The document MUST list:
  - the admitted and rejected reads;
  - every segment, with its reference interval and FASTA record name;
  - every site whose call differs from the reference, that was not called outright, that had disagreeing clean observations, or that holds a run of length evidence `phase_loss` or `reference_frame`, with its supporting, opposing, and uninformative reads and, for a stretch decided by its run structure, each run's length and length evidence (SRS-CONS-008);
  - decision counts, including the runs whose length is `phase_loss` or `reference_frame`.

  It MUST NOT carry a genotype, heteroplasmy fraction, or sample verdict beyond the consensus sequence (SRS-VAR-009).
- **SRS-CONS-008:** A run whose length no read shows MUST take, in order:
  - `confirmed_end`: the count of two reads that agree on it, each entering the run in phase from an opposite end, reading every call of the run unmasked, and losing phase at the other end on a call that does not read the run's base; every such read MUST agree;
  - `reference_frame`: the reference positions between its ends, when every read that shows an end places it alike and every position between is read cleanly by some read;
  - `phase_loss`: otherwise, the longest run read in phase by a read that shows one end and loses phase at the other (its next call masked, unresolved, or unclean and not reading another base);
  - otherwise it stays undecided.

  Each run's length evidence (`in_phase`, `anchored_end`, `confirmed_end`, `phase_loss`, `reference_frame`) MUST be published. The rule MUST NOT use target positions, regions, or conventions.
- **SRS-CONS-009:** A decided site that changes the reference length MUST stay undecided while a position beside it is undecided or uncovered, unless it was decided by its run structure with every run's length shown by a read (`in_phase` or `anchored_end`), since those reads show the run's bounds: otherwise the bases beside an insertion or deletion place it, and inside a repeat they decide which unit changes. Undoing one site can leave another one's neighbour undecided, so the rule MUST apply until nothing changes.
- **SRS-CONS-010:** A reference-only call (ADR-0080) whose base is the reference's at its position, scored and lying beyond `variant_calling.read_end_margin` of every unvouched read end (its modality's own edge verdict does not stand against it), MUST make that position's base observation clean as a confirmation of the reference, and nothing else: its run role MUST stay its call's, and no stretch, by its run structure or position by position, MUST count it. A confirmation MUST be withdrawn where any read's call at the position shows a second peak its modality does not explain, at or beside a position where any read observes something other than the reference, cleanly or not (another base, a deletion, or an insertion at a junction beside it), where the position forms a reference run of at least `consensus.seeded_run` bases with a neighbour across a junction no read observes cleanly and no read observes the position cleanly otherwise, where a read of the other strand shows the confirming call's own second base standing out there, and within or beside a decided site that differs from the reference (a site written partly as `N` included), after which the positions MUST be decided again until no confirmation is withdrawn. The positions that only confirmations observe cleanly MUST be published as `reference_only`, and calling the decisions (SRS-VAR-021) MUST treat the call at such a position as showing its base only, so that it bounds no neighbouring run.

## Implemented in

<!-- docs:begin backlinks -->
- **Implemented in:** [design/consensus](../design/consensus.md), [reference/configuration](../reference/configuration.md), [requirements/output](output.md), [requirements/variants](variants.md)
<!-- docs:end backlinks -->
