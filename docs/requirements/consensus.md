# Sample consensus requirements

**Requirement namespace:** `SRS-CONS-*`

Design: [consensus method](../design/consensus.md). Contract:
[consensus result](../reference/consensus.md). Decision:
[ADR-0073](../decisions/adr/0073-sample-consensus.md), [ADR-0074](../decisions/adr/0074-run-structure-consensus.md).

- **SRS-CONS-001:** `dna consensus <sample-id> <trace>... --reference <fasta>` MUST process the traces as `dna sample` does, with rejected reads (SRS-SAMPLE-027, SRS-SAMPLE-029). It MUST write `results/<sample-id>.consensus.json` (`dna.consensus/v1`) and `results/<sample-id>.consensus.fasta`. It MUST fail without output when every read is rejected or when either target exists.
- **SRS-CONS-002:** The consensus MUST be modality-neutral. It MUST read admitted reads only as `CalledRead` records (`ReadEvidence`, alignment) and MUST NOT use filenames, primer labels, or profile windows.
- **SRS-CONS-003:** Only clean observations MUST decide. A clean observation's calls are unmasked, trusted under `variant_calling.read_end_margin`, and raise no support veto. Indels longer than `variant_calling.max_indel_length` are never clean, and neither is a deletion whose neighbouring calls are not clean. An insertion or deletion that changes the length of a run whose ends the read does not resolve (ADR-0071) MUST NOT be clean.
- **SRS-CONS-004:** Every place where a read cleanly differs from the reference MUST be grown to whole reference runs, together with each neighbouring reference run that a clean substituted base or a clean deletion joins. Overlapping or touching places MUST merge into one stretch, and a stretch MUST NOT cross the reference origin. A stretch MUST be decided by its run structure (ADR-0074):
  - its composition (the ordered run bases) from reads that read every run cleanly and lose phase at neither edge; otherwise from the single composition joining one read's clean start to another read's clean end where both reach a common reference position; otherwise the reference's, when a read shows its run structure; otherwise position by position;
  - its run lengths first as whole sequences, from reads that show every run's length, a tie leaving the stretch to the position-by-position decision; otherwise run by run, from reads that show both ends of a run with unmasked run calls;
  - position by position from reads without a clean indel inside it.

  A run end is shown only by a call beyond the run that reads another base: unmasked (`in_phase`), or masked but anchoring, possibly beyond one unresolved call (`anchored_end`). An anchored end placed elsewhere than where reads show that end in phase MUST NOT count. A partly clean run at the inner edge of a read's clean start or end MUST contribute only its base.
- **SRS-CONS-005:** Every decision MUST follow rule version 1:
  - all clean observations agree: that sequence;
  - a strict majority wins;
  - one read against one read: the reference sequence, `contested`;
  - any other tie: no sequence, `contested`;
  - no clean observation: `unresolved`.

  A change to the rule MUST raise the `consensus` plugin's method version.
- **SRS-CONS-006:** Segments MUST be maximal runs of decided positions. Undecided positions MUST be written as `N`. A gap of more than `alignment.minimum_callable_bases` undecided positions MUST split a segment, and a segment with fewer resolved bases than that minimum MUST be dropped.
- **SRS-CONS-007:** The document MUST list:
  - the admitted and rejected reads;
  - every segment, with its reference interval and FASTA record name;
  - every site whose call differs from the reference, that was not called outright, or that had disagreeing clean observations, with its supporting, opposing, and uninformative reads and, for a stretch decided by its run structure, each run's length and length evidence (SRS-CONS-008);
  - decision counts, including the runs whose length is `phase_loss` or `reference_frame`.

  It MUST NOT carry a genotype, heteroplasmy fraction, or sample verdict beyond the consensus sequence (SRS-VAR-009).
- **SRS-CONS-008:** A run whose length no read shows MUST take, in order:
  - `reference_frame`: the reference positions between its ends, when every read that shows an end places it alike and every position between is read cleanly by some read;
  - `phase_loss`: otherwise, the longest run read in phase by a read that shows one end and loses phase at the other (its next call masked, unresolved, or unclean and not reading another base);
  - otherwise it stays undecided.

  Each run's length evidence (`in_phase`, `anchored_end`, `phase_loss`, `reference_frame`) MUST be published. The rule MUST NOT use target positions, regions, or conventions.
