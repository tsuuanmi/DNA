# Sample consensus requirements

**Requirement namespace:** `SRS-CONS-*`

Design: [consensus method](../design/consensus.md). Contract:
[consensus result](../reference/consensus.md). Decision:
[ADR-0073](../decisions/adr/0073-sample-consensus.md).

- **SRS-CONS-001:** `dna consensus <sample-id> <trace>... --reference <fasta>` MUST process the traces as `dna sample` does, with rejected reads (SRS-SAMPLE-027, SRS-SAMPLE-029). It MUST write `results/<sample-id>.consensus.json` (`dna.consensus/v1`) and `results/<sample-id>.consensus.fasta`. It MUST fail without output when every read is rejected or when either target exists.
- **SRS-CONS-002:** The consensus MUST be modality-neutral. It MUST read admitted reads only as `CalledRead` records (`ReadEvidence`, alignment) and MUST NOT use filenames, primer labels, or profile windows.
- **SRS-CONS-003:** Only clean observations MUST decide. A clean observation's calls are unmasked, trusted under `variant_calling.read_end_margin`, and raise no support veto. Indels longer than `variant_calling.max_indel_length` are never clean, and neither is a deletion whose neighbouring calls are not clean.
- **SRS-CONS-004:** Every place where a read cleanly differs from the reference MUST be grown to whole reference runs. Overlapping or touching places MUST merge into one stretch, and a stretch MUST NOT cross the reference origin. A stretch MUST be decided:
  - as a whole, when a read cleanly observes all of it with another run structure than the reference's;
  - otherwise run by run, from reads with the reference run structure whose run calls are clean and whose bases on both sides of the run are canonical, informative, and unmasked or anchoring, with no unresolved call next to them;
  - otherwise position by position, from reads without a clean indel inside it.
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
  - every site whose call differs from the reference, that was not called outright, or that had disagreeing clean observations, with its supporting, opposing, and uninformative reads;
  - decision counts.

  It MUST NOT carry a genotype, heteroplasmy fraction, or sample verdict beyond the consensus sequence (SRS-VAR-009).
