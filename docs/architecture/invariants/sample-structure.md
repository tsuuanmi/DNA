# Sample Structure and Topology Invariants

Part of the canonical [read/sample invariant family](sample-boundaries.md).

- **INV-SAMPLE-001:** A future consensus sequence is a downstream projection, not the authoritative source of sample variants or discordance.
- **INV-SAMPLE-002:** Missing coverage is distinct from reference support. For each read, mapped reference segments define coverage; inside coverage, omission from sparse `locus_differences[]` means canonical reference match, while positions outside coverage remain uncovered.
- **INV-SAMPLE-003:** Sample evidence defines read identity, unique reviewer-facing filename stem, orientation, and post-trim coverage once in a SHA-sorted read registry. Public locus and normalized-variant evidence reference reads by that unique name; internal aggregation MAY use deterministic indexes but MUST NOT expose them as reviewer-facing identifiers.
- **INV-SAMPLE-004:** Duplicate trace content cannot be counted twice within one sample evidence result.
- **INV-SAMPLE-005:** A locus record exists only when at least one covering read is alternate, unresolved, or deleted. Every covering read at that retained locus remains explicit so reference support and its focused quality evidence are not lost.
- **INV-SAMPLE-006:** Normalized variant support preserves configured eligibility, exclusion reasons, and reviewer-facing reference-oriented base/peak/quality evidence; internal call mappings remain authoritative for scientific traceability, and read-level filtering MUST NOT erase a normalized observation.
- **INV-SAMPLE-007:** Pairwise overlap is discovered only after independent read placement from shared reference coordinates. Non-overlapping reads have no edge; missing an overlapping or canonical F/R partner does not invalidate a read.
- **INV-SAMPLE-008:** Pairwise nucleotide agreement uses only coordinates where both reads carry canonical A/C/G/T query bases. Unresolved symbols and deletions remain outside that denominator; gap/indel evidence is never converted into fabricated nucleotide agreement.
- **INV-SAMPLE-009:** Overlap eligibility is downstream evidence for future consensus and cannot rewrite read placement, read-level observations, or variant eligibility.
- **INV-SAMPLE-010:** Sample coverage topology derives only from selected mapped reference segments and counts all independently placed reads. It cannot inherit pairwise overlap eligibility as read rejection, and orientation depth is not equivalent to nucleotide agreement, consensus confidence, or biological strand independence.
- **INV-SAMPLE-011:** Normalized-variant support topology is a lossless summary of the existing per-read variant observations across eligibility and selected-orientation dimensions. It cannot add supporting reads, erase ineligible observations, count reference/unresolved/competing-event coverage as support for that variant, or become a confidence/independence verdict.
