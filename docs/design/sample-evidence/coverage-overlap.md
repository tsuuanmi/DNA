# Coverage and Pairwise Overlap

Part of the canonical [sample evidence aggregation method](README.md).

Before pairwise/locus aggregation, DNA derives a run-length reference coverage topology from every selected post-trim read segment. Each maximal interval records total read depth plus forward/reverse orientation depth. This counts all independently placed reads regardless of later pairwise eligibility and does not imply nucleotide agreement or consensus admission.

DNA then builds a deterministic pairwise overlap graph from the SHA-sorted read registry. Every unordered pair is compared only at shared selected-alignment reference coordinates. `shared_positions` counts all shared coordinates, while the agreement denominator includes only positions where both query observations are canonical A/C/G/T. Equal canonical observations are agreements; unequal canonical observations are conflicts. Unresolved symbols and deletions do not enter that nucleotide denominator, so gap/indel evidence remains separate. A pair is eligible for later consensus reconciliation only when the comparable-base count reaches `sample_reconciliation.minimum_comparable_bases` and the agreement fraction reaches `sample_reconciliation.minimum_overlap_agreement`. Non-overlapping reads produce no pair edge and remain valid sample evidence.
