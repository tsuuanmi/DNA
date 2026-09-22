# Primary-sequence Difference Requirements

**Requirement namespace:** `SRS-VAR-*`

These requirements are part of the canonical [DNA SRS](README.md).

- **SRS-VAR-001:** Variants MUST derive only from the selected primary-sequence alignment.
- **SRS-VAR-002:** Canonical primary A/C/G/T mismatches MAY produce normalized SNV observations; unresolved primary N differences MUST be excluded and counted. A mixed call may still produce a primary-sequence difference from its strongest channel, but mixed supporting signal MUST remain explicit and MUST NOT be presented as an ordinary clean SNV.
- **SRS-VAR-003:** Contiguous gaps MAY produce insertions/deletions no longer than configured changed length, excluding the anchor.
- **SRS-VAR-004:** Indel representation MUST preserve the canonical gap placement selected by alignment under `SRS-ALN-012` in the [alignment requirements](alignment.md); variant construction MUST NOT independently left-shift, rotate, or otherwise move a repeat-equivalent indel. When no aligned left flank exists, the actual reference predecessor MUST be derived when possible; a true leading linear insertion/deletion MUST use a real right-anchor representation. Circular origin-spanning representation MUST preserve the alignment-selected side of the rCRS seam.
- **SRS-VAR-005:** Reported variants MUST be normalized and include only 1-based `position`, reference/alternate alleles, kind, and a direct `calls` array. Contig, classification, and normalization labels MUST NOT be duplicated in compact JSON.
- **SRS-VAR-006:** Every public variant-associated call MUST preserve its role, reference-oriented called `base`, co-located reference-oriented A/C/G/T primary-event channel heights, and uncalibrated `quality`. Original call index, ABIF PLOC, mapped call position, trace-strand primary/ambiguity, peak positions/sources, and vendor evidence MUST remain internal. Deletions MUST contain real flanks only and MUST NOT fabricate deleted-base signal/quality.
- **SRS-VAR-007:** Variant construction MUST preserve observed call/reference mappings and the alignment-selected canonical indel topology. Reported allele anchoring MAY add the required adjacent reference base but MUST NOT relocate the event through a repeat or across the rCRS origin seam.
- **SRS-VAR-008:** Variant alleles MUST use the reference strand. Public call bases and A/C/G/T peak labels MUST be projected to the reference strand, including reverse alignments. Only the four co-located primary-event heights are emitted; selected peak positions/sources and broader per-call peak objects remain internal.
- **SRS-VAR-009:** Genotype, zygosity, homoplasmy, heteroplasmy fraction, phase, PHFinder, pathogenicity, and clinical significance are prohibited in DNA output.
- **SRS-VAR-010:** Every reported variant's normalized 1-based anchor position MUST lie in at least one configured inclusive biological region.
- **SRS-VAR-011:** Every SNV supporting call and every inserted-base supporting call MUST meet the configured maximum-channel peak floor and strictly exceed the configured relative-quality threshold. Insertion flanks and deletion flanks MUST NOT be used for this supporting-evidence gate.
- **SRS-VAR-012:** A normalized SNV whose supporting call retains more than one co-localized qualifying channel under the active basecalling method MUST remain in observed/sample evidence but MUST be ineligible for clean SNV reporting with exclusion reason `mixed_supporting_dna`. This rule MUST NOT be generalized to insertions/deletions without a separately specified mixed-length method.
