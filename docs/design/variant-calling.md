# Variant Calling Method

Part of the canonical [DNA pipeline](pipeline.md).

Extracts evidence-backed primary-sequence differences from the selected
alignment and constructs validated anchored REF/ALT alleles. Only differences
in the primary sequence are considered; no allele-frequency, genotype,
heteroplasmy, post-calling haplotype canonicalization, or target nomenclature
inference is performed.

### Substep 6.1 — Difference extraction

Walking the alignment columns:

- a column where the query is `-` is a **deletion** of the reference bases;
- a column where the reference is `-` is an **insertion** of the query bases;
- a column with unequal canonical query and reference bases is an **SNV**.

Differences whose allele contains a non-canonical base, or whose indel length
exceeds `max_indel_length`, increment the excluded-candidate warning count rather
than being reported.

### Substep 6.2 — Allele construction and anchoring

The caller converts each selected-alignment difference into a validated
reference-oriented allele representation.

For indels, the current production implementation preserves the canonical gap
placement selected by alignment under `SRS-ALN-012`; it does not independently
left-shift, rotate, or right-shift the event through a repeat. When an aligned
left flank is unavailable, the caller derives the real reference predecessor
where possible. A true leading linear indel uses a real right anchor. Circular
origin-spanning representation preserves the alignment-selected side of the
rCRS seam.

Internally each variant retains its contig, 1-based position, reference/alternate
alleles, kind, and direct call mappings. Compact v7 emits only `position`,
`reference`, `alternate`, `kind`, and `calls`. Every public call contains
only its supporting/flanking `role`, reference-oriented called `base`,
co-located reference-oriented A/C/G/T primary-event channel heights in `peaks`,
and uncalibrated `quality`. Original call index, PLOC, mapped call position,
trace-strand symbols, selected-peak positions/sources, penalties, and vendor
evidence remain internal. Deletions carry real aligned flanks only and never
fabricate deleted-base dna. A flanking call whose strongest channels tie or
carry no signal is an unresolved `N` with no primary event, so it is omitted
from the public calls instead of receiving fabricated heights. The emitted reference allele is validated against
the supplied reference.

The implemented post-calling `variant_normalization` capability defined by
ADR-0060 may move sequence-equivalent events only after reconstructing and
preserving the complete represented haplotype. That capability is intentionally
separate from this caller.

### Substep 6.3 — Configured eligibility

A called candidate is retained only when its 1-based anchor `position` lies
inside at least one configured inclusive region. SNV supporting calls and every
inserted-base supporting call must each have a highest A/C/G/T peak greater than
or equal to `minimum_peak_height` and an uncalibrated relative score strictly
greater than `relative_quality_threshold`. Insertion flanks are not evaluated.
Deletions have no supporting trace base, so their flanks are not subjected to
peak or quality thresholds; their caller anchor must still be in a region.
Vendor PCON is not used by this filter. For SNVs, a supporting call with more than one co-localized qualifying channel is retained as an observed called difference but is ineligible for clean-SNV reporting with `mixed_supporting_dna`. Insertions and deletions are not subjected to this point-mixed-signal gate; persistent mixed-length evidence is a separate method boundary.

Each removed candidate increments `excluded_variant_candidates` once, even when
it fails more than one eligibility condition. The pure variant stage also returns
a concise exclusion diagnostic containing kind, contig, caller position when
available, and all failed rules. Pipeline orchestration writes one WARN record per
diagnostic without reference/alternate alleles. Sample aggregation logs aggregate
counts of differential-locus observations and variant-associated calls that
retain a basecall-independent profile; those operational counts do not alter the
scientific result.

### Substep 6.4 — Ordering

Reported variants are sorted by `(contig, position, reference, alternate)` and
deduplicated.
