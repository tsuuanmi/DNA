# Variant Calling Method

Part of the canonical [DNA pipeline](pipeline.md).

Extracts evidence-backed primary-sequence differences from the selected
alignment and constructs validated anchored REF/ALT alleles. Only differences
in the primary sequence are considered; no allele-frequency, genotype,
heteroplasmy, post-calling haplotype canonicalization, or target nomenclature
inference is performed.

### Substep 7.1 — Difference extraction

Walking the alignment columns:

- a column where the query is `-` is a **deletion** of the reference bases;
- a column where the reference is `-` is an **insertion** of the query bases;
- a column with unequal canonical query and reference bases is an **SNV**.

A column on a masked call ([callability](callability.md)) is no SNV candidate
and produces no record. An insertion with a masked inserted call is excluded
with the mask reasons of those calls. Other differences whose allele contains a
non-canonical base, or whose indel length exceeds `max_indel_length`, increment
the excluded-candidate warning count rather than being reported.

### Substep 7.2 — Allele construction and anchoring

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
alleles, kind, and direct call mappings. Compact v8 emits only `position`,
`reference`, `alternate`, `kind`, and `calls`. Every public call contains
only its supporting/flanking `role`, reference-oriented called `base`,
co-located reference-oriented A/C/G/T primary-event channel heights in `peaks`,
and uncalibrated `quality`. Original call index, PLOC, mapped call position,
trace-strand symbols, selected-peak positions/sources, penalties, and vendor
evidence remain internal. Deletions carry real aligned flanks only and never
fabricate deleted-base dna. A flanking call whose primary base is unresolved
(`N`: its strongest channels tie, carry no signal, or all qualify) is omitted
from the public calls, which admit only A/C/G/T bases. The emitted reference allele is validated against
the supplied reference.

The implemented post-calling `variant_normalization` capability defined by
ADR-0060 may move sequence-equivalent events only after reconstructing and
preserving the complete represented haplotype. That capability is intentionally
separate from this caller.

### Substep 7.3 — Configured eligibility

A called candidate is retained only when its 1-based anchor `position` lies
inside at least one configured inclusive region. SNV supporting calls and every
inserted-base supporting call must each have a highest A/C/G/T peak greater than
or equal to `sanger_evidence.minimum_peak_height` and an uncalibrated relative
score strictly greater than `sanger_evidence.relative_quality_threshold`. Insertion flanks are not evaluated.
Deletions have no supporting trace base, so their flanks are not subjected to
peak or quality thresholds; their caller anchor must still be in a region.
Read eligibility (ADR-0062, ADR-0067, SRS-VAR-013) then marks a variant
ineligible when any of its mapped calls lies within `read_end_margin` calls of
an uninformative call (`read_end`). A call is informative when it lies inside
the trim interval and is unmasked or dephased; calls beyond the trim interval
and masked calls aligned as unresolved are not. A difference therefore needs
that many informative calls between it and the nearest place where the
alignment has no information, at the read's ends and around internal masked
segments alike. A masked
evidence call — a supporting call, or a flanking call of a deletion — adds the
reason of its segment: `post_homopolymer` when the segment starts in the window
after a repeat run, otherwise `dephased_signal`, `mixed_signal`, `weak_signal`,
or `irregular_spacing`. Trace order is the sequencing direction for both
strands, so forward and reverse reads are handled alike without primer
knowledge.

Vendor PCON is not used by this filter. For SNVs, a supporting call with more than one co-localized qualifying channel is retained as an observed called difference but is ineligible for clean-SNV reporting with `mixed_supporting_dna`. Insertions and deletions are not subjected to this point-mixed-signal gate; persistent mixed-length evidence is a separate method boundary.

The caller reads each read only as `ReadEvidence`
([ADR-0069](../decisions/adr/0069-plugin-first-modality-core-post-calling.md)).
The region gate, `indel_length_exceeded`, `non_canonical_allele`, and
`read_end` are core rules. Everything else arrives from the modality as data:
- the Sanger adapter (`read_processing::evidence`) raises the support vetoes
  `peak_below_minimum` and `relative_quality_not_above_threshold`, which apply
  to SNVs and insertions, and `mixed_supporting_dna`, which applies to SNVs
  only;
- it labels masked calls with their segment reason.

A variant's reasons are reported in this order:
1. `outside_target_region`;
2. the union of its supporting calls' vetoes, in the adapter's vocabulary
   order (the order listed above);
3. `read_end`;
4. the mask reasons of its evidence calls, in mapping order and without
   duplicates.

Each removed candidate increments `excluded_variant_candidates` once, even when
it fails more than one eligibility condition. The pure variant stage also returns
a concise exclusion diagnostic containing kind, contig, caller position when
available, and all failed rules. Pipeline orchestration writes one WARN record per
diagnostic without reference/alternate alleles. Sample aggregation logs aggregate
counts of differential-locus observations and variant-associated calls that
retain a basecall-independent profile; those operational counts do not alter the
scientific result.

### Substep 7.4 — Ordering

Reported variants are sorted by `(contig, position, reference, alternate)` and
deduplicated.
