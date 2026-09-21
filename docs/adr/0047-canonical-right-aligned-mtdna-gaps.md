# ADR-0047: Canonical right-aligned mtDNA gap placement

- **Status:** Accepted
- **Date:** 2026-09-19

## Context

Affine-gap alignment can produce multiple equally scoring tracebacks that represent the
same biological sequence difference but place an insertion/deletion at different coordinates
inside a homopolymer or tandem repeat. Treating those paths as distinct placements makes
alignment topology depend on arbitrary traceback order and can move downstream sample-locus
states or normalized indel coordinates even though the observed haplotype is unchanged.

This ambiguity is especially important for human mitochondrial DNA. Forensic mtDNA
alignment/nomenclature practice has traditionally placed indels at the 3' end of
homopolymeric and repeat regions relative to the light strand of the rCRS. The Wilson
rules and later forensic guidance use this convention to prevent identical mtDNA sequences
from receiving different descriptions. ISFG guidance additionally recommends phylogenetic
alignment, and EMPOP applies region-specific conventions in several unstable mtDNA regions.

Generic VCF tooling commonly uses the opposite representational convention: left-aligning
and trimming indels. That convention is useful for VCF interoperability but must not silently
define DNA's internal mtDNA alignment topology.

ADR-0029 currently treats modulo-distinct equally scoring placements as an error. That is
too strict when the placements differ only by repeat-equivalent gap shifting and therefore
have one canonical mtDNA representation.

## Decision

DNA alignment MUST have one deterministic canonical traceback for repeat-equivalent
optimal indel placements.

### 1. Score optimality comes first

Right alignment is a canonicalization rule, not an alternate scoring objective.

DNA MUST first retain only tracebacks with the maximum scientific alignment score under
the configured profile-aware Gotoh model. Canonicalization MUST NOT choose a lower-scoring
path merely because its gap lies farther to the right.

### 2. Canonicalization applies only to repeat-equivalent optimal paths

Two optimal paths may be canonicalized together only when their difference is a shift of the
same insertion/deletion event through sequence-equivalent homopolymer or tandem-repeat
context and does not change the reconstructed query haplotype.

Canonicalization MUST NOT silently collapse genuinely different substitution/indel
explanations merely because their numeric scores tie.

### 3. Prefer contiguous indel structure

When equivalent paths differ only by splitting versus combining the same gap content,
DNA SHOULD prefer the representation with fewer gap runs before positional
canonicalization. This follows established mtDNA alignment practice of keeping equivalent
indels contiguous where possible.

### 4. Right-align equivalent indels

After the equivalent indel structure is identified, each gap run MUST be shifted as far
3' as sequence equivalence permits with respect to the light strand of the rCRS.

For the ordinary linearized rCRS coordinate system this means the highest valid reference
coordinate: the **right-most** equivalent placement.

Examples:

~~~text
reference:  C A A A A G
query:      C A A A - G

equivalent gap placements inside the A run
    -> canonical gap is the furthest-right / 3'-most placement
~~~

For a tandem repeat, one whole repeat unit may be shifted only when the resulting ungapped
query/reference relationship is sequence-equivalent.

### 5. Deterministic ordering for multiple gap runs

If more than one repeat-equivalent gap run remains, canonicalization MUST be deterministic.
The ordered vector of gap anchors is compared from the 3'-most event toward the 5' end;
the candidate with the lexicographically greater right-shifted anchor vector wins after
gap-run count equivalence.

An implementation MAY use another mathematically equivalent total ordering if it produces
the same canonical 3'-most alignment.

### 6. Circular rCRS has a canonical seam

Human mtDNA is circular biologically, but alignment/reporting coordinates use the rCRS
origin. Canonical right shifting MUST treat the rCRS 16569|1 origin as a fixed seam rather
than permitting an equivalent gap to rotate indefinitely across the doubled working
reference.

Origin-spanning alignment may still use the existing doubled-reference algorithm. Gap
canonicalization is applied in the projected rCRS coordinate system without shifting an
indel across the canonical origin solely to obtain a numerically larger coordinate.

Region-specific phylogenetic exceptions, if required for full EMPOP-compatible mtDNA
nomenclature, need a separate explicit method/ADR rather than hidden special cases.

ADR-0057 defines the complementary comparison/reporting boundary: different event
representations that reconstruct the same resolved haplotype are biologically equivalent
for validation, while DNA's deterministic right-aligned internal representation remains
unchanged.

### 7. Canonical alignment and variant representation are separate contracts

The selected alignment columns MUST use the canonical right-aligned topology.

Variant extraction/representation MUST NOT move an indel back to a contradictory placement without explicitly declaring a separate output-representation convention. DNA's production variant builder therefore preserves the alignment-selected canonical gap anchor and performs no independent repeat shifting or circular rotation.

VCF-style left normalization may be offered only at an explicit serialization/interchange
boundary in the future; it is not the scientific source of truth for DNA's internal
mtDNA alignment.

### 8. Primary sequence and profile evidence remain distinct

This decision does not change ADR-0029 substitution scoring. `EvidenceProfile` continues
to determine profile-aware substitution scores, while the retained primary sequence remains
the source for callable identity and primary-difference extraction.

Canonical right alignment resolves equivalent gap topology only after optimal-score
placement has been established.

## Consequences

- Repeat-equivalent optimal alignments become stable instead of failing or depending on
  traceback ordering.
- Homopolymer/repeat indels use the mtDNA forensic 3'/right-most convention internally.
- Sample locus states and call-to-reference mappings become deterministic across equivalent
  gap placements.
- Alignment canonicalization no longer inherits VCF left-normalization by accident.
- Variant representation preserves canonical alignment placement; positional repeat canonicalization has one authoritative implementation in `alignment::canonical`.
- Phylogenetic/EMPOP special-region notation remains a separate future layer.

## Validation requirements

Implementation MUST add fixtures covering at least:

- one-base deletion at every equivalent position of a homopolymer;
- one-base insertion at every equivalent position of a homopolymer;
- multi-base tandem-repeat insertion/deletion;
- multiple equivalent gap runs where contiguity and right-most ordering are tested;
- a tie that is not repeat-equivalent and MUST remain ambiguous/error;
- forward and reverse trace orientations producing the same reference-oriented canonical
  alignment;
- circular origin-adjacent repeats where canonicalization MUST NOT rotate across the rCRS
  seam;
- the local m.16431-m.16435-style indel neighborhood as a validation pattern, without
  encoding private sample data in repository fixtures.

## References

- Wilson MR et al. Recommendations for consistent treatment of length variants in the
  human mitochondrial DNA control region. Forensic Sci Int. 2002.
  https://doi.org/10.1016/S0379-0738(02)00206-2
- Parson W et al. DNA Commission of the ISFG: revised and extended guidelines for
  mitochondrial DNA typing. Forensic Sci Int Genet. 2014.
  https://doi.org/10.1016/j.fsigen.2014.07.010
- EMPOP documentation, alignment/nomenclature conventions:
  https://empop-use.readthedocs.io/en/latest/using_empop/
- GATK LeftAlignAndTrimVariants documents the contrasting VCF left-alignment convention:
  https://gatk.broadinstitute.org/hc/en-us/articles/360042913511-LeftAlignAndTrimVariants

## Non-goals

This ADR does not define phylogenetic haplogroup-aware realignment, EMPOP special-region
exceptions, a VCF compatibility format, heteroplasmy interpretation, indel confidence, or
thresholds.
