# `src/variant_calling/normalize.rs`

## Purpose

Builds minimal SNV/indel allele representations while preserving the canonical gap placement already chosen by alignment and the mapped original-call evidence.

## Responsibilities

- Build SNV, insertion, and deletion `Variant` records with one-based positions.
- Preserve the observed canonical alignment anchor for linear and circular indels; do not shift repeat-equivalent events a second time.
- Derive a real predecessor anchor when one is available and use the existing right-anchor representation only for a true leading linear indel.
- Preserve circular origin-spanning anchors from canonical alignment without rotating the event around the rCRS seam.
- Validate call mappings and emitted reference alleles.

## Non-responsibilities

No extraction, gap canonicalization, repeat shifting, phylogenetic alignment, mapping construction, or serialization.

## Key types and functions

- `snv(reference, position, alternate, calls) -> Result<Variant>`.
- `insertion(reference, previous, next, inserted, calls) -> Result<Variant>`.
- `deletion(reference, previous, first_deleted, next, deleted, calls) -> Result<Variant>`.
- `observed_anchor`, `build_insertion`, `build_deletion`, `validated`, and `reference_base`: representation helpers.

## Invariants and errors

- Indel alleles must be non-empty; otherwise `Error::Variant`.
- Reference positions and anchors are bounds-checked; out-of-range access returns `Error::Variant`.
- Variant representation cannot move a canonical alignment gap through a homopolymer or tandem repeat.
- When no aligned left flank exists on a linear reference, the actual reference predecessor is derived when possible; a true linear-origin insertion/deletion uses a real right anchor.
- Circular alignment placement remains on the same side of the canonical rCRS seam selected by alignment.
- Emitted reference alleles are validated against the supplied reference; disagreement returns `Error::Variant`.
- Call mappings are validated by `mapping` and preserved unchanged.

## Dependencies

- `model::coordinate` for one-based reference positions.
- `model::reference` for `Reference` and `ReferenceTopology`.
- `model::variant` for `Variant`, `VariantCallMapping`, and `VariantKind`.
- `mapping` for call validation.
- `error` for typed failure.

## Biological semantics

Repeat/homopolymer positional canonicality belongs to alignment under ADR-0047. Variant representation therefore reflects that scientific topology instead of applying an independent VCF-style left shift or circular rotation. A future interchange exporter may choose another explicit representation convention, but it cannot redefine DNA's authoritative mtDNA alignment.

## Tests

- Canonical homopolymer insertion/deletion anchors are preserved.
- Mapped call positions survive representation unchanged.
- Circular origin-seam anchors remain observed rather than being rotated independently.
- Leading linear/circular predecessor derivation and reference-allele integrity remain covered.

## Traceability

ADR-0047; `SRS-VAR-004`, `SRS-VAR-007`; `INV-ALN-004`.

## Status

Implemented.
