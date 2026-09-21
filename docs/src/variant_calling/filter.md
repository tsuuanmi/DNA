# `src/variant_calling/filter.rs`

## Purpose

Applies configured report eligibility to normalized primary-sequence variant
candidates.

## Responsibilities

- Keep candidates whose normalized 1-based anchor lies in the union of configured
  inclusive regions.
- Require every SNV or inserted-base supporting call to meet the configured
  maximum-channel peak floor and strictly exceed the relative-quality threshold.
- Require an SNV supporting call to retain exactly one co-localized qualifying channel; mixed supporting calls remain observed but receive `mixed_supporting_dna` and are not clean-report eligible.
- Do not apply the point-mixed-signal rule to insertion/deletion evidence.
- Exempt insertion flanks and deletion flanks from supporting-signal thresholds.
- Retain every normalized canonical candidate in the observed stream with all failed region/peak/relative-quality rules, even when it is not reportable.
- Add each non-reportable candidate once to the concise exclusion diagnostics used for warnings/logging.

## Non-responsibilities

No extraction, normalization, vendor-quality filtering, genotype inference,
report projection, report-only label assignment, or logging.

## Key functions

- `apply(extracted, calls, quality, config) -> Result<VariantCallingResult>`:
  evaluates normalized candidates, returning both configured-eligible reported variants and the complete normalized observed stream.
- `supporting_evidence_reasons`: selects relevant supporting mappings by kind and
  returns each failed evidence rule once.
- `assess_call`: joins one mapping to its base-call peaks, relative quality, and authoritative qualifying-channel state, then evaluates peak, quality, and mixed-signal evidence independently.

## Invariants and errors

Region positions and normalized variant anchors are 1-based and inclusive.
Call/PLOC indexes remain 0-based. Filtering eligibility never erases a normalized observation. Every non-reportable candidate creates exactly one allele-free diagnostic, even when multiple rules fail. Missing or mismatched
call/quality mappings return `Error::Variant`. Vendor PCON does not affect
eligibility.

## Tests

Unit tests cover threshold boundaries, region endpoints, mixed-SNV exclusion without applying the rule to insertions, multi-base insertions, deletion exemption, mapping errors, exclusion accounting, and retained observed-candidate eligibility.

## Status

Implemented.
