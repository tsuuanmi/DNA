# Variant Calling

Owns extraction, direct evidence mapping, allele anchoring, and configured
eligibility of primary-sequence SNVs and supported small indels. Reportable
regions come from the target profile and are passed in by the caller.

Key children: `extract.rs`, `mapping.rs`, `anchor.rs`, `eligibility.rs`, and
`filter.rs`. The module reads each read only as `ReadEvidence`
([ADR-0069](../../../../docs/decisions/adr/0069-plugin-first-modality-core-post-calling.md))
and imports no Sanger type.
- `eligibility.rs` decides which calls of a read can support a variant. Calls
  the modality places at a signal edge (`CallEvidence::edge`), or within the
  margin of a read end it does not vouch for, give `read_end` (ADR-0062). A
  call vouched for by its base only supports a substitution but no indel and
  shows no run's end; one that may sit one position off gives
  `unresolved_flank` (SRS-VAR-018). An edit that changes the length of
  a run of one base whose end the read does not resolve gives
  `run_boundary` (ADR-0071); the sample consensus applies the same run-end
  test to the indels it weighs. An insertion into a repeat the read bounds is
  judged for support vetoes on its placements along the repeat
  (`repeat_unvetoed`, SRS-VAR-020), which the consensus shares. Masked
  evidence calls give their mask reason, verbatim.
- `filter.rs` applies the region gate and reports the union of the supporting
  calls' modality vetoes, in vocabulary order and filtered by scope, except
  for an insertion `repeat_unvetoed` clears; it reads the reference around a
  substitution for the one-position-shift test.
- A masked call is no SNV candidate in `extract.rs`.

The core never interprets a modality label; the Sanger labels are produced in
[read_processing](../../../dna-sanger/src/read_processing/README.md).
`anchor.rs` constructs and validates anchored REF/ALT alleles while preserving
the alignment-selected event placement; post-calling haplotype normalization
and target nomenclature belong to `variant_normalization` and
`variant_nomenclature`.

This module does not infer genotype, heteroplasmy, phase, pathogenicity, or
clinical significance.

`config.rs` owns the `[variant_calling]`, with the indel-length cap configuration section: its raw record, validated
record, and rules (ADR-0069).

See [variant requirements](../../../../docs/requirements/variants.md),
[variant-calling method](../../../../docs/design/variant-calling.md), and
[scientific-state invariants](../../../../docs/architecture/invariants/scientific-state.md).
