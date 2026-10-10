# Consensus

Owns the sample consensus
([ADR-0073](../../../../docs/decisions/adr/0073-sample-consensus.md),
[ADR-0074](../../../../docs/decisions/adr/0074-run-structure-consensus.md),
[method](../../../../docs/design/consensus.md)). It adjudicates between a
sample's placed reads and assembles consensus segments, with the evidence
behind every non-trivial decision.

- `observe.rs`: each read's observations per reference position and junction,
  with each base's trust level (`Trust`) and reference position. It reuses
  `variant_calling::eligibility` for the trusted calls and for the run-end
  test of insertions and deletions (ADR-0071).
- `decide.rs`: decision rule version 1 over clean observations.
- `mod.rs`: the entry point `build`, which handles stretches of whole
  reference runs, segment assembly, and the `consensus` plugin descriptor.
  Each stretch is decided by its run structure (ADR-0074): composition, then
  lengths with their evidence (`in_phase`, `anchored_end`, `phase_loss`,
  `reference_frame`), or position by position.

The module reads only `CalledRead` records, so it is modality-neutral
(ADR-0069). It uses no profile window or filename. It never calls variants:
`dna call` calls the consensus FASTA.

See [dna-core](../README.md).
