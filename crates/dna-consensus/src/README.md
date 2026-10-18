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
- `decided.rs`: the consensus's decided edits as called records
  (SRS-VAR-021, ADR-0078, plugin `decided`): each segment's columns from its
  sites, each site's call aligned against its own interval (an `N` stands for
  its position and is not callable), gaps moved 3'-most, then the record rules
  of `variant_calling`; a length change of a segment-edge run, and a change
  at the 3' end of a run of length evidence `phase_loss` or
  `reference_frame`, excluded as `run_boundary`, with that end left out of
  the callable segments.
- `mod.rs`: the entry point `build`, which handles stretches of whole
  reference runs (seeded by clean differences and by disputed runs no read
  bounds), segment assembly, and the `consensus` plugin descriptor.
  Each stretch is decided by its run structure (ADR-0074): composition, then
  lengths with their evidence (`in_phase`, `anchored_end`, `phase_loss`,
  `reference_frame`), or position by position.

The module reads only `CalledRead` records, so it is modality-neutral
(ADR-0069). It uses no profile window or filename. `decided.rs` turns its
decisions into called records; it never places a segment again.

See [dna-core](../README.md).
