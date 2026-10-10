# Consensus

Owns the sample consensus
([ADR-0073](../../../../docs/decisions/adr/0073-sample-consensus.md),
[method](../../../../docs/design/consensus.md)). It adjudicates between a
sample's placed reads and assembles consensus segments, with the evidence
behind every non-trivial decision.

- `observe.rs`: each read's observations per reference position and junction,
  with their cleanliness and run-end visibility. It reuses
  `variant_calling::eligibility` for the trusted calls.
- `decide.rs`: decision rule version 1 over clean observations.
- `mod.rs`: the entry point `build`, which handles stretches of whole
  reference runs (decided as a whole, run by run, or position by position),
  segment assembly, and the `consensus` plugin descriptor.

The module reads only `CalledRead` records, so it is modality-neutral
(ADR-0069). It uses no profile window or filename. It never calls variants:
`dna call` calls the consensus FASTA.

See [dna-core](../README.md).
