# Consensus unit tests

Submodules of the crate's unit tests ([tests.rs](../tests.rs)), which own the
shared fixtures.

- `n_reads.rs` — the passes over more than one read per strand (ADR-0090):
  set-aside, the one-strand insertion contest, confirmed ends, the
  two-strand rule and its cross-strand pairs with two forward reads and one
  or two reverse reads; ties of three or more sequences; and a seeded
  property test that a set-aside keeps every decision's value and that the
  set-aside passes end.

See [dna-consensus](../README.md).
