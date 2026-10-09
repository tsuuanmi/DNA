# Sample aggregation

Owns the sample-level entry point: it validates that reads share one reference
and configuration identity and contain no duplicate trace content, orders reads
deterministically, and assembles coverage, overlap, locus, and variant evidence
from the sibling `sample` modules.

- `mod.rs` — validation, ordering, and evidence assembly.
- `tests.rs` — end-to-end aggregation over synthetic `CalledRead` records.

See [sample methods](../../../docs/design/sample-evidence/README.md) and
[sample invariants](../../../docs/architecture/invariants/sample-boundaries.md).
