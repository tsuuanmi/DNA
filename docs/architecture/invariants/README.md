# System Invariants

This is the canonical index for cross-cutting DNA invariants. Invariant IDs remain stable even when their physical document changes.

| Family | Scope | Document |
|---|---|---|
| `INV-EVID-*` | immutable source and derived signal evidence | [Evidence](evidence.md) |
| `INV-COORD-*`, `INV-ID-*` | coordinate domains and identity preservation | [Coordinates and identity](coordinates-identity.md) |
| `INV-ALN-*` | alignment optimality and canonical topology | [Alignment canonicality](alignment.md) |
| `INV-BIO-*` | scientific interpretation boundaries | [Scientific state](scientific-state.md) |
| `INV-READ-*`, `INV-SAMPLE-*` | independent read placement and sample evidence | [Read/sample boundaries](sample-boundaries.md) |
| `INV-PIPE-*` | stage ordering, determinism, and failure behavior | [Pipeline](pipeline.md) |
| `INV-OUT-*` | result publication, schemas, and operational separation | [Output and operations](output-operations.md) |
| `INV-RUST-*` | Rust/source-policy correctness constraints | [Rust implementation](rust.md) |

## Use

Read only the invariant families relevant to the change, then follow [traceability](../../validation/traceability.md) to the owning requirements, methods, contracts, implementation, and tests.

Do not duplicate an invariant into another architecture document. Other documents should reference the invariant ID and this canonical family instead.
