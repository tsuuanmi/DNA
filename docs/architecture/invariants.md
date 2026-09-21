# System Invariants

This is the canonical index for cross-cutting DNA invariants. Invariant IDs remain stable even when their physical document changes.

| Family | Scope | Document |
|---|---|---|
| `INV-EVID-*` | immutable source and derived signal evidence | [Evidence](invariants/evidence.md) |
| `INV-COORD-*`, `INV-ID-*` | coordinate domains and identity preservation | [Coordinates and identity](invariants/coordinates-identity.md) |
| `INV-ALN-*` | alignment optimality and canonical topology | [Alignment canonicality](invariants/alignment.md) |
| `INV-BIO-*` | scientific interpretation boundaries | [Scientific state](invariants/scientific-state.md) |
| `INV-READ-*`, `INV-SAMPLE-*` | independent read placement and sample evidence | [Read/sample boundaries](invariants/sample-boundaries.md) |
| `INV-PIPE-*` | stage ordering, determinism, and failure behavior | [Pipeline](invariants/pipeline.md) |
| `INV-OUT-*` | result publication, schemas, and operational separation | [Output and operations](invariants/output-operations.md) |
| `INV-RUST-*` | Rust/source-policy correctness constraints | [Rust implementation](invariants/rust.md) |

## Use

Read only the invariant families relevant to the change, then follow [traceability](../traceability.md) to the owning requirements, methods, contracts, implementation, and tests.

Do not duplicate an invariant into another architecture document. Other documents should reference the invariant ID and this canonical family instead.
