# Source Mirror

`docs/src/` mirrors first-party Rust source ownership. Every `src/**/*.rs` file covered by the repository mirror policy has a same-relative-path manual under `docs/src/`.

## Module manual template

A manual should answer:

1. **Responsibility** — what this module owns.
2. **Inputs** — validated values it consumes.
3. **Outputs** — typed values it produces.
4. **Invariants** — cross-cutting `INV-*` rules and local invariants it preserves.
5. **Dependencies** — modules it may depend on and important boundaries it must not cross.
6. **Failure modes** — expected typed failures and unsupported states.
7. **Algorithm boundary** — what scientific/technical decision is implemented here, without translating code line by line.
8. **Traceability** — relevant requirements, decisions, public contracts, and tests.

## Synchronization rule

- Source responsibility change -> update the mirror in the same change.
- New mirrored source file -> add the matching manual.
- Removed source file -> remove or relocate the manual.
- Pure implementation detail that does not change responsibility does not require rewriting the manual.

The mirror is descriptive implementation documentation. It does not override normative requirements, public contracts, accepted decisions, or current method documentation.
