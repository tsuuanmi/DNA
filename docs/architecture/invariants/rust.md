# Rust implementation Invariants

These invariants are part of the canonical [system invariant set](../invariants.md).

- **INV-RUST-001:** First-party production code forbids unsafe Rust and denies deprecated API use; first-party source cannot suppress that diagnostic under the source-policy gate.
- **INV-RUST-002:** Production paths do not use `unwrap` or `expect` for recoverable external conditions.
- **INV-RUST-003:** Types and module boundaries should encode coordinate, topology, strand, and validated-state distinctions when doing so removes a concrete failure mode.
- **INV-RUST-004:** Production source must not hide obsolete code behind deprecated declarations, legacy/backward-compatibility feature paths, compatibility-named declarations, or warning suppressions for deprecated/dead/unreachable/unused code; CI enforces this explicit-source policy in addition to compiler and Clippy diagnostics.
