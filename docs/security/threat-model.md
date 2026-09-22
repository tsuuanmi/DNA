# Threat Model

DNA is a local scientific CLI. Its primary security boundary is processing untrusted files and paths while protecting sensitive biological data and host resources.

## Assets

- integrity of scientific results;
- source chromatogram/reference/configuration identity;
- sensitive biological sequence/signal data;
- host filesystem integrity;
- reproducible release artifacts.

## Trust boundaries and threats

### Untrusted ABIF/FASTA/TOML bytes

Threats include malformed offsets/lengths, oversized allocation requests, parser ambiguity, and invalid cardinality.

Controls: input-size caps, checked arithmetic/slices, strict parsing, typed errors, no unsafe Rust.

### Filesystem paths and publication

Threats include path confusion, overwrite, partial-result publication, and unsafe cleanup.

Controls: validated paths, atomic no-overwrite result publication, explicit cleanup targets, typed I/O errors.

### Scientific over-interpretation

Threats include presenting unresolved/mixed evidence as stronger biological claims.

Controls: evidence hierarchy, explicit unresolved states, requirements/design invariants, independent validation.

### Sensitive data leakage

Threats include full sequences or dense signal appearing in logs, tests, issue reports, or unapproved fixtures.

Controls: [data governance](../governance/data.md), minimal operational logging, synthetic committed fixtures by default.

### Dependency/supply-chain compromise

Threats include vulnerable or unreviewed third-party packages and non-reproducible release inputs.

Controls: lockfiles, dependency review, advisory/license checks, artifact identity, [supply-chain policy](supply-chain.md).

Any future network service, unsafe/FFI boundary, secret-bearing integration, or privileged deployment requires threat-model revision and an explicit architectural/security decision.
