# Reference Comparison, Quality, and Validation Requirements

**Requirement namespaces:** `SRS-COMPAT-*`, `SRS-NFR-*`, `SRS-VAL-*`

These requirements are part of the canonical [DNA SRS](README.md).

- **SRS-COMPAT-001:** Apollo comparisons MUST follow [reference validation policy](../governance/reference-validation.md); known defects are intentional divergences, not parity failures or backward-compatibility obligations.
- **SRS-COMPAT-002:** Approved differential evidence MUST compare exact decoded arrays and unaffected deterministic results; normalized variants compare by full tuple without ignoring extras/missing calls.
- **SRS-NFR-001:** Production code MUST forbid unsafe Rust, deny deprecated API use at crate roots, and avoid production `unwrap`/`expect`. First-party source MUST NOT suppress the deprecated-use diagnostic.
- **SRS-NFR-002:** Scientific stage functions MUST be side-effect-free and return typed results; only pipeline-level operational logging and report publication write files.
- **SRS-NFR-003:** Representative 500–1,000 base release analysis SHOULD complete within 30 seconds and 512 MiB on a documented host.
- **SRS-NFR-004:** Every directory under `src/` MUST contain an up-to-date `README.md` describing implementation ownership and routing to canonical documentation; file-only modules MAY rely on rustdoc and source comments.
- **SRS-NFR-005:** First-party production Rust MUST NOT retain deprecated compatibility APIs, explicit legacy/backward-compatibility feature paths, or suppress `deprecated`, dead-code, unreachable-code, or unused-code diagnostics to preserve obsolete source. The required Rust source-policy validator MUST fail such scaffolding.
- **SRS-VAL-001:** Parser, calling, signal processing, QC, alignment, normalization, JSON, and atomic publication MUST have focused boundary/adversarial tests.
- **SRS-VAL-002:** A synthetic canonical ABIF MUST exercise end-to-end forward/reverse and variant behavior without identifying data.
- **SRS-VAL-003:** Real-trace release evidence MUST follow [data policy](../governance/data.md); ignored local data MUST never be a build/test prerequisite.
- **SRS-VAL-004:** Rust source-policy validation, documentation-structure validation, format, check, Clippy warnings-denied, all tests, rustdoc, schema/example validation, TOML validation, source-directory README coverage, and rCRS identity gates MUST pass.
