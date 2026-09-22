# Quality Attributes

**Requirement namespace:** `SRS-NFR-*`

These requirements are part of the canonical [DNA SRS](SRS.md).

- **SRS-NFR-001:** Production code MUST forbid unsafe Rust, deny deprecated API use at crate roots, and avoid production `unwrap`/`expect`. First-party source MUST NOT suppress the deprecated-use diagnostic.
- **SRS-NFR-002:** Scientific stage functions MUST be side-effect-free and return typed results; only pipeline-level operational logging and report publication write files.
- **SRS-NFR-003:** Representative 500–1,000 base release analysis SHOULD complete within 30 seconds and 512 MiB on a documented host.
- **SRS-NFR-004:** Every directory under `src/` MUST contain an up-to-date `README.md` describing implementation ownership and routing to canonical documentation; file-only modules MAY rely on rustdoc and source comments.
- **SRS-NFR-005:** First-party production Rust MUST NOT retain deprecated compatibility APIs, explicit legacy/backward-compatibility feature paths, or suppress `deprecated`, dead-code, unreachable-code, or unused-code diagnostics to preserve obsolete source. The required Rust source-policy validator MUST fail such scaffolding.
