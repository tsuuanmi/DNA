# Quality Attributes

**Requirement namespace:** `SRS-NFR-*`

These requirements are part of the canonical [DNA SRS](SRS.md).

- **SRS-NFR-001:** Production code MUST forbid unsafe Rust, deny deprecated API use at crate roots, and avoid production `unwrap`/`expect`. First-party source MUST NOT suppress the deprecated-use diagnostic.
- **SRS-NFR-002:** Scientific stage functions MUST be side-effect-free and return typed results; only pipeline-level operational logging and report publication write files.
- **SRS-NFR-003:** Representative 500–1,000 base release analysis SHOULD complete within 30 seconds and 512 MiB on a documented host.
- **SRS-NFR-004:** Every directory under `src/` and `crates/*/src/` MUST contain an up-to-date `README.md` describing implementation ownership and routing to canonical documentation; file-only modules MAY rely on rustdoc and source comments.
- **SRS-NFR-005:** First-party production Rust MUST NOT retain deprecated compatibility APIs, explicit legacy/backward-compatibility feature paths, or suppress `deprecated`, dead-code, unreachable-code, or unused-code diagnostics to preserve obsolete source. The required Rust source-policy validator MUST fail such scaffolding.

- **SRS-NFR-006:** Before introducing a custom parser, algorithm, data structure, or infrastructure component, implementation work MUST evaluate maintained ecosystem implementations and prefer reuse or adaptation when they satisfy DNA's semantic, validation, operational, licensing, security, and platform requirements. A custom implementation MUST have a documented unmet requirement or constraint.
- **SRS-NFR-007:** Production support for standard bioinformatics formats MUST use a maintained standards-aware implementation when one satisfies the required contract. DNA MUST NOT introduce bespoke SAM/BAM/CRAM or VCF/BCF parsers under that condition; FASTA/FASTQ parsing SHOULD likewise reuse maintained ecosystem support unless a documented DNA-specific requirement justifies otherwise.
