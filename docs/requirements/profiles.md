# Target Profile Requirements

**Requirement namespace:** `SRS-PRF-*`

These requirements are part of the canonical [DNA SRS](SRS.md). The rationale is
[ADR-0063](../decisions/adr/0063-target-profiles.md); the file format is the
[profile reference](../reference/profiles.md).

- **SRS-PRF-001:** Knowledge about a sequencing target (reference identity and topology, reportable regions, and the representation chain) MUST come from one target profile named by the configuration. Production code MUST NOT recognise a target by its reference sequence or hard-code target windows, bases, or regions.
- **SRS-PRF-002:** A profile MUST be strict TOML of at most 1 MiB with `schema_version = 1`. Unknown keys, missing required keys, unsupported enums or versions, and every violated constraint MUST fail at load, before any trace is decoded.
- **SRS-PRF-003:** A profile `id` MUST be 1-64 lowercase ASCII letters, digits, `.`, `_` or `-`, starting with a letter or digit. A declared `reference.sequence_sha256` MUST be 64 lowercase hexadecimal digits.
- **SRS-PRF-004:** `variant_calling.regions` MUST be a non-empty list of inclusive 1-based ranges within `1..=50000`; their union defines where reported variants may lie (SRS-VAR-010).
- **SRS-PRF-005:** `analyze`, `sample`, and `variant_nomenclature::apply` MUST load the reference with the profile topology and MUST fail closed for a reference that does not carry every nomenclature window sequence at its position or, when the profile declares `reference.sequence_sha256`, whose normalized-sequence SHA-256 differs.
- **SRS-PRF-006:** `normalization` and `notation` MUST be declared together. Nomenclature windows MUST be in reference order without overlap, have unique non-empty names and A/C/G/T sequences, and declare at least one rule without repetition. Every rule MUST be compatible with the window structure and its preconditions (SRS-NOM-004) MUST hold for the declared sequence; `canonical_haplotype` data MUST be given exactly when that rule is used and list unique in-window substitutions that change the reference base.
- **SRS-PRF-007:** Analysis and sample-evidence results MUST record the profile `id` and the SHA-256 of the profile file bytes. Basecalling MUST NOT read the profile.
