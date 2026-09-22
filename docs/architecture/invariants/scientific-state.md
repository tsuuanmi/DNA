# Scientific state Invariants

These invariants are part of the canonical [system invariant set](README.md).

- **INV-BIO-001:** Unresolved evidence is not equivalent to a reference call, absence of variation, or absence of coverage.
- **INV-BIO-002:** A mixed/secondary signal is an observation, not automatically heteroplasmy, genotype, contamination, or mixture.
- **INV-BIO-003:** A single chromatogram produces read-level evidence, not a sample-level biological conclusion.
- **INV-BIO-004:** A derived confidence value is not an error probability or Phred score unless separately calibrated and validated.
- **INV-BIO-005:** A strongest canonical base with more than one co-localized qualifying channel is mixed signal evidence, not an ordinary clean substitution. If it yields a normalized SNV observation, the observation remains preserved but is not clean-report eligible.
