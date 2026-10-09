# Core model

Owns the records the core caller produces and consumes:
- `alignment`: placement, orientation, metrics;
- `variant`: internal variants with call mappings and eligibility;
- `coordinate`;
- `called_read`: one read after the core;
- `reference_call`: resolves a variant's public calls and their reference-strand bases;
- `sample_evidence`: sample aggregation.

These records hold no modality evidence. Reports join that evidence by read
identity and call index.

See [dna-core](../README.md).
