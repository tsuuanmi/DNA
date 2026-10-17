# Core model

Owns the records the core caller produces and consumes:
- `alignment`: placement, orientation, metrics;
- `variant`: internal variants with call mappings and eligibility;
- `coordinate`;
- `called_read`: one read after the core, and the call workflow's cut of its
  excluded variants' positions from its callable segments (SRS-VAR-019);
- `reference_call`: resolves a variant's public calls and their reference-strand bases;
- `sample_evidence`: sample aggregation;
- `trust`: the core's per-call trust verdict (`CallTrust`, `Withholder`),
  diagnostic evidence for evaluation, not a calling contract;
- `consensus`: consensus segments, decided sites, and decision counts.

These records hold no modality evidence. Reports join that evidence by read
identity and call index.

See [dna-core](../README.md).
