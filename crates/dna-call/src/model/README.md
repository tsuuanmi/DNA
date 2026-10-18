# Call model

Owns the records of one called read:
- `variant`: internal variants with call mappings and eligibility, and the
  candidates excluded before they were observed;
- `coordinate`;
- `called_read`: one read after the core, and the call workflow's cut of its
  excluded variants' positions from its callable segments (SRS-VAR-019);
- `reference_call`: resolves a variant's public calls and their
  reference-strand bases;
- `trust`: the core's per-call trust verdict (`CallTrust`, `Withholder`),
  diagnostic evidence for evaluation, not a calling contract.

These records hold no modality evidence. Reports join that evidence by read
identity and call index.

See [dna-call](../README.md).
