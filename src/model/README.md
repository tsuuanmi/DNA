# Model

Owns validated domain vocabulary shared across scientific stages.

Model modules contain typed state and serialization records; they do not own
filesystem access, CLI parsing, configuration loading, or algorithm orchestration.

Important children cover alignment, basecalls, callability, coordinates, locus
evidence, nucleotides, quality, references, read observations, sample
evidence/results, signal, canonical Sanger evidence, and variants. `nucleotide` is the single
authority for canonical-base parsing and A/C/G/T channel order, and
`reference_call` resolves a variant call mapping to its reference-strand called
base for both sample aggregation and reporting.

`EvidenceProfile` belongs to the modality-neutral `read_evidence` contract, not
to `model`. The modality-neutral children are `alignment`, `called_read` (the
core's per-read products), `coordinate`, `nucleotide`, `reference`,
`reference_call`, `sample_evidence`, and `variant`. `read_observation` pairs a
`CalledRead` with its `SangerAttachment`. The other children are Sanger
vocabulary, and the module validator keeps neutral modules from importing them
([ADR-0069](../../docs/decisions/adr/0069-plugin-first-modality-core-post-calling.md)).

See [system architecture](../../docs/architecture/overview.md) and
[system invariants](../../docs/architecture/invariants/README.md).
