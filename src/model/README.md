# Model

Owns validated domain vocabulary shared across scientific stages.

Model modules contain typed state and serialization records; they do not own
filesystem access, CLI parsing, configuration loading, or algorithm orchestration.

Important children cover alignment, basecalls, callability, coordinates, locus
evidence, nucleotides, quality, references, read observations, sample
evidence/results, signal, canonical Sanger evidence, and variants. `nucleotide` is the single
authority for canonical-base parsing and A/C/G/T channel order, and
`reference_call` resolves a variant call mapping to its reference-strand call
evidence for both sample aggregation and reporting.

See [system architecture](../../docs/architecture/overview.md) and
[system invariants](../../docs/architecture/invariants/README.md).
