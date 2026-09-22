# Model

Owns validated domain vocabulary shared across scientific stages.

Model modules contain typed state and serialization records; they do not own
filesystem access, CLI parsing, configuration loading, or algorithm orchestration.

Important children cover alignment, basecalls, coordinates, locus evidence,
nucleotides, quality, references, read observations, sample evidence/results,
signal, traces, and variants.

See [system architecture](../../docs/architecture/system.md) and
[system invariants](../../docs/architecture/invariants/README.md).
