# Reference

Owns the DNA reference contract: strict single-record FASTA validation, topology,
normalization, and stable sequence identity.

Commodity FASTA record parsing is delegated to `noodles-fasta`; dependency-specific
record types do not cross this module boundary. DNA remains responsible for the
single-record rule, accepted A/C/G/T/N alphabet, size and length limits, identifier
semantics, topology, and SHA-256 identity.

The module does not perform alignment or reference selection across multiple
contigs.

See [input requirements](../../docs/requirements/input.md) and
[system architecture](../../docs/architecture/overview.md).
