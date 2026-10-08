# Dependency Policy

Dependencies and lockfiles are part of the production codebase.

## Reuse-first policy

DNA owns scientific semantics and canonical contracts; maintained ecosystem
libraries should provide commodity bioinformatics machinery whenever they meet
those contracts.

Before implementing a parser, algorithm, data structure, or infrastructure
component:

1. define the smallest DNA-owned semantic contract;
2. inventory maintained ecosystem implementations;
3. evaluate scientific semantics, maintenance, security, licensing, platform,
   performance, and supply-chain cost;
4. prefer direct reuse or a narrow adapter over a new first-party
   implementation;
5. document the concrete unmet requirement before retaining or adding custom
   machinery.

This is governed by
[ADR-0059](../decisions/adr/0059-reuse-ecosystem-machinery-behind-dna-contracts.md)
and requirements `SRS-NFR-006` / `SRS-NFR-007`.

## Bioinformatics candidates

These are evaluation candidates, not dependencies that should be added
speculatively.

| Need | First candidates | Notes |
|---|---|---|
| pairwise/sequence algorithms | `rust-bio` | evaluate before extending custom alignment machinery |
| FASTA / FASTQ | `rust-bio`, `noodles` | prefer maintained parsing; DNA retains reference/evidence semantics |
| SAM / BAM / CRAM | `noodles`, `rust-htslib` | do not implement bespoke format parsers |
| VCF / BCF | `noodles`, `rust-htslib` | external serialization/transport, not canonical DNA domain types |
| pileup / indexed HTS traversal | `rust-htslib`, `noodles` where applicable | choose by scientific and operational requirements |

Current upstream documentation:

- rust-bio: <https://docs.rs/bio/latest/bio/>
- noodles: <https://docs.rs/noodles/latest/noodles/>
- rust-htslib: <https://docs.rs/rust-htslib/latest/rust_htslib/>

## Current implementation evaluations

Evaluation is per DNA contract, not a permanent ranking of libraries.

| Current need | Evaluation | Decision |
|---|---|---|
| reference FASTA parsing | `noodles-fasta` satisfies the record-parsing need behind DNA's single-reference validation contract; the production line is constrained by the repository MSRV | reuse `noodles-fasta` behind `Reference` |
| Sanger evidence-profile alignment | `rust-bio` pairwise alignment was evaluated, but its current affine-gap convention and byte-pair substitution callback do not directly express DNA's required profile-weighted, fixed-point scoring contract; DNA additionally owns deterministic repeat placement and circular-placement semantics | retain the current custom Gotoh implementation |
| operational logging | `tracing` with a `tracing-subscriber` registry expresses stage spans and structured events; a small DNA layer renders the SRS-OUT-008 record format | reuse `tracing` behind `operation_log` ([ADR-0061](../decisions/adr/0061-tracing-for-operational-logging.md)) |
| Sanger ABIF decoding | the Applied Biosystems ABIF specification and available Rust implementations were reviewed; the evaluated `bio_files::ab1` implementation is Biopython-derived and currently documents unresolved offset handling, while DNA requires strict bounded untrusted-input parsing and an exact Sanger tag contract | retain the current bounded ABIF implementation; reevaluate only against equivalent contract and corpus evidence |
| future FASTQ / BAM / CRAM / VCF / BCF | no current production consumer | do not add dependencies speculatively |

A retained custom implementation is not exempt from reuse-first review. Its
owning design must state the concrete semantic or operational gap, and the
decision should be revisited when either DNA's contract or the ecosystem
implementation materially changes.

## Adding a dependency

Before adding or materially changing a dependency:

1. confirm there is a current production consumer;
2. confirm the capability is not already available through an existing
   dependency;
3. justify maintenance, security, licensing, platform, and supply-chain cost;
4. keep dependency-specific representations behind DNA-owned boundaries unless
   the external type is itself the intentional public contract;
5. use the repository's package manager and lockfile workflow;
6. review the resulting lockfile diff;
7. update security/release evidence when the dependency changes the trust
   boundary.

Production dependencies must not be added speculatively. Release dependency
review is tracked by [supply-chain security](../security/supply-chain.md).

## Custom implementation threshold

A custom implementation is justified only when a concrete requirement cannot be
met cleanly through a maintained implementation or adapter. Examples include a
DNA-specific scientific invariant, unsupported scoring semantics, incompatible
platform/dependency constraints, or validation evidence that an available
implementation is unsuitable.

When custom code is retained despite an ecosystem alternative, the owning design
or decision documentation should state the gap explicitly.
