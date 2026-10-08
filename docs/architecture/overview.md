# Architecture Overview

## Principles

- Biological target, sequencing modality, external file format, and analysis
  capability are separate architectural dimensions.
- Canonicalize formats at adapters; canonicalize biological evidence only where
  semantics actually converge.
- DNA owns scientific semantics and canonical contracts; maintained ecosystem
  implementations provide commodity bioinformatics machinery when they satisfy
  those contracts.
- One strict configuration and one command-specific JSON result per CLI
  invocation.
- Untrusted binary/text input is checked before slicing, conversion, and large
  allocation.
- Models enforce cardinality and coordinate invariants; scientific stages avoid
  filesystem side effects.
- The CLI and operating-system boundary remain thin.
- Algorithms are deterministic and biologically explicit.
- Every source directory has a colocated README that routes implementation
  ownership.

## Module boundaries

| Module | Owns | Excludes |
|---|---|---|
| `cli` | command syntax | I/O and algorithms |
| `config` | strict parsing, validation, caps | per-value environment overrides |
| `error` | typed cross-stage failures | logging and recovery policy |
| `operation_log` | append-only operational records rendered from `tracing` events | scientific decisions |
| `checksum` | stable SHA-256 byte identity | file I/O and policy |
| `model` | validated domain vocabulary including canonical Sanger evidence | filesystem and algorithms |
| `input` | source/modality adapters; currently Sanger ABIF | CLI publication paths and scientific algorithms |
| `reference` | validated reference identity/model | alignment |
| `read_processing` | shared reference-free Sanger read processing | reference interpretation |
| `basecalling` | signal-derived calls | trimming/reference knowledge |
| `signal_processing` | observation-only Sanger signal evidence | reference interpretation |
| `quality_control` | relative quality/end trimming | variant filtering |
| `alignment` | current profile-aware pairwise placement/orientation | variant extraction and input-format parsing |
| `variant_calling` | evidence-backed differences, mapping, allele anchoring, eligibility | target nomenclature, genotype/clinical interpretation |
| `variant_analysis` | reference-guided one-read scientific composition and public Variant Analysis capability | CLI logging/JSON publication |
| `sample` | multi-read evidence aggregation | input discovery/consensus |
| `report` | contract projection, serialization, atomic publish | scientific decisions |
| `pipeline` | CLI/sample orchestration, path/log/publication lifecycle | scientific implementation ownership |

Dependencies point toward shared model/config/error and capability boundaries;
cycles are forbidden.

External libraries may implement commodity algorithms or formats behind these
boundaries, but dependency-specific types must not become accidental DNA domain
contracts. See [dependency policy](../engineering/dependencies.md) and
[ADR-0059](../decisions/adr/0059-reuse-ecosystem-machinery-behind-dna-contracts.md).

## Resource bounds

Config/FASTA source files, ABIF input, normalized reference length, changed indel
length, and alignment cells are explicitly bounded. Exact current limits and
requirements are owned by [requirements](../requirements/README.md) and
[configuration reference](../reference/configuration.md).
