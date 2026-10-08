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

| Module | Layer | Owns | Excludes |
|---|---:|---|---|
| `cli` | 4 | command syntax | I/O and algorithms |
| `config` | 1 | strict parsing, validation, caps; names the target profile | per-value environment overrides, target knowledge |
| `profile` | 1 | target profiles: reference identity/topology, regions, nomenclature windows, notation | running rules, loading references |
| `error` | 0 | typed cross-stage failures | logging and recovery policy |
| `operation_log` | 4 | append-only operational records rendered from `tracing` events | scientific decisions |
| `checksum` | 0 | stable SHA-256 byte identity | file I/O and policy |
| `locus` | 0 | Sanger locus window geometry | signal interpretation |
| `model` | 0 | validated domain vocabulary including canonical Sanger evidence | filesystem and algorithms |
| `variant` | 0 | public canonical called-variant contracts | evidence and call mappings |
| `input` | 3 | source/modality adapters; currently Sanger ABIF | CLI publication paths and scientific algorithms |
| `reference` | 1 | validated reference identity/model | alignment |
| `read_processing` | 2 | shared reference-free Sanger read processing | reference interpretation |
| `basecalling` | 2 | signal-derived calls | trimming/reference knowledge |
| `signal_processing` | 2 | observation-only Sanger signal evidence | reference interpretation |
| `quality_control` | 2 | relative quality/end trimming | variant filtering |
| `alignment` | 2 | current profile-aware pairwise placement/orientation | variant extraction and input-format parsing |
| `variant_calling` | 2 | evidence-backed differences, mapping, allele anchoring, eligibility | target nomenclature, genotype/clinical interpretation |
| `variant_analysis` | 3 | reference-guided one-read scientific composition and public Variant Analysis capability | CLI logging/JSON publication |
| `sample` | 2 | multi-read evidence aggregation | input discovery/consensus |
| `variant_representation` | 2 | haplotype-preserving edit conversion, application, rendering | policy choices |
| `variant_normalization` | 2 | optional sequence-equivalent normalization policies | nomenclature windows |
| `variant_nomenclature` | 2 | profile-driven window representation engine | target knowledge, notation rendering |
| `report` | 4 | contract projection, serialization, atomic publish | scientific decisions |
| `pipeline` | 4 | CLI/sample orchestration, path/log/publication lifecycle | scientific implementation ownership |

Layers are core (0), target data (1), science (2), adapters and capabilities (3),
and delivery (4). A module depends only on the same or a lower layer, without
cycles, so a layer can later become a crate ([ADR-0064](../decisions/adr/0064-crate-ready-module-layering.md));
CI enforces this with `validate_module_layers.py`.

Dependencies point toward shared model/config/error and capability boundaries;
cycles are forbidden.

External libraries may implement commodity algorithms or formats behind these
boundaries, but dependency-specific types must not become accidental DNA domain
contracts. See [dependency policy](../engineering/dependencies.md) and
[ADR-0059](../decisions/adr/0059-reuse-ecosystem-machinery-behind-dna-contracts.md).

## Resource bounds

Config/profile/FASTA source files, ABIF input, normalized reference length, changed indel
length, and alignment cells are explicitly bounded. Exact current limits and
requirements are owned by [requirements](../requirements/README.md) and
[configuration reference](../reference/configuration.md).
