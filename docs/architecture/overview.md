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
| `plugin` | 0 | static plugin registry and workflow compositions, validated at compile time; plugin provenance identities | running stages, configuration values |
| `read_evidence` | 0 | modality → core per-read evidence contract: bases, profiles, labelled masks, support vetoes, informative interval | modality algorithms, interpreting reason labels |
| `input` | 3 | source/modality adapters; currently Sanger ABIF | CLI publication paths and scientific algorithms |
| `reference` | 1 | validated reference identity/model | alignment |
| `read_processing` | 2 | shared reference-free Sanger read processing; the Sanger adapter that builds `ReadEvidence` (support vetoes, mask reasons) | reference interpretation |
| `basecalling` | 2 | signal-derived calls | trimming/reference knowledge |
| `signal_processing` | 2 | observation-only Sanger signal evidence | reference interpretation |
| `callability` | 2 | signal-derived per-read phase state, typed mask, callable span; modality-generic core plus one Sanger adapter | channel/call mutation, basecall rescue, reference or profile knowledge, variant decisions |
| `quality_control` | 2 | relative quality; trim interval from the callable span | deciding callability, variant filtering |
| `alignment` | 2 | current profile-aware pairwise placement/orientation of `ReadEvidence` | variant extraction, input-format parsing, modality types |
| `variant_calling` | 2 | evidence-backed differences, mapping, allele anchoring, core eligibility gates, reporting modality vetoes and mask reasons | target nomenclature, genotype/clinical interpretation, modality types |
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

DNA composes three plugin families
([ADR-0069](../decisions/adr/0069-plugin-first-modality-core-post-calling.md)):
modality plugins (Sanger: `input`, `read_processing`, `basecalling`,
`signal_processing`, `callability`, `quality_control`), the core caller
(`read_evidence`, `alignment`, `variant_calling`, `sample`), and post-calling
plugins (`variant_representation`, `variant_normalization`,
`variant_nomenclature`). The same validator rejects any dependency from a
modality-neutral module (`plugin`, `read_evidence`, `variant`, `alignment`,
`variant_calling`, and the post-calling modules) on a Sanger module or a Sanger
child of `model`. `sample` still reads Sanger read observations until the
attachment split in [PROP-0002](../proposals/0002-plugin-first-architecture.md).

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
