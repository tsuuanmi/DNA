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

| Module | Crate | Owns | Excludes |
|---|---|---|---|
| `cli` | dna | command syntax | I/O and algorithms |
| `config` | dna | strict loading of the configuration envelope that composes the plugin-owned sections; names the target profile | section rules (owned by each plugin), per-value environment overrides, target knowledge |
| `profile` | kernel | target profiles: reference identity/topology, regions, nomenclature windows, notation | running rules, loading references |
| `error` | kernel | typed cross-stage failures | logging and recovery policy |
| `operation_log` | dna | append-only operational records rendered from `tracing` events | scientific decisions |
| `checksum` | kernel | stable SHA-256 byte identity | file I/O and policy |
| `locus` | sanger | Sanger locus window geometry | signal interpretation |
| `model` | kernel, core, sanger, dna | validated domain vocabulary; each child belongs to the crate of its family (Sanger evidence and attachments to `sanger`) | filesystem and algorithms |
| `variant` | kernel | public canonical called-variant contracts | evidence and call mappings |
| `plugin` | kernel | plugin descriptor types and their compile-time validation; each plugin declares its own descriptor, and `pipeline::plugins` lists the registry and workflow compositions | running stages, configuration values |
| `read_evidence` | kernel | modality → core per-read evidence contract: bases, profiles, labelled masks, support vetoes, informative interval | modality algorithms, interpreting reason labels |
| `input` | dna, sanger | loading orchestration (configuration, profile, reference) and source adapters: Sanger ABIF decoding (`sanger`), reviewed consensus sequences with their `ReadEvidence` adapter, and variants documents | CLI publication paths and scientific algorithms |
| `reference` | kernel | validated reference identity/model | alignment |
| `read_processing` | sanger | shared reference-free Sanger read processing; the Sanger adapter that builds `ReadEvidence` (support vetoes, mask reasons) | reference interpretation |
| `basecalling` | sanger | signal-derived calls | trimming/reference knowledge |
| `signal_processing` | sanger | observation-only Sanger signal evidence | reference interpretation |
| `callability` | sanger | signal-derived per-read phase state, typed mask, callable span; modality-generic core plus one Sanger adapter | channel/call mutation, basecall rescue, reference or profile knowledge, variant decisions |
| `quality_control` | sanger | relative quality; trim interval from the callable span | deciding callability, variant filtering |
| `alignment` | core | current profile-aware pairwise placement/orientation of `ReadEvidence` | variant extraction, input-format parsing, modality types |
| `variant_calling` | core | evidence-backed differences, mapping, allele anchoring, core eligibility gates, reporting modality vetoes and mask reasons | target nomenclature, genotype/clinical interpretation, modality types |
| `variant_analysis` | dna | Sanger one-read composition (modality plus core) and the public Variant Analysis capability | CLI logging/JSON publication |
| `read_call` | core | the core's one-read path from `ReadEvidence` to a `CalledRead`, and the core-owned configuration | modality types |
| `sample` | core | multi-read evidence aggregation over modality-neutral `CalledRead` records | input discovery/consensus, modality types |
| `variant_representation` | post | haplotype-preserving edit conversion, application, rendering | policy choices |
| `variant_normalization` | post | optional sequence-equivalent normalization policies | nomenclature windows |
| `variant_nomenclature` | post | profile-driven window representation engine | target knowledge, notation rendering |
| `conformance` | post | reporting represented calls that break the profile's notation conventions | rewriting calls, target knowledge |
| `report` | dna | contract projection, serialization, atomic publish | scientific decisions |
| `pipeline` | dna | CLI/sample orchestration, path/log/publication lifecycle | scientific implementation ownership |

Every module belongs to one crate of the plugin-first workspace
([ADR-0069](../decisions/adr/0069-plugin-first-modality-core-post-calling.md),
[PROP-0002](../proposals/0002-plugin-first-architecture.md) phase 5):

- `kernel` holds the shared contracts;
- `core` is the core caller;
- `sanger` is the Sanger modality;
- `post` holds the post-calling plugins;
- `dna` is the facade that composes them.

The plugin crates depend only on `kernel`, and only `dna` composes them. Each
plugin owns its configuration sections and its descriptor. CI enforces the
crate map and an acyclic module graph with `validate_module_layers.py`.

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
