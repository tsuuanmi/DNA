# Architecture Overview

## Principles

- One strict configuration and one command-specific JSON result per invocation.
- Untrusted binary/text input is checked before slicing, conversion, and large allocation.
- Models enforce cardinality and coordinate invariants; scientific stages avoid filesystem side effects.
- The CLI and operating-system boundary remain thin.
- Algorithms are deterministic and biologically explicit.
- Every source directory has a colocated README that routes implementation ownership.

## Module boundaries

| Module | Owns | Excludes |
|---|---|---|
| `cli` | command syntax | I/O and algorithms |
| `config` | strict parsing, validation, caps | per-value environment overrides |
| `error` | typed cross-stage failures | logging and recovery policy |
| `logger` | append-only operational records | scientific decisions |
| `checksum` | stable SHA-256 byte identity | file I/O and policy |
| `model` | validated domain vocabulary | filesystem and algorithms |
| `trace` | canonical ABIF decode | base calling |
| `reference` | one-record FASTA and identity | alignment |
| `basecalling` | signal-derived calls | trimming/reference knowledge |
| `signal_processing` | observation-only signal evidence | reference interpretation |
| `quality_control` | relative quality/end trimming | variant filtering |
| `alignment` | profile-aware placement/orientation | variant extraction |
| `variant_calling` | differences, mapping, normalization, eligibility | genotype/clinical interpretation |
| `sample` | multi-read evidence aggregation | input discovery/consensus |
| `report` | contract projection, serialization, atomic publish | scientific decisions |
| `pipeline` | command sequencing | algorithm internals |

Dependencies point toward shared model/config/error boundaries; cycles are forbidden.

## Resource bounds

Config/FASTA source files, AB1 input, normalized reference length, changed indel length, and alignment cells are explicitly bounded. Exact current limits and requirements are owned by [requirements](../requirements/README.md) and [configuration reference](../reference/configuration.md).
