# DNA Architecture

## Principles

- One strict config and one command-specific JSON result per invocation. Single-read commands accept one AB1; sample evidence accepts one or more AB1 traces. Reference-guided commands additionally require one FASTA.
- Untrusted binary input is checked before every slice, conversion, and allocation.
- Models enforce cardinality and coordinate invariants; scientific functions have no filesystem side effects.
- The CLI and operating-system boundary remain thin.
- Algorithms are deterministic and biologically explicit; known Apollo defects are not parity requirements.
- Every source directory has an up-to-date colocated `README.md` that routes readers to its code and canonical documentation.

## Flow

```text
AB1 -> decode -> basecalling -> signal_processing -> quality_control
                                                   |
                                                   +-> basecall report v2
                                                   |
FASTA -----------------------------------------> alignment -> variant_calling -> ReadObservation
                                                                                   |
                                                            +----------------------+
                                                            |                      |
                                                            v                      v
                                                   analysis report v7       sample aggregation
                                                                                   |
                                                                                   v
                                                                     sample_evidence/v8
```

`pipeline::observation` is the one authoritative reference-guided read path.
`analyze` consumes one observation; `sample` independently creates one
observation per trace and only then calls `sample::aggregate`.

The shared `checksum` module provides the stable SHA-256 identities used by
`config`, `trace`, and `reference` loading.

## Boundaries

| Module | Owns | Excludes |
|---|---|---|
| `cli` | command syntax | I/O and algorithms |
| `config` | path resolution, strict parsing, validation, caps | per-value environment overrides |
| `error` | typed cross-stage failures | logging and recovery policy |
| `logger` | append-only timestamped operational records | scientific decisions and JSON output |
| `checksum` | shared stable SHA-256 byte identity | file I/O and policy |
| `model` | validated vocabulary and JSON result records | filesystem and algorithms |
| `trace` | canonical ABIF decode | base calling |
| `reference` | one-record FASTA and identity | alignment |
| `basecalling` | basecall peak selection and primary/ambiguity calls inside shared PLOC locus geometry | trimming and reference knowledge |
| `signal_processing` | rolling sample-domain SNR, basecall-independent `LocusEvidence`/`EvidenceProfile`, trace-integrity observations, and merged candidate-noisy regions | channel mutation, calibrated quality, artifact reclassification, basecall classification, reference interpretation, and variant eligibility |
| `quality_control` | penalties, relative scores, end trimming | Phred calibration and variant filtering |
| `alignment` | fixed-point evidence-profile Gotoh scoring, traceback, orientation, circular projection | variant extraction and evidence mutation |
| `variant_calling` | SNV/indel extraction, call/reference mapping, normalization, configured region/supporting-evidence filters | genotype and clinical interpretation |
| `sample` | deterministic read ordering, run-length reference coverage topology, pairwise reference-coordinate overlap admission, sparse differential-locus evidence, reference-oriented preservation of basecall-independent call profiles, and normalized-variant aggregation | input loading, filename/HV pairing, consensus and interpretation |
| `report` | analysis-v7/basecalls-v2/sample-evidence-v8 projection, shared serialization, atomic publish | scientific decisions and alternate/legacy output projection |
| `pipeline` | command sequencing plus shared reference-independent `read` and reference-guided `observation` paths | algorithm internals |

Dependencies point toward `model`, `config`, and `error`; cycles are forbidden. Shared `locus` geometry is reference-free and classification-free. `signal_processing` derives locus profiles from `Chromatogram` channel evidence directly; alignment consumes those immutable profiles for placement without mutating them or the upstream base calls. Existing rolling noisy-window analysis still consumes basecall window records. No algorithm module depends back on signal processing.

## Coordinates and strand

Trace samples, rolling signal-window call indexes, and original call indexes are 0-based. Internal reference intervals are 0-based half-open. Variant positions are 1-based. Reverse alignments retain an explicit oriented-query to original-call mapping. Circular alignments may contain two reference segments when they cross the origin.

## Output projection and transaction

Completed typed scientific results are projected by the report layer and serialized before filesystem publication. Exact public fields, schema versions, and result semantics are owned by the [production contracts](../contracts/README.md), not by architecture.

The core publication boundary is atomic and no-overwrite: a failed core invocation leaves no scientific result and never replaces an existing target. Operational logs remain separate from deterministic scientific result contracts. The detailed publication and logging invariants are canonical under [output/operations invariants](invariants/output-operations.md).

## External batch orchestration

Local corpus discovery, selected cleanup, per-trace execution, and aggregate publication are outside the core scientific architecture. Their authoritative operational contract is [local batch orchestration](../operations/batch.md).

## Resource bounds

Config/FASTA source files are capped at 1/4 MiB before reading, AB1 input at 64 MiB, normalized references at 50,000 bases, indels at 50 changed bases, and Gotoh traceback at 100 million cells. Checked arithmetic rejects an over-limit job before allocation.
