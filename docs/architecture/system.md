# DNA Architecture

## Principles

- One strict config and one command-specific JSON result per invocation. Single-read commands accept one AB1; sample evidence accepts one or more AB1 traces. Reference-guided commands additionally require one FASTA.
- Untrusted binary input is checked before every slice, conversion, and allocation.
- Models enforce cardinality and coordinate invariants; scientific functions have no filesystem side effects.
- The CLI and operating-system boundary remain thin.
- Algorithms are deterministic and biologically explicit; known Apollo defects are not compatibility requirements.
- Every Rust source has a same-path manual under `docs/src/`.

## Flow

```text
AB1 -> decode -> basecalling -> signal_processing -> quality_control
                                                   |
                                                   +-> basecall report v2
                                                   |
FASTA -----------------------------------------> alignment -+-> variant_calling -+
                                                           +-> phase -----------+-> ReadObservation
                                                         /          \
                                                        v            v
                                               analysis report v7   SampleEvidence
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
| `report` | analysis-v7/basecalls-v2/sample-evidence-v8 projection, shared serialization, atomic publish | scientific decisions and compatibility output |
| `pipeline` | command sequencing plus shared reference-independent `read` and reference-guided `observation` paths | algorithm internals |

Dependencies point toward `model`, `config`, and `error`; cycles are forbidden. Shared `locus` geometry is reference-free and classification-free. `signal_processing` derives locus profiles from `Chromatogram` channel evidence directly; alignment consumes those immutable profiles for placement without mutating them or the upstream base calls. Existing rolling noisy-window analysis still consumes basecall window records. No algorithm module depends back on signal processing.

## Coordinates and strand

Trace samples, rolling signal-window call indexes, and original call indexes are 0-based. Internal reference intervals are 0-based half-open. Variant positions are 1-based. Reverse alignments retain an explicit oriented-query to original-call mapping. Circular alignments may contain two reference segments when they cross the origin.

## Output projection and transaction

Compact `dna.analysis/v7` projects one completed read observation including concise trace-integrity evidence. `dna.basecalls/v2` projects the shared reference-independent stages. `dna.sample_evidence/v8` projects independently placed reads into one deterministic read registry, run-length coverage/orientation topology, pairwise overlap/admission evidence, sparse differential-locus observations, and normalized variants with factorized read/eligibility/orientation support topology. Read identity/orientation/coverage are factored once at top level; public locus and variant records use unique reviewer-facing read names while internal aggregation remains deterministically SHA-ordered. It omits dense all-reference loci and contains no consensus or sample-level variant verdict. Configuration remains schema version 5, and no compatibility result is assembled.

The completed typed result is serialized before filesystem publication. The core CLI writes a sibling temporary file, flushes and synchronizes it, creates the final path without overwrite, removes the temporary link, and synchronizes the directory. A failed core invocation leaves no command result and never replaces an existing file. Operational logs are deliberately separate, timestamped, run-correlated, escaped to one physical line, and append-only. Pipeline orchestration records aggregate metrics and elapsed time at every stage boundary, each removed variant's kind/position/reasons without alleles, the final warning categories, and stage-aware terminal failures. Mandatory pre-publication records are synchronized before the result transaction begins; no required record is written after a successful publication.

## External batch orchestration

`scripts/analyze_samples.py` is outside the core CLI boundary. It validates the manifest, selected traces, identities, destinations, and cleanup targets; rejects ambiguous matches, trace-stem collisions, and symlinks; then builds or validates the binary before deleting anything. Cleanup removes only selected sample directories plus matching per-trace and sample logs, preserving unselected artifacts. Each selected trace first runs through the one-file no-overwrite CLI in isolation. When every trace for a sample succeeds, the script invokes `signal sample` with that complete trace set and atomically publishes the aggregate as `results/<sample>/<sample>.json` beside the per-trace JSON files. If any trace fails, the aggregate is skipped so no incomplete sample result is presented as complete. Because cleanup is intentionally destructive and execution is sequential, a later failure may leave partial new outputs from earlier successful traces; it does not restore the removed prior batch.

## Resource bounds

Config/FASTA source files are capped at 1/4 MiB before reading, AB1 input at 64 MiB, normalized references at 50,000 bases, indels at 50 changed bases, and Gotoh traceback at 100 million cells. Checked arithmetic rejects an over-limit job before allocation.
