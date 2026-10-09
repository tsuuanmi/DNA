# DNA Pipeline

This document is the canonical orchestration map for the production scientific pipeline. Stage-specific algorithms live in dedicated method documents so a change can load only the relevant method.

## Overview

```text
AB1 + TOML ──► decode ──► basecalling ──► signal_processing ──► callability ──► quality_control
                                                                                      ├─► basecalls/v3
FASTA reference ──────────────────────────────────────────────────────────────────────┴─► alignment ─► variant_calling
                                                                                                               │
                                                                                                        ReadObservation
                                                                                                          ├─► analysis/v9
                                                                                                          └─► sample aggregation

consensus FASTA + TOML ──► input::sequence ──► ReadEvidence ──► alignment ─► variant_calling ──► CalledRead ──► variants/v1

variants/v1 + TOML + FASTA ──► input::variants ──► normalization ─► nomenclature ─► conformance ──► notation/v1
```

`analyze` and `basecall` consume exactly one AB1 trace. `sample` consumes one or more AB1 traces and processes each independently through the same reference-guided observation path. `analyze` and `sample` additionally consume one single-record FASTA reference; `basecall` performs no reference I/O. `call` consumes one or more FASTA files of reviewed consensus sequences plus the reference and runs the core alone: each record becomes `ReadEvidence` through the sequence modality, then follows the same `dna_core::read_call` path that Sanger reads take after their evidence adapter ([ADR-0069](../decisions/adr/0069-plugin-first-modality-core-post-calling.md)). `notation` reads one variants document and runs only the post-calling plugins: the per-read representation of [variant nomenclature](variant-nomenclature.md#notation-command-and-conformance), then the profile's conformance checks. Each
stage consumes the validated output of the previous stage and produces a new
typed result; no stage mutates shared state.

## Inputs

- **Trace:** one regular ABIF/AB1 file, decoded into four A/C/G/T signal
  channels, basecall positions (`PLOC.2`), and optional vendor evidence
  (`PBAS.2`, `PCON.2`). `P2BA.1` is ignored. Vendor base strings retain uppercase
  IUPAC symbols, and PCON accepts the ABIF one-byte byte or char representation.
- **Reference (`analyze`, `sample`, `call`, and `notation`):** one plain FASTA record of A/C/G/T/N bases, up to 50,000 bases, loaded with the target profile's linear or circular topology. `basecall` does not accept or load a reference.
- **Configuration:** one strict TOML file selected by `DNA_CONFIG` or
  `config/dna.toml`. Unknown keys, missing sections, and out-of-range values
  are errors.

## Stage methods

| Stage | Method |
|---|---|
| 1 | [ABIF decoding](abif-decoding.md) |
| 2 | [Basecalling](basecalling.md) |
| 3 | [Signal processing](signal-processing/README.md) |
| 4 | [Callability](callability.md) |
| 5 | [Quality control](quality-control.md) |
| 6 | [Alignment](alignment.md) |
| 7 | [Variant calling](variant-calling.md) |
| sample aggregation | [Sample evidence aggregation](sample-evidence/README.md) |

## Ownership boundaries

The stage sequence is shared across multiple application surfaces, but ownership
does not follow the CLI command tree:

- `input` owns source validation and loading of the configuration, target
  profile, and reference; `input::sanger` decodes traces through
  `dna_sanger::abif`;
- `dna_sanger::read_processing` owns the shared reference-free
  basecalling/signal/callability/QC path used by basecall and Variant Analysis;
- `variant_analysis` owns the reference-guided one-read observation path,
  including alignment and variant calling;
- `pipeline` owns CLI/sample orchestration, application filesystem naming and
  overwrite checks, operational logging lifecycle, report projection, and
  publication.

This direction keeps reusable scientific capabilities independent of the
command-line orchestration layer while preserving one authoritative scientific
implementation.

## One-read observation boundary

After selected alignment and variant calling, DNA materializes a `ReadObservation` for exactly one trace. It has two parts ([ADR-0069](../decisions/adr/0069-plugin-first-modality-core-post-calling.md)):

- `CalledRead` is the core's modality-neutral part: input identity, the `ReadEvidence` the core consumed, the selected alignment, and the read-level variant result. Sample aggregation reads only this part.
- `SangerAttachment` holds the base calls, basecall-independent locus/signal observations, read callability, and quality-control result. Reports join it to core records by read identity and call index.

The read has already located itself at this boundary. Its orientation and covered reference segments come from evidence-driven semi-global alignment and circular projection; filenames or nominal HV/F/R labels are not placement inputs. This same one-read product feeds both the current analysis report and implemented sample-level reconciliation.

## Output boundary

Public serialization is not defined by the method layer. See [production contracts](../reference/README.md) for `dna.basecalls/v3`, `dna.analysis/v9`, `dna.sample_evidence/v10`, `dna.variants/v1`, and `dna.notation/v1`, including schemas, examples, coordinate semantics, and publication-visible fields.

## Interpretation boundary

The pipeline produces signal evidence, read interpretation, and primary-sequence differences under the [scientific-state invariants](../architecture/invariants/scientific-state.md). It does not turn a single chromatogram into genotype, quantitative heteroplasmy, pathogenicity, or clinical significance.
