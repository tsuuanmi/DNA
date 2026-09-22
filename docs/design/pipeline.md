# DNA Pipeline

This document is the canonical orchestration map for the production scientific pipeline. Stage-specific algorithms live in dedicated method documents so a change can load only the relevant method.

## Overview

```text
AB1 + TOML ──► decode ──► basecalling ──► signal_processing ──► quality_control
                                                                     ├─► basecalls/v2
FASTA reference ─────────────────────────────────────────────────────┴─► alignment ─► variant_calling
                                                                                              │
                                                                                       ReadObservation
                                                                                         ├─► analysis/v7
                                                                                         └─► sample aggregation
```

`analyze` and `basecall` consume exactly one AB1 trace. `sample` consumes one or more AB1 traces and processes each independently through the same reference-guided observation path. `analyze` and `sample` additionally consume one single-record FASTA reference; `basecall` performs no reference I/O. Each
stage consumes the validated output of the previous stage and produces a new
typed result; no stage mutates shared state.

## Inputs

- **Trace:** one regular ABIF/AB1 file, decoded into four A/C/G/T signal
  channels, basecall positions (`PLOC.2`), and optional vendor evidence
  (`PBAS.2`, `PCON.2`). `P2BA.1` is ignored. Vendor base strings retain uppercase
  IUPAC symbols, and PCON accepts the ABIF one-byte byte or char representation.
- **Reference (`analyze` and `sample`):** one plain FASTA record of A/C/G/T/N bases, up to 50,000 bases, interpreted as linear or circular per configuration. `basecall` does not accept or load a reference.
- **Configuration:** one strict TOML file selected by `DNA_CONFIG` or
  `config/dna.toml`. Unknown keys, missing sections, and out-of-range values
  are errors.

## Stage methods

| Stage | Method |
|---|---|
| 1 | [ABIF decoding](abif-decoding.md) |
| 2 | [Basecalling](basecalling.md) |
| 3 | [Signal processing](signal-processing/README.md) |
| 4 | [Quality control](quality-control.md) |
| 5 | [Alignment](alignment.md) |
| 6 | [Variant calling](variant-calling.md) |
| sample aggregation | [Sample evidence aggregation](sample-evidence/README.md) |

## One-read observation boundary

After selected alignment and variant calling, DNA materializes a `ReadObservation` that owns the input identity, base calls, basecall-independent locus/signal observations, quality-control result, selected alignment, and read-level variant result for exactly one trace.

The read has already located itself at this boundary. Its orientation and covered reference segments come from evidence-driven semi-global alignment and circular projection; filenames or nominal HV/F/R labels are not placement inputs. This same one-read product feeds both the current analysis report and implemented sample-level reconciliation.

## Output boundary

Public serialization is not defined by the method layer. See [production contracts](../contracts/README.md) for `dna.basecalls/v2`, `dna.analysis/v7`, and `dna.sample_evidence/v8`, including schemas, examples, coordinate semantics, and publication-visible fields.

## Interpretation boundary

The pipeline produces signal evidence, read interpretation, and primary-sequence differences under the [scientific-state invariants](../architecture/invariants/scientific-state.md). It does not turn a single chromatogram into genotype, quantitative heteroplasmy, pathogenicity, or clinical significance.
