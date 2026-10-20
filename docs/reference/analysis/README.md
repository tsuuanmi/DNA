---
decided-by: [ADR-0026, ADR-0053]
---

# DNA Compact JSON Output

This document covers reference-guided analysis. Reference-free output is the separate [`dna.basecalls/v3` contract](../basecalls.md).

`dna analyze <trace.ab1> --reference <reference.fasta>` writes one deterministic file named `results/<trace-stem>.json`. The `results/` directory is created when publication begins. The core CLI never overwrites an existing result. After validating a non-empty UTF-8 trace stem, Rust separately appends nondeterministic operational records to `$DNA_LOG_DIR/<trace-stem>.log` (default `logs/`); that sidecar is outside the JSON contract.

`dna analyze --evidence <trace-stem>.evidence.json --reference <reference.fasta>` writes the same document from the read's [`dna.read_evidence/v1`](../read-evidence.md#as-an-input) alone, without the trace: byte-identical to the trace's analysis under the same configuration, under the same `results/<trace-stem>.json` and log name (the document's file stem without `.evidence`). The trace and `--evidence` exclude each other. Every Sanger field of the analysis document (`read.call_count`, `read.trim`, `read.callability`, `signal_quality`, the Sanger `warnings` counts, and each variant call's `peaks` and `quality`) is read from the evidence document; nothing is defaulted.

The authoritative contract is [`schemas/analysis-v10.schema.json`](../schemas/analysis-v10.schema.json); a synthetic example is [`examples/analysis-v10.example.json`](../examples/analysis-v10.example.json). DNA emits `dna.analysis/v10` only; v10 adds `calls`, each aligned call's reference column and the core's trust in it ([PROP-0011](../../proposals/0011-step-contracts.md)), so one document shows the chain sanger → align → call for its read. Earlier result versions, compatibility documents, and duplicate legacy fields are not emitted. The strict scientific configuration remains schema version 17.

## Top-level fields

| Field | Meaning |
|---|---|
| `schema_version` | Always `dna.analysis/v10`. |
| `provenance` | Input, reference, configuration, target-profile, and plugin identities. |
| `read` | Original call count, the retained 0-based half-open trim interval, and the read callability view. |
| `signal_quality` | Merged candidate-noisy call/sample regions only. |
| `alignment` | Selected-orientation alignment summary and reference segments. |
| `calls` | Every call of the selected alignment, in call order, with its reference column and the core's trust in it. |
| `variants` | Normalized primary-sequence differences with concise mapped calls. |
| `warnings` | Counts of unresolved primary calls, multi-channel unresolved calls, and excluded variant candidates. |

All objects are closed by the schema. Compact v10 deliberately omits trace filenames, full primary/ambiguity/retained sequences, individual rolling windows, per-position callability features and masks, gapped alignment rows (each call's column is in `calls`), operation runs, alignment score and redundant match counts, method constants, selected per-channel peak position/source objects, vendor PBAS/PCON data, variant contig/classification/normalization labels, warning totals, and duplicated origin-wrap or vendor-disagreement fields.

## Contract sections

- [Provenance, read, and signal quality](provenance-read-signal.md)
- [Alignment, aligned calls, and variants](alignment-variants.md)
- [Coordinate contract](../coordinates.md)

The shared coordinate contract is authoritative for coordinate bases, interval semantics, reverse orientation, inserted/deleted evidence, normalization, and circular projection. Analysis-specific fields follow that contract rather than redefining it here.

## Interpretation boundary

Variant alleles and reviewer-facing base/peak labels are written on the supplied reference strand. Public `quality` remains an uncalibrated relative score; neither channel height nor quality implies genotype, zygosity, allele fraction, heteroplasmy, or clinical significance.
