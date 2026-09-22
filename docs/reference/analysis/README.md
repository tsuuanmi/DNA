# DNA Compact JSON Output

This document covers reference-guided analysis. Reference-free output is the separate [`dna.basecalls/v2` contract](../basecalls.md).

`signal analyze <trace.ab1> --reference <reference.fasta>` writes one deterministic file named `results/<trace-stem>.json`. The `results/` directory is created when publication begins. The core CLI never overwrites an existing result. After validating a non-empty UTF-8 trace stem, Rust separately appends nondeterministic operational records to `$DNA_LOG_DIR/<trace-stem>.log` (default `logs/`); that sidecar is outside the JSON contract.

The authoritative contract is [`schemas/analysis-v7.schema.json`](../schemas/analysis-v7.schema.json); a synthetic example is [`examples/analysis-v7.example.json`](../examples/analysis-v7.example.json). DNA emits `dna.analysis/v7` only. Earlier result versions, compatibility documents, and duplicate legacy fields are not emitted. The strict scientific configuration remains schema version 5.

## Top-level fields

| Field | Meaning |
|---|---|
| `schema_version` | Always `dna.analysis/v7`. |
| `provenance` | Input, reference, and configuration identities. |
| `read` | Original call count and the retained 0-based half-open trim interval. |
| `signal_quality` | Merged candidate-noisy call/sample regions only. |
| `alignment` | Selected-orientation alignment summary and reference segments. |
| `variants` | Normalized primary-sequence differences with concise mapped calls. |
| `warnings` | Counts of unresolved primary calls, multi-channel unresolved calls, and excluded variant candidates. |

All objects are closed by the schema. Compact v7 deliberately omits trace filenames, full primary/ambiguity/retained sequences, individual rolling windows, gapped alignment rows, operation runs, alignment score and redundant match counts, method constants, selected per-channel peak position/source objects, vendor PBAS/PCON data, variant contig/classification/normalization labels, warning totals, and duplicated origin-wrap or vendor-disagreement fields.

## Contract sections

- [Provenance, read, and signal quality](provenance-read-signal.md)
- [Alignment and variants](alignment-variants.md)
- [Coordinate contract](../coordinates.md)

The shared coordinate contract is authoritative for coordinate bases, interval semantics, reverse orientation, inserted/deleted evidence, normalization, and circular projection. Analysis-specific fields follow that contract rather than redefining it here.

## Interpretation boundary

Variant alleles and reviewer-facing base/peak labels are written on the supplied reference strand. Public `quality` remains an uncalibrated relative score; neither channel height nor quality implies genotype, zygosity, allele fraction, heteroplasmy, or clinical significance.
