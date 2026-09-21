# DNA Compact JSON Output

This document covers reference-guided analysis. Reference-free output is the separate [`dna.basecalls/v2` contract](basecalls.md).

`signal analyze <trace.ab1> --reference <reference.fasta>` writes one deterministic file named `results/<trace-stem>.json`. The `results/` directory is created when publication begins. The core CLI never overwrites an existing result. After validating a non-empty UTF-8 trace stem, Rust separately appends nondeterministic operational records to `$DNA_LOG_DIR/<trace-stem>.log` (default `logs/`); that sidecar is outside the JSON contract.

The authoritative contract is [`schemas/analysis-v7.schema.json`](schemas/analysis-v7.schema.json); a synthetic example is [`examples/analysis-v7.example.json`](examples/analysis-v7.example.json). Output v7 is intentionally incompatible with earlier result versions, and DNA emits no compatibility document or duplicate legacy fields. The strict scientific configuration remains schema version 5.

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

## Provenance

`provenance` retains the information needed to identify a deterministic run without exposing the trace filename:

- input AB1 `sha256`;
- reference `name`, `topology`, and sequence `sha256`;
- `configuration_sha256`.

Software/build identity, local input/configuration paths, expanded configuration, program constants, timestamps, host data, and method identifiers are not serialized. Software/build provenance is deferred until a stable versioning strategy is defined. Effective scientific settings remain in the strict configuration selected for the run.

## Read and signal-quality summary

`read.call_count` is the number of decoded PLOC call loci. `read.trim.start` and `read.trim.end` delimit the retained calls as a 0-based half-open interval. No sequence string is emitted.

`signal_quality.noisy_regions` contains only merged candidate-noisy regions. Each region has 0-based half-open `calls` and `samples` intervals plus `minimum_primary_snr`. Full-width stride-one windows are still calculated internally by `dna.windowed_snr/v1`, but v7 does not serialize them. The regions remain observational and do not alter trimming, alignment, warning counts, or variant eligibility.

## Alignment summary

The selected alignment reports:

- `orientation` (`forward` or `reverse`);
- `callable_bases` and callable `identity`;
- `unresolved_bases` and `gap_opens`;
- one or two 0-based half-open `reference_segments`;
- `wraps_origin`.

Gapped query/reference rows, operation runs, fixed-point profile score, exact-match/mismatch redundancy, and traceback columns remain internal. Placement consumes the retained post-trim `EvidenceProfile` sequence, while the retained primary sequence remains attached to traceback columns for callable/identity metrics and downstream primary-sequence variant extraction. `reference_segments` therefore describe post-trim mapped coverage.

## Variant calls

Each normalized variant contains only `position`, `reference`, `alternate`,
`kind`, and direct `calls`. The variant position is 1-based on the supplied
reference strand.

Each call is optimized for scientific review rather than implementation
traceability:

- `role`: `supporting` or `flanking`;
- `base`: called base projected onto the reference strand;
- `peaks`: raw analyzed A/C/G/T channel heights sampled together at the unique
  primary-event coordinate and projected to reference orientation;
- `quality`: the existing uncalibrated relative score under a concise public
  field name.

Original call indexes, ABIF PLOC coordinates, mapped call positions, selected-peak
positions/sources, penalties, calibration flags, and vendor scores remain internal.

For reverse reads, both `base` and `peaks` are reference-oriented, so the
reviewer can compare the normalized alternate allele directly with the signal
channels.

### SNVs

An SNV contains one or more supporting calls. The top-level variant `position`
already identifies the biological coordinate, so that coordinate is not repeated
inside each call.

### Insertions

An insertion contains one supporting call per inserted base plus available
flanking calls. Supporting-call bases and peaks show the observed inserted
sequence directly.

### Deletions

A deletion has no trace signal at the deleted reference base. It therefore
contains available flanking calls with their real peak/quality evidence and never
fabricates deleted-base dna.

For repeat-associated indels, normalization can move the reported allele
representation away from the observed alignment gap. The normalized
`position/reference/alternate` remain authoritative for the biological variant.


## Coordinate conventions

| Field | Coordinate system |
|---|---|
| variant `position` | 1-based biological reference coordinate |
| trim, segment, noisy-region `start`/`end` | 0-based half-open interval `[start, end)` |

Variant alleles, call `base`, and peak labels are written on the supplied reference strand. `quality` remains an uncalibrated relative score; neither channel height nor quality implies genotype, zygosity, allele fraction, heteroplasmy, or clinical significance.

## Trace integrity

`signal_quality.integrity` preserves concise evidence about the trace foundation:

- PLOC count;
- optional PBAS/PCON counts;
- minimum/median/maximum adjacent PLOC spacing when at least two loci exist;
- exact signed-16-bit clipped channel-sample count;
- optional maximum-to-median corrected event-signal ratio.

A PBAS/PCON length mismatch is non-fatal and does not create/remove calls:
DNA still processes exactly the valid PLOC-defined loci. Exact clipping and
event-signal imbalance are observations only and do not change calls, trim,
alignment, or variants. The ratio is not an artifact probability or dye-blob
classification.
