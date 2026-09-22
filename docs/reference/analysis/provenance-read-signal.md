# Analysis Provenance, Read, and Signal Quality

Part of the canonical [analysis contract](README.md).

## Provenance

`provenance` retains the information needed to identify a deterministic run without exposing the trace filename:

- input AB1 `sha256`;
- reference `name`, `topology`, and sequence `sha256`;
- `configuration_sha256`.

Software/build identity, local input/configuration paths, expanded configuration, program constants, timestamps, host data, and method identifiers are not serialized. Software/build provenance is deferred until a stable versioning strategy is defined. Effective scientific settings remain in the strict configuration selected for the run.

## Read and signal-quality summary

`read.call_count` is the number of decoded PLOC call loci. `read.trim.start` and `read.trim.end` delimit the retained calls as a 0-based half-open interval. No sequence string is emitted.

`signal_quality.noisy_regions` contains only merged candidate-noisy regions. Each region has 0-based half-open `calls` and `samples` intervals plus `minimum_primary_snr`. Full-width stride-one windows are still calculated internally by `dna.windowed_snr/v1`, but v7 does not serialize them. The regions remain observational and do not alter trimming, alignment, warning counts, or variant eligibility.

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
