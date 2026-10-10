# Analysis Provenance, Read, and Signal Quality

Part of the canonical [analysis contract](README.md).

## Provenance

`provenance` retains the information needed to identify a deterministic run without exposing the trace filename:

- input AB1 `sha256`;
- reference `name`, `topology`, and sequence `sha256`;
- `configuration_sha256`;
- `profile`: the [target profile](../profiles.md) `id` and the `sha256` of the profile file;
- `plugins`: the plugins that produced the document ([ADR-0069](../../decisions/adr/0069-plugin-first-modality-core-post-calling.md)), in execution order, each with its `id`, `family` (`modality`, `core`, or `post_calling`), and method `version` ([versioning](../../governance/versioning.md#plugins)). `analyze` runs `sanger`, then `core`.

Software/build identity, local input/configuration paths, expanded configuration, program constants, timestamps, and host data are not serialized; plugin method versions are the only method identifiers. Software/build provenance is deferred until a stable versioning strategy is defined. Effective scientific settings remain in the strict configuration selected for the run.

## Read and signal-quality summary

The v9 serialized integrity object retains `ploc_*` field names because those
names are part of the current closed JSON contract. Internally, `PLOC.2` is
projected to canonical Sanger locus positions at the ABIF boundary.

`read.call_count` is the number of canonical Sanger call loci decoded from `PLOC.2`, after merging repeated positions. `read.trim.start` and `read.trim.end` delimit the retained calls as a 0-based half-open interval: the callable span plus at most `read_end_margin` adjacent dephased calls per side ([quality control](../../design/quality-control.md)). No sequence string is emitted.

`read.callability` is the `dna.read_callability/v1` view ([method](../../design/callability.md)): `callable_span` (0-based half-open calls from the first to the last unmasked call; empty when every call is masked), `segments` (ordered 0-based half-open call intervals that partition the read, each with `state` `in_phase`, `dephased`, `mixed`, `weak`, or `irregular`, `after_repeat`, true when the segment starts in the window after a long repeat run, and, exactly on `dephased` segments, `shadow_offsets`: the ascending call offsets whose shadow the segment carries, always including `-1` or `1`), and `masked_calls`. The view sets the trim interval, the masked alignment query, and variant eligibility (ADR-0067). Phase states are read-level signal states, not error probabilities or artifact classes.

`signal_quality.noisy_regions` contains only merged candidate-noisy regions. Each region has 0-based half-open `calls` and `samples` intervals plus `minimum_primary_snr`. Full-width stride-one windows are still calculated internally by `dna.windowed_snr/v1`, but v9 does not serialize them. The regions remain observational and do not alter trimming, alignment, warning counts, or variant eligibility.

## Trace integrity

`signal_quality.integrity` preserves concise evidence about the trace foundation:

- serialized `ploc_count`, representing canonical locus count;
- `duplicate_ploc_positions`: `PLOC.2` positions equal to their predecessor,
  merged into that locus with their vendor entries;
- optional PBAS/PCON counts;
- serialized `*_ploc_spacing`, representing adjacent canonical locus spacing when at least two loci exist;
- exact signed-16-bit clipped channel-sample count;
- optional maximum-to-median corrected event-signal ratio.

A PBAS/PCON length mismatch is non-fatal and does not create/remove calls:
DNA still processes exactly the valid canonical loci decoded from `PLOC.2`. Exact clipping and
event-signal imbalance are observations only and do not change calls, trim,
alignment, or variants. The ratio is not an artifact probability or dye-blob
classification.
