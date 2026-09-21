# `src/basecalling/peak.rs`

## Purpose

Selects the strongest positive local peak per channel inside shared validated PLOC locus windows.

## Responsibilities

- Reuse shared symmetric neighboring-midpoint half-open locus windows.
- For each channel, find the strongest positive local maximum inside the window,
  falling back to the PLOC sample when no positive local maximum exists.
- Return each channel's base, selected height, 0-based sample position, and
  selection source.

## Non-responsibilities

No ratio thresholding, IUPAC mapping, or call orchestration.

## Key types and functions

- `LocusWindow`: shared half-open `[start, end)` sample geometry from `src/locus.rs`.
- `windows(trace) -> Result<Vec<LocusWindow>>`: maps shared geometry errors into the basecalling stage.
- `peaks(trace, window, ploc) -> [ChannelPeak; 4]`: returns the selected peak for
  each of the four channels.
- `midpoint(left, right) -> Result<usize>`: checked midpoint arithmetic.

## Invariants and errors

- At least two PLOC positions are required; otherwise `Error::Basecalling`.
- Every window must satisfy `start < end`, `end <= sample_count`, and contain its
  PLOC position; violations return `Error::Basecalling`.
- Shared locus-window geometry errors are mapped to `Error::Basecalling`.
- A peak is a positive local maximum; when none exists, the PLOC sample is used
  and recorded as `PeakSource::PlocFallback` at the PLOC position.
- Every selected peak position remains within its call window; call orchestration
  rejects a violated invariant.

## Dependencies

- `locus` for authoritative PLOC window geometry.
- `model::basecalls` for `ChannelPeak` and `PeakSource`.
- `model::nucleotide` for `Nucleotide`.
- `model::trace` for `Chromatogram`.
- `error` for `Error`/`Result`.

## Biological semantics

Each PLOC locus is a vendor-identified base position. The window around it
captures the local signal for that base. Selecting the strongest positive local
maximum per channel recovers the peak height used to rank channels and apply the
secondary-peak threshold. Call orchestration also uses the uniquely strongest
selected peak's coordinate to sample every channel at the same primary event;
each channel's selected coordinate remains retained as internal evidence.

## Tests

- `selects_each_channel_independently_and_falls_back_to_ploc`: verifies independent local maxima and per-channel PLOC fallback.

## Status

Implemented.
