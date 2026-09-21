# `src/model/alignment.rs`

## Purpose

Defines selected alignment records, orientation, reference segments, and alignment metrics.

## Responsibilities

- Represent forward/reverse orientation relative to the supplied reference and retain the selected fixed-point profile-alignment score internally.
- Project a trace-strand canonical base to the reference strand.
- Project A/C/G/T channel-height order and floating-point signal-channel arrays to reference orientation.
- Represent 0-based half-open mapped reference segments and circular-origin wrap state.

## Orientation projection

Forward evidence is unchanged. Reverse evidence complements base labels and reorders A/C/G/T channel heights and signal-value arrays as raw T/G/C/A into reference-oriented A/C/G/T. Alignment score is an internal fixed-point value in ADR-0029 units and remains omitted from public JSON.

## Status

Implemented.
