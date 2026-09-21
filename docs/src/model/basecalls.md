# `src/model/basecalls.rs`

## Purpose

Defines signal-derived base calls and peak evidence.

## Responsibilities

- Retain selected per-channel peaks, optional primary-event evidence, primary/ambiguity calls, qualifying channels, and vendor-agreement state.
- Retain one shared-position `PrimaryPeakEvidence` containing raw analyzed A/C/G/T channel heights at the uniquely strongest primary-event coordinate.

## Evidence distinction

Selected per-channel peaks are the strongest events for each channel within the call window and may occur at different sample positions. `PrimaryPeakEvidence.channel_heights` samples all four channels at one shared primary-event coordinate.

Reviewer-facing variant output uses the latter because it answers the direct question: what were A/C/G/T signals at the call that produced this base? Peak positions, sources, and broader per-window objects remain internal.

## Status

Implemented.
