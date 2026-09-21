# `src/model/result.rs`

## Purpose

Defines compact serializable `dna.analysis/v7` records and shared public result types.

## Responsibilities

- Represent provenance, read/trim summary, signal-quality summary, trace-integrity evidence, alignment summary, normalized variants, and warnings.
- Represent reviewer-facing variant-call evidence as `role + base + peaks + quality`.
- Represent co-located A/C/G/T primary-event heights with stable `A/C/G/T` JSON keys.

## Variant call contract

`VariantCallResult` intentionally omits original call index, mapped call position, ABIF PLOC, trace-strand symbols, and peak positions/sources.

`base` and `peaks` are projected to the reference strand. `quality` is the existing uncalibrated relative score under a concise public name.

Deletion calls contain real flank evidence only; no deleted-base signal is fabricated.

## Coordinates

Variant `position` is 1-based. Trim/reference/noisy-region intervals remain 0-based half-open.

## Status

Implemented.
