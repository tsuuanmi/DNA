# `src/model/trace.rs`

## Purpose

Defines the validated four-channel chromatogram and the optional vendor evidence
retained from the AB1 file.

## Responsibilities

- Represent the decoded A/C/G/T channel samples, base locations, and source
  identity.
- Represent optional vendor primary calls and qualities as evidence only.

## Non-responsibilities

No parsing, decoding, or algorithm logic.

## Key types and functions

- `VendorEvidence`: optional `primary` call string and `qualities`.
- `Chromatogram`: source name, source SHA-256, the four A/C/G/T channel arrays,
  base locations, and vendor evidence.
- `Chromatogram::sample_count()`: number of samples in every channel.
- `Chromatogram::call_count()`: number of vendor-defined base loci.

## Invariants and errors

- All four channels have equal length (enforced at decode time).
- Vendor evidence is optional and is never used as a final algorithm result.

## Dependencies

None beyond the standard library.

## Biological semantics

A chromatogram is the raw signal output of a Sanger sequencing run: four
fluorescence channels sampled over time, with base locations marking where the
ABI basecaller placed each base. Vendor calls are retained for comparison but do
not drive DNA's own re-calling.

## Tests

No dedicated unit tests; behavior is exercised through `trace::decode` and the
integration tests.

## Status

Implemented.
