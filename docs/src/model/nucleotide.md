# `src/model/nucleotide.rs`

## Purpose

Defines the canonical DNA bases and the complement/reverse-complement operations
used across the pipeline.

## Responsibilities

- Define the `Nucleotide` enum (`A`, `C`, `G`, `T`) with its canonical channel
  order and channel index.
- Provide IUPAC complement and reverse-complement helpers for sequence strings.

## Non-responsibilities

No parsing, validation, or algorithm logic.

## Key types and functions

- `Nucleotide` enum with `ALL`, `as_char`, and `channel_index`.
- `complement_iupac(value: char) -> char`: complements an uppercase IUPAC symbol.
- `reverse_complement(sequence: &str) -> String`: reverses and complements a
  sequence.

## Invariants and errors

- `Nucleotide::ALL` is the canonical `A/C/G/T` order used throughout DNA.
- `complement_iupac` is total: any unrecognized character maps to `N`.
- `reverse_complement` preserves length and operates on uppercase symbols.

## Dependencies

- `serde` for serialization/deserialization of `Nucleotide`.

## Biological semantics

The four canonical bases and their IUPAC complements model DNA strand
complementarity. Reverse complement is used to align reads sequenced on the
opposite strand to the reference.

## Tests

No dedicated unit tests; behavior is exercised through `basecalling` and
`alignment`.

## Status

Implemented.
