# `src/basecalling/iupac.rs`

## Purpose

Maps a set of qualifying channels to the canonical or two-base IUPAC ambiguity
symbol.

## Responsibilities

- Convert a sorted set of `Nucleotide` bases into a single IUPAC code.
- Return the canonical base for a single channel and the two-base code for a
  pair; return `N` for three or more distinct channels.

## Non-responsibilities

No peak detection, ratio thresholding, or call orchestration.

## Key types and functions

- `code(bases: &[Nucleotide]) -> char`: the IUPAC mapping. It marks which of the
  four channels are present and matches the resulting bit pattern.

## Invariants and errors

- The function is total: any subset of the four channels maps to a valid symbol.
- A single channel maps to its canonical base; exactly two channels map to the
  standard two-base code (`M`, `R`, `W`, `S`, `Y`, `K`); three or more map to `N`.

## Dependencies

- `model::nucleotide` for `Nucleotide` and `channel_index`.

## Biological semantics

IUPAC ambiguity codes represent positions where more than one base is supported
by the dna. Two-base codes (`R`, `Y`, `S`, `W`, `K`, `M`) are the standard
representation of a heterozygous or mixed position in Sanger data.

## Tests

- `maps_two_base_codes`: verifies `A+G -> R`, `C+T -> Y`, and three channels
  map to `N`.

## Status

Implemented.
