# `src/reference/fasta.rs`

## Purpose

Loads and validates one non-empty plain FASTA record.

## Responsibilities

- Reject an empty or over-4-MiB reference source before reading it, then require UTF-8.
- Parse exactly one record: a `>` header with an identifier and a sequence body.
- Normalize supported DNA symbols (`A/C/G/T/N`) and enforce the length cap.

## Non-responsibilities

No gzip, multi-contig selection, indexing, or alignment.

## Key types and functions

- `load(path, topology) -> Result<Reference>`: the module entry point,
  re-exported from `mod.rs`.

## Invariants and errors

- The file must contain exactly one record; a second `>` line returns
  `Error::Fasta`.
- The sequence must be non-empty and at most `MAX_REFERENCE_LENGTH`; otherwise
  `Error::Fasta`.
- Unsupported symbols return `Error::Fasta`.
- The first line must contain a FASTA identifier; otherwise `Error::Fasta`.
- `sequence_sha256` identifies the normalized, validated sequence.

## Dependencies

- `config` for `MAX_REFERENCE_BYTES` and `MAX_REFERENCE_LENGTH`.
- `checksum` for `hex_sha256`.
- `model::reference` for `Reference` and `ReferenceTopology`.
- `error` for `Error`/`Result`.

## Biological semantics

The reference is the sequence the read is aligned against. Topology (linear or
circular) is supplied by configuration and recorded on the reference.

## Tests

- `rejects_multiple_records`: verifies a second record causes a typed error.

## Status

Implemented.
