# `src/alignment/mod.rs`

## Purpose

Aligns the retained basecall-independent evidence-profile sequence to the reference using bounded deterministic fixed-point affine-gap Gotoh scoring, while retaining the primary sequence for traceback metrics and downstream variant extraction.

## Responsibilities

- Re-export `align_best` as the module boundary.
- Coordinate fixed-point profile scoring, traceback, score-preserving canonical right-gap placement, and forward/reverse orientation selection.
- Return the selected internal orientation, score, reference segments, metrics,
  and per-column coordinates without duplicate gapped-row strings.

## Non-responsibilities

No FM-index, multi-contig search, hardcoded HV rescue, variant extraction, or
output formatting.

## Key types and functions

- `align_best(qc, signal, reference, config) -> Result<Alignment>`: the public entry point, re-exported from `orient`.
- Child modules: `scoring` (fixed-point profile substitution scores and state ordering), `gotoh` (DP matrices), `traceback` (primary-sequence aligned columns and metrics), `canonical` (score-verified repeat-equivalent 3'/right-most gap placement), and `orient` (profile orientation, strand selection, and coordinate projection).

## Invariants and errors

The complete retained query is consumed; reference flanks may be free. The selected model retains one ordered column representation rather than duplicate row strings. DP tie order and fixed-point scores are explicit and deterministic. Exact forward/reverse profile-score ties fail rather than falling back to primary-basecall match counts. Ambiguous placement or sub-threshold primary-sequence admission metrics return `Error::Alignment`.

## Dependencies

- `config` for `AlignmentConfig` and `MAX_ALIGNMENT_CELLS`.
- `model::alignment`, `model::locus_evidence`, `model::nucleotide`, `model::quality`, `model::reference`, and `model::signal`.
- `error` for `Error`/`Result`.

## Apollo mapping

`apollo/include/apollo/alignment/gotoh.h` and the semi-global configuration used
by Apollo alignment commands.

## Requirements and decisions

ADR-0004, ADR-0029, and ADR-0047; `SRS-ALN-001` through `SRS-ALN-012`.

## Tests

Unit tests in `gotoh` cover free reference flanks and affine gap scoring. The
end-to-end `tests/analyze.rs` integration tests exercise strand selection and
coordinate projection.

## Status

Implemented.
