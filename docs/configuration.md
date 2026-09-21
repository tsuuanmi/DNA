# Configuration and Environment

All core commands (`analyze`, `basecall`, and `sample`) load exactly one strict TOML file selected by `DNA_CONFIG` or `config/dna.toml`. The binary does not parse `.env`, and environment variables do not override individual scientific values.

Unknown keys, missing sections, duplicate TOML keys, unsupported schema versions, non-finite numbers, invalid ranges, and config sources above 1 MiB are errors. No scientific value falls back silently.

## Schema version 5

| Section | Key | Value | Validation |
|---|---|---:|---|
| root | `schema_version` | `5` | exactly 5 |
| `reference` | `topology` | `circular` | `linear` or `circular` |
| `basecalling` | `secondary_peak_ratio` | `0.33` | finite `(0,1]` |
| `signal_processing` | `window_size_bases` | `10` | integer `5..=10` |
| | `minimum_primary_snr` | `3.0` | finite and positive |
| | `minimum_noisy_windows` | `2` | integer at least `2` |
| `quality_control` | `trim_window_size` | `10` | positive |
| | `best_section_fraction` | `0.10` | finite `(0,1]` |
| | `max_relative_quality_score` | `60` | positive `u8` |
| | `trim_stringency` | `7.0` | finite `[0,9]` |
| | `minimum_retained_bases` | `20` | positive |
| `alignment` | `match_score` | `3` | positive |
| | `mismatch_score` | `-5` | negative |
| | `ambiguous_score` | `0` | integer |
| | `gap_open_score` | `-10` | negative |
| | `gap_extension_score` | `-4` | negative |
| | `minimum_callable_bases` | `20` | positive |
| | `minimum_identity` | `0.80` | finite `(0,1]` |
| `sample_reconciliation` | `minimum_comparable_bases` | `25` | positive |
| | `minimum_overlap_agreement` | `0.50` | finite `(0,1]` |
| `variant_calling` | `max_indel_length` | `50` | `1..=50` |
| | `minimum_peak_height` | `150` | `1..=32767` |
| | `relative_quality_threshold` | `30` | less than `max_relative_quality_score`; comparison is strict `>` |
| | `regions` | `[[16024, 16365], [73, 340], [438, 576]]` | non-empty inclusive 1-based ranges within `1..=50000` |

For a uniquely strongest basecalling peak, `secondary_peak_ratio` applies both to each channel's selected peak relative to that primary peak and to the channel signal sampled at the primary peak position. Both comparisons are inclusive; this prevents a remote maximum elsewhere in the same PLOC window from qualifying as ambiguity evidence.

`basecall` consumes the basecalling, signal-processing, and quality-control settings; it still validates the complete schema and records the complete configuration checksum. Reference, alignment, and variant-calling settings are used by reference-guided operations; sample-reconciliation settings are consumed only by `sample` after every read has completed independent placement. DNA-processing values control observation-only annotations and never change calls, trim bounds, alignments, or variants. The region list is treated as a union in the supplied reference coordinate system. Region order and overlap do not change eligibility. Compact result contracts record the raw configuration checksum but omit method constants and expanded effective values. The sample-reconciliation defaults are Tracy-derived pre-consensus admission controls: they require at least 25 comparable canonical-base positions and at least 0.50 canonical-base agreement for an overlapping read pair to be eligible for later consensus. Gaps, deletions, and unresolved symbols do not enter this nucleotide denominator. Effective values and configuration schema version 5 remain in the strict TOML selected for the run; the local path is omitted.

## `.env`

`.env.example` is a shell-tooling template for `DNA_CONFIG=config/dna.toml` and the operational `DNA_LOG_DIR=logs`. `DNA_LOG_DIR` changes only the append-only log destination; it does not alter scientific settings or their checksum. Local `.env` remains ignored. Shells, IDEs, and containers may export these values; DNA itself never reads dotenv files.
