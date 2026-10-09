# Configuration and Environment

All core commands (`analyze`, `basecall`, and `sample`) load exactly one strict TOML file selected by `DNA_CONFIG` or `config/dna.toml`. `analyze` and `sample` also load the [target profile](profiles.md) it names; `basecall` does not read the profile. The binary does not parse `.env`, and environment variables do not override individual scientific values.

Unknown keys, missing sections, duplicate TOML keys, unsupported schema versions, non-finite numbers, invalid ranges, and config sources above 1 MiB are errors. No scientific value falls back silently.

## Schema version 7

| Section | Key | Value | Validation |
|---|---|---:|---|
| root | `schema_version` | `7` | exactly 7 |
| | `profile` | `profiles/human-mtdna-rcrs.toml` | non-empty path to a [target profile](profiles.md), resolved against this file's directory when relative |
| `basecalling` | `secondary_peak_ratio` | `0.33` | finite `(0,1]` |
| `signal_processing` | `window_size_bases` | `10` | integer `5..=10` |
| | `minimum_primary_snr` | `3.0` | finite and positive |
| | `minimum_noisy_windows` | `2` | integer at least `2` |
| `callability` | `window_calls` | `16` | integer `8..=64` |
| | `onset_defect_fraction` | `0.375` | finite `(0,1]` |
| | `exit_defect_fraction` | `0.125` | finite `[0,1]`, below `onset_defect_fraction` |
| | `minimum_main_share` | `0.35` | finite `(0,1]` |
| | `maximum_far_share` | `0.12` | finite `[0,1]` |
| | `minimum_shadow_share` | `0.1` | finite `(0,1]` |
| | `weak_amplitude_fraction` | `0.1` | finite `[0,0.5]`; `0` disables the weak defect |
| `quality_control` | `penalty_window_size` | `10` | positive |
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
| | `read_end_margin` | `8` | non-negative integer; `0` disables |
| | `homopolymer_min_length` | `8` | integer at least `2` |
| | `post_homopolymer_window` | `7` | non-negative integer; `0` disables |

For a uniquely strongest basecalling peak, `secondary_peak_ratio` applies both to each channel's selected peak relative to that primary peak and to the channel signal sampled at the primary peak position. Both comparisons are inclusive; this prevents a remote maximum elsewhere in the same PLOC window from qualifying as ambiguity evidence.

`basecall` consumes the basecalling, signal-processing, callability, and quality-control settings; it still validates the complete schema and records the complete configuration checksum. Alignment and variant-calling settings are used by reference-guided operations; sample-reconciliation settings are consumed only by `sample` after every read has completed independent placement. Signal-processing values control observation-only annotations and never change calls, trim bounds, alignments, or variants. Callability values control the read callability view (ADR-0067); while that view is observation-only they change no trim bound, alignment, or variant either, and the double-peak threshold and repeat length reuse `basecalling.secondary_peak_ratio` and `variant_calling.homopolymer_min_length`. Reference topology and reportable regions belong to the target profile. Compact result contracts record the raw configuration checksum and the profile identity but omit method constants and expanded effective values. The sample-reconciliation defaults are Tracy-derived pre-consensus admission controls: they require at least 25 comparable canonical-base positions and at least 0.50 canonical-base agreement for an overlapping read pair to be eligible for later consensus. Gaps, deletions, and unresolved symbols do not enter this nucleotide denominator. Effective values and configuration schema version 7 remain in the strict TOML selected for the run; the local path is omitted.

The read-callability keys (ADR-0062) make a variant ineligible when any of its mapped calls lies within `read_end_margin` calls of either end of the retained interval (`read_end`), or within the `post_homopolymer_window` calls that start at the last call of a run of at least `homopolymer_min_length` identical canonical primary calls in sequencing order (`post_homopolymer`). They use only the read's own calls; no primer or file-name knowledge is involved.

## `.env`

`.env.example` is a shell-tooling template for `DNA_CONFIG=config/dna.toml` and the operational `DNA_LOG_DIR=logs`. `DNA_LOG_DIR` changes only the append-only log destination; it does not alter scientific settings or their checksum. Local `.env` remains ignored. Shells, IDEs, and containers may export these values; DNA itself never reads dotenv files.
