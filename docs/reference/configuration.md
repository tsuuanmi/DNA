# Configuration and Environment

All core commands (`analyze`, `basecall`, `sample`, `call`, and `notation`) load exactly one strict TOML file selected by `DNA_CONFIG` or `config/dna.toml`. Every command except `basecall` also loads the [target profile](profiles.md) it names; `basecall` does not read the profile. The binary does not parse `.env`, and environment variables do not override individual scientific values.

Unknown keys, missing sections, duplicate TOML keys, unsupported schema versions, non-finite numbers, invalid ranges, and config sources above 1 MiB are errors. No scientific value falls back silently.

## Schema version 7

| Section | Key | Value | Validation |
|---|---|---:|---|
| root | `schema_version` | `14` | exactly 14 |
| | `profile` | `profiles/human-mtdna-rcrs.toml` | non-empty path to a [target profile](profiles.md), resolved against this file's directory when relative |
| `basecalling` | `secondary_peak_ratio` | `0.33` | finite `(0,1]` |
| `signal_processing` | `window_size_bases` | `10` | integer `5..=10` |
| | `minimum_primary_snr` | `3.0` | finite and positive |
| | `minimum_noisy_windows` | `2` | integer at least `2` |
| | `minimum_noise_sigma` | `1.0` | finite and positive: floor of a channel's noise scale (SRS-SIG-004) |
| | `second_floor` | `0.20` | finite fraction: lowest second-channel threshold on the peak scale, the threshold without a usable run level (SRS-CALL-020) |
| | `run_signal_minimum_reads` | `40` | integer at least `1`: fewest reads a run signal needs to give a run level (PROP-0009) |
| | `run_signal_minimum_calls` | `40` | integer at least `1`: fewest calls of a key, the sample's own taken out, to give a run level |
| | `run_signal_margin` | `0.05` | finite fraction: smallest margin above a key's median |
| | `run_signal_mad_factor` | `4.0` | finite and positive: median absolute deviations above a key's median |
| `callability` | `window_calls` | `16` | integer `8..=64` |
| | `onset_defect_fraction` | `0.375` | finite `(0,1]` |
| | `exit_defect_fraction` | `0.125` | finite `[0,1]`, below `onset_defect_fraction` |
| | `minimum_main_share` | `0.35` | finite `(0,1]` |
| | `maximum_far_share` | `0.12` | finite `[0,1]` |
| | `minimum_shadow_share` | `0.1` | finite `(0,1]` |
| | `weak_amplitude_fraction` | `0.2` | finite `[0,0.5]`; `0` disables the weak defect |
| | `repeat_min_length` | `8` | integer at least `2` |
| | `minimum_callable_calls` | `20` | positive |
| | `repeat_window_calls` | `16` | positive: calls after a repeat run that form its prior window |
| | `double_peak_ratio` | `0.33` | finite `(0,1]`: secondary ratio at or above which a position carries a double peak |
| | `spacing_defect` | `0.5` | finite and positive: relative spacing deviation of a spacing defect |
| | `spacing_half_window` | `8` | positive: spacings on each side whose median anchors a position's spacing deviation |
| `quality_control` | `penalty_window_size` | `10` | positive |
| | `max_relative_quality_score` | `60` | positive `u8` |
| | `trim_context_margin` | `12` | non-negative integer: dephased context calls the trim keeps past the callable span (SRS-QC-004) |
| `alignment` | `match_score` | `3` | positive |
| | `mismatch_score` | `-5` | negative |
| | `ambiguous_score` | `0` | integer |
| | `gap_open_score` | `-10` | negative |
| | `gap_extension_score` | `-4` | negative |
| | `minimum_callable_bases` | `20` | positive |
| | `minimum_identity` | `0.80` | finite `(0,1]` |
| | `max_stretch` | `64` | positive: longest reference stretch a run-structure re-expression spans (SRS-ALN-016) |
| `sample_reconciliation` | `minimum_comparable_bases` | `25` | positive |
| | `minimum_overlap_agreement` | `0.50` | finite `(0,1]` |
| `sanger_evidence` | `minimum_peak_height` | `150` | `1..=32767` |
| `sanger_evidence` | `peak_rescue_quality` | `50` | above `relative_quality_threshold`: a call below the peak floor whose relative quality is at least this still shows its base and indel, when its second base holds less than `peak_rescue_second_share` |
| `sanger_evidence` | `peak_rescue_second_share` | `0.20` | `[0, 1]`: the peak share a rescued call's second base must stay below |
| | `relative_quality_threshold` | `30` | less than `max_relative_quality_score`; comparison is strict `>` |
| | `signal_edge_margin` | `12` | non-negative integer: calls from an end of the usable signal at a signal edge (SRS-CALL-013) |
| | `leading_height_fraction` | `1.0` | finite and positive: smallest primary peak of a leading call, against the read's median (SRS-CALL-014) |
| | `post_repeat_margin` | `12` | non-negative integer (SRS-CALL-016) |
| | `post_repeat_minimum_calls` | `2` | positive (SRS-CALL-016) |
| | `edge_identity_margin` | `12` | non-negative integer: reach of the edge identity rule from an informative end (SRS-CALL-017) |
| | `edge_height_fraction` | `1.0` | finite and positive (SRS-CALL-017) |
| | `start_identity_margin` | `12` | non-negative integer: reach of the start identity rule past the callable span's start (SRS-CALL-018) |
| | `start_norm_calls` | `30` | positive (SRS-CALL-018) |
| | `start_height_fraction` | `0.5` | finite and positive (SRS-CALL-018) |
| | `start_second_share` | `0.15` | finite `[0,1]`: second-channel share below which a start call's signal is clean (SRS-CALL-018) |
| | `edge_second_share` | `0.10` | finite `[0,1]`: second-channel share below which an edge or start call's signal is clean |
| | `edge_width_factor` | `1.5` | finite and positive: widest primary peak, against the median width |
| | `edge_spacing_minimum` | `0.5` | finite and positive, not above `edge_spacing_maximum` |
| | `edge_spacing_maximum` | `1.5` | finite and positive |
| `variant_calling` | `max_indel_length` | `50` | `1..=50` |
| | `read_end_margin` | `12` | non-negative integer; `0` disables |
| `consensus` | `seeded_run` | `4` | integer at least `2`: shortest disputed, unbounded homopolymer grown as one stretch (SRS-CONS-004) |
| `review` | `run_growth` | `3` | positive: a decided run at least this many bases longer than its reference run is a check when its length rests on an anchored end or a single read |
| `review` | `edge_zone` | `20` | positions from a reportable region's first or last position within which a deletion resting on a single read is a check |
| `review` | `phase_run` | `4` | at least 2: equal bases in a row after which a read's in-phase end of an item is no longer compared with the decision for `read_disagrees` |
| | `maximum_undecided_gap` | `20` | non-negative integer: undecided positions kept inside a segment as `N` (SRS-CONS-006) |
| | `minimum_resolved_bases` | `20` | positive: fewest resolved bases of a segment (SRS-CONS-006) |
| | `site_mismatch_cost` | `2` | positive; `site_mismatch_cost + site_gap_extension_cost > site_gap_open_cost` and `site_mismatch_cost < 2 × site_gap_open_cost` |
| | `site_gap_open_cost` | `3` | positive: opening an insertion or deletion when a decided site's call is aligned to its interval |
| | `site_gap_extension_cost` | `2` | positive: extending one |

Each section is owned by one plugin of the static registry
([ADR-0069](../decisions/adr/0069-plugin-first-modality-core-post-calling.md)):
- the Sanger modality plugin owns `basecalling`, `signal_processing`,
  `callability`, `quality_control`, and `sanger_evidence`;
- the core plugin owns `alignment`, `sample_reconciliation`, and
  `variant_calling`;
- the consensus plugin owns `consensus`, which its decided-edits call reads
  too.

The post-calling plugins take their policy from the target profile. Core
sections hold no Sanger key. `sanger_evidence` holds the thresholds that the
Sanger evidence adapter turns into support vetoes, and the parameters of its
read-edge rules. The Sanger plugin's margins (`signal_edge_margin`,
`post_repeat_margin`, `edge_identity_margin`, `start_identity_margin`,
`trim_context_margin`) are its own; the core's `read_end_margin` is read by the
core alone (SRS-CFG-009).

For a uniquely strongest basecalling peak, `secondary_peak_ratio` applies both to each channel's selected peak relative to that primary peak and to the channel signal sampled at the primary peak position. Both comparisons are inclusive; this prevents a remote maximum elsewhere in the same PLOC window from qualifying as ambiguity evidence.

`basecall` consumes the basecalling, signal-processing, callability, and quality-control settings; it still validates the complete schema and records the complete configuration checksum. Alignment and variant-calling settings are used by reference-guided operations; sample-reconciliation settings are consumed only by `sample` after every read has completed independent placement. Signal-processing values control observation-only annotations and never change calls, trim bounds, alignments, or variants. Callability values control the read callability view (ADR-0067), which sets the trim interval, the masked alignment query, variant eligibility, and masked sample observations; the double-peak threshold is `callability.double_peak_ratio`, set apart from `basecalling.secondary_peak_ratio`. Reference topology and reportable regions belong to the target profile. Compact result contracts record the raw configuration checksum and the profile identity but omit method constants and expanded effective values. The sample-reconciliation defaults are Tracy-derived pre-consensus admission controls: they require at least 25 comparable canonical-base positions and at least 0.50 canonical-base agreement for an overlapping read pair to be eligible for later consensus. Gaps, deletions, and unresolved symbols do not enter this nucleotide denominator. Effective values and configuration schema version 14 remain in the strict TOML selected for the run; the local path is omitted.

`read_end_margin` (ADR-0062) makes a variant ineligible when any of its mapped calls lies within that many calls of an uninformative call (`read_end`): beyond the trim interval, which is the callable span plus at most `trim_context_margin` dephased context calls per side, or inside a masked segment that aligns as unresolved. Masked evidence calls add the reason of their phase segment (ADR-0067). These rules use only the read's own calls and signal; no primer or file-name knowledge is involved.

## `.env`

`.env.example` is a shell-tooling template for `DNA_CONFIG=config/dna.toml` and the operational `DNA_LOG_DIR=logs`. `DNA_LOG_DIR` changes only the append-only log destination; it does not alter scientific settings or their checksum. Local `.env` remains ignored. Shells, IDEs, and containers may export these values; DNA itself never reads dotenv files.
