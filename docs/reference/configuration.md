---
implements: [SRS-ALN, SRS-CALL, SRS-CFG, SRS-CONS, SRS-QC, SRS-SIG]
decided-by: [ADR-0062, ADR-0067, ADR-0069, ADR-0084]
covers:
  - config/dna.toml
  - crates/dna-kernel/src/bounds.rs
  - crates/dna-sanger/src/read_processing/config.rs
  - crates/dna-sanger/src/tract/config.rs
  - src/config/
---

# Configuration and Environment

Every command loads exactly one strict TOML file selected by `DNA_CONFIG` or `config/dna.toml`. Every command except `basecall`, `evidence`, and `signal` also loads the [[reference/profiles|target profile]] it names; those three read no reference and no profile. The binary does not parse `.env`, and environment variables do not override individual scientific values.

Unknown keys, missing sections, duplicate TOML keys, unsupported schema versions, non-finite numbers, invalid ranges, and config sources above 1 MiB are errors. No scientific value falls back silently.

## Schema version 24

| Section | Key | Value | Validation |
|---|---|---:|---|
| root | `schema_version` | `25` | exactly 25 |
| | `profile` | `profiles/human-mtdna-rcrs.toml` | non-empty path to a [target profile](profiles.md), resolved against this file's directory when relative |
| `basecalling` | `secondary_peak_ratio` | `0.33` | finite `(0,1]` |
| `signal_processing` | `window_size_bases` | `10` | integer `5..=10` |
| | `minimum_primary_snr` | `3.0` | finite and positive |
| | `minimum_noisy_windows` | `2` | integer at least `2` |
| | `minimum_noise_sigma` | `1.0` | finite and positive: floor of a channel's noise scale (SRS-SIG-004) |
| | `second_floor` | `0.20` | finite fraction: lowest second-channel threshold on the peak scale, the threshold without a usable run level (SRS-CALL-020) |
| | `run_signal_flank` | `2` | integer in `1..=4`: calls on each side of a call that form its run signal key's context (a 5-mer at `2`); a run signal made with another flank is refused |
| | `run_signal_resolution` | `1000` | integer in `100..=10000`: steps of a second-channel share from 0 to 1 that a run signal counts; a run signal made with another resolution is refused |
| | `run_signal_minimum_reads` | `40` | integer at least `1`: fewest reads a run signal needs to give a run level (PROP-0009) |
| | `run_signal_minimum_calls` | `40` | integer at least `1`: fewest calls of a key, the sample's own taken out, to give a run level |
| | `run_signal_margin` | `0.05` | finite fraction: smallest margin above a key's median |
| | `run_signal_mad_factor` | `4.0` | finite and positive: median absolute deviations above a key's median |
| | `second_cause_noise_ratio` | `3.0` | finite and positive: below this multiple of the call's lowest channel a second peak is `noise` (SRS-CALL-021) |
| | `second_cause_tail_offset` | `0.4` | finite and positive: offset from the call's centre, in call spacings, from which a second peak is a `neighbour_tail` or `off_centre` |
| | `second_cause_blob_width` | `1.8` | finite and positive: multiple of the neighbouring calls' median primary width from which a centred second peak is a `dye_blob` |
| | `decompose_unexplained_quantile` | `0.99` | finite, in `(0, 1)`: nearest-rank quantile of a read's own in-phase residuals above which a call of a phase mixture is unexplained (SRS-CALL-023) |
| | `decompose_minimum_calls` | `20` | positive: in-phase single-channel calls a read needs to measure its own noise before any phase mixture is fitted (SRS-CALL-023) |
| | `decompose_onset_window` | `16` | positive: calls from a dephased segment's start within which the run that opens its phase mixture ends (SRS-CALL-023) |
| | `second_cause_width_calls` | `6` | integer at least `1`: calls on each side whose primary widths give that median |
| | `spike_height` | `3.0` | finite, at least `1`: multiple of the local level, the median primary peak of the neighbouring calls, that every channel's apex of an electrophoretic spike reaches (SRS-SIG-012) |
| | `spike_span` | `1.5` | finite and positive: median call spacings within which a spike's four apexes lie |
| | `spike_level_calls` | `6` | integer at least `1`: calls on each side whose primary peaks give the local level |
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
| | `edge_identity_margin` | `12` | non-negative integer: reach of the edge identity rule from a stop of the read's in-phase signal (SRS-CALL-017) |
| | `edge_height_fraction` | `1.0` | finite and positive: smallest primary peak at an informative end, against the read's median (SRS-CALL-017) |
| | `end_norm_calls` | `30` | positive: last scored single-channel calls before a stop (the informative end, or a post-repeat phase loss inside it) whose medians judge the calls before the run that closes it (SRS-CALL-017) |
| | `end_height_fraction` | `0.3` | finite and positive: smallest primary peak of such a call, against those medians (SRS-CALL-017) |
| | `start_identity_margin` | `12` | non-negative integer: reach of the start identity rule past the callable span's start (SRS-CALL-018) |
| | `poor_start_calls` | `12` | positive: consecutive sound calls that end a read's poor start (SRS-CALL-024) |
| | `start_norm_calls` | `30` | positive (SRS-CALL-018) |
| | `start_height_fraction` | `0.5` | finite and positive (SRS-CALL-018) |
| | `start_second_share` | `0.15` | finite `[0,1]`: second-channel share below which a start call's signal is clean (SRS-CALL-018) |
| | `edge_second_share` | `0.10` | finite `[0,1]`: second-channel share below which an edge or start call's signal is clean |
| | `edge_width_factor` | `1.5` | finite and positive: widest primary peak, against the median width |
| | `edge_spacing_minimum` | `0.5` | finite and positive, not above `edge_spacing_maximum` |
| | `edge_spacing_maximum` | `1.5` | finite and positive |
| `trace` | `locate_min_run` | `4` | positive: fewest equal primary calls that open a tract before a call where the read leaves phase (SRS-CALL-025) |
| | `locate_min_before` | `2` | non-negative integer: fewest run-base calls before a single-call mark that join the run before it to the tract |
| | `span_below` | `10` | non-negative integer: a candidate total length starts this many calls below the run extent |
| | `span_above` | `3` | non-negative integer: and ends this many above it |
| | `min_length` | `3` | positive: shortest candidate total length |
| | `downstream_calls` | `60` | positive: calls after the longest candidate the decode fits |
| | `joint_lengths` | `6` | integer `1..=8`, a computation cap: most consecutive lengths whose downstream sequence is decoded jointly (4^K windows per call); the `tracts_decoded` log counts the tracts where it binds |
| | `iterations` | `6` | positive, a computation cap: most rounds of fit and decode per start; the `tracts_decoded` log counts the tracts where it binds |
| | `signature_calls` | `200` | positive: most calls before the tract that the read's base signatures and leakage null read |
| | `signature_gap` | `3` | non-negative integer: calls just before the tract that they leave out |
| | `leak_min_run` | `3` | positive: the leakage null leaves out runs of at least this many equal primaries |
| | `leak_pad` | `2` | non-negative integer: and this many calls on each side of them |
| | `leak_min_calls` | `10` | positive: fewest calls the leakage null needs |
| | `near_equal_z` | `2.0` | finite and positive: a marginal's top two lengths are near-equal when their margin, in standard errors deflated by the read's residual autocorrelation, is below it; a margin that is not finite is undefined, never clear |
| `variant_calling` | `max_indel_length` | `50` | `1..=50` |
| | `read_end_margin` | `12` | non-negative integer, at most the longest supported reference (50,000); `0` disables |
| `consensus` | `seeded_run` | `4` | integer at least `2`: shortest disputed, unbounded homopolymer grown as one stretch (SRS-CONS-004) |
| `review` | `run_growth` | `3` | positive: a decided run at least this many bases longer than its reference run is a check when its length rests on an anchored end or on one strand, however many reads |
| `review` | `edge_zone` | `40` | positions from a reportable region's end that the call does not analyse (SRS-CONS-012) within which any difference is a check (`beside_not_analysed`) |
| `review` | `deletion_zone` | `20` | positions from a reportable region's first or last position within which a deletion resting on one strand is a check (`edge_deletion`) |
| `review` | `phase_run` | `4` | at least 2: equal bases in a row after which a read's in-phase end of an item is no longer compared with the decision for `read_disagrees` |
| `consensus` | `maximum_undecided_gap` | `20` | non-negative integer, at most the longest supported reference (50,000): undecided positions kept inside a segment as `N` (SRS-CONS-006) |
| | `minimum_resolved_bases` | `20` | positive: fewest resolved bases of a segment (SRS-CONS-006) |
| | `site_mismatch_cost` | `2` | positive; `site_mismatch_cost + site_gap_extension_cost > site_gap_open_cost` and `site_mismatch_cost < 2 × site_gap_open_cost` |
| | `site_gap_open_cost` | `3` | positive: opening an insertion or deletion when a decided site's call is aligned to its interval |
| | `site_gap_extension_cost` | `2` | positive: extending one |
| | `dissent_edge` | `40` | positive: calls from a read's ends within which its insertion or deletion against a reference decision is set aside, when a read of the other strand reads the site cleanly at least this far from its own ends (SRS-CONS-011) |
| | `length_confidence` | `2.326` | finite and positive: standard errors by which every read's share of a window run's length must clear the strands' spread and noise for the two-strand rule to decide it (SRS-CONS-013; one-sided 1 %) |
| | `length_noise_forward` | `[0.065, 0.111]` | two values in `0..=0.5` (½ is the most a share can clear, `length::MAXIMUM_NOISE`), runs 1 and 2: the forward strand's replicate noise τ in a read's share (SRS-CONS-013) |
| | `length_noise_reverse` | `[0.111, 0.128]` | two values in `0..=0.5`, runs 1 and 2: the reverse strand's (SRS-CONS-013) |
| | `length_context_minimum_pairs` | `10` | at least 1: fewest cross-strand pairs of a run in a batch's run context, the sample's own out, for it to rescale τ; below it τ is the configured one (SRS-CONS-014; batches with fewer pairs carry too little evidence, run-calibrated-length research 2026-10-21) |
| | `length_context_prior_pairs` | `[4.3, 64.0]` | two finite positive values, runs 1 and 2: the corpus prior's weight in pairs when a batch's dispersion is shrunk towards it (SRS-CONS-014). Method of moments over the 1,000 dev v5 batch contexts (every read; batches with at least 10 pairs), from their `sum_z2` and `sum_z4`: n0 = per-pair relative variance of z² / between-batch variance of the relative dispersion. Run 1: 0.452 / 0.104 = 4.3 (739 batches). Run 2: the between-batch variance, 0.004 against a sampling variance of 0.707, is not distinguishable from zero, so 64 is kept from the research's HV3R instrument-run grouping |
| | `length_context_reference` | `[0.4577, 0.0542]` | two finite positive values, runs 1 and 2: the corpus mean squared normalised strand difference, the dispersion at which τ keeps its configured value (SRS-CONS-014; 13,529 and 14,542 pairs over every read of the 1,000 dev v5 batches' full folders, as `dna run-context` gathers them, re-sequencing reads included; HV2F and HV3R alone give 0.4643 and 0.0563) |
| | `length_context_minimum_scale` | `1.0` | finite and positive: the smallest scale of τ, so a quiet batch never lowers τ below the replicate estimate (SRS-CONS-014) |
| | `length_context_edges` | `[-3.0, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 3.0]` | finite, strictly ascending, at least one: the bin edges of a run context's histogram of normalised strand differences, recorded in the context so a consumer can take a robust dispersion without new contexts (SRS-CONS-014) |
| `normalization` | `shortest_dominant_ratio` | `0.6` | finite, in `(0, 0.85]`; the section is optional: the length convention reads the next-shorter length of a read's tract decode when its share is at least this fraction of the dominant's (ISFG 2014 3.3.2), a labelled suggestion in the review, never a call ([[reference/review-lengths]]). Without the section no convention applies and a near tie stays a check. `0.6` from the 2026-10-19 decoder spike (53 confident HVS-II run-1 disagreements, all with the truth at the mode − 1 where the shorter share was at least 0.6 of the mode's); above about 0.85 a decode no longer resolves the two shares (trace-simulator benchmark, 2026-10-21) |

Each section is owned by one plugin of the static registry
([[ADR-0069]]):
- the Sanger modality plugin owns `basecalling`, `signal_processing`,
  `callability`, `quality_control`, `sanger_evidence`, and `trace`, the
  constants of the tract decode `dna-trace` computes for it ([[ADR-0084]]);
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

`basecall` consumes the basecalling, signal-processing, callability, and quality-control settings; it still validates the complete schema and records the complete configuration checksum. Alignment and variant-calling settings are used by reference-guided operations; sample-reconciliation settings are consumed only by `sample` after every read has completed independent placement. Signal-processing values control observation-only annotations and never change calls, trim bounds, alignments, or variants. Callability values control the read callability view (ADR-0067), which sets the trim interval, the masked alignment query, variant eligibility, and masked sample observations; the double-peak threshold is `callability.double_peak_ratio`, set apart from `basecalling.secondary_peak_ratio`. Reference topology and reportable regions belong to the target profile. Compact result contracts record the raw configuration checksum and the profile identity but omit method constants and expanded effective values. The sample-reconciliation defaults are Tracy-derived pre-consensus admission controls: they require at least 25 comparable canonical-base positions and at least 0.50 canonical-base agreement for an overlapping read pair to be eligible for later consensus. Gaps, deletions, and unresolved symbols do not enter this nucleotide denominator. Effective values and the configuration schema version remain in the strict TOML selected for the run; the local path is omitted.

`read_end_margin` (ADR-0062) makes a variant ineligible when any of its mapped calls lies within that many calls of an uninformative call (`read_end`): beyond the trim interval, which is the callable span plus at most `trim_context_margin` dephased context calls per side, or inside a masked segment that aligns as unresolved. Masked evidence calls add the reason of their phase segment (ADR-0067). These rules use only the read's own calls and signal; no primer or file-name knowledge is involved.

## `.env`

`.env.example` is a shell-tooling template for `DNA_CONFIG=config/dna.toml` and the operational `DNA_LOG_DIR=logs`. `DNA_LOG_DIR` changes only the append-only log destination; it does not alter scientific settings or their checksum. Local `.env` remains ignored. Shells, IDEs, and containers may export these values; DNA itself never reads dotenv files.
