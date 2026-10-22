---
implements: [SRS-CALL, SRS-PRF]
verified-by: [src/config]
---

# Configuration Requirements

**Requirement namespace:** `SRS-CFG-*`

These requirements are part of the canonical [[requirements/SRS|DNA SRS]].

- **SRS-CFG-001:** DNA MUST load `DNA_CONFIG` or `config/dna.toml`; it MUST NOT parse `.env`.
- **SRS-CFG-002:** Schema version MUST be the one integer the build accepts, which the [[reference/configuration|configuration reference]] states; all documented sections/keys are required, including a non-empty root `profile` naming the target profile (SRS-PRF-001), resolved against the configuration file's directory when relative.
- **SRS-CFG-003:** Unknown/duplicate keys, missing fields, non-finite values, unsupported enums/versions, and invalid ranges MUST fail.
- **SRS-CFG-004:** Environment MUST NOT override individual scientific settings. `DNA_LOG_DIR` MAY select only the operational log directory.
- **SRS-CFG-005:** JSON MUST include the configuration checksum but MUST omit method constants, the local config path, and the effective-value expansion. Effective values and the configuration schema version remain in the selected strict TOML; target knowledge is identified separately by the profile identity (SRS-PRF-007).
- **SRS-CFG-006:** Hard caps MUST bound AB1 bytes, reference length, alignment cells, indel length, and reachable peak thresholds before unsafe allocation.
- **SRS-CFG-007:** Configuration MUST require a positive `sample_reconciliation.minimum_comparable_bases` and finite `sample_reconciliation.minimum_overlap_agreement` in `(0,1]`; these settings apply only after independently placed reads enter sample reconciliation.
- **SRS-CFG-008:** Configuration MUST require a `[callability]` section with `window_calls` in `8..=64`, finite `onset_defect_fraction` in `(0,1]`, finite `exit_defect_fraction` in `[0,1]` below the onset fraction, finite `minimum_main_share` in `(0,1]`, finite `maximum_far_share` in `[0,1]`, finite `minimum_shadow_share` in `(0,1]`, finite `weak_amplitude_fraction` in `[0,0.5]`, `repeat_min_length` of at least `2`, a positive `minimum_callable_calls`, a positive `repeat_window_calls`, finite `double_peak_ratio` in `(0,1]`, a finite positive `spacing_defect`, and a positive `spacing_half_window`; these settings control read callability, which the trim interval, alignment, eligibility, and sample aggregation consume (SRS-CALL-009..012).
- **SRS-CFG-009:** Every threshold, factor, window, and minimum count a scientific rule applies MUST be configuration, validated and identified by the configuration checksum, except where a requirement names it a method constant of a model's definition. Configuration MUST require the Sanger read-edge parameters in `[sanger_evidence]` (`signal_edge_margin`, `leading_height_fraction`, `post_repeat_margin`, a positive `post_repeat_minimum_calls`, `edge_identity_margin`, `edge_height_fraction`, a positive `end_norm_calls`, `end_height_fraction`, `edge_second_share` in `[0,1]`, `edge_width_factor`, `edge_spacing_minimum` not above `edge_spacing_maximum`, `start_identity_margin`, a positive `start_norm_calls`, `start_height_fraction`, and `start_second_share` in `[0,1]`, every factor and fraction finite and positive), `quality_control.trim_context_margin`, a finite positive `signal_processing.minimum_noise_sigma`, a positive `alignment.max_stretch`, and the `[consensus]` section (`seeded_run` of at least 2, `maximum_undecided_gap` at most the longest supported reference, a positive `minimum_resolved_bases`, and positive site costs with `site_mismatch_cost + site_gap_extension_cost > site_gap_open_cost` and `site_mismatch_cost < 2 × site_gap_open_cost`, and a positive `dissent_edge`). `variant_calling.read_end_margin` MUST be at most the longest supported reference. The Sanger margins MUST be read only by the Sanger plugin and `variant_calling.read_end_margin` only by the core, so one key never moves the other's rules.

## Implemented in

<!-- docs:begin backlinks -->
- **Implemented in:** [reference/configuration](../reference/configuration.md), [requirements/callability](callability.md)
<!-- docs:end backlinks -->
