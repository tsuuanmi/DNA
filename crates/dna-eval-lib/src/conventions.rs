//! The conventions the harness reproduces, each with its source (the
//! laboratory's comparison rules and names, and the measurement definitions
//! the user set), and the defaults of a run's options. The geometry comes
//! from the target profile ([`crate::geometry`]); every option is given or
//! defaulted on the command line and recorded in `run.json`.

/// Positions whose insertions and deletions the laboratory's comparison
/// leaves out of the primary measure ("special" positions, with their `.n`
/// insertion positions): `mtdna_raw` `src/core/comparison.py`, the rules the
/// primary measure adapts (ADR-0075).
pub(crate) const SPECIAL_POSITIONS: [&str; 4] = ["455", "463", "573", "309"];

/// The laboratory's name of the profile's reportable region at `index`:
/// `HV1`, `HV2`, `HV3` for HVS-I, HVS-II, and HVS-III, in the profile's
/// order. `scores.tsv` (`callable_yield`), the review windows, and the zone
/// names (`after@HV2`) use it.
pub(crate) fn region_name(index: usize) -> String {
    format!("HV{}", index + 1)
}

/// The limit of detection of a second base, in thousandths of a read's
/// normalised A/C/G/T profile: 10 % (user decision 2026-10-16). The
/// `haplotype_lod` measure compares a truth heteroplasmy no read shows at it
/// compatibly; `run.json` records it as `constants.lod_share`.
pub(crate) const LOD: u16 = 100;

/// The smallest second channel, in thousandths of a read's profile, that the
/// discrepancy causes read as a faint second base: 3 %. `run.json` records
/// it as `constants.faint_share`.
pub(crate) const FAINT: u16 = 30;

/// The default window of the review workload, in positions
/// (`--review-window`).
pub(crate) const REVIEW_WINDOW: i64 = 10;

/// The default widths that place the review zones, in positions
/// (`--zone-after`, `--zone-edge`, `--zone-overrun`, set with the review
/// zones on 2026-10-19): after a repeat window, at a region's ends, and past
/// a window an item still counts for it.
pub(crate) const ZONE_AFTER: i64 = 20;
pub(crate) const ZONE_EDGE: i64 = 40;
pub(crate) const ZONE_OVERRUN: i64 = 2;

/// The default fewest reads of a batch for its own run signal
/// (`--run-signal-batch-reads`, PROP-0009).
pub(crate) const RUN_SIGNAL_BATCH_READS: usize = 40;

/// The default share of the gated group's samples the clear samples may move
/// by before `diff` warns (`--clear-tolerance`): one percentage point.
pub(crate) const CLEAR_TOLERANCE: f64 = 0.01;

/// The default share of the base's positions to review they may move by
/// before `diff` warns (`--positions-tolerance`): 5 %.
pub(crate) const POSITIONS_TOLERANCE: f64 = 0.05;

/// The default points a gated measure of the simulated bench may worsen by
/// against its baseline (`sim-bench --tolerance`): none, since the frozen
/// cases render the same reads on every run, so any worsening is the
/// decoder's.
pub(crate) const SIM_BENCH_TOLERANCE: f64 = 0.0;
