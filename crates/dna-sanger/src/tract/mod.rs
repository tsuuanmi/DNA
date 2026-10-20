//! A read's homopolymer tracts, decoded by `dna-trace` (PROP-0012
//! increment 2, ADR-0084): where the read leaves phase, the run that opens
//! the tract before it, its structure (no mark, or the mark between two
//! runs), and its length marginals, from the read's own signal.

mod config;

pub use config::{RawTraceConfig, TraceConfig};

use dna_kernel::read_evidence::{
    CleanedLengths, LeakageShare, LengthMarginal, Margin, StructureFit, TractLengths,
};
use dna_trace::tract::{Marginal, Tract, analyse};

use crate::model::basecalls::BaseCalls;
use crate::model::callability::{PhaseState, ReadCallability};
use crate::model::sanger::Chromatogram;

/// The bases in channel order.
const BASES: [char; 4] = ['A', 'C', 'G', 'T'];

/// The read's tracts (SRS-CALL-025): one before each call where callability
/// sees the read leave phase (the start of a segment not in phase), each decoded from the
/// read's analysed channels at its calls' locus samples, normalised by the
/// row sum. A tract whose fit holds no signal, or that the kernel would
/// refuse, is left out.
pub(crate) fn tracts(
    trace: &Chromatogram,
    calls: &BaseCalls,
    callability: &ReadCallability,
    config: &TraceConfig,
) -> Vec<TractLengths> {
    let profiles = calls
        .calls
        .iter()
        .map(|call| {
            let heights: [f64; 4] = std::array::from_fn(|channel| {
                trace.channels[channel]
                    .get(call.locus_position_0based)
                    .copied()
                    .map_or(0.0, f64::from)
            });
            let sum = heights.iter().sum::<f64>().max(1.0);
            heights.map(|height| height / sum)
        })
        .collect::<Vec<_>>();
    let onsets = callability
        .segments
        .iter()
        .filter(|segment| segment.state != PhaseState::InPhase)
        .map(|segment| segment.call_start_0based)
        .collect::<Vec<_>>();
    let tracts = analyse(&profiles, &onsets, &config.analysis);
    let iterations_capped = tracts
        .iter()
        .filter(|tract| tract.choice.decode.iterations_capped)
        .count();
    let joint_capped = tracts
        .iter()
        .filter(|tract| tract.choice.decode.joint_capped)
        .count();
    tracing::info!(
        event = "tracts_decoded",
        onsets = onsets.len(),
        tracts = tracts.len(),
        iterations_capped,
        joint_capped,
    );
    let asserted = tracts
        .iter()
        .filter_map(|tract| assertion(tract, config.near_equal_z))
        .collect::<Vec<_>>();
    // A decode the kernel would refuse is left out and counted, never
    // failing the read: the decodes are evidence only.
    let (kept, refused): (Vec<_>, Vec<_>) = asserted
        .into_iter()
        .partition(|tract| tract.valid(calls.calls.len()));
    if !refused.is_empty() {
        tracing::warn!(event = "tracts_refused", refused = refused.len());
    }
    kept
}

/// A tract as the kernel asserts it; `None` when its fit holds no signal.
fn assertion(tract: &Tract, near_equal_z: f64) -> Option<TractLengths> {
    let decode = &tract.choice.decode;
    let cleaned = tract.cleaned.as_ref().and_then(|cleaned| {
        Some(CleanedLengths {
            leakage: tract.leakage.as_ref().map(|shares| {
                shares
                    .iter()
                    .map(|&(offset, share)| LeakageShare {
                        offset: i32::try_from(offset).unwrap_or_default(),
                        share,
                    })
                    .collect()
            }),
            residual: cleaned.residual,
            first: marginal(&cleaned.first, near_equal_z)?,
        })
    });
    Some(TractLengths {
        onset: tract.onset,
        call_start: tract.start,
        extent: decode.extent,
        window: decode.window,
        run: BASES[tract.run],
        mark: tract.choice.mark.map(|mark| BASES[mark]),
        candidates: tract
            .choice
            .residuals
            .iter()
            .map(|&(mark, residual)| StructureFit {
                mark: mark.map(|mark| BASES[mark]),
                residual,
            })
            .collect(),
        residual: decode.residual,
        sigma: decode.sigma,
        autocorrelation: decode.autocorrelation,
        first: marginal(&decode.first, near_equal_z)?,
        total: marginal(&decode.total, near_equal_z)?,
        cleaned,
    })
}

/// A marginal as the kernel asserts it: near-equal when its deflated margin
/// is below `near_equal_z`, and undefined, never clear, when two lengths show
/// without a finite deflated margin; `None` when no length holds signal.
fn marginal(marginal: &Marginal, near_equal_z: f64) -> Option<LengthMarginal> {
    marginal.dominant?;
    let finite = |value: f64| value.is_finite().then_some(value);
    let z_adjusted = finite(marginal.z_adjusted);
    // The shares are normalised by their sum, so rounding can leave one a
    // unit in the last place outside [0, 1].
    let share = |share: f64| if share <= 0.0 { 0.0 } else { share.min(1.0) };
    Some(LengthMarginal {
        lengths: marginal.lengths.clone(),
        shares: marginal.shares.iter().copied().map(share).collect(),
        errors: marginal.errors.clone(),
        dominant: marginal.dominant,
        next: marginal.next,
        z: finite(marginal.z),
        z_adjusted,
        margin: match (marginal.next, z_adjusted) {
            (Some(_), None) => Margin::Undefined,
            (Some(_), Some(z)) if z < near_equal_z => Margin::NearEqual,
            _ => Margin::Clear,
        },
    })
}
