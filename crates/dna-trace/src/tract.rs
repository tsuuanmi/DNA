//! A homopolymer tract's length mixture, decoded from the read's own signal
//! (PROP-0012 increment 2).
//!
//! Within a tract that a modality flags, the read is a mixture of molecules
//! whose runs differ in length. Each molecule's sequence over the window is
//! `run^a mark run^b D`, or `run^L D` without a mark, where `D` is the shared
//! downstream sequence, unknown, and the same for every molecule. The decode
//! fits the window's normalised profiles as non-negative fractions of every
//! molecule (Lawson–Hanson NNLS) and the downstream sequence exactly given
//! the fractions (dynamic programming over the strongest lengths), in
//! alternation, from several starts. It reads no reference: `e`, the run
//! base, and the mark come from the modality.
//!
//! The decode reports each length's share of the signal, marginal over the
//! first run and over the total (the only identifiable ones), with a standard error from the
//! fit, and the margin of the dominant length over the next, `z`, deflated by
//! the read's own residual autocorrelation (`z_adj`). Shares are signal
//! shares, never heteroplasmy levels.

use std::collections::HashMap;

use crate::locate::{LocateParams, locate};

/// The method's constants: the span of candidate lengths around the read's
/// own run extent, the downstream calls fitted, and the computation caps.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TractParams {
    /// Candidate total lengths start this many calls below the run extent.
    pub span_below: usize,
    /// And end this many calls above it.
    pub span_above: usize,
    /// The shortest candidate total length.
    pub min_length: usize,
    /// Calls after the longest candidate that the window fits.
    pub downstream: usize,
    /// Most consecutive lengths whose downstream sequence the dynamic
    /// programme decodes jointly (a computation cap: 4^K windows).
    pub joint_lengths: usize,
    /// Most alternations of fit and decode per start (a convergence cap).
    pub iterations: usize,
}

/// One length marginal: each length's share of the fitted signal and its
/// standard error, ascending, and the margin of the dominant length.
#[derive(Debug, Clone, PartialEq)]
pub struct Marginal {
    /// The candidate lengths, ascending.
    pub lengths: Vec<usize>,
    /// Each length's share of the fitted signal; the shares sum to one.
    pub shares: Vec<f64>,
    /// Each share's standard error from the fit.
    pub errors: Vec<f64>,
    /// The length with the largest share.
    pub dominant: Option<usize>,
    /// The length with the next largest positive share.
    pub next: Option<usize>,
    /// The dominant share's margin over the next, in standard errors of
    /// their difference; `f64::INFINITY` when only one length shows.
    pub z: f64,
    /// `z` deflated by the read's residual autocorrelation.
    pub z_adjusted: f64,
}

/// A tract's decode.
#[derive(Debug, Clone, PartialEq)]
pub struct TractDecode {
    /// Calls from the tract's first one that read the run or mark base.
    pub extent: usize,
    /// Calls fitted.
    pub window: usize,
    /// The decoded downstream sequence, as channel indices.
    pub downstream: Vec<usize>,
    /// Root mean square of the fit's residual per channel row.
    pub residual: f64,
    /// The residual scaled to the fit's degrees of freedom.
    pub sigma: f64,
    /// Lag-1 autocorrelation of the residual along the calls.
    pub autocorrelation: f64,
    /// The first run's lengths. A molecule's split between its runs is not
    /// identifiable: each column is its total length's run plus its mark's
    /// row, so only the first run's and the total's marginals are.
    pub first: Marginal,
    /// The total lengths, run, mark, and run.
    pub total: Marginal,
    /// Whether the kept start ran out of iterations before its downstream
    /// sequence settled (`TractParams::iterations` bound).
    pub iterations_capped: bool,
    /// Whether more lengths held signal than the dynamic programme decodes
    /// jointly (`TractParams::joint_lengths` bound) in the kept start's last
    /// decode.
    pub joint_capped: bool,
}

/// A start's fit: residual, fractions, downstream sequence, and whether its
/// iterations and its joint lengths were capped.
type Fit = (f64, Vec<f64>, Vec<usize>, bool, bool);

/// Gradient below which no column enters the NNLS active set.
const TOLERANCE: f64 = 1e-10;
/// The cost of an excluded window in the dynamic programme.
const EXCLUDED: f64 = 1e18;
/// The margin reported when a single length shows.
const SINGLE: f64 = f64::INFINITY;

/// One-hot profiles: row `b` is base `b`'s expected profile.
pub const ONE_HOT: [[f64; 4]; 4] = [
    [1.0, 0.0, 0.0, 0.0],
    [0.0, 1.0, 0.0, 0.0],
    [0.0, 0.0, 1.0, 0.0],
    [0.0, 0.0, 0.0, 1.0],
];

/// The window a tract's fit reads: the run extent, the candidate lengths,
/// the calls fitted and their profiles, and the molecules.
struct Setup {
    primaries: Vec<usize>,
    extent: usize,
    lengths: Vec<usize>,
    window: usize,
    observed: Vec<f64>,
    molecules: Vec<Molecule>,
}

impl Setup {
    fn new(
        profiles: &[[f64; 4]],
        start: usize,
        run: usize,
        mark: Option<usize>,
        window: Option<usize>,
        params: &TractParams,
    ) -> Option<Self> {
        let primaries = profiles.iter().map(argmax).collect::<Vec<_>>();
        let extent = extent(&primaries, start, run, mark);
        let lowest = extent
            .saturating_sub(params.span_below)
            .max(params.min_length);
        let mut lengths = (lowest..=extent + params.span_above).collect::<Vec<_>>();
        let longest = *lengths.last()?;
        let window = profiles
            .len()
            .checked_sub(start)?
            .min(window.unwrap_or(longest + params.downstream));
        // A read that ends within the candidates fits only the lengths its
        // window holds: a longer molecule's column would be its run alone,
        // the same as every other one past the read's end.
        lengths.retain(|&length| length <= window);
        if lengths.is_empty() {
            return None;
        }
        let observed = profiles[start..start + window]
            .iter()
            .flatten()
            .copied()
            .collect::<Vec<_>>();
        let molecules = lengths
            .iter()
            .flat_map(|&length| -> Vec<Molecule> {
                match mark {
                    Some(_) => (1..length.saturating_sub(1))
                        .map(|first| Molecule {
                            first,
                            second: Some(length - first - 1),
                        })
                        .collect(),
                    None => vec![Molecule {
                        first: length,
                        second: None,
                    }],
                }
            })
            .collect::<Vec<_>>();
        Some(Self {
            primaries,
            extent,
            lengths,
            window,
            observed,
            molecules,
        })
    }
}

/// Decodes the tract starting at call `start` of `profiles` (normalised
/// A/C/G/T rows), whose run reads channel `run` and whose mark, if any,
/// channel `mark`, with `signature` as each base's expected profile, over
/// `window` calls (at most; by default its longest length's and the
/// downstream calls). `None` when the window holds no call.
#[must_use]
pub fn decode(
    profiles: &[[f64; 4]],
    start: usize,
    run: usize,
    mark: Option<usize>,
    window: Option<usize>,
    signature: &[[f64; 4]; 4],
    params: &TractParams,
) -> Option<TractDecode> {
    let setup = Setup::new(profiles, start, run, mark, window, params)?;
    let model = Model {
        molecules: &setup.molecules,
        window: setup.window,
        run,
        mark,
        signature,
        leakage: &[],
        before: &[],
    };
    let observed = &setup.observed;
    // Each step is a function of the downstream sequence alone: the fit of
    // its design, then the decode given that fit. Starts often meet on a
    // sequence, so each sequence's fit and decode are computed once.
    let mut fits: HashMap<Vec<usize>, Vec<f64>> = HashMap::new();
    let mut steps: HashMap<Vec<usize>, (Vec<usize>, bool)> = HashMap::new();
    let mut fit = |downstream: &[usize]| -> Vec<f64> {
        fits.entry(downstream.to_vec())
            .or_insert_with(|| nnls(&model.design(downstream), observed))
            .clone()
    };
    let mut best: Option<Fit> = None;
    for &initial in &setup.lengths {
        let mut downstream = (0..setup.window)
            .map(|m| {
                setup
                    .primaries
                    .get(start + initial + m)
                    .copied()
                    .unwrap_or(0)
            })
            .collect::<Vec<_>>();
        if downstream.first() == Some(&run) {
            downstream[0] = (run + 1) % 4;
        }
        let mut previous: Option<Vec<usize>> = None;
        let mut settled = false;
        let mut joint_capped = false;
        for _ in 0..params.iterations {
            let (next, capped) = if let Some(step) = steps.get(&downstream) {
                step.clone()
            } else {
                let fractions = fit(&downstream);
                let step = model.decode_downstream(
                    observed,
                    &fractions,
                    &downstream,
                    params.joint_lengths,
                );
                steps.insert(downstream.clone(), step.clone());
                step
            };
            downstream = next;
            joint_capped = capped;
            if previous.as_ref() == Some(&downstream) {
                settled = true;
                break;
            }
            previous = Some(downstream.clone());
        }
        let fractions = fit(&downstream);
        let residual = rms(&residuals(&model.design(&downstream), &fractions, observed));
        if best.as_ref().is_none_or(|(kept, ..)| residual < *kept) {
            best = Some((residual, fractions, downstream, !settled, joint_capped));
        }
    }
    let (residual, fractions, downstream, iterations_capped, joint_capped) = best?;
    let design = model.design(&downstream);
    let mut decoded = summarise(
        &setup.molecules,
        &design,
        &fractions,
        observed,
        residual,
        downstream,
        setup.extent,
        setup.window,
    );
    decoded.iterations_capped = iterations_capped;
    decoded.joint_capped = joint_capped;
    Some(decoded)
}

/// A tract's structure: no mark, or the mark that fits best, with its decode.
#[derive(Debug, Clone, PartialEq)]
pub struct StructureChoice {
    /// The chosen mark's channel; `None` for no mark.
    pub mark: Option<usize>,
    /// Each candidate's residual over the shared window, in evaluation
    /// order: no mark, then the other bases in A, C, G, T order.
    pub residuals: Vec<(Option<usize>, f64)>,
    /// The chosen structure's decode over the shared window.
    pub decode: TractDecode,
}

/// Chooses the tract's structure: no mark (`run^L D`), or a mark of each
/// base but the run (`run^a mark run^b D`), whichever decode fits with the
/// lowest residual. Each candidate keeps its own run extent and lengths; all
/// fit the same calls, the longest extent's window, so their residuals
/// compare. An exact tie keeps the earlier candidate. `None` when the window
/// holds no call.
#[must_use]
pub fn choose(
    profiles: &[[f64; 4]],
    start: usize,
    run: usize,
    signature: &[[f64; 4]; 4],
    params: &TractParams,
) -> Option<StructureChoice> {
    let primaries = profiles.iter().map(argmax).collect::<Vec<_>>();
    let candidates = std::iter::once(None)
        .chain((0..4).filter(|&base| base != run).map(Some))
        .collect::<Vec<_>>();
    let longest = candidates
        .iter()
        .map(|&mark| extent(&primaries, start, run, mark))
        .max()?;
    let window = longest + params.span_above + params.downstream;
    let decodes = std::thread::scope(|scope| {
        let handles = candidates
            .iter()
            .map(|&mark| {
                scope.spawn(move || {
                    decode(profiles, start, run, mark, Some(window), signature, params)
                })
            })
            .collect::<Vec<_>>();
        handles
            .into_iter()
            .map(|handle| handle.join().ok().flatten())
            .collect::<Option<Vec<_>>>()
    })?;
    let residuals = candidates
        .iter()
        .zip(&decodes)
        .map(|(&mark, decoded)| (mark, decoded.residual))
        .collect::<Vec<_>>();
    let chosen = (1..decodes.len()).fold(0, |kept, index| {
        if decodes[index].residual < decodes[kept].residual {
            index
        } else {
            kept
        }
    });
    let decode = decodes.into_iter().nth(chosen)?;
    Some(StructureChoice {
        mark: candidates[chosen],
        residuals,
        decode,
    })
}

/// Every constant of a read's tract analysis.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AnalysisParams {
    /// Where a tract opens.
    pub locate: LocateParams,
    /// The decode.
    pub tract: TractParams,
    /// The read's signature and leakage null.
    pub leakage: LeakageParams,
}

/// One tract of a read: where it opens, its chosen structure and decode,
/// and the decode cleaned by the read's own signature and leakage (SL).
#[derive(Debug, Clone, PartialEq)]
pub struct Tract {
    /// The call where the modality sees the read leave phase.
    pub onset: usize,
    /// The tract's first call.
    pub start: usize,
    /// The run base's channel.
    pub run: usize,
    /// The chosen structure and its one-hot decode over the shared window.
    pub choice: StructureChoice,
    /// Each leakage offset's share, from the read's calls before the tract;
    /// `None` when too few calls measure it.
    pub leakage: Option<Vec<(isize, f64)>>,
    /// The chosen structure refitted with the read's signature and leakage
    /// over the signature decode's downstream sequence, on the same window.
    pub cleaned: Option<TractDecode>,
}

/// The tracts that open before each of `onsets` (calls where the modality
/// sees the read leave phase), one per distinct first call, in onset order:
/// located from the read's own primaries, their structure chosen, decoded,
/// and cleaned. An onset with no run before it yields none.
#[must_use]
pub fn analyse(profiles: &[[f64; 4]], onsets: &[usize], params: &AnalysisParams) -> Vec<Tract> {
    let primaries = profiles.iter().map(argmax).collect::<Vec<_>>();
    let mut tracts: Vec<Tract> = Vec::new();
    for &onset in onsets {
        let Some((start, run)) = locate(&primaries, onset, &params.locate) else {
            continue;
        };
        if tracts.iter().any(|tract| tract.start == start) {
            continue;
        }
        let Some(choice) = choose(profiles, start, run, &ONE_HOT, &params.tract) else {
            continue;
        };
        let rows = signature(profiles, start, &params.leakage);
        let shares = leakage(profiles, start, &rows, &params.leakage);
        let mark = choice.mark;
        let window = Some(choice.decode.window);
        let cleaned =
            decode(profiles, start, run, mark, window, &rows, &params.tract).and_then(|signed| {
                refit(
                    profiles,
                    start,
                    run,
                    mark,
                    &rows,
                    shares.as_deref().unwrap_or(&[]),
                    &signed.downstream,
                    &params.tract,
                )
            });
        tracts.push(Tract {
            onset,
            start,
            run,
            choice,
            leakage: shares,
            cleaned,
        });
    }
    tracts
}

/// The leakage model: each base's mean profile over the read's in-phase calls
/// before the tract, and the share of the signal that each offset call leaks
/// into a call (a neighbour's peak tail).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LeakageParams {
    /// Calls before the tract left out of the signature and the null.
    pub gap: usize,
    /// The leakage null leaves out the calls within `pad` calls of a run of
    /// at least `min_run` equal primaries, where the tails stack.
    pub min_run: usize,
    /// See `min_run`.
    pub pad: usize,
    /// Fewest calls the null needs.
    pub min_rows: usize,
    /// Calls before the tract that the signature and the null read, at most:
    /// the read's in-phase context near the tract.
    pub context: usize,
}

/// The calls before the tract a refit reads for the leakage of its first
/// calls: as far as the farthest offset reaches, and one more.
const BEFORE_CALLS: usize = 3;

/// The offsets the leakage null measures, in the order they are added.
pub const LEAKAGE_OFFSETS: [isize; 4] = [-2, -1, 1, 2];

/// Each base's mean profile over the context calls before `start − gap`
/// whose primary it is; a base no such call reads keeps its one-hot profile.
#[must_use]
pub fn signature(profiles: &[[f64; 4]], start: usize, params: &LeakageParams) -> [[f64; 4]; 4] {
    let before = context(profiles, start, params);
    let mut rows = ONE_HOT;
    for (base, row) in rows.iter_mut().enumerate() {
        let calls = before
            .iter()
            .filter(|profile| argmax(profile) == base)
            .collect::<Vec<_>>();
        if calls.is_empty() {
            continue;
        }
        let mut mean = [0.0; 4];
        for profile in &calls {
            for channel in 0..4 {
                mean[channel] += profile[channel];
            }
        }
        for value in &mut mean {
            *value /= calls.len() as f64;
        }
        *row = mean;
    }
    rows
}

/// The calls before `start − gap`, at most `context` calls before `start`.
fn context<'a>(profiles: &'a [[f64; 4]], start: usize, params: &LeakageParams) -> &'a [[f64; 4]] {
    let end = start.saturating_sub(params.gap).min(profiles.len());
    &profiles[start.saturating_sub(params.context).min(end)..end]
}

/// The leakage null: over the read's context calls before `start − gap`, away from
/// runs, the share each offset call in [`LEAKAGE_OFFSETS`] adds to a call's
/// profile relative to the call's own base, by NNLS. `None` with fewer than
/// `min_rows` calls, or when the call's own base takes no share.
#[must_use]
pub fn leakage(
    profiles: &[[f64; 4]],
    start: usize,
    signature: &[[f64; 4]; 4],
    params: &LeakageParams,
) -> Option<Vec<(isize, f64)>> {
    let before = context(profiles, start, params);
    let bases = before.iter().map(argmax).collect::<Vec<_>>();
    let count = bases.len();
    let mut masked = vec![false; count];
    let mut index = 0;
    while index < count {
        let mut end = index;
        while end + 1 < count && bases[end + 1] == bases[index] {
            end += 1;
        }
        if end - index + 1 >= params.min_run {
            for cell in masked
                .iter_mut()
                .take((end + params.pad + 1).min(count))
                .skip(index.saturating_sub(params.pad))
            {
                *cell = true;
            }
        }
        index = end + 1;
    }
    let reach_back = LEAKAGE_OFFSETS
        .iter()
        .copied()
        .max()
        .unwrap_or(0)
        .unsigned_abs();
    let reach_ahead = LEAKAGE_OFFSETS
        .iter()
        .copied()
        .min()
        .unwrap_or(0)
        .min(0)
        .unsigned_abs();
    let rows = (reach_back..count.saturating_sub(reach_ahead))
        .filter(|&row| !masked[row])
        .collect::<Vec<_>>();
    if rows.len() < params.min_rows {
        return None;
    }
    let offsets = std::iter::once(0)
        .chain(LEAKAGE_OFFSETS)
        .collect::<Vec<isize>>();
    let design = offsets
        .iter()
        .map(|&offset| {
            rows.iter()
                .flat_map(|&row| {
                    let source = row.checked_add_signed(-offset).unwrap_or(row);
                    signature[bases[source]]
                })
                .collect::<Vec<f64>>()
        })
        .collect::<Vec<_>>();
    let observed = rows
        .iter()
        .flat_map(|&row| before[row])
        .collect::<Vec<f64>>();
    let shares = nnls(&design, &observed);
    if shares[0] <= 0.0 {
        return None;
    }
    Some(
        offsets
            .iter()
            .zip(&shares)
            .skip(1)
            .map(|(&offset, &share)| (offset, share / shares[0]))
            .collect(),
    )
}

/// Refits the tract with `signature` and `leakage` over the downstream
/// sequence `downstream` (the signature decode's), keeping the decode's
/// window: the leakage-corrected shares. The calls before the tract extend
/// each sequence to the left for the leakage of its first calls.
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn refit(
    profiles: &[[f64; 4]],
    start: usize,
    run: usize,
    mark: Option<usize>,
    signature: &[[f64; 4]; 4],
    leakage: &[(isize, f64)],
    downstream: &[usize],
    params: &TractParams,
) -> Option<TractDecode> {
    let setup = Setup::new(profiles, start, run, mark, Some(downstream.len()), params)?;
    let before = setup.primaries[start.saturating_sub(BEFORE_CALLS)..start].to_vec();
    let model = Model {
        molecules: &setup.molecules,
        window: setup.window,
        run,
        mark,
        signature,
        leakage,
        before: &before,
    };
    let design = model.design(downstream);
    let fractions = nnls(&design, &setup.observed);
    let residual = rms(&residuals(&design, &fractions, &setup.observed));
    Some(summarise(
        &setup.molecules,
        &design,
        &fractions,
        &setup.observed,
        residual,
        downstream.to_vec(),
        setup.extent,
        setup.window,
    ))
}

#[derive(Debug, Clone, Copy)]
struct Molecule {
    first: usize,
    second: Option<usize>,
}

impl Molecule {
    fn length(self) -> usize {
        self.second
            .map_or(self.first, |second| self.first + second + 1)
    }
}

struct Model<'a> {
    molecules: &'a [Molecule],
    window: usize,
    run: usize,
    mark: Option<usize>,
    signature: &'a [[f64; 4]; 4],
    /// Offset shares added to each call, in order; none in a decode.
    leakage: &'a [(isize, f64)],
    /// The calls before the tract, for the leakage of its first calls.
    before: &'a [usize],
}

impl Model<'_> {
    /// The base of `molecule`'s tract at row `row`, inside its tract.
    fn tract_base(&self, molecule: Molecule, row: usize) -> usize {
        match (self.mark, molecule.second) {
            (Some(mark), Some(_)) if row == molecule.first => mark,
            _ => self.run,
        }
    }

    /// The design's columns, one per molecule, each `4 × window` rows: the
    /// profile of the molecule's sequence, plus each leakage offset's share of
    /// the profile of the call that offset away.
    fn design(&self, downstream: &[usize]) -> Vec<Vec<f64>> {
        self.molecules
            .iter()
            .map(|&molecule| {
                let length = molecule.length();
                let sequence = (0..self.window)
                    .map(|row| {
                        if row < length {
                            self.tract_base(molecule, row)
                        } else {
                            downstream[row - length]
                        }
                    })
                    .collect::<Vec<_>>();
                let mut column = vec![0.0; self.window * 4];
                for (row, &base) in sequence.iter().enumerate() {
                    column[row * 4..row * 4 + 4].copy_from_slice(&self.signature[base]);
                }
                let lead = self.before.len();
                for &(offset, share) in self.leakage {
                    for row in 0..self.window {
                        let Some(index) = (row + lead).checked_add_signed(-offset) else {
                            continue;
                        };
                        let base = if index < lead {
                            self.before[index]
                        } else if let Some(&base) = sequence.get(index - lead) {
                            base
                        } else {
                            continue;
                        };
                        for channel in 0..4 {
                            column[row * 4 + channel] += share * self.signature[base][channel];
                        }
                    }
                }
                column
            })
            .collect()
    }

    /// The downstream sequence that best fits `observed` given `fractions`,
    /// decoded jointly over the window of at most `joint` consecutive total
    /// lengths holding the most signal; the other lengths keep `current`.
    /// Also whether lengths outside that window held signal.
    #[allow(clippy::too_many_lines)]
    fn decode_downstream(
        &self,
        observed: &[f64],
        fractions: &[f64],
        current: &[usize],
        joint: usize,
    ) -> (Vec<usize>, bool) {
        let mut weight = std::collections::BTreeMap::<usize, f64>::new();
        for (molecule, &fraction) in self.molecules.iter().zip(fractions) {
            if fraction > 0.0 {
                *weight.entry(molecule.length()).or_default() += fraction;
            }
        }
        let Some(window_start) = weight
            .keys()
            .copied()
            .map(|from| {
                let held = weight
                    .range(from..from + joint)
                    .map(|(_, value)| *value)
                    .sum::<f64>();
                (from, held)
            })
            .fold(None::<(usize, f64)>, |best, (from, held)| match best {
                Some((_, kept)) if held <= kept => best,
                _ => Some((from, held)),
            })
            .map(|(from, _)| from)
        else {
            return (current.to_vec(), false);
        };
        let joined = weight
            .range(window_start..window_start + joint)
            .map(|(&length, _)| length)
            .collect::<Vec<_>>();
        let capped = joined.len() < weight.len();
        let (Some(&shortest), Some(&longest)) = (joined.first(), joined.last()) else {
            return (current.to_vec(), capped);
        };
        let width = longest - shortest + 1;
        let rows = self.window;
        // Every molecule's tract, and the lengths outside the window with
        // the current downstream sequence.
        let mut base = vec![[0.0; 4]; rows];
        for (&molecule, &fraction) in self.molecules.iter().zip(fractions) {
            if fraction > 0.0 {
                for (row, cell) in base.iter_mut().enumerate().take(molecule.length()) {
                    let profile = &self.signature[self.tract_base(molecule, row)];
                    for channel in 0..4 {
                        cell[channel] += fraction * profile[channel];
                    }
                }
            }
        }
        for (&length, &value) in &weight {
            if joined.contains(&length) {
                continue;
            }
            for row in length..rows {
                let profile = &self.signature[current[row - length]];
                for channel in 0..4 {
                    base[row][channel] += value * profile[channel];
                }
            }
        }
        // Window position i holds D[m − K + 1 + i], the downstream base of
        // the molecules of length `longest − i`.
        let weights = (0..width)
            .map(|i| weight.get(&(longest - i)).copied().unwrap_or(0.0))
            .collect::<Vec<_>>();
        let windows = 4usize.pow(u32::try_from(width).unwrap_or(0));
        let states = windows / 4;
        let digit = |window: usize, i: usize| (window >> (2 * (width - 1 - i))) & 3;
        let prefix = |window: usize| window / 4;
        let suffix = |window: usize| window % states.max(1);
        let contribution = |window: usize, from: usize| {
            let mut sum = [0.0; 4];
            for (i, &value) in weights.iter().enumerate().skip(from) {
                let profile = &self.signature[digit(window, i)];
                for channel in 0..4 {
                    sum[channel] += value * profile[channel];
                }
            }
            sum
        };
        let full = (0..windows)
            .map(|window| contribution(window, 0))
            .collect::<Vec<_>>();
        let steps = rows.saturating_sub(shortest);
        let mut cost = vec![EXCLUDED; states.max(1)];
        cost[0] = 0.0;
        let mut back: Vec<Vec<usize>> = Vec::with_capacity(steps);
        for m in 0..steps {
            let row = m + shortest;
            // Positions before the downstream sequence starts.
            let first_real = (width - 1).saturating_sub(m);
            let mut next = vec![f64::INFINITY; states.max(1)];
            let mut chosen = vec![usize::MAX; states.max(1)];
            if m > 0 && first_real == 0 {
                // Every position is downstream and the run is no longer
                // excluded: the same arithmetic without the per-window
                // checks.
                let target: [f64; 4] =
                    std::array::from_fn(|channel| observed[row * 4 + channel] - base[row][channel]);
                let mask = states.max(1) - 1;
                for (window, added) in full.iter().enumerate() {
                    let d0 = target[0] - added[0];
                    let d1 = target[1] - added[1];
                    let d2 = target[2] - added[2];
                    let d3 = target[3] - added[3];
                    let squared = d0 * d0 + d1 * d1 + d2 * d2 + d3 * d3;
                    let mut total = cost[window >> 2];
                    if total < EXCLUDED {
                        total += squared;
                    }
                    let state = window & mask;
                    if total < next[state] {
                        next[state] = total;
                        chosen[state] = window;
                    }
                }
                back.push(chosen);
                cost = next;
                continue;
            }
            for window in 0..windows {
                let mut total = cost[prefix(window)];
                if m == 0 && digit(window, width - 1) == self.run {
                    total = EXCLUDED;
                }
                if (0..first_real).any(|i| digit(window, i) != 0) {
                    total = EXCLUDED;
                }
                let added = if first_real == 0 {
                    full[window]
                } else {
                    contribution(window, first_real)
                };
                let mut squared = 0.0;
                for channel in 0..4 {
                    let difference =
                        observed[row * 4 + channel] - base[row][channel] - added[channel];
                    squared += difference * difference;
                }
                if total < EXCLUDED {
                    total += squared;
                }
                let state = suffix(window);
                if total < next[state] {
                    next[state] = total;
                    chosen[state] = window;
                }
            }
            back.push(chosen);
            cost = next;
        }
        let mut decoded = current.to_vec();
        let mut state = cost
            .iter()
            .enumerate()
            .fold((0, f64::INFINITY), |(kept, low), (index, &value)| {
                if value < low {
                    (index, value)
                } else {
                    (kept, low)
                }
            })
            .0;
        for m in (0..steps).rev() {
            let window = back[m][state];
            if window == usize::MAX {
                break;
            }
            decoded[m] = digit(window, width - 1);
            state = prefix(window);
        }
        (decoded, capped)
    }
}

/// The calls from `start` whose primary reads the run or the mark.
fn extent(primaries: &[usize], start: usize, run: usize, mark: Option<usize>) -> usize {
    primaries
        .iter()
        .skip(start)
        .take_while(|&&base| base == run || Some(base) == mark)
        .count()
}

/// The base with the largest profile entry, the lower channel on a tie.
fn argmax(profile: &[f64; 4]) -> usize {
    (1..4).fold(0, |best, channel| {
        if profile[channel] > profile[best] {
            channel
        } else {
            best
        }
    })
}

fn residuals(design: &[Vec<f64>], fractions: &[f64], observed: &[f64]) -> Vec<f64> {
    let mut fitted = vec![0.0; observed.len()];
    for (column, &fraction) in design.iter().zip(fractions) {
        if fraction != 0.0 {
            for (cell, value) in fitted.iter_mut().zip(column) {
                *cell += fraction * value;
            }
        }
    }
    fitted
        .iter()
        .zip(observed)
        .map(|(fit, y)| fit - y)
        .collect()
}

fn rms(values: &[f64]) -> f64 {
    (values.iter().map(|value| value * value).sum::<f64>() / values.len().max(1) as f64).sqrt()
}

/// Non-negative least squares by the Lawson–Hanson active-set method over
/// `design`'s columns; the passive least-squares problems take the
/// minimum-norm solution.
pub(crate) fn nnls(design: &[Vec<f64>], observed: &[f64]) -> Vec<f64> {
    let columns = design.len();
    let mut solution = vec![0.0; columns];
    let mut passive = vec![false; columns];
    // Each column's non-zero rows: a one-hot design is three-quarters zeros,
    // and a zero adds nothing to a fit or a gradient.
    let sparse = design
        .iter()
        .map(|column| {
            column
                .iter()
                .enumerate()
                .filter(|(_, value)| **value != 0.0)
                .map(|(row, &value)| (row, value))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let gradient = |solution: &[f64]| -> Vec<f64> {
        // The fit, column by column in column order, over the positive
        // fractions only.
        let mut fitted = vec![-0.0; observed.len()];
        for (column, &x) in sparse.iter().zip(solution) {
            if x != 0.0 {
                for &(row, value) in column {
                    fitted[row] += value * x;
                }
            }
        }
        let residual = observed
            .iter()
            .zip(&fitted)
            .map(|(y, fit)| y - fit)
            .collect::<Vec<_>>();
        sparse
            .iter()
            .map(|column| {
                column
                    .iter()
                    .map(|&(row, value)| value * residual[row])
                    .sum()
            })
            .collect()
    };
    let mut factor = Factor::default();
    let mut descent = gradient(&solution);
    let mut outer = 0;
    while outer < 3 * columns {
        let Some(entering) =
            (0..columns)
                .filter(|&column| !passive[column])
                .fold(None::<usize>, |best, column| match best {
                    Some(kept) if descent[column] <= descent[kept] => Some(kept),
                    _ => Some(column),
                })
        else {
            break;
        };
        if descent[entering] <= TOLERANCE {
            break;
        }
        outer += 1;
        passive[entering] = true;
        loop {
            let chosen = (0..columns)
                .filter(|&column| passive[column])
                .collect::<Vec<_>>();
            let least = least_squares(design, observed, &chosen, &mut factor);
            let mut trial = vec![0.0; columns];
            for (&column, &value) in chosen.iter().zip(&least) {
                trial[column] = value;
            }
            if chosen.iter().all(|&column| trial[column] > TOLERANCE) {
                solution = trial;
                break;
            }
            let step = chosen
                .iter()
                .filter(|&&column| trial[column] <= TOLERANCE)
                .map(|&column| solution[column] / (solution[column] - trial[column]))
                .fold(f64::INFINITY, f64::min);
            for column in 0..columns {
                solution[column] += step * (trial[column] - solution[column]);
            }
            for column in 0..columns {
                passive[column] = passive[column] && solution[column] > TOLERANCE;
                if !passive[column] {
                    solution[column] = 0.0;
                }
            }
            if !passive.iter().any(|&flag| flag) {
                break;
            }
        }
        descent = gradient(&solution);
    }
    solution
}

/// Minimum-norm least squares over `chosen` columns of `design`, by a
/// one-sided Jacobi singular value decomposition: singular values below
/// machine epsilon times the larger dimension times the largest one count as
/// zero, as LAPACK's minimum-norm solver treats them. The decode's columns
/// can be linearly dependent, so the minimum-norm solution is the one that
/// fixes the active set.
fn least_squares(
    design: &[Vec<f64>],
    observed: &[f64],
    chosen: &[usize],
    factor: &mut Factor,
) -> Vec<f64> {
    let size = chosen.len();
    if size == 0 {
        return Vec::new();
    }
    // A full-rank passive set has one least-squares solution, which QR finds
    // fast; only a rank-deficient one needs the minimum-norm solve.
    if let Some(solution) = full_rank_solution(design, observed, chosen, factor) {
        return solution;
    }
    let mut left = chosen
        .iter()
        .map(|&column| design[column].clone())
        .collect::<Vec<_>>();
    let mut right = (0..size)
        .map(|i| {
            (0..size)
                .map(|j| if i == j { 1.0 } else { 0.0 })
                .collect::<Vec<f64>>()
        })
        .collect::<Vec<_>>();
    for _ in 0..80 {
        let mut rotated = false;
        for p in 0..size {
            for q in p + 1..size {
                let alpha = left[p].iter().map(|value| value * value).sum::<f64>();
                let beta = left[q].iter().map(|value| value * value).sum::<f64>();
                let gamma = left[p]
                    .iter()
                    .zip(&left[q])
                    .map(|(a, b)| a * b)
                    .sum::<f64>();
                if gamma.abs() <= f64::EPSILON * (alpha * beta).sqrt() || gamma == 0.0 {
                    continue;
                }
                rotated = true;
                let zeta = (beta - alpha) / (2.0 * gamma);
                let tangent = zeta.signum() / (zeta.abs() + (1.0 + zeta * zeta).sqrt());
                let cos_theta = 1.0 / (1.0 + tangent * tangent).sqrt();
                let sin_theta = cos_theta * tangent;
                for matrix in [&mut left, &mut right] {
                    let (first, second) = matrix.split_at_mut(q);
                    for (x, y) in first[p].iter_mut().zip(second[0].iter_mut()) {
                        let (old_x, old_y) = (*x, *y);
                        *x = cos_theta * old_x - sin_theta * old_y;
                        *y = sin_theta * old_x + cos_theta * old_y;
                    }
                }
            }
        }
        if !rotated {
            break;
        }
    }
    let values = left
        .iter()
        .map(|column| column.iter().map(|value| value * value).sum::<f64>().sqrt())
        .collect::<Vec<_>>();
    let largest = values.iter().copied().fold(0.0, f64::max);
    let cutoff = f64::EPSILON * observed.len().max(size) as f64 * largest;
    let mut solution = vec![0.0; size];
    for (k, &value) in values.iter().enumerate() {
        if value <= cutoff || value == 0.0 {
            continue;
        }
        let projection = left[k]
            .iter()
            .zip(observed)
            .map(|(a, y)| a * y)
            .sum::<f64>()
            / (value * value);
        for (cell, v) in solution.iter_mut().zip(&right[k]) {
            *cell += projection * v;
        }
    }
    solution
}

/// A modified Gram–Schmidt QR of some columns, in order, kept between
/// solves: the factor of a column depends only on the columns before it, so
/// a solve whose columns share a prefix with the last one's reuses it.
#[derive(Default)]
struct Factor {
    /// The columns factored, in order.
    columns: Vec<usize>,
    /// The orthonormal basis, one vector per column.
    basis: Vec<Vec<f64>>,
    /// Column `k` of the upper triangle: its `k + 1` entries.
    upper: Vec<Vec<f64>>,
}

/// The least-squares solution over `chosen` columns by modified Gram–Schmidt
/// QR, `None` when a column is, within the rank cutoff of [`least_squares`],
/// a combination of the ones before it.
fn full_rank_solution(
    design: &[Vec<f64>],
    observed: &[f64],
    chosen: &[usize],
    factor: &mut Factor,
) -> Option<Vec<f64>> {
    let size = chosen.len();
    let largest = chosen
        .iter()
        .map(|&column| {
            design[column]
                .iter()
                .map(|value| value * value)
                .sum::<f64>()
                .sqrt()
        })
        .fold(0.0, f64::max);
    let cutoff = 1e-8 * largest.max(f64::MIN_POSITIVE);
    let shared = factor
        .columns
        .iter()
        .zip(chosen)
        .take_while(|(kept, wanted)| kept == wanted)
        .count();
    factor.columns.truncate(shared);
    factor.basis.truncate(shared);
    factor.upper.truncate(shared);
    if factor
        .upper
        .iter()
        .enumerate()
        .any(|(k, column)| column[k] <= cutoff)
    {
        return None;
    }
    for &column in &chosen[shared..] {
        let mut vector = design[column].clone();
        let mut entries = Vec::with_capacity(factor.basis.len() + 1);
        for q in &factor.basis {
            let dot = q.iter().zip(&vector).map(|(a, b)| a * b).sum::<f64>();
            entries.push(dot);
            for (cell, value) in vector.iter_mut().zip(q) {
                *cell -= dot * value;
            }
        }
        let norm = vector.iter().map(|value| value * value).sum::<f64>().sqrt();
        entries.push(norm);
        if norm <= cutoff {
            return None;
        }
        for cell in &mut vector {
            *cell /= norm;
        }
        factor.columns.push(column);
        factor.basis.push(vector);
        factor.upper.push(entries);
    }
    let projected = factor
        .basis
        .iter()
        .map(|q| q.iter().zip(observed).map(|(a, y)| a * y).sum::<f64>())
        .collect::<Vec<_>>();
    let mut solution = vec![0.0; size];
    for k in (0..size).rev() {
        let mut sum = projected[k];
        for (column, value) in factor.upper.iter().zip(&solution).skip(k + 1) {
            sum -= column[k] * value;
        }
        solution[k] = sum / factor.upper[k][k];
    }
    Some(solution)
}

/// The pseudo-inverse of a symmetric positive semi-definite matrix by the
/// cyclic Jacobi eigenvalue method: eigenvalues below 1e-15 times the
/// largest count as zero, as numpy's `pinv` treats them.
fn inverse(matrix: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let size = matrix.len();
    let mut a = matrix.to_vec();
    let mut vectors = (0..size)
        .map(|i| {
            (0..size)
                .map(|j| if i == j { 1.0 } else { 0.0 })
                .collect::<Vec<f64>>()
        })
        .collect::<Vec<_>>();
    for _ in 0..100 {
        let off = (0..size)
            .flat_map(|i| (0..size).filter(move |&j| j != i).map(move |j| (i, j)))
            .map(|(i, j)| a[i][j] * a[i][j])
            .sum::<f64>();
        if off == 0.0 {
            break;
        }
        for p in 0..size {
            for q in p + 1..size {
                if a[p][q] == 0.0 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                let tangent = theta.signum() / (theta.abs() + (1.0 + theta * theta).sqrt());
                let tangent = if theta == 0.0 { 1.0 } else { tangent };
                let cos_theta = 1.0 / (1.0 + tangent * tangent).sqrt();
                let sin_theta = tangent * cos_theta;
                for row in &mut a {
                    let (kp, kq) = (row[p], row[q]);
                    row[p] = cos_theta * kp - sin_theta * kq;
                    row[q] = sin_theta * kp + cos_theta * kq;
                }
                let (upper, lower) = a.split_at_mut(q);
                for (pk, qk) in upper[p].iter_mut().zip(lower[0].iter_mut()) {
                    let (old_p, old_q) = (*pk, *qk);
                    *pk = cos_theta * old_p - sin_theta * old_q;
                    *qk = sin_theta * old_p + cos_theta * old_q;
                }
                for row in &mut vectors {
                    let (vp, vq) = (row[p], row[q]);
                    row[p] = cos_theta * vp - sin_theta * vq;
                    row[q] = sin_theta * vp + cos_theta * vq;
                }
            }
        }
    }
    let values = (0..size).map(|i| a[i][i]).collect::<Vec<_>>();
    let largest = values.iter().map(|value| value.abs()).fold(0.0, f64::max);
    let cutoff = 1e-15 * largest;
    let mut result = vec![vec![0.0; size]; size];
    for (k, &value) in values.iter().enumerate() {
        if value.abs() <= cutoff {
            continue;
        }
        for i in 0..size {
            for j in 0..size {
                result[i][j] += vectors[i][k] * vectors[j][k] / value;
            }
        }
    }
    result
}

#[allow(clippy::too_many_arguments)]
fn summarise(
    molecules: &[Molecule],
    design: &[Vec<f64>],
    fractions: &[f64],
    observed: &[f64],
    residual: f64,
    downstream: Vec<usize>,
    extent: usize,
    window: usize,
) -> TractDecode {
    let sum = fractions.iter().sum::<f64>();
    let scale = if sum == 0.0 { 1.0 } else { sum };
    let positive = (0..fractions.len())
        .filter(|&k| fractions[k] > 0.0)
        .collect::<Vec<_>>();
    let rows = observed.len();
    let sigma =
        residual * (rows as f64 / (rows.saturating_sub(positive.len()).max(1)) as f64).sqrt();
    // Covariance of the normalised fractions over the positive columns.
    let gram = positive
        .iter()
        .map(|&i| {
            positive
                .iter()
                .map(|&j| design[i].iter().zip(&design[j]).map(|(a, b)| a * b).sum())
                .collect()
        })
        .collect::<Vec<Vec<f64>>>();
    let inverted = inverse(&gram);
    let factor = sigma * sigma / (scale * scale);
    let covariance = |i: usize, j: usize| inverted[i][j] * factor;
    // Lag-1 autocorrelation of the residual along the calls.
    let errors = residuals(design, fractions, observed);
    let calls = rows / 4;
    let mut means = [0.0; 4];
    for row in 0..calls {
        for channel in 0..4 {
            means[channel] += errors[row * 4 + channel] / calls.max(1) as f64;
        }
    }
    let centred = |row: usize, channel: usize| errors[row * 4 + channel] - means[channel];
    let (mut lagged, mut variance) = (0.0, 0.0);
    for row in 0..calls.saturating_sub(1) {
        for channel in 0..4 {
            lagged += centred(row, channel) * centred(row + 1, channel);
            variance += centred(row, channel) * centred(row, channel);
        }
    }
    let autocorrelation = if variance > 0.0 {
        lagged / variance
    } else {
        0.0
    };
    let inflation = if autocorrelation > 0.0 {
        ((1.0 + autocorrelation) / (1.0 - autocorrelation)).sqrt()
    } else {
        1.0
    };
    let marginal = |key: &dyn Fn(Molecule) -> Option<usize>| -> Option<Marginal> {
        let mut lengths = molecules.iter().filter_map(|&m| key(m)).collect::<Vec<_>>();
        lengths.sort_unstable();
        lengths.dedup();
        if lengths.is_empty() {
            return None;
        }
        // Positive columns grouped by the key's length.
        let group = |length: usize| {
            positive
                .iter()
                .enumerate()
                .filter(|&(_, &k)| key(molecules[k]) == Some(length))
                .map(|(index, _)| index)
                .collect::<Vec<_>>()
        };
        let shares = lengths
            .iter()
            .map(|&length| {
                group(length)
                    .iter()
                    .map(|&index| fractions[positive[index]] / scale)
                    .sum::<f64>()
            })
            .collect::<Vec<_>>();
        let variance_of = |weights: &[(usize, f64)]| {
            let mut total = 0.0;
            for &(i, wi) in weights {
                for &(j, wj) in weights {
                    total += wi * wj * covariance(i, j);
                }
            }
            total.max(0.0)
        };
        let errors = lengths
            .iter()
            .map(|&length| {
                let weights = group(length)
                    .into_iter()
                    .map(|i| (i, 1.0))
                    .collect::<Vec<_>>();
                variance_of(&weights).sqrt()
            })
            .collect::<Vec<_>>();
        let mut order = (0..lengths.len())
            .filter(|&index| shares[index] > 0.0)
            .collect::<Vec<_>>();
        order.sort_by(|&left, &right| {
            shares[right]
                .total_cmp(&shares[left])
                .then(lengths[left].cmp(&lengths[right]))
        });
        let dominant = order.first().map(|&index| lengths[index]);
        let next = order.get(1).map(|&index| lengths[index]);
        let z = match (order.first(), order.get(1)) {
            (Some(&top), Some(&second)) => {
                let mut weights = group(lengths[top])
                    .into_iter()
                    .map(|i| (i, 1.0))
                    .collect::<Vec<_>>();
                weights.extend(group(lengths[second]).into_iter().map(|i| (i, -1.0)));
                (shares[top] - shares[second]) / variance_of(&weights).sqrt()
            }
            _ => SINGLE,
        };
        Some(Marginal {
            lengths,
            shares,
            errors,
            dominant,
            next,
            z,
            z_adjusted: z / inflation,
        })
    };
    let empty = || Marginal {
        lengths: Vec::new(),
        shares: Vec::new(),
        errors: Vec::new(),
        dominant: None,
        next: None,
        z: SINGLE,
        z_adjusted: SINGLE,
    };
    TractDecode {
        extent,
        window,
        downstream,
        residual,
        sigma,
        autocorrelation,
        first: marginal(&|m| Some(m.first)).unwrap_or_else(empty),
        total: marginal(&|m| Some(m.length())).unwrap_or_else(empty),
        iterations_capped: false,
        joint_capped: false,
    }
}
