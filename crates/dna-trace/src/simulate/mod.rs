//! A forward model of the Sanger electropherogram and its renderer
//! (PROP-0015 increment 1): test support for the tract decoder's simulated
//! bench, behind the `simulate` feature and never on a runtime path.
//!
//! It ports the research simulator (`tools/python/research/trace_simulator`,
//! `model.py`). A read is a mixture of templates; every template's fragment
//! that ends `n` bases after the primer migrates to the same grid point, so a
//! length mixture shares one peak grid and differs only in the base each
//! template shows there. Per peak:
//!
//! - centre: the grid (spacing and its drift), a per-dye offset, a
//!   trinucleotide context offset, a sequence-determined pattern smoothed
//!   along the template, and a per-read residual;
//! - log height: an envelope along the read, a per-dye gain, a trinucleotide
//!   context effect, a sequence-determined pattern, and a per-read residual,
//!   plus the log of the template's fraction;
//! - shape: a Gaussian whose sigma is a line in the call index.
//!
//! The analysed channels add minor peaks on the grid, spectral cross-talk, a
//! slowly varying floor, white noise and a noisy read start, and are clipped
//! at zero. Calls come from an oracle-grid basecaller: one call per grid
//! point, at the strongest apex within a third of a spacing. A read's
//! parameters are drawn from an aggregate table ([`Table`]); the crate's
//! boundary holds: plain data in, no instrument format, reference, or truth.

mod random;
mod table;

pub use random::{Rng, hash_normal};
pub use table::{PrimerTable, Table, TableError};

/// The channel order of the rendered data.
const BASES: [u8; 4] = *b"ACGT";

/// One read's physical parameters (samples, RFU, natural-log heights).
#[derive(Debug, Clone, PartialEq)]
pub struct Params {
    /// The sample of call 0.
    pub first_sample: f64,
    /// Samples per call at call 0.
    pub spacing0: f64,
    /// The spacing's change per call.
    pub spacing_slope: f64,
    /// Gaussian sigma (samples) at call 0.
    pub sigma0: f64,
    /// Sigma's change per call.
    pub sigma_slope: f64,
    /// Per-dye position offsets (samples), A/C/G/T.
    pub dye_offset: [f64; 4],
    /// The log envelope's knots, `(call, log height)`, by call.
    pub envelope: Vec<(f64, f64)>,
    /// Per-dye log gains, A/C/G/T.
    pub dye_gain: [f64; 4],
    /// Log-height effect of a trinucleotide (previous, own, next base).
    pub height_context: std::collections::BTreeMap<String, f64>,
    /// Position effect (samples) of a trinucleotide.
    pub position_context: std::collections::BTreeMap<String, f64>,
    /// SD of the sequence-determined log-height pattern.
    pub height_det_sd: f64,
    /// SD of the per-read log-height residual.
    pub height_iid_sd: f64,
    /// SD of the sequence-determined position pattern (samples).
    pub position_det_sd: f64,
    /// SD of the per-read position residual (samples).
    pub position_iid_sd: f64,
    /// The position pattern's correlation length (calls).
    pub position_smooth: f64,
    /// Log of a minor peak's height relative to the envelope.
    pub minor_mu: f64,
    /// Its spread.
    pub minor_sd: f64,
    /// Spectral cross-talk, `crosstalk[from][to]`.
    pub crosstalk: [[f64; 4]; 4],
    /// Each channel's floor (RFU).
    pub floor: [f64; 4],
    /// The slowly varying floor's SD (RFU).
    pub floor_sd: f64,
    /// Its Gaussian smoothing (samples).
    pub floor_len: f64,
    /// White noise SD (RFU).
    pub white_sd: f64,
    /// The seed of the read's own random terms.
    pub seed: u64,
    /// The seed of the sequence-determined patterns, shared across reads.
    pub pattern_seed: u64,
    /// Calls of unresolved signal at the read start.
    pub start_calls: usize,
    /// Its level, relative to the envelope at call 0.
    pub start_noise: f64,
    /// The vendor quality table: primary shares, ascending.
    pub qv_share: Vec<f64>,
    /// The quality at each share.
    pub qv_value: Vec<f64>,
}

impl Params {
    /// The same read with its noise (white, floor, minor peaks) times
    /// `factor`.
    #[must_use]
    pub fn scaled_noise(&self, factor: f64) -> Self {
        Self {
            white_sd: self.white_sd * factor,
            floor_sd: self.floor_sd * factor,
            minor_mu: if factor > 0.0 {
                self.minor_mu + factor.ln()
            } else {
                -50.0
            },
            ..self.clone()
        }
    }

    /// The sample of call `n`.
    #[must_use]
    pub fn grid(&self, n: f64) -> f64 {
        (0.5 * self.spacing_slope * n).mul_add(n, self.spacing0.mul_add(n, self.first_sample))
    }

    fn spacing(&self, n: f64) -> f64 {
        self.spacing_slope.mul_add(n, self.spacing0)
    }

    fn sigma(&self, n: f64) -> f64 {
        self.sigma_slope.mul_add(n, self.sigma0).max(0.8)
    }

    /// The log envelope at call `n`: linear between knots, constant before
    /// the first, and extended past the last by the last segment's slope
    /// when it falls.
    #[must_use]
    pub fn log_envelope(&self, n: f64) -> f64 {
        let knots = &self.envelope;
        let Some(&(first_call, first_value)) = knots.first() else {
            return 0.0;
        };
        if n <= first_call {
            return first_value;
        }
        for pair in knots.windows(2) {
            let ((x0, y0), (x1, y1)) = (pair[0], pair[1]);
            if n <= x1 {
                return y0 + (y1 - y0) * (n - x0) / (x1 - x0);
            }
        }
        let (x1, y1) = knots[knots.len() - 1];
        let slope = if knots.len() >= 2 {
            let (x0, y0) = knots[knots.len() - 2];
            ((y1 - y0) / (x1 - x0)).min(0.0)
        } else {
            0.0
        };
        slope.mul_add(n - x1, y1)
    }
}

/// An optional artefact.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Artefact {
    /// Four narrow channels at once, `height` times the local envelope.
    Spike {
        /// Where (grid call index).
        call: f64,
        /// Relative height.
        height: f64,
    },
    /// One wide channel.
    Blob {
        /// Where.
        call: f64,
        /// The channel, A/C/G/T as 0-3.
        channel: usize,
        /// Relative height.
        height: f64,
        /// Full width at half maximum (calls).
        width: f64,
    },
    /// A burst of noise in every channel over `width` calls.
    Noise {
        /// Where it starts.
        call: f64,
        /// Relative level.
        height: f64,
        /// Its length (calls).
        width: f64,
    },
}

/// A rendered read: its analysed channels and the oracle basecaller's calls.
#[derive(Debug, Clone, PartialEq)]
pub struct Rendered {
    /// The analysed channels, A/C/G/T, clipped to `0..=32767`.
    pub channels: [Vec<f64>; 4],
    /// Each call's sample.
    pub ploc: Vec<usize>,
    /// Each call's base.
    pub pbas: Vec<u8>,
    /// Each call's vendor quality.
    pub pcon: Vec<u8>,
    /// Each call's grid position.
    pub grid: Vec<f64>,
}

/// The whole part of a non-negative sample or call position, 0 below 0, as
/// the research simulator's `int()` takes it; positions stay far below
/// `usize::MAX`.
#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn whole(value: f64) -> usize {
    value.max(0.0) as usize
}

/// A vendor quality, `0..=255`, from its interpolated value.
#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn quality(value: f64) -> u8 {
    value.round_ties_even().clamp(0.0, 255.0) as u8
}

/// Adds a Gaussian peak to a channel.
fn add(channel: &mut [f64], centre: f64, height: f64, width: f64) {
    let lo = whole(centre - 5.0 * width);
    let hi = (whole(centre + 5.0 * width) + 2).min(channel.len());
    for (sample, value) in channel.iter_mut().enumerate().take(hi).skip(lo) {
        let z = (sample as f64 - centre) / width;
        *value += height * (-0.5 * z * z).exp();
    }
}

/// `values` convolved with a symmetric `kernel` (odd length), the output
/// centred on the input and as long, zero beyond its ends.
fn convolve_same(values: &[f64], kernel: &[f64]) -> Vec<f64> {
    let half = kernel.len() / 2;
    (0..values.len())
        .map(|index| {
            kernel
                .iter()
                .enumerate()
                .filter_map(|(k, weight)| {
                    (index + k)
                        .checked_sub(half)
                        .and_then(|at| values.get(at))
                        .map(|value| value * weight)
                })
                .sum()
        })
        .collect()
}

/// A Gaussian kernel of `width` over `-reach..=reach`, scaled to unit
/// sum of squares, so white input keeps unit variance.
fn kernel(width: f64, reach: usize) -> Vec<f64> {
    let raw = (0..=2 * reach)
        .map(|k| {
            let z = (k as f64 - reach as f64) / width;
            (-0.5 * z * z).exp()
        })
        .collect::<Vec<_>>();
    let norm = raw.iter().map(|value| value * value).sum::<f64>().sqrt();
    raw.into_iter().map(|value| value / norm).collect()
}

/// The template's bases `from..to`, `N` outside it.
fn slice(padded: &[u8], from: usize, to: usize) -> &[u8] {
    &padded[from.min(padded.len())..to.min(padded.len())]
}

/// Renders the analysed channels of a mixture of `templates` (read-order
/// sequence, fraction) over `calls` calls (the shortest template's length
/// when `None`), with `artefacts`.
#[must_use]
pub fn render(
    templates: &[(&[u8], f64)],
    params: &Params,
    calls: Option<usize>,
    artefacts: &[Artefact],
) -> Rendered {
    let mut rng = Rng::new(params.seed);
    let total = templates.iter().map(|(_, fraction)| fraction).sum::<f64>();
    let count = calls.unwrap_or_else(|| {
        templates
            .iter()
            .map(|(template, _)| template.len())
            .min()
            .unwrap_or(0)
    });
    let grid = (0..count)
        .map(|n| params.grid(n as f64))
        .collect::<Vec<_>>();
    let sigma = (0..count)
        .map(|n| params.sigma(n as f64))
        .collect::<Vec<_>>();
    let envelope = (0..count)
        .map(|n| params.log_envelope(n as f64))
        .collect::<Vec<_>>();
    let last = grid.last().copied().unwrap_or(params.first_sample);
    let samples = whole(last + 6.0 * params.spacing0).max(1);
    let mut y: [Vec<f64>; 4] = std::array::from_fn(|_| vec![0.0; samples]);

    for &(template, fraction) in templates {
        if fraction <= 0.0 || total <= 0.0 {
            continue;
        }
        let share = fraction / total;
        let mut padded = Vec::with_capacity(template.len() + 9);
        padded.push(b'N');
        padded.extend_from_slice(template);
        padded.extend_from_slice(b"NNNNNNNN");
        // Sequence-determined mobility: a hashed value per context, smoothed
        // along the template so neighbouring peaks drift together.
        let hashed = (0..count)
            .map(|i| {
                let mut context = vec![b'p'];
                context.extend_from_slice(slice(&padded, i.saturating_sub(4), i + 3));
                hash_normal(&context, params.pattern_seed)
            })
            .collect::<Vec<_>>();
        let smooth = params.position_smooth;
        let drift = if smooth > 0.0 {
            convolve_same(&hashed, &kernel(smooth, whole(3.0 * smooth) + 1))
        } else {
            hashed
        };
        for i in 0..count {
            let Some(b) = template
                .get(i)
                .and_then(|base| BASES.iter().position(|known| known == base))
            else {
                continue;
            };
            let tri = String::from_utf8_lossy(slice(&padded, i, i + 3)).into_owned();
            let centre = grid[i]
                + params.dye_offset[b]
                + params.position_context.get(&tri).copied().unwrap_or(0.0)
                + params.position_det_sd * drift[i]
                + params.position_iid_sd * rng.normal();
            let mut context = vec![b'h'];
            context.extend_from_slice(slice(&padded, i.saturating_sub(4), i + 3));
            let log_height = envelope[i]
                + params.dye_gain[b]
                + params.height_context.get(&tri).copied().unwrap_or(0.0)
                + params.height_det_sd * hash_normal(&context, params.pattern_seed)
                + params.height_iid_sd * rng.normal();
            let height = share * log_height.exp();
            add(&mut y[b], centre, height, sigma[i]);
            for (c, channel) in y.iter_mut().enumerate() {
                let leak = params.crosstalk[b][c];
                if c != b && leak > 0.0 {
                    add(channel, centre, height * leak, sigma[i]);
                }
            }
        }
    }
    // Background of minor peaks on the grid, every channel.
    let minor = (0..count)
        .map(|_| {
            std::array::from_fn::<f64, 4, _>(|_| {
                (params.minor_mu + params.minor_sd * rng.normal()).exp()
            })
        })
        .collect::<Vec<_>>();
    for i in 0..count {
        let level = envelope[i].exp();
        for (c, channel) in y.iter_mut().enumerate() {
            add(
                channel,
                grid[i] + params.dye_offset[c],
                level * minor[i][c],
                sigma[i],
            );
        }
    }
    let mut all = artefacts.to_vec();
    if params.start_calls > 0 && params.start_noise > 0.0 {
        all.push(Artefact::Noise {
            call: -0.5,
            height: params.start_noise,
            width: params.start_calls as f64 + 0.5,
        });
    }
    for artefact in all {
        let at = |call: f64| whole(call).min(count.saturating_sub(1));
        match artefact {
            Artefact::Spike { call, height } => {
                let level = envelope.get(at(call)).copied().unwrap_or(0.0).exp();
                for channel in &mut y {
                    add(channel, params.grid(call), level * height, 0.9);
                }
            }
            Artefact::Blob {
                call,
                channel,
                height,
                width,
            } => {
                let level = envelope.get(at(call)).copied().unwrap_or(0.0).exp();
                if let Some(channel) = y.get_mut(channel) {
                    add(
                        channel,
                        params.grid(call),
                        level * height,
                        width * params.spacing0 / 2.355,
                    );
                }
            }
            Artefact::Noise {
                call,
                height,
                width,
            } => {
                let level = envelope.get(at(call)).copied().unwrap_or(0.0).exp();
                let lo = whole(params.grid(call));
                let hi = whole(params.grid(call + width)).min(samples);
                for channel in &mut y {
                    for value in channel.iter_mut().take(hi).skip(lo) {
                        *value += level * height * rng.normal().abs();
                    }
                }
            }
        }
    }
    // Floor: a slowly varying component plus white noise, clipped at zero.
    if params.floor_sd > 0.0 {
        let floor_kernel = kernel(params.floor_len, whole(4.0 * params.floor_len));
        for (c, channel) in y.iter_mut().enumerate() {
            let noise = (0..samples).map(|_| rng.normal()).collect::<Vec<_>>();
            let slow = convolve_same(&noise, &floor_kernel);
            for (value, slow) in channel.iter_mut().zip(slow) {
                *value += params.floor_sd.mul_add(slow, params.floor[c]);
            }
        }
    }
    for channel in &mut y {
        for value in channel.iter_mut() {
            *value = params
                .white_sd
                .mul_add(rng.normal(), *value)
                .clamp(0.0, 32767.0);
        }
    }
    let spacing = (0..count)
        .map(|n| params.spacing(n as f64))
        .collect::<Vec<_>>();
    let (ploc, pbas, pcon) = basecall(&y, &grid, &spacing, &params.qv_share, &params.qv_value);
    Rendered {
        channels: y,
        ploc,
        pbas,
        pcon,
        grid,
    }
}

/// `x` interpolated in a table of ascending `xs` and their `ys`, clamped at
/// both ends.
fn interpolate(x: f64, xs: &[f64], ys: &[f64]) -> f64 {
    let (Some(&first), Some(&last)) = (xs.first(), xs.last()) else {
        return 0.0;
    };
    if x <= first {
        return ys[0];
    }
    if x >= last {
        return ys[ys.len() - 1];
    }
    let upper = xs.partition_point(|&value| value < x);
    let (x0, x1, y0, y1) = (xs[upper - 1], xs[upper], ys[upper - 1], ys[upper]);
    y0 + (y1 - y0) * (x - x0) / (x1 - x0)
}

/// The oracle-grid basecaller: one call per grid point, at the strongest
/// apex within a third of a spacing of it (never before the previous call),
/// its base the strongest channel there and its quality the vendor's for
/// that primary share.
#[must_use]
pub fn basecall(
    y: &[Vec<f64>; 4],
    grid: &[f64],
    spacing: &[f64],
    qv_share: &[f64],
    qv_value: &[f64],
) -> (Vec<usize>, Vec<u8>, Vec<u8>) {
    let samples = y[0].len();
    let top = (0..samples)
        .map(|sample| {
            y.iter()
                .map(|channel| channel[sample])
                .fold(f64::MIN, f64::max)
        })
        .collect::<Vec<_>>();
    let mut ploc = Vec::with_capacity(grid.len());
    let mut previous = 0;
    for (&g, &s) in grid.iter().zip(spacing) {
        let lo = whole((g - s / 3.0).floor());
        let hi = (whole((g + s / 3.0).ceil()) + 1).min(samples);
        let mut best = lo.min(samples.saturating_sub(1));
        for sample in lo..hi {
            if top[sample] > top[best] {
                best = sample;
            }
        }
        previous = previous.max(best);
        ploc.push(previous);
    }
    let mut pbas = Vec::with_capacity(ploc.len());
    let mut pcon = Vec::with_capacity(ploc.len());
    for &sample in &ploc {
        let at: [f64; 4] = std::array::from_fn(|c| y[c].get(sample).copied().unwrap_or(0.0));
        let mut best = 0;
        for c in 1..4 {
            if at[c] > at[best] {
                best = c;
            }
        }
        let max = at[best];
        let share = max / at.iter().sum::<f64>().max(1e-9);
        pbas.push(BASES[best]);
        pcon.push(if max <= 0.0 {
            1
        } else {
            quality(interpolate(share, qv_share, qv_value))
        });
    }
    (ploc, pbas, pcon)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn params(table: &Table, seed: u64) -> Params {
        let mut rng = Rng::new(seed);
        table
            .params("HV2F", &mut rng)
            .unwrap_or_else(|| panic!("the table holds HV2F"))
    }

    // PROP-0015 increment 1
    #[test]
    fn renders_a_template_the_oracle_basecaller_reads_back() {
        let table = Table::pooled().unwrap_or_else(|error| panic!("{error}"));
        let template = b"ACGTTGCAACGGTACCATGGCTAGCTAGGATCCATGCAAGTCGATCGTACGATCGATCGGATC".repeat(4);
        let rendered = render(&[(&template, 1.0)], &params(&table, 3), None, &[]);
        assert_eq!(rendered.pbas.len(), template.len());
        // The read start is noise; past it the calls are the template's.
        let agree = rendered
            .pbas
            .iter()
            .zip(&template)
            .skip(20)
            .filter(|(call, base)| call == base)
            .count();
        assert!(
            agree >= template.len() - 20 - 2,
            "{agree} of {}",
            template.len() - 20
        );
        assert!(rendered.ploc.windows(2).all(|pair| pair[0] <= pair[1]));
        // The same parameters render the same read.
        assert_eq!(
            render(&[(&template, 1.0)], &params(&table, 3), None, &[]),
            rendered
        );
    }

    #[test]
    fn convolution_keeps_length_and_centre() {
        let out = convolve_same(&[0.0, 0.0, 1.0, 0.0, 0.0], &[0.25, 0.5, 0.25]);
        assert_eq!(out, vec![0.0, 0.25, 0.5, 0.25, 0.0]);
    }
}
