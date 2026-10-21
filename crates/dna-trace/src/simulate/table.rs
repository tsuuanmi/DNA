//! The simulator's aggregate parameter table: per primer, quantiles over
//! the fitted reads of each read parameter, never a read's own values; the
//! pooled context tables; and the calibration constants. Built by the
//! research script `table.py` from the research fit (PROP-0015 increment 1).

use std::collections::BTreeMap;

use serde::Deserialize;

use super::{Params, Rng};

/// The committed table.
const POOLED: &str = include_str!("../../tables/simulator-params.json");

/// A table that cannot be read.
#[derive(Debug, thiserror::Error)]
pub enum TableError {
    /// The JSON is not a table.
    #[error("invalid simulator parameter table: {0}")]
    Parse(#[from] serde_json::Error),
    /// A quantile list does not match the levels.
    #[error("simulator parameter table: {0}")]
    Shape(&'static str),
}

/// One primer's quantiles over reads.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PrimerTable {
    /// Reads the quantiles are over.
    pub reads: usize,
    /// Samples per call at call 0.
    pub spacing0: Vec<f64>,
    /// The spacing's change per call.
    pub spacing_slope: Vec<f64>,
    /// Sigma at call 0.
    pub sigma0: Vec<f64>,
    /// Sigma's change per call.
    pub sigma_slope: Vec<f64>,
    /// SD of the sequence-determined log-height pattern.
    pub height_det_sd: Vec<f64>,
    /// SD of the per-read log-height residual.
    pub height_iid_sd: Vec<f64>,
    /// SD of the sequence-determined position pattern.
    pub position_det_sd: Vec<f64>,
    /// SD of the per-read position residual.
    pub position_iid_sd: Vec<f64>,
    /// The slowly varying floor's SD.
    pub floor_sd: Vec<f64>,
    /// White noise SD.
    pub white_sd: Vec<f64>,
    /// Per-dye position offsets.
    pub dye_offset: [Vec<f64>; 4],
    /// Per-dye log gains.
    pub dye_gain: [Vec<f64>; 4],
    /// Per-channel floors.
    pub floor: [Vec<f64>; 4],
    /// The log envelope's quantiles at each of the table's envelope calls.
    pub envelope: Vec<Vec<f64>>,
    /// Spectral cross-talk from one dye into another, by `"<from><to>"`
    /// (A/C/G/T as 0-3), off the diagonal.
    pub crosstalk: BTreeMap<String, Vec<f64>>,
    /// The dependence between a read's parameters.
    pub copula: Copula,
}

/// A Gaussian copula: the correlation of the parameters' normal scores
/// (ranks) over reads.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Copula {
    /// The parameters it ties, by name: `spacing0`, `spacing_slope`,
    /// `sigma0`, `sigma_slope`, `height_det_sd`, `position_det_sd`,
    /// `floor_sd`, `white_sd`, `dye_offset.<c>`, `dye_gain.<c>`,
    /// `floor.<c>`, and `envelope.<call>`.
    pub parameters: Vec<String>,
    /// Their correlation matrix.
    pub correlation: Vec<Vec<f64>>,
}

/// The weight on the identity that keeps a rounded correlation matrix
/// positive definite for its Cholesky factor.
const RIDGE: f64 = 1e-3;

/// The lower Cholesky factor of `(matrix + ridge * I) / (1 + ridge)`;
/// `None` when it is not positive definite even so.
fn cholesky(matrix: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let n = matrix.len();
    let mut lower = vec![vec![0.0; n]; n];
    for i in 0..n {
        for j in 0..=i {
            let target = (matrix[i][j] + if i == j { RIDGE } else { 0.0 }) / (1.0 + RIDGE);
            let sum = (0..j).map(|k| lower[i][k] * lower[j][k]).sum::<f64>();
            if i == j {
                let diagonal = target - sum;
                if diagonal <= 0.0 {
                    return None;
                }
                lower[i][j] = diagonal.sqrt();
            } else {
                lower[i][j] = (target - sum) / lower[j][j];
            }
        }
    }
    Some(lower)
}

/// The standard normal distribution function (complementary error function
/// by Chebyshev fitting, relative error below 1.2e-7).
fn normal_cdf(x: f64) -> f64 {
    let z = x.abs() / std::f64::consts::SQRT_2;
    let t = 1.0 / (1.0 + 0.5 * z);
    let poly = [
        -1.265_512_23,
        1.000_023_68,
        0.374_091_96,
        0.096_784_18,
        -0.186_288_06,
        0.278_868_07,
        -1.135_203_98,
        1.488_515_87,
        -0.822_152_23,
        0.170_872_77,
    ]
    .iter()
    .rev()
    .fold(0.0_f64, |acc, coefficient| acc.mul_add(t, *coefficient));
    let erfc = t * (-z * z + poly).exp();
    if x >= 0.0 {
        1.0 - 0.5 * erfc
    } else {
        0.5 * erfc
    }
}

/// The simulator's parameter table.
#[derive(Debug, Clone, Deserialize)]
pub struct Table {
    /// The table's format version.
    pub version: u32,
    /// Where it comes from.
    pub provenance: BTreeMap<String, serde_json::Value>,
    /// The quantile levels, ascending in `(0, 1)`.
    pub levels: Vec<f64>,
    /// The calls of the envelope's knots.
    pub envelope_calls: Vec<f64>,
    /// The sample of call 0.
    pub first_sample: f64,
    /// Log of a minor peak's height relative to the envelope.
    pub minor_mu: f64,
    /// Its spread.
    pub minor_sd: f64,
    /// The floor's smoothing (samples).
    pub floor_len: f64,
    /// The position pattern's correlation length (calls).
    pub position_smooth: f64,
    /// Calls of unresolved signal at the read start.
    pub start_calls: usize,
    /// Its level relative to the envelope.
    pub start_noise: f64,
    /// The vendor quality table's shares.
    pub qv_share: Vec<f64>,
    /// Its qualities.
    pub qv_value: Vec<f64>,
    /// Pooled trinucleotide log-height effects, centred per own base.
    pub height_context: BTreeMap<String, f64>,
    /// Pooled trinucleotide position effects, centred per own base.
    pub position_context: BTreeMap<String, f64>,
    /// Per primer.
    pub primers: BTreeMap<String, PrimerTable>,
    /// Each primer's copula factor, computed on reading.
    #[serde(skip)]
    factors: BTreeMap<String, Vec<Vec<f64>>>,
}

impl Table {
    /// The committed table.
    ///
    /// # Errors
    ///
    /// When the committed data are not a valid table.
    pub fn pooled() -> Result<Self, TableError> {
        Self::parse(POOLED)
    }

    /// A table from its JSON.
    ///
    /// # Errors
    ///
    /// When `text` is not a table, or a quantile list does not match the
    /// levels.
    pub fn parse(text: &str) -> Result<Self, TableError> {
        let mut table: Self = serde_json::from_str(text)?;
        let levels = table.levels.len();
        if levels < 2 || table.levels.windows(2).any(|pair| pair[0] >= pair[1]) {
            return Err(TableError::Shape("levels must ascend"));
        }
        for primer in table.primers.values() {
            let lists = [
                &primer.spacing0,
                &primer.spacing_slope,
                &primer.sigma0,
                &primer.sigma_slope,
                &primer.height_det_sd,
                &primer.height_iid_sd,
                &primer.position_det_sd,
                &primer.position_iid_sd,
                &primer.floor_sd,
                &primer.white_sd,
            ];
            let vectors = primer
                .dye_offset
                .iter()
                .chain(&primer.dye_gain)
                .chain(&primer.floor)
                .chain(&primer.envelope);
            if lists
                .into_iter()
                .chain(vectors)
                .chain(primer.crosstalk.values())
                .any(|list| list.len() != levels)
                || primer.envelope.len() != table.envelope_calls.len()
                || primer.crosstalk.keys().cloned().collect::<Vec<_>>() != crosstalk_keys()
            {
                return Err(TableError::Shape(
                    "a quantile list does not match the levels",
                ));
            }
            let copula = &primer.copula;
            let expected = copula_parameters(&table.envelope_calls);
            if copula.parameters != expected
                || copula.correlation.len() != expected.len()
                || copula
                    .correlation
                    .iter()
                    .any(|row| row.len() != expected.len())
            {
                return Err(TableError::Shape(
                    "the copula does not tie the expected parameters",
                ));
            }
        }
        let mut factors = BTreeMap::new();
        for (name, primer) in &table.primers {
            let factor = cholesky(&primer.copula.correlation).ok_or(TableError::Shape(
                "a copula correlation is not positive definite",
            ))?;
            factors.insert(name.clone(), factor);
        }
        table.factors = factors;
        Ok(table)
    }

    /// The value at quantile `u` of a quantile list, linear between levels
    /// and clamped beyond the outer ones.
    fn at(&self, quantiles: &[f64], u: f64) -> f64 {
        super::interpolate(u, &self.levels, quantiles)
    }

    /// A read of `primer`, its parameters drawn from the table by `rng`;
    /// `None` for a primer the table lacks. The parameters are drawn
    /// together through the primer's Gaussian copula, so their dependence
    /// over real reads is kept (a read with a low envelope has a low floor
    /// too); a pattern's and a residual's SDs share one quantile, since
    /// both scale one fitted residual.
    #[must_use]
    pub fn params(&self, primer: &str, rng: &mut Rng) -> Option<Params> {
        let table = self.primers.get(primer)?;
        let factor = self.factors.get(primer)?;
        let normal = (0..factor.len()).map(|_| rng.normal()).collect::<Vec<_>>();
        let quantile = table
            .copula
            .parameters
            .iter()
            .zip(factor)
            .map(|(name, row)| {
                let z = row
                    .iter()
                    .zip(&normal)
                    .map(|(weight, value)| weight * value)
                    .sum::<f64>();
                (name.as_str(), normal_cdf(z))
            })
            .collect::<BTreeMap<_, _>>();
        let u = |name: &str| quantile.get(name).copied().unwrap_or(0.5);
        let indexed = |prefix: &str, c: usize| u(&format!("{prefix}.{c}"));
        Some(Params {
            first_sample: self.first_sample,
            spacing0: self.at(&table.spacing0, u("spacing0")),
            spacing_slope: self.at(&table.spacing_slope, u("spacing_slope")),
            sigma0: self.at(&table.sigma0, u("sigma0")),
            sigma_slope: self.at(&table.sigma_slope, u("sigma_slope")),
            dye_offset: std::array::from_fn(|c| {
                self.at(&table.dye_offset[c], indexed("dye_offset", c))
            }),
            envelope: self
                .envelope_calls
                .iter()
                .zip(&table.envelope)
                .map(|(&call, quantiles)| {
                    (call, self.at(quantiles, u(&format!("envelope.{call}"))))
                })
                .collect(),
            dye_gain: std::array::from_fn(|c| self.at(&table.dye_gain[c], indexed("dye_gain", c))),
            height_context: self.height_context.clone(),
            position_context: self.position_context.clone(),
            height_det_sd: self.at(&table.height_det_sd, u("height_det_sd")),
            height_iid_sd: self.at(&table.height_iid_sd, u("height_det_sd")),
            position_det_sd: self.at(&table.position_det_sd, u("position_det_sd")),
            position_iid_sd: self.at(&table.position_iid_sd, u("position_det_sd")),
            position_smooth: self.position_smooth,
            minor_mu: self.minor_mu,
            minor_sd: self.minor_sd,
            crosstalk: std::array::from_fn(|from| {
                std::array::from_fn(|to| {
                    let key = format!("{from}{to}");
                    table.crosstalk.get(&key).map_or(0.0, |quantiles| {
                        self.at(quantiles, u(&format!("crosstalk.{key}")))
                    })
                })
            }),
            floor: std::array::from_fn(|c| self.at(&table.floor[c], indexed("floor", c))),
            floor_sd: self.at(&table.floor_sd, u("floor_sd")),
            floor_len: self.floor_len,
            white_sd: self.at(&table.white_sd, u("white_sd")),
            seed: rng.next_u64(),
            pattern_seed: 7,
            start_calls: self.start_calls,
            start_noise: self.start_noise,
            qv_share: self.qv_share.clone(),
            qv_value: self.qv_value.clone(),
        })
    }
}

/// The parameters a copula ties, in the table's order.
fn copula_parameters(envelope_calls: &[f64]) -> Vec<String> {
    let mut names = [
        "spacing0",
        "spacing_slope",
        "sigma0",
        "sigma_slope",
        "height_det_sd",
        "position_det_sd",
        "floor_sd",
        "white_sd",
    ]
    .map(String::from)
    .to_vec();
    for prefix in ["dye_offset", "dye_gain", "floor"] {
        names.extend((0..4).map(|c| format!("{prefix}.{c}")));
    }
    names.extend(envelope_calls.iter().map(|call| format!("envelope.{call}")));
    names.extend(
        crosstalk_keys()
            .into_iter()
            .map(|key| format!("crosstalk.{key}")),
    );
    names
}

/// The cross-talk keys, `"<from><to>"` off the diagonal, in order.
fn crosstalk_keys() -> Vec<String> {
    (0..4)
        .flat_map(|from| {
            (0..4)
                .filter(move |&to| to != from)
                .map(move |to| format!("{from}{to}"))
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_committed_table_holds_every_bench_primer() {
        let table = Table::pooled().unwrap_or_else(|error| panic!("{error}"));
        for primer in ["HV1F", "HV1R", "HV2F", "HV3R"] {
            let mut rng = Rng::new(1);
            let params = table
                .params(primer, &mut rng)
                .unwrap_or_else(|| panic!("{primer} is missing"));
            assert!((11.0..13.0).contains(&params.spacing0), "{params:?}");
            assert!((2.0..4.0).contains(&params.sigma0));
        }
        assert!(table.params("HV9X", &mut Rng::new(1)).is_none());
    }

    #[test]
    fn the_normal_distribution_function_is_accurate() {
        for (x, expected) in [
            (0.0, 0.5),
            (1.0, 0.841_344_746),
            (-1.96, 0.024_997_895),
            (3.0, 0.998_650_102),
        ] {
            assert!((normal_cdf(x) - expected).abs() < 1e-6, "{x}");
        }
        let factor = cholesky(&[vec![1.0, 0.5], vec![0.5, 1.0]]).unwrap_or_default();
        assert!((factor[1][0] - 0.5 / (1.0 + RIDGE)).abs() < 1e-3);
    }
}
