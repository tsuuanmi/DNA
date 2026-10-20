//! The `[trace]` configuration section: the constants of a read's tract
//! decode (PROP-0012 increment 2), which `dna-trace` computes.

use serde::Deserialize;

use dna_kernel::bounds;
use dna_kernel::error::{ConfigError, Result};
use dna_trace::locate::LocateParams;
use dna_trace::tract::{AnalysisParams, LeakageParams, TractParams};

/// Most consecutive lengths the downstream decode may join: it weighs
/// 4^K windows per call.
const MAX_JOINT_LENGTHS: usize = 8;

/// The tract decode's settings.
#[derive(Debug, Clone, Copy)]
pub struct TraceConfig {
    /// The decode's constants, as `dna-trace` takes them.
    pub(crate) analysis: AnalysisParams,
    /// Margin, in standard errors deflated by the read's residual
    /// autocorrelation, below which a marginal's top two lengths are
    /// near-equal.
    pub(crate) near_equal_z: f64,
}

/// The `[trace]` section as written in the configuration.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RawTraceConfig {
    pub(crate) locate_min_run: usize,
    pub(crate) locate_min_before: usize,
    pub(crate) span_below: usize,
    pub(crate) span_above: usize,
    pub(crate) min_length: usize,
    pub(crate) downstream_calls: usize,
    pub(crate) joint_lengths: usize,
    pub(crate) iterations: usize,
    pub(crate) signature_gap: usize,
    pub(crate) leak_min_run: usize,
    pub(crate) leak_pad: usize,
    pub(crate) leak_min_calls: usize,
    pub(crate) signature_calls: usize,
    pub(crate) near_equal_z: f64,
}

impl RawTraceConfig {
    /// Validates the section.
    pub(crate) fn validate(self) -> Result<TraceConfig> {
        for (value, message) in [
            (self.locate_min_run, "trace.locate_min_run must be positive"),
            (self.min_length, "trace.min_length must be positive"),
            (
                self.downstream_calls,
                "trace.downstream_calls must be positive",
            ),
            (self.iterations, "trace.iterations must be positive"),
            (self.leak_min_run, "trace.leak_min_run must be positive"),
            (self.leak_min_calls, "trace.leak_min_calls must be positive"),
            (
                self.signature_calls,
                "trace.signature_calls must be positive",
            ),
        ] {
            if value == 0 {
                return Err(ConfigError::Constraint(message).into());
            }
        }
        if !(1..=MAX_JOINT_LENGTHS).contains(&self.joint_lengths) {
            return Err(
                ConfigError::Constraint("trace.joint_lengths must be between 1 and 8").into(),
            );
        }
        bounds::positive_finite("trace.near_equal_z", self.near_equal_z)?;
        Ok(TraceConfig {
            analysis: AnalysisParams {
                locate: LocateParams {
                    min_run: self.locate_min_run,
                    min_before: self.locate_min_before,
                },
                tract: TractParams {
                    span_below: self.span_below,
                    span_above: self.span_above,
                    min_length: self.min_length,
                    downstream: self.downstream_calls,
                    joint_lengths: self.joint_lengths,
                    iterations: self.iterations,
                },
                leakage: LeakageParams {
                    gap: self.signature_gap,
                    min_run: self.leak_min_run,
                    pad: self.leak_pad,
                    min_rows: self.leak_min_calls,
                    context: self.signature_calls,
                },
            },
            near_equal_z: self.near_equal_z,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::RawTraceConfig;

    fn raw() -> RawTraceConfig {
        RawTraceConfig {
            locate_min_run: 4,
            locate_min_before: 2,
            span_below: 10,
            span_above: 3,
            min_length: 3,
            downstream_calls: 60,
            joint_lengths: 6,
            iterations: 6,
            signature_calls: 200,
            signature_gap: 3,
            leak_min_run: 3,
            leak_pad: 2,
            leak_min_calls: 10,
            near_equal_z: 2.0,
        }
    }

    // SRS-CALL-025
    #[test]
    fn validates_the_trace_section() {
        assert!(raw().validate().is_ok());
        let invalid: [fn(&mut RawTraceConfig); 6] = [
            |raw| raw.locate_min_run = 0,
            |raw| raw.downstream_calls = 0,
            |raw| raw.joint_lengths = 0,
            |raw| raw.joint_lengths = 9,
            |raw| raw.near_equal_z = 0.0,
            |raw| raw.near_equal_z = f64::NAN,
        ];
        for change in invalid {
            let mut config = raw();
            change(&mut config);
            assert!(config.validate().is_err());
        }
    }
}
