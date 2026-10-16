//! Draft review triage for phase 0 of PROP-0004: research-grade stage-1
//! detectors over the evidence DNA already produces, and their score against
//! the truth. It is not production triage and not a DNA contract; it runs
//! only when `dna-eval run --triage` asks for it.
//!
//! The runner takes each sample's compact triage evidence
//! ([`crate::evidence`], stage E); then
//!
//! - [`detect`] runs the detectors and decides each region, trace, and
//!   sample under a policy (stage R), with the minor-signal rules of
//!   [`minor`];
//! - [`score`] measures A/B/C, leaks, catches, and a policy sweep against the
//!   truth (stage C), and writes the `triage*` outputs.

pub(crate) mod detect;
pub(crate) mod minor;
pub(crate) mod repeats;
pub(crate) mod score;
