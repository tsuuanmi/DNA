//! Public `variant_analysis` capability reporting stage progress through `tracing`.
//!
//! This test lives in its own binary: `tracing` caches callsite interest
//! process-wide, so a concurrent test that runs the same capability without a
//! subscriber can race the registration of this test's scoped dispatcher.

mod support;

use std::sync::{Arc, Mutex};

use tempfile::tempdir;
use tracing::Subscriber;
use tracing::span::{Attributes, Id};
use tracing_subscriber::Registry;
use tracing_subscriber::layer::{Context, Layer, SubscriberExt};

use dna::variant_analysis;
use support::{write_abif, write_config, write_reference};

const QUERY: &str = "ACGTCAGTACGATCGTACCTGAGTACGA";

/// Records the names of the spans a capability opens.
struct StageNames(Arc<Mutex<Vec<&'static str>>>);

impl<S: Subscriber> Layer<S> for StageNames {
    fn on_new_span(&self, attributes: &Attributes<'_>, _id: &Id, _context: Context<'_, S>) {
        if let Ok(mut names) = self.0.lock() {
            names.push(attributes.metadata().name());
        }
    }
}

#[test]
fn sanger_analysis_reports_stage_spans_to_the_callers_subscriber_without_writing_logs()
-> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let trace = directory.path().join("trace.ab1");
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    write_abif(&trace, QUERY)?;
    write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
    write_config(&config, "linear")?;
    let names = Arc::new(Mutex::new(Vec::new()));
    let subscriber = Registry::default().with(StageNames(Arc::clone(&names)));

    tracing::subscriber::with_default(subscriber, || {
        variant_analysis::analyze_sanger(&trace, &reference, &config)
    })?;

    let names = names
        .lock()
        .map_err(|_| "stage-name lock poisoned")?
        .clone();
    assert_eq!(
        names,
        [
            "basecalling",
            "signal_processing",
            "quality_control",
            "alignment",
            "variant_calling"
        ]
    );
    assert!(!directory.path().join("logs").exists());
    Ok(())
}
