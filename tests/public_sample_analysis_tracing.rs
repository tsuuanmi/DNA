//! Public `sample_analysis` capability reporting each read's stages to the
//! caller's `tracing` dispatcher, although reads run on worker threads.
//!
//! This test lives in its own binary: `tracing` caches callsite interest
//! process-wide, so a concurrent test that runs the same capability without a
//! subscriber can race the registration of this test's scoped dispatcher.

pub mod support;

use std::sync::{Arc, Mutex};

use tempfile::tempdir;
use tracing::Subscriber;
use tracing::span::{Attributes, Id};
use tracing_subscriber::Registry;
use tracing_subscriber::layer::{Context, Layer, SubscriberExt};
use tracing_subscriber::registry::LookupSpan;

use dna::sample_analysis;
use support::{write_abif, write_config, write_reference};

const QUERY: &str = "ACGTCAGTACGATCGTACCTGAGTACGA";

/// Records the name of each span a capability opens.
struct SpanNames(Arc<Mutex<Vec<&'static str>>>);

impl<S> Layer<S> for SpanNames
where
    S: Subscriber + for<'lookup> LookupSpan<'lookup>,
{
    fn on_new_span(&self, attributes: &Attributes<'_>, _id: &Id, _context: Context<'_, S>) {
        if let Ok(mut names) = self.0.lock() {
            names.push(attributes.metadata().name());
        }
    }
}

#[test]
fn sample_analysis_reports_every_reads_stages_to_the_callers_subscriber()
-> Result<(), Box<dyn std::error::Error>> {
    let directory = tempdir()?;
    let reference = directory.path().join("reference.fa");
    let config = directory.path().join("dna.toml");
    let reverse = QUERY
        .chars()
        .rev()
        .map(|base| match base {
            'A' => 'T',
            'C' => 'G',
            'G' => 'C',
            _ => 'A',
        })
        .collect::<String>();
    let alternate = QUERY.replacen("GAT", "GAA", 1);
    let mut traces = Vec::new();
    for (index, sequence) in [QUERY, &reverse, &alternate].into_iter().enumerate() {
        let trace = directory.path().join(format!("read-{index}.ab1"));
        write_abif(&trace, sequence)?;
        traces.push(trace);
    }
    write_reference(&reference, &format!("TTTT{QUERY}CCCC"))?;
    write_config(&config, "linear")?;
    let names = Arc::new(Mutex::new(Vec::new()));
    let subscriber = Registry::default().with(SpanNames(Arc::clone(&names)));

    tracing::subscriber::with_default(subscriber, || {
        sample_analysis::analyze_sanger(&traces, &reference, &config)
    })?;

    let names = names.lock().map_err(|_| "span-name lock poisoned")?.clone();
    for stage in ["basecalling", "signal_processing", "alignment"] {
        assert_eq!(
            names.iter().filter(|name| **name == stage).count(),
            traces.len(),
            "every read reports its {stage} stage"
        );
    }
    assert!(!directory.path().join("logs").exists());
    Ok(())
}
