//! Public `sample_analysis` capability reporting each read's stages to the
//! caller's `tracing` dispatcher, inside the caller's span, although reads run
//! on worker threads.
//!
//! This test lives in its own binary: `tracing` caches callsite interest
//! process-wide, so a concurrent test that runs the same capability without a
//! subscriber can race the registration of this test's scoped dispatcher.

pub mod support;

use std::num::NonZeroUsize;
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

/// A span name and the name of the span it was opened inside, if any.
type OpenedSpan = (&'static str, Option<&'static str>);

/// Records each span a capability opens with the span it was opened inside.
struct SpanParents(Arc<Mutex<Vec<OpenedSpan>>>);

impl<S> Layer<S> for SpanParents
where
    S: Subscriber + for<'lookup> LookupSpan<'lookup>,
{
    fn on_new_span(&self, attributes: &Attributes<'_>, _id: &Id, context: Context<'_, S>) {
        let parent = context
            .current_span()
            .metadata()
            .map(tracing::Metadata::name);
        if let Ok(mut spans) = self.0.lock() {
            spans.push((attributes.metadata().name(), parent));
        }
    }
}

// SRS-API-011
#[test]
fn sample_analysis_reports_every_reads_stages_inside_the_callers_span()
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
    let workers = NonZeroUsize::new(traces.len()).ok_or("at least one trace")?;
    let spans = Arc::new(Mutex::new(Vec::new()));
    let subscriber = Registry::default().with(SpanParents(Arc::clone(&spans)));

    tracing::subscriber::with_default(subscriber, || {
        let _caller = tracing::info_span!("caller").entered();
        sample_analysis::analyze_sanger(&traces, &reference, &config, workers)
    })?;

    let spans = spans.lock().map_err(|_| "span lock poisoned")?.clone();
    for stage in ["basecalling", "signal_processing", "alignment"] {
        let parents = spans
            .iter()
            .filter(|(name, _)| *name == stage)
            .map(|(_, parent)| *parent)
            .collect::<Vec<_>>();
        assert_eq!(
            parents,
            vec![Some("caller"); traces.len()],
            "every read opens its {stage} stage inside the caller's span"
        );
    }
    assert!(!directory.path().join("logs").exists());
    Ok(())
}
