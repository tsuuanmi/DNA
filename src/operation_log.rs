//! Per-operation append-only operational log.
//!
//! Scientific stages emit structured [`tracing`] events and enter one span per
//! stage; they never choose a destination. A CLI operation opens an
//! [`OperationLog`] for its trace or sample, runs inside [`OperationLog::in_scope`],
//! and the [`RecordLayer`] renders each crate event at `INFO` or above as one
//! physical line:
//!
//! ```text
//! 2026-10-08 12:00:00.000 | INFO     | dna::pipeline::analyze:20 - run_id=<id> event=analysis_started …
//! ```
//!
//! Fields render as `name=value` in macro order, `?`-formatted values keep
//! their `Debug` quoting, and a free-text message follows the fields. Record
//! delimiters and control characters are escaped so one event is one line.
//! Write failures are captured and surfaced by [`OperationLog::check`] or
//! [`OperationLog::sync`], because a `tracing` layer cannot return errors.

use std::env;
use std::fmt::{self, Write as _};
use std::fs::{self, File, OpenOptions};
use std::io::{self, Write as _};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex, MutexGuard, PoisonError};

use chrono::Local;
use tracing::field::{Field, Visit};
use tracing::span::Id;
use tracing::{Dispatch, Event, Level, Metadata, Subscriber};
use tracing_subscriber::layer::{Context, Layer, SubscriberExt};
use tracing_subscriber::registry::{LookupSpan, Registry};

use crate::error::{Error, Result};

const DEFAULT_LOG_DIRECTORY: &str = "logs";
const CRATE_TARGET: &str = env!("CARGO_CRATE_NAME");
static RUN_SEQUENCE: AtomicU64 = AtomicU64::new(0);

/// One append-only log file and the dispatcher that writes to it.
pub(crate) struct OperationLog {
    path: PathBuf,
    state: Arc<Mutex<LogState>>,
    dispatch: Dispatch,
}

/// Mutable state shared between the operation and its layer.
struct LogState {
    file: File,
    write_error: Option<io::Error>,
    stage: Option<&'static str>,
    failed_stage: Option<&'static str>,
}

impl OperationLog {
    /// Opens `<stem>.log` under `DNA_LOG_DIR`, or `logs/` when it is unset.
    pub(crate) fn open(stem: &str) -> Result<Self> {
        let directory = env::var_os("DNA_LOG_DIR")
            .map_or_else(|| PathBuf::from(DEFAULT_LOG_DIRECTORY), PathBuf::from);
        Self::open_in(&directory, stem)
    }

    fn open_in(directory: &Path, stem: &str) -> Result<Self> {
        if directory.as_os_str().is_empty() {
            return Err(Error::Path {
                kind: "log directory",
                path: directory.to_path_buf(),
                reason: "path must be non-empty",
            });
        }
        fs::create_dir_all(directory).map_err(|source| Error::Log {
            path: directory.to_path_buf(),
            source,
        })?;
        let path = directory.join(format!("{stem}.log"));
        let file = OpenOptions::new()
            .create(true)
            .append(true)
            .open(&path)
            .map_err(|source| Error::Log {
                path: path.clone(),
                source,
            })?;
        let state = Arc::new(Mutex::new(LogState {
            file,
            write_error: None,
            stage: None,
            failed_stage: None,
        }));
        let layer = RecordLayer {
            run_id: new_run_id(),
            state: Arc::clone(&state),
        };
        let dispatch = Dispatch::new(Registry::default().with(layer));
        Ok(Self {
            path,
            state,
            dispatch,
        })
    }

    /// Runs `operation` with this log as the current thread's subscriber.
    pub(crate) fn in_scope<T>(&self, operation: impl FnOnce() -> T) -> T {
        tracing::dispatcher::with_default(&self.dispatch, operation)
    }

    /// The stage an operation failure belongs to: the stage whose record first
    /// failed to write, otherwise the most recently entered stage span.
    pub(crate) fn stage(&self) -> Option<&'static str> {
        let state = lock(&self.state);
        state.failed_stage.or(state.stage)
    }

    /// Surfaces, once, the first record write that failed since the last check.
    pub(crate) fn check(&self) -> Result<()> {
        match lock(&self.state).write_error.take() {
            Some(source) => Err(self.log_error(source)),
            None => Ok(()),
        }
    }

    /// Surfaces any write failure, then durably synchronizes all records.
    pub(crate) fn sync(&self) -> Result<()> {
        self.check()?;
        let state = lock(&self.state);
        state
            .file
            .sync_all()
            .map_err(|source| self.log_error(source))
    }

    fn log_error(&self, source: io::Error) -> Error {
        Error::Log {
            path: self.path.clone(),
            source,
        }
    }
}

/// `tracing` layer rendering crate events into the operation's log file.
struct RecordLayer {
    run_id: String,
    state: Arc<Mutex<LogState>>,
}

impl<S> Layer<S> for RecordLayer
where
    S: Subscriber + for<'lookup> LookupSpan<'lookup>,
{
    fn enabled(&self, metadata: &Metadata<'_>, _context: Context<'_, S>) -> bool {
        is_operational(metadata)
    }

    fn on_enter(&self, id: &Id, context: Context<'_, S>) {
        if let Some(span) = context.span(id) {
            lock(&self.state).stage = Some(span.name());
        }
    }

    fn on_event(&self, event: &Event<'_>, _context: Context<'_, S>) {
        let record = self.render(event);
        let mut state = lock(&self.state);
        if let Err(error) = state.file.write_all(record.as_bytes())
            && state.write_error.is_none()
        {
            state.write_error = Some(error);
            state.failed_stage = state.failed_stage.or(state.stage);
        }
    }
}

impl RecordLayer {
    fn render(&self, event: &Event<'_>) -> String {
        let metadata = event.metadata();
        let mut fields = RecordFields::default();
        event.record(&mut fields);
        format!(
            "{} | {:<8} | {}:{} - run_id={} {}\n",
            Local::now().format("%Y-%m-%d %H:%M:%S%.3f"),
            metadata.level().as_str(),
            metadata.module_path().unwrap_or_else(|| metadata.target()),
            metadata.line().unwrap_or_default(),
            self.run_id,
            fields.finish()
        )
    }
}

/// Selects this crate's spans and events at `INFO` or above.
fn is_operational(metadata: &Metadata<'_>) -> bool {
    let target = metadata.target();
    *metadata.level() <= Level::INFO
        && target
            .strip_prefix(CRATE_TARGET)
            .is_some_and(|rest| rest.is_empty() || rest.starts_with("::"))
}

/// Accumulates `name=value` pairs and an optional trailing message.
#[derive(Default)]
struct RecordFields {
    pairs: String,
    message: Option<String>,
}

impl RecordFields {
    fn push(&mut self, field: &Field, value: fmt::Arguments<'_>) {
        let value = escape_record_field(&value.to_string());
        if field.name() == "message" {
            self.message = Some(value);
            return;
        }
        if !self.pairs.is_empty() {
            self.pairs.push(' ');
        }
        // Writing into a `String` cannot fail.
        let _ = write!(self.pairs, "{}={value}", field.name());
    }

    fn finish(self) -> String {
        match self.message {
            Some(message) if self.pairs.is_empty() => message,
            Some(message) => format!("{} {message}", self.pairs),
            None => self.pairs,
        }
    }
}

impl Visit for RecordFields {
    fn record_f64(&mut self, field: &Field, value: f64) {
        self.push(field, format_args!("{value}"));
    }

    fn record_i64(&mut self, field: &Field, value: i64) {
        self.push(field, format_args!("{value}"));
    }

    fn record_u64(&mut self, field: &Field, value: u64) {
        self.push(field, format_args!("{value}"));
    }

    fn record_i128(&mut self, field: &Field, value: i128) {
        self.push(field, format_args!("{value}"));
    }

    fn record_u128(&mut self, field: &Field, value: u128) {
        self.push(field, format_args!("{value}"));
    }

    fn record_bool(&mut self, field: &Field, value: bool) {
        self.push(field, format_args!("{value}"));
    }

    fn record_str(&mut self, field: &Field, value: &str) {
        self.push(field, format_args!("{value}"));
    }

    fn record_debug(&mut self, field: &Field, value: &dyn fmt::Debug) {
        self.push(field, format_args!("{value:?}"));
    }
}

fn new_run_id() -> String {
    let opened_at = Local::now();
    let sequence = RUN_SEQUENCE.fetch_add(1, Ordering::Relaxed);
    format!(
        "{}-{}-{sequence}",
        opened_at.format("%Y%m%dT%H%M%S%.3f"),
        std::process::id()
    )
}

fn escape_record_field(value: &str) -> String {
    let mut escaped = String::with_capacity(value.len());
    for character in value.chars() {
        match character {
            '\n' => escaped.push_str("\\n"),
            '\r' => escaped.push_str("\\r"),
            '\t' => escaped.push_str("\\t"),
            '|' => escaped.push_str("\\|"),
            character if character.is_control() => {
                escaped.extend(character.escape_unicode());
            }
            character => escaped.push(character),
        }
    }
    escaped
}

/// Locks shared log state; no code path panics while holding the lock, so a
/// poisoned lock still guards consistent state.
fn lock(state: &Mutex<LogState>) -> MutexGuard<'_, LogState> {
    state.lock().unwrap_or_else(PoisonError::into_inner)
}

#[cfg(test)]
mod tests {
    use std::fs;

    use tempfile::tempdir;

    use crate::error::Error;

    use super::*;

    type TestResult = std::result::Result<(), Box<dyn std::error::Error>>;

    fn records(path: &Path) -> std::result::Result<Vec<String>, io::Error> {
        Ok(fs::read_to_string(path)?
            .lines()
            .map(str::to_owned)
            .collect())
    }

    #[test]
    fn appends_one_line_record_per_event_across_operations() -> TestResult {
        let directory = tempdir()?;
        let first = OperationLog::open_in(directory.path(), "trace")?;
        first.in_scope(|| {
            tracing::info!(event = "first", count = 3_usize, name = ?"a.ab1");
            tracing::warn!(event = "second", ratio = %format_args!("{:.4}", 0.5), ok = true);
        });
        first.sync()?;
        let second = OperationLog::open_in(directory.path(), "trace")?;
        second.in_scope(|| tracing::error!(event = "third"));
        second.sync()?;

        let lines = records(&directory.path().join("trace.log"))?;
        assert_eq!(lines.len(), 3);
        let module = module_path!();
        assert!(lines[0].contains(&format!(" | INFO     | {module}:")));
        assert!(lines[0].contains(" - run_id="));
        assert!(lines[0].ends_with(r#" event=first count=3 name="a.ab1""#));
        assert!(lines[1].contains(" | WARN     | "));
        assert!(lines[1].ends_with(" event=second ratio=0.5000 ok=true"));
        assert!(lines[2].contains(" | ERROR    | "));
        assert!(lines[2].ends_with(" event=third"));
        assert_eq!(lines[0].as_bytes().get(4), Some(&b'-'));
        assert_eq!(lines[0].as_bytes().get(23), Some(&b' '));
        let run_id = |line: &str| {
            line.split("run_id=")
                .nth(1)
                .and_then(|rest| rest.split(' ').next())
                .map(str::to_owned)
        };
        assert_eq!(run_id(&lines[0]), run_id(&lines[1]));
        assert_ne!(run_id(&lines[0]), run_id(&lines[2]));
        Ok(())
    }

    #[test]
    fn escapes_record_delimiters_and_control_characters() -> TestResult {
        let directory = tempdir()?;
        let log = OperationLog::open_in(directory.path(), "trace")?;
        log.in_scope(|| tracing::info!(event = "test", value = "left|right\nforged\tfield\r"));
        log.sync()?;

        let text = fs::read_to_string(directory.path().join("trace.log"))?;
        assert_eq!(text.lines().count(), 1);
        assert!(text.contains("value=left\\|right\\nforged\\tfield\\r"));
        Ok(())
    }

    #[test]
    fn renders_message_text_after_structured_fields() -> TestResult {
        let directory = tempdir()?;
        let log = OperationLog::open_in(directory.path(), "trace")?;
        log.in_scope(|| tracing::info!(event = "summary", reads = 2_usize, "a=1 b=2"));
        log.sync()?;

        let lines = records(&directory.path().join("trace.log"))?;
        assert!(lines[0].ends_with(" event=summary reads=2 a=1 b=2"));
        Ok(())
    }

    #[test]
    fn tracks_the_last_entered_stage_span() -> TestResult {
        let directory = tempdir()?;
        let log = OperationLog::open_in(directory.path(), "trace")?;
        assert_eq!(log.stage(), None);
        log.in_scope(|| {
            let _loading = tracing::info_span!("input_loading").entered();
            let _alignment = tracing::info_span!("alignment").entered();
        });
        assert_eq!(log.stage(), Some("alignment"));
        Ok(())
    }

    #[test]
    fn records_only_crate_events_at_info_or_above_inside_the_scope() -> TestResult {
        let directory = tempdir()?;
        let log = OperationLog::open_in(directory.path(), "trace")?;
        tracing::info!(event = "outside_scope");
        log.in_scope(|| {
            tracing::debug!(event = "too_verbose");
            tracing::info!(target: "dependency", event = "foreign");
            let _span = tracing::info_span!(target: "dependency", "foreign_stage").entered();
            tracing::info!(event = "kept");
        });
        log.sync()?;

        let lines = records(&directory.path().join("trace.log"))?;
        assert_eq!(lines.len(), 1);
        assert!(lines[0].ends_with(" event=kept"));
        assert_eq!(log.stage(), None);
        Ok(())
    }

    #[cfg(target_os = "linux")]
    #[test]
    fn surfaces_the_first_write_failure_once() -> TestResult {
        let directory = tempdir()?;
        std::os::unix::fs::symlink("/dev/full", directory.path().join("trace.log"))?;
        let log = OperationLog::open_in(directory.path(), "trace")?;
        log.in_scope(|| tracing::info!(event = "lost"));

        assert!(matches!(log.check(), Err(Error::Log { path, .. }) if path.ends_with("trace.log")));
        assert!(log.check().is_ok());
        Ok(())
    }

    #[cfg(target_os = "linux")]
    #[test]
    fn attributes_a_failed_record_to_the_stage_it_was_written_in() -> TestResult {
        let directory = tempdir()?;
        std::os::unix::fs::symlink("/dev/full", directory.path().join("trace.log"))?;
        let log = OperationLog::open_in(directory.path(), "trace")?;
        log.in_scope(|| {
            let alignment = tracing::info_span!("alignment").entered();
            tracing::info!(event = "alignment_completed");
            drop(alignment);
            let _reporting = tracing::info_span!("reporting").entered();
        });

        assert!(log.check().is_err());
        assert_eq!(log.stage(), Some("alignment"));
        Ok(())
    }

    #[test]
    fn rejects_an_empty_log_directory() {
        assert!(matches!(
            OperationLog::open_in(Path::new(""), "trace"),
            Err(Error::Path {
                kind: "log directory",
                ..
            })
        ));
    }
}
