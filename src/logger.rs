//! Append-only, per-operation operational logging with Apollo-style records.

use std::env;
use std::fmt;
use std::fs::{self, File, OpenOptions};
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use chrono::Local;

use crate::error::{Error, Result};

const DEFAULT_LOG_DIRECTORY: &str = "logs";
static RUN_SEQUENCE: AtomicU64 = AtomicU64::new(0);

#[derive(Clone, Copy)]
enum LogLevel {
    Info,
    Warn,
    Error,
}

impl LogLevel {
    const fn label(self) -> &'static str {
        match self {
            Self::Info => "INFO",
            Self::Warn => "WARN",
            Self::Error => "ERROR",
        }
    }
}

/// One append-only log file for a command operation.
pub(crate) struct Logger {
    path: PathBuf,
    file: File,
    run_id: String,
}

impl Logger {
    /// Opens `logs/<operation-stem>.log`, honoring `DNA_LOG_DIR` when set.
    pub(crate) fn open(operation_stem: &str) -> Result<Self> {
        let directory = env::var_os("DNA_LOG_DIR")
            .map(PathBuf::from)
            .unwrap_or_else(|| PathBuf::from(DEFAULT_LOG_DIRECTORY));
        if directory.as_os_str().is_empty() {
            return Err(Error::Path {
                kind: "log directory",
                path: directory,
                reason: "path must be non-empty".into(),
            });
        }
        Self::open_in(&directory, operation_stem)
    }

    fn open_in(directory: &Path, operation_stem: &str) -> Result<Self> {
        fs::create_dir_all(directory).map_err(|source| Error::Log {
            path: directory.to_path_buf(),
            source,
        })?;
        let path = directory.join(format!("{operation_stem}.log"));
        let file = OpenOptions::new()
            .create(true)
            .append(true)
            .open(&path)
            .map_err(|source| Error::Log {
                path: path.clone(),
                source,
            })?;
        let opened_at = Local::now();
        let sequence = RUN_SEQUENCE.fetch_add(1, Ordering::Relaxed);
        let run_id = format!(
            "{}-{}-{sequence}",
            opened_at.format("%Y%m%dT%H%M%S%.3f"),
            std::process::id()
        );
        Ok(Self { path, file, run_id })
    }

    /// Appends an informational record.
    pub(crate) fn info(
        &mut self,
        module: &str,
        line: u32,
        message: fmt::Arguments<'_>,
    ) -> Result<()> {
        self.write(LogLevel::Info, module, line, message)
    }

    /// Appends a non-fatal warning record.
    pub(crate) fn warn(
        &mut self,
        module: &str,
        line: u32,
        message: fmt::Arguments<'_>,
    ) -> Result<()> {
        self.write(LogLevel::Warn, module, line, message)
    }

    /// Appends an error record.
    pub(crate) fn error(
        &mut self,
        module: &str,
        line: u32,
        message: fmt::Arguments<'_>,
    ) -> Result<()> {
        self.write(LogLevel::Error, module, line, message)
    }

    /// Flushes and synchronizes all records written so far.
    pub(crate) fn sync(&mut self) -> Result<()> {
        self.file.flush().map_err(|source| Error::Log {
            path: self.path.clone(),
            source,
        })?;
        self.file.sync_all().map_err(|source| Error::Log {
            path: self.path.clone(),
            source,
        })
    }

    fn write(
        &mut self,
        level: LogLevel,
        module: &str,
        line: u32,
        message: fmt::Arguments<'_>,
    ) -> Result<()> {
        let timestamp = Local::now().format("%Y-%m-%d %H:%M:%S%.3f");
        let message = escape_record_field(&message.to_string());
        let record = format!(
            "{timestamp} | {:<8} | {module}:{line} - run_id={} {message}\n",
            level.label(),
            self.run_id
        );
        self.file
            .write_all(record.as_bytes())
            .map_err(|source| Error::Log {
                path: self.path.clone(),
                source,
            })
    }
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

#[cfg(test)]
mod tests {
    use std::fs;

    use tempfile::tempdir;

    use super::*;

    #[test]
    fn creates_and_appends_apollo_style_records()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        let directory = tempdir()?;
        {
            let mut logger = Logger::open_in(directory.path(), "trace")?;
            logger.info("signal::test", 12, format_args!("first"))?;
            logger.warn("signal::test", 13, format_args!("second"))?;
            logger.sync()?;
        }
        {
            let mut logger = Logger::open_in(directory.path(), "trace")?;
            logger.error("signal::test", 14, format_args!("third"))?;
            logger.sync()?;
        }

        let text = fs::read_to_string(directory.path().join("trace.log"))?;
        let lines = text.lines().collect::<Vec<_>>();
        assert_eq!(lines.len(), 3);
        assert!(lines[0].contains(" | INFO     | signal::test:12 - run_id="));
        assert!(lines[0].ends_with(" first"));
        assert!(lines[1].contains(" | WARN     | signal::test:13 - run_id="));
        assert!(lines[1].ends_with(" second"));
        assert!(lines[2].contains(" | ERROR    | signal::test:14 - run_id="));
        assert!(lines[2].ends_with(" third"));
        assert_eq!(lines[0].as_bytes().get(4), Some(&b'-'));
        assert_eq!(lines[0].as_bytes().get(23), Some(&b' '));
        Ok(())
    }

    #[test]
    fn escapes_record_delimiters_and_control_characters()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        let directory = tempdir()?;
        let mut logger = Logger::open_in(directory.path(), "trace")?;
        logger.info(
            "signal::test",
            20,
            format_args!("event=test value=left|right\nforged\tfield\r"),
        )?;
        logger.sync()?;

        let text = fs::read_to_string(directory.path().join("trace.log"))?;
        assert_eq!(text.lines().count(), 1);
        assert!(text.contains("value=left\\|right\\nforged\\tfield\\r"));
        Ok(())
    }
}
