//! The read evidence cache (`--evidence-cache`): each trace's
//! `dna.read_evidence/v1` document, stored once per trace and read stage, so
//! a run whose read stage is unchanged skips it (PROP-0011's step contracts).
//!
//! A document depends on the trace, the configuration, and the read stage's
//! code alone: the read stage reads no profile and no reference (ADR-0069),
//! and nothing of the core, the consensus, or the review. So an entry is
//! keyed by the trace's content SHA-256 under a generation directory named
//! for the configuration file's SHA-256 (which the document's provenance
//! records and the reader checks), the `sanger` plugin version, and the read
//! stage's code identity: the git object ids of every source on the path
//! from a trace to its document, recorded at build time. A build from a
//! checkout whose read-stage sources differ from their commit has no code
//! identity, and the cache is refused: a plugin version does not move on a
//! work-in-progress commit, and a stale document would fake a measurement.
//!
//! The same directory holds the placement cache: each read evidence
//! document's `dna.placed_read/v1` document, keyed by the evidence
//! document's SHA-256 under a generation named for the configuration, the
//! placement stage's code identity (the sources from an evidence document to
//! its placed-read document), and the reference. A run whose placement stage
//! is unchanged skips the aligner's search; `dna-align` re-checks every
//! document it reads back against the read's evidence and the reference.
//!
//! The same directory also holds the index of the truth's written forms
//! (`forms.tsv`, [`crate::cause::forms`]), one per generation named for the
//! truth database file's SHA-256, the reference sequence's SHA-256, the
//! target profile's repeat windows, and the index's code identity: a run on
//! the same truth, reference, and windows reads it instead of indexing the
//! whole database, and a renamed or moved window is a new generation. A
//! build whose index sources differ from their commit indexes the database
//! and stores nothing.
//!
//! Entries are gzip-compressed and written to a temporary file then renamed,
//! so a concurrent run never reads half an entry. An entry that does not
//! decompress, is not JSON, or names another trace is a miss: it is
//! recomputed and replaced, never an error.

use std::fs;
use std::io::{self, Read as _, Write as _};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};

use flate2::Compression;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};

/// The read stage's code identity this harness was built with:
/// `path=object;` per read-stage source, empty when a source differed from
/// its commit or without git.
const CODE: &str = env!("DNA_EVAL_READ_STAGE");

/// The placement stage's code identity, as [`CODE`] for the read stage.
const PLACEMENT_CODE: &str = env!("DNA_EVAL_PLACEMENT_STAGE");

/// The forms index's code identity, as [`CODE`] for the read stage.
const FORMS_CODE: &str = env!("DNA_EVAL_FORMS_STAGE");

/// An entry's file name suffix.
const SUFFIX: &str = ".json.gz";

/// The forms index's entry in its generation.
const FORMS_ENTRY: &str = "forms.tsv.gz";

/// Largest forms index read back, decompressed.
const MAX_FORMS_BYTES: u64 = 256 * 1024 * 1024;

/// Why the cache cannot be used.
#[derive(Debug, thiserror::Error)]
pub(crate) enum CacheError {
    /// The harness was built from read-stage sources that differ from their
    /// commit, or without git.
    #[error(
        "--evidence-cache needs a build from committed read-stage and placement-stage sources; this build has no code identity for one of them"
    )]
    NoCodeIdentity,
    /// The cache directory cannot be created or read.
    #[error("evidence cache {path}: {source}")]
    Io {
        /// The directory or entry.
        path: PathBuf,
        /// The failure.
        source: io::Error,
    },
}

/// The parts of a generation's key, as `run.json` records them.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub(crate) struct Key {
    /// The `sanger` plugin version.
    pub(crate) sanger: u32,
    /// The configuration file's SHA-256.
    pub(crate) configuration_sha256: String,
    /// The SHA-256 of the read stage's code identity.
    pub(crate) code_sha256: String,
}

impl Key {
    /// The generation directory's name: the SHA-256 of the key's parts.
    fn generation(&self) -> String {
        let mut hasher = Sha256::new();
        hasher.update(format!(
            "sanger={}\nconfiguration={}\ncode={}\n",
            self.sanger, self.configuration_sha256, self.code_sha256
        ));
        format!("{:x}", hasher.finalize())
    }
}

/// A run's use of the cache, as `run.json` records it.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub(crate) struct Usage {
    /// The cache directory, absolute.
    pub(crate) directory: String,
    /// The generation's key.
    pub(crate) key: Key,
    /// Traces read from the cache.
    pub(crate) hits: usize,
    /// Traces processed and stored: absent or unreadable entries.
    pub(crate) misses: usize,
    /// Of the misses, entries that were present but unreadable.
    pub(crate) replaced: usize,
    /// Bytes of the generation's entries after the run.
    pub(crate) bytes: u64,
}

/// What the index of the truth's written forms depends on, as `run.json`
/// records it.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub(crate) struct FormsKey {
    /// The truth database file's SHA-256.
    pub(crate) truth: String,
    /// The reference sequence's SHA-256.
    pub(crate) reference: String,
    /// The SHA-256 of the target profile's repeat windows: one
    /// `<name>\t<first>\t<last>` line per window, in the profile's order.
    pub(crate) windows: String,
    /// The SHA-256 of the index's code identity.
    pub(crate) code: String,
}

impl FormsKey {
    /// The key of the index of the truth database `truth_sha256` on the
    /// reference sequence `reference_sha256` over `windows` (name, first,
    /// last), or `None` when this build has no code identity for the index.
    pub(crate) fn new<'a>(
        truth_sha256: &str,
        reference_sha256: &str,
        windows: impl IntoIterator<Item = (&'a str, i64, i64)>,
    ) -> Option<Self> {
        Self::with_code(truth_sha256, reference_sha256, windows, FORMS_CODE)
    }

    fn with_code<'a>(
        truth_sha256: &str,
        reference_sha256: &str,
        windows: impl IntoIterator<Item = (&'a str, i64, i64)>,
        code: &str,
    ) -> Option<Self> {
        if code.is_empty() {
            return None;
        }
        let mut hasher = Sha256::new();
        for (name, first, last) in windows {
            hasher.update(format!("{name}\t{first}\t{last}\n"));
        }
        Some(Self {
            truth: truth_sha256.to_owned(),
            reference: reference_sha256.to_owned(),
            windows: format!("{:x}", hasher.finalize()),
            code: format!("{:x}", Sha256::digest(code.as_bytes())),
        })
    }

    /// The generation directory's name: the SHA-256 of the key's parts.
    fn generation(&self) -> String {
        format!(
            "{:x}",
            Sha256::digest(format!(
                "forms\ntruth={}\nreference={}\nwindows={}\ncode={}\n",
                self.truth, self.reference, self.windows, self.code
            ))
        )
    }
}

/// A run's use of the forms index cache, as `run.json` records it.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub(crate) struct FormsUsage {
    /// The cache directory, absolute.
    pub(crate) directory: String,
    /// The generation's key.
    pub(crate) key: FormsKey,
    /// Whether the index was read from the cache rather than built.
    pub(crate) hit: bool,
}

/// The generation of a cache directory that holds one forms index.
pub(crate) struct FormsCache {
    root: PathBuf,
    generation: PathBuf,
    key: FormsKey,
}

impl FormsCache {
    /// Opens the generation of `root` for `key`.
    ///
    /// # Errors
    ///
    /// Returns [`CacheError::Io`] when the directory cannot be created.
    pub(crate) fn open(root: &Path, key: FormsKey) -> Result<Self, CacheError> {
        let io = |path: &Path| {
            let path = path.to_path_buf();
            move |source| CacheError::Io { path, source }
        };
        let root = std::path::absolute(root).map_err(io(root))?;
        let generation = root.join(key.generation());
        fs::create_dir_all(&generation).map_err(io(&generation))?;
        Ok(Self {
            root,
            generation,
            key,
        })
    }

    /// The stored index's text, or `None` when it is absent, does not
    /// decompress, is not UTF-8, or is larger than its bound: a miss.
    pub(crate) fn read(&self) -> Option<String> {
        let file = fs::File::open(self.generation.join(FORMS_ENTRY)).ok()?;
        let mut text = Vec::new();
        GzDecoder::new(file)
            .take(MAX_FORMS_BYTES + 1)
            .read_to_end(&mut text)
            .ok()?;
        if text.len() as u64 > MAX_FORMS_BYTES {
            return None;
        }
        String::from_utf8(text).ok()
    }

    /// Stores the index's text, replacing an entry that was a miss. A
    /// failure to store is not an error: the next run builds it again.
    pub(crate) fn store(&self, table: &str) {
        let temporary = self
            .generation
            .join(format!(".{FORMS_ENTRY}.{}.tmp", std::process::id()));
        if write_entry(&temporary, table.as_bytes()).is_err()
            || fs::rename(&temporary, self.generation.join(FORMS_ENTRY)).is_err()
        {
            let _ = fs::remove_file(&temporary);
        }
    }

    /// The run's use of the cache.
    pub(crate) fn usage(&self, hit: bool) -> FormsUsage {
        FormsUsage {
            directory: self.root.display().to_string(),
            key: self.key.clone(),
            hit,
        }
    }
}

/// A run's caches: read evidence documents, and the placed-read documents
/// of those documents.
pub(crate) struct Caches {
    pub(crate) evidence: EvidenceCache,
    pub(crate) placed: EvidenceCache,
}

/// One generation of the cache, filled and read by a run.
pub(crate) struct EvidenceCache {
    root: PathBuf,
    generation: PathBuf,
    key: Key,
    hits: AtomicUsize,
    misses: AtomicUsize,
    replaced: AtomicUsize,
    sequence: AtomicU64,
}

impl EvidenceCache {
    /// Opens the generation of `root` for the configuration at `config`
    /// under this build's read stage.
    ///
    /// # Errors
    ///
    /// Returns [`CacheError`] without a code identity, or when the
    /// directory cannot be created or the configuration read.
    pub(crate) fn open(root: &Path, config: &Path) -> Result<Self, CacheError> {
        Self::with_code(root, config, CODE)
    }

    /// Opens the placement cache's generation of `root` for the
    /// configuration at `config` and the reference whose sequence SHA-256 is
    /// `reference_sha256`, under this build's placement stage.
    ///
    /// # Errors
    ///
    /// As [`open`](Self::open), without a placement-stage code identity.
    pub(crate) fn open_placed(
        root: &Path,
        config: &Path,
        reference_sha256: &str,
    ) -> Result<Self, CacheError> {
        if PLACEMENT_CODE.is_empty() {
            return Err(CacheError::NoCodeIdentity);
        }
        Self::with_code(
            root,
            config,
            &format!("placement;{PLACEMENT_CODE}reference={reference_sha256};"),
        )
    }

    fn with_code(root: &Path, config: &Path, code: &str) -> Result<Self, CacheError> {
        if code.is_empty() {
            return Err(CacheError::NoCodeIdentity);
        }
        let io = |path: &Path| {
            let path = path.to_path_buf();
            move |source| CacheError::Io { path, source }
        };
        let configuration = fs::read(config).map_err(io(config))?;
        let key = Key {
            sanger: dna_sanger::read_processing::PLUGIN.version,
            configuration_sha256: format!("{:x}", Sha256::digest(&configuration)),
            code_sha256: format!("{:x}", Sha256::digest(code.as_bytes())),
        };
        let root = std::path::absolute(root).map_err(io(root))?;
        let generation = root.join(key.generation());
        fs::create_dir_all(&generation).map_err(io(&generation))?;
        Ok(Self {
            root,
            generation,
            key,
            hits: AtomicUsize::new(0),
            misses: AtomicUsize::new(0),
            replaced: AtomicUsize::new(0),
            sequence: AtomicU64::new(0),
        })
    }

    /// The read evidence document of the trace at `trace`, whose content
    /// SHA-256 is `sha256`: read from the cache, or made by `make` and
    /// stored.
    ///
    /// # Errors
    ///
    /// Returns `make`'s failure. A failure to store the document is not an
    /// error: the next run makes it again.
    pub(crate) fn document<E>(
        &self,
        sha256: &str,
        make: impl FnOnce() -> Result<Vec<u8>, E>,
    ) -> Result<Vec<u8>, E> {
        let path = self.generation.join(format!("{sha256}{SUFFIX}"));
        match read_entry(&path, sha256) {
            Entry::Hit(document) => {
                self.hits.fetch_add(1, Ordering::Relaxed);
                return Ok(document);
            }
            Entry::Unreadable => {
                self.replaced.fetch_add(1, Ordering::Relaxed);
            }
            Entry::Absent => {}
        }
        self.misses.fetch_add(1, Ordering::Relaxed);
        let document = make()?;
        let temporary = self.generation.join(format!(
            ".{sha256}.{}.{}.tmp",
            std::process::id(),
            self.sequence.fetch_add(1, Ordering::Relaxed)
        ));
        if write_entry(&temporary, &document).is_ok() {
            if fs::rename(&temporary, &path).is_err() {
                let _ = fs::remove_file(&temporary);
            }
        } else {
            let _ = fs::remove_file(&temporary);
        }
        Ok(document)
    }

    /// The run's use of the cache so far.
    pub(crate) fn usage(&self) -> Usage {
        Usage {
            directory: self.root.display().to_string(),
            key: self.key.clone(),
            hits: self.hits.load(Ordering::Relaxed),
            misses: self.misses.load(Ordering::Relaxed),
            replaced: self.replaced.load(Ordering::Relaxed),
            bytes: entries(&self.generation)
                .iter()
                .map(|(_, bytes, _)| bytes)
                .sum(),
        }
    }
}

/// What a cache file holds.
enum Entry {
    Absent,
    Unreadable,
    Hit(Vec<u8>),
}

/// Reads the entry at `path`: a gzip stream of a JSON document whose
/// provenance names the trace `sha256`.
fn read_entry(path: &Path, sha256: &str) -> Entry {
    let Ok(file) = fs::File::open(path) else {
        return Entry::Absent;
    };
    let mut document = Vec::new();
    if GzDecoder::new(file).read_to_end(&mut document).is_err() {
        return Entry::Unreadable;
    }
    let names = serde_json::from_slice::<serde_json::Value>(&document)
        .is_ok_and(|value| value["provenance"]["input"]["sha256"] == sha256);
    if names {
        Entry::Hit(document)
    } else {
        Entry::Unreadable
    }
}

fn write_entry(path: &Path, document: &[u8]) -> io::Result<()> {
    let mut encoder = GzEncoder::new(fs::File::create(path)?, Compression::default());
    encoder.write_all(document)?;
    encoder.finish()?.sync_all()
}

/// Each entry of a generation: its path, size, and modification time.
fn entries(generation: &Path) -> Vec<(PathBuf, u64, std::time::SystemTime)> {
    fs::read_dir(generation)
        .into_iter()
        .flatten()
        .filter_map(Result::ok)
        .filter_map(|entry| {
            let metadata = entry.metadata().ok()?;
            Some((entry.path(), metadata.len(), metadata.modified().ok()?))
        })
        .collect()
}

/// Removes every generation of `root` but the `keep` most recently used
/// (by their newest entry), and returns the removed directories with their
/// sizes.
///
/// # Errors
///
/// Returns [`CacheError::Io`] when the root cannot be read or a generation
/// removed.
pub(crate) fn prune(root: &Path, keep: usize) -> Result<Vec<(PathBuf, u64)>, CacheError> {
    let io = |path: &Path| {
        let path = path.to_path_buf();
        move |source| CacheError::Io { path, source }
    };
    let mut generations = fs::read_dir(root)
        .map_err(io(root))?
        .filter_map(Result::ok)
        .filter(|entry| entry.file_type().is_ok_and(|kind| kind.is_dir()))
        .map(|entry| {
            let path = entry.path();
            let files = entries(&path);
            let newest = files.iter().map(|(_, _, modified)| *modified).max();
            let bytes = files.iter().map(|(_, bytes, _)| bytes).sum::<u64>();
            (newest, path, bytes)
        })
        .collect::<Vec<_>>();
    generations.sort_by(|left, right| right.0.cmp(&left.0).then(left.1.cmp(&right.1)));
    let mut removed = Vec::new();
    for (_, path, bytes) in generations.into_iter().skip(keep) {
        fs::remove_dir_all(&path).map_err(io(&path))?;
        removed.push((path, bytes));
    }
    Ok(removed)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn document(sha256: &str) -> Vec<u8> {
        format!("{{\"provenance\":{{\"input\":{{\"sha256\":\"{sha256}\"}}}}}}").into_bytes()
    }

    fn cache(root: &Path, code: &str) -> Result<EvidenceCache, CacheError> {
        let config = root.join("dna.toml");
        if !config.exists() {
            fs::write(&config, "schema_version = 19\n").map_err(|source| CacheError::Io {
                path: config.clone(),
                source,
            })?;
        }
        EvidenceCache::with_code(&root.join("cache"), &config, code)
    }

    /// A miss stores the document, a hit reads it back, and a build without
    /// a code identity is refused.
    #[test]
    fn stores_then_reads_a_document() -> Result<(), Box<dyn std::error::Error>> {
        let root = tempfile::tempdir()?;
        let cache = cache(root.path(), "crates/dna-sanger=abc;")?;
        let made = cache.document("aa", || Ok::<_, CacheError>(document("aa")))?;
        let read = cache.document("aa", || -> Result<Vec<u8>, CacheError> {
            panic!("a stored document is made again")
        })?;
        assert_eq!(made, read);
        let usage = cache.usage();
        assert_eq!((usage.hits, usage.misses, usage.replaced), (1, 1, 0));
        assert!(usage.bytes > 0);
        assert!(matches!(
            self::cache(root.path(), ""),
            Err(CacheError::NoCodeIdentity)
        ));
        Ok(())
    }

    /// An entry stored under another code identity is never read, and an
    /// unreadable or forged entry is replaced, not an error.
    #[test]
    fn ignores_other_generations_and_replaces_bad_entries() -> Result<(), Box<dyn std::error::Error>>
    {
        let root = tempfile::tempdir()?;
        let old = cache(root.path(), "crates/dna-sanger=old;")?;
        old.document("aa", || Ok::<_, CacheError>(document("aa")))?;
        let new = cache(root.path(), "crates/dna-sanger=new;")?;
        let mut made = 0;
        new.document("aa", || {
            made += 1;
            Ok::<_, CacheError>(document("aa"))
        })?;
        assert_eq!(made, 1);
        let path = new.generation.join(format!("bb{SUFFIX}"));
        fs::write(&path, b"not gzip")?;
        new.document("bb", || Ok::<_, CacheError>(document("bb")))?;
        write_entry(&new.generation.join(format!("cc{SUFFIX}")), &document("dd"))?;
        let replaced = new.document("cc", || Ok::<_, CacheError>(document("cc")))?;
        assert_eq!(replaced, document("cc"));
        let usage = new.usage();
        assert_eq!((usage.hits, usage.misses, usage.replaced), (0, 3, 2));
        assert_eq!(
            new.document("bb", || -> Result<Vec<u8>, CacheError> {
                panic!("a replaced entry is made again")
            })?,
            document("bb")
        );
        Ok(())
    }

    /// A stored forms index reads back under its key alone; another
    /// window, reference, or truth is another generation, and a build without
    /// a code identity has no key.
    #[test]
    fn stores_a_forms_index_per_key() -> Result<(), Box<dyn std::error::Error>> {
        let root = tempfile::tempdir()?;
        let key = |truth: &str, reference: &str, window: (&'static str, i64, i64)| {
            FormsKey::with_code(truth, reference, [window], "forms.rs=abc;")
                .ok_or("a code identity gives a key")
        };
        let hv1 = ("hv1_c_tract", 16184, 16193);
        let cache = FormsCache::open(root.path(), key("t", "r", hv1)?)?;
        assert_eq!(cache.read(), None);
        cache.store("table\n");
        assert_eq!(cache.read().as_deref(), Some("table\n"));
        assert!(cache.usage(true).hit);
        for other in [
            key("t", "r", ("ac_repeat", 16184, 16193))?,
            key("t", "r", ("hv1_c_tract", 16184, 16194))?,
            key("t", "other", hv1)?,
            key("other", "r", hv1)?,
        ] {
            assert_eq!(FormsCache::open(root.path(), other)?.read(), None);
        }
        assert_eq!(FormsKey::with_code("t", "r", [hv1], ""), None);
        fs::write(cache.generation.join(FORMS_ENTRY), b"not gzip")?;
        assert_eq!(cache.read(), None);
        Ok(())
    }

    /// Pruning keeps the most recently used generations.
    #[test]
    fn prunes_all_but_the_newest_generations() -> Result<(), Box<dyn std::error::Error>> {
        let root = tempfile::tempdir()?;
        let old = cache(root.path(), "crates/dna-sanger=old;")?;
        old.document("aa", || Ok::<_, CacheError>(document("aa")))?;
        std::thread::sleep(std::time::Duration::from_millis(20));
        let new = cache(root.path(), "crates/dna-sanger=new;")?;
        new.document("aa", || Ok::<_, CacheError>(document("aa")))?;
        let removed = prune(&root.path().join("cache"), 1)?;
        assert_eq!(removed.len(), 1);
        assert_eq!(removed[0].0, old.generation);
        assert!(new.generation.exists());
        Ok(())
    }
}
