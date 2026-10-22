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
//!
//! Each generation labels itself: `key.json` records its kind (`evidence`,
//! `placement`, or `forms`) and key when it is first opened, and
//! `last_used` the time of the latest run that opened it, rewritten by every
//! run. `cache list` reads them; `cache prune` keeps the most recently used
//! generations of each kind, and removes the temporary files of runs no
//! longer alive. A generation written before the labels is `unlabelled`
//! until a run opens it again, and its recency is its newest entry's.

use std::fmt::Write as _;
use std::fs;
use std::io::{self, Read as _, Write as _};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::time::{Duration, SystemTime, UNIX_EPOCH};

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

/// A generation's label: its kind and key.
const KEY_JSON: &str = "key.json";

/// When a run last opened a generation, in seconds since the Unix epoch.
const LAST_USED: &str = "last_used";

/// The kinds of generation a cache directory holds.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub(crate) enum Kind {
    /// Read evidence documents.
    Evidence,
    /// Placed-read documents.
    Placement,
    /// The index of the truth's written forms.
    Forms,
}

impl Kind {
    fn label(kind: Option<Self>) -> &'static str {
        match kind {
            Some(Self::Evidence) => "evidence",
            Some(Self::Placement) => "placement",
            Some(Self::Forms) => "forms",
            None => "unlabelled",
        }
    }
}

/// `key.json`: a generation's kind and key.
#[derive(Debug, Serialize, Deserialize)]
struct Label {
    kind: Kind,
    key: serde_json::Value,
}

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
        let (root, generation) = open_generation(root, Kind::Forms, &key.generation(), &key)?;
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
            .join(format!(".{FORMS_ENTRY}.{}.0.tmp", std::process::id()));
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
        Self::with_code(root, config, (Kind::Evidence, CODE))
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
            (
                Kind::Placement,
                &format!("placement;{PLACEMENT_CODE}reference={reference_sha256};"),
            ),
        )
    }

    fn with_code(
        root: &Path,
        config: &Path,
        (kind, code): (Kind, &str),
    ) -> Result<Self, CacheError> {
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
        let (root, generation) = open_generation(root, kind, &key.generation(), &key)?;
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

/// Creates the generation `name` of `root` and labels it: `key.json` with
/// `kind` and `key` when it has none, and `last_used` now. Returns the
/// absolute root and the generation. A label that cannot be written is not
/// an error: the generation then lists as unlabelled or by its entries.
fn open_generation(
    root: &Path,
    kind: Kind,
    name: &str,
    key: &impl Serialize,
) -> Result<(PathBuf, PathBuf), CacheError> {
    let io = |path: &Path| {
        let path = path.to_path_buf();
        move |source| CacheError::Io { path, source }
    };
    let root = std::path::absolute(root).map_err(io(root))?;
    let generation = root.join(name);
    fs::create_dir_all(&generation).map_err(io(&generation))?;
    if !generation.join(KEY_JSON).exists()
        && let Ok(key) = serde_json::to_value(key)
        && let Ok(text) = serde_json::to_string_pretty(&Label { kind, key })
    {
        replace(&generation, KEY_JSON, format!("{text}\n").as_bytes());
    }
    let now = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_or(0, |since| since.as_secs());
    replace(&generation, LAST_USED, format!("{now}\n").as_bytes());
    Ok((root, generation))
}

/// Writes `contents` as the file `name` of `generation` through a temporary
/// file, so a reader sees the old file or the new one; a failure leaves the
/// old one.
fn replace(generation: &Path, name: &str, contents: &[u8]) {
    let temporary = generation.join(format!(".{name}.{}.0.tmp", std::process::id()));
    if fs::write(&temporary, contents).is_err()
        || fs::rename(&temporary, generation.join(name)).is_err()
    {
        let _ = fs::remove_file(&temporary);
    }
}

/// Each entry of a generation, without its label and temporary files: its
/// path, size, and modification time.
fn entries(generation: &Path) -> Vec<(PathBuf, u64, SystemTime)> {
    fs::read_dir(generation)
        .into_iter()
        .flatten()
        .filter_map(Result::ok)
        .filter(|entry| {
            let name = entry.file_name();
            let name = name.to_string_lossy();
            !name.starts_with('.') && name != KEY_JSON && name != LAST_USED
        })
        .filter_map(|entry| {
            let metadata = entry.metadata().ok()?;
            Some((entry.path(), metadata.len(), metadata.modified().ok()?))
        })
        .collect()
}

/// One generation of a cache directory, as `cache list` shows it.
#[derive(Debug)]
pub(crate) struct Generation {
    /// The generation directory.
    pub(crate) path: PathBuf,
    /// Its kind, `None` when it has no readable `key.json`.
    pub(crate) kind: Option<Kind>,
    /// Its key, as `key.json` records it.
    pub(crate) key: Option<serde_json::Value>,
    /// When a run last opened it (`last_used`), else its newest entry's
    /// modification time; `None` when it has neither.
    pub(crate) last_used: Option<SystemTime>,
    /// Its entries.
    pub(crate) entries: usize,
    /// Their bytes.
    pub(crate) bytes: u64,
}

/// A cache directory's generations, and the temporary files left by runs
/// that are no longer alive.
#[derive(Debug)]
pub(crate) struct Listing {
    pub(crate) generations: Vec<Generation>,
    /// Each orphaned temporary file and its size.
    pub(crate) orphans: Vec<(PathBuf, u64)>,
}

/// Lists the generations of `root`, by kind and then most recently used
/// first.
///
/// # Errors
///
/// Returns [`CacheError::Io`] when the root cannot be read.
pub(crate) fn list(root: &Path) -> Result<Listing, CacheError> {
    let io = |path: &Path| {
        let path = path.to_path_buf();
        move |source| CacheError::Io { path, source }
    };
    let mut generations = Vec::new();
    let mut orphans = Vec::new();
    for entry in fs::read_dir(root).map_err(io(root))? {
        let entry = entry.map_err(io(root))?;
        if !entry.file_type().is_ok_and(|kind| kind.is_dir()) {
            continue;
        }
        let path = entry.path();
        let files = entries(&path);
        let label = fs::read(path.join(KEY_JSON))
            .ok()
            .and_then(|text| serde_json::from_slice::<Label>(&text).ok());
        let last_used = fs::read_to_string(path.join(LAST_USED))
            .ok()
            .and_then(|text| text.trim().parse::<u64>().ok())
            .map(|seconds| UNIX_EPOCH + Duration::from_secs(seconds))
            .or_else(|| files.iter().map(|(_, _, modified)| *modified).max());
        orphans.extend(orphaned(&path));
        generations.push(Generation {
            kind: label.as_ref().map(|label| label.kind),
            key: label.map(|label| label.key),
            last_used,
            entries: files.len(),
            bytes: files.iter().map(|(_, bytes, _)| bytes).sum(),
            path,
        });
    }
    generations.sort_by(|left, right| {
        (left.kind.is_none(), left.kind)
            .cmp(&(right.kind.is_none(), right.kind))
            .then(right.last_used.cmp(&left.last_used))
            .then(left.path.cmp(&right.path))
    });
    orphans.sort();
    Ok(Listing {
        generations,
        orphans,
    })
}

/// The temporary files of `generation` whose writing process is no longer
/// alive on this host: a temporary file is named `.<name>.<pid>.<n>.tmp`.
/// Without `/proc` no process can be checked, and none is orphaned.
fn orphaned(generation: &Path) -> Vec<(PathBuf, u64)> {
    let processes = Path::new("/proc");
    if !processes.is_dir() {
        return Vec::new();
    }
    fs::read_dir(generation)
        .into_iter()
        .flatten()
        .filter_map(Result::ok)
        .filter_map(|entry| {
            let name = entry.file_name().into_string().ok()?;
            let stem = name.strip_prefix('.')?.strip_suffix(".tmp")?;
            let mut parts = stem.rsplit('.');
            parts.next()?;
            let pid = parts.next()?.parse::<u32>().ok()?;
            let alive = processes.join(pid.to_string()).exists();
            (!alive).then(|| (entry.path(), entry.metadata().map_or(0, |data| data.len())))
        })
        .collect()
}

/// What `prune` removed.
#[derive(Debug)]
pub(crate) struct Pruned {
    /// Each removed generation, with its kind and bytes.
    pub(crate) generations: Vec<(PathBuf, Option<Kind>, u64)>,
    /// Each removed temporary file and its size.
    pub(crate) orphans: Vec<(PathBuf, u64)>,
}

/// Removes every generation of `root` but the `keep` most recently used of
/// each kind (unlabelled generations being a kind of their own), and the
/// temporary files of runs no longer alive; with `dry_run`, removes nothing
/// and returns what it would remove.
///
/// # Errors
///
/// Returns [`CacheError::Io`] when the root cannot be read or a generation
/// or file removed.
pub(crate) fn prune(root: &Path, keep: usize, dry_run: bool) -> Result<Pruned, CacheError> {
    let io = |path: &Path| {
        let path = path.to_path_buf();
        move |source| CacheError::Io { path, source }
    };
    let listing = list(root)?;
    let mut kept = std::collections::BTreeMap::<Option<Kind>, usize>::new();
    let mut pruned = Pruned {
        generations: Vec::new(),
        orphans: Vec::new(),
    };
    for generation in listing.generations {
        let count = kept.entry(generation.kind).or_default();
        if *count < keep {
            *count += 1;
            continue;
        }
        if !dry_run {
            fs::remove_dir_all(&generation.path).map_err(io(&generation.path))?;
        }
        pruned
            .generations
            .push((generation.path, generation.kind, generation.bytes));
    }
    for (path, bytes) in listing.orphans {
        if pruned
            .generations
            .iter()
            .any(|(generation, _, _)| path.starts_with(generation))
        {
            continue;
        }
        if !dry_run {
            fs::remove_file(&path).map_err(io(&path))?;
        }
        pruned.orphans.push((path, bytes));
    }
    Ok(pruned)
}

/// `cache list`'s table: one line per generation, then the orphaned
/// temporary files.
pub(crate) fn listing_table(listing: &Listing) -> String {
    let mut text = String::from("generation\tkind\tlast_used\tentries\tbytes\tkey\n");
    for generation in &listing.generations {
        let name = generation
            .path
            .file_name()
            .map_or_else(String::new, |name| name.to_string_lossy().into_owned());
        let _ = writeln!(
            text,
            "{name}\t{}\t{}\t{}\t{}\t{}",
            Kind::label(generation.kind),
            generation.last_used.map_or_else(|| "-".to_owned(), utc),
            generation.entries,
            generation.bytes,
            generation
                .key
                .as_ref()
                .map_or_else(|| "-".to_owned(), ToString::to_string),
        );
    }
    let bytes = listing.orphans.iter().map(|(_, bytes)| bytes).sum::<u64>();
    let _ = writeln!(
        text,
        "# {} orphaned temporary file(s), {bytes} bytes",
        listing.orphans.len()
    );
    text
}

/// `cache prune`'s summary.
pub(crate) fn pruned_summary(pruned: &Pruned, dry_run: bool) -> String {
    let mut kinds = std::collections::BTreeMap::<&str, (usize, u64)>::new();
    for (_, kind, bytes) in &pruned.generations {
        let entry = kinds.entry(Kind::label(*kind)).or_default();
        entry.0 += 1;
        entry.1 += bytes;
    }
    let generations = kinds
        .iter()
        .map(|(kind, (count, bytes))| format!("{kind} {count} ({bytes} bytes)"))
        .collect::<Vec<_>>();
    let orphans = pruned.orphans.iter().map(|(_, bytes)| bytes).sum::<u64>();
    format!(
        "dna-eval: {} {} cache generation(s){}{}; {} orphaned temporary file(s), {orphans} bytes",
        if dry_run { "would remove" } else { "removed" },
        pruned.generations.len(),
        if generations.is_empty() { "" } else { ": " },
        generations.join(", "),
        pruned.orphans.len(),
    )
}

/// A time as `YYYY-MM-DDTHH:MM:SSZ`.
fn utc(time: SystemTime) -> String {
    let seconds = time
        .duration_since(UNIX_EPOCH)
        .map_or(0, |since| since.as_secs());
    let (days, rest) = (seconds / 86_400, seconds % 86_400);
    // Civil date from days since 1970-01-01 (proleptic Gregorian).
    let shifted = days + 719_468;
    let era = shifted / 146_097;
    let day_of_era = shifted % 146_097;
    let year_of_era =
        (day_of_era - day_of_era / 1460 + day_of_era / 36_524 - day_of_era / 146_096) / 365;
    let day_of_year = day_of_era - (365 * year_of_era + year_of_era / 4 - year_of_era / 100);
    let month_index = (5 * day_of_year + 2) / 153;
    let day = day_of_year - (153 * month_index + 2) / 5 + 1;
    let month = if month_index < 10 {
        month_index + 3
    } else {
        month_index - 9
    };
    let year = year_of_era + era * 400 + u64::from(month <= 2);
    format!(
        "{year:04}-{month:02}-{day:02}T{:02}:{:02}:{:02}Z",
        rest / 3600,
        rest / 60 % 60,
        rest % 60
    )
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
        EvidenceCache::with_code(&root.join("cache"), &config, (Kind::Evidence, code))
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

    /// Every generation labels itself with its kind and key, and records
    /// when a run last opened it; a generation without a label lists as
    /// unlabelled, used when its newest entry was written.
    #[test]
    fn lists_generations_by_kind_and_last_use() -> Result<(), Box<dyn std::error::Error>> {
        let root = tempfile::tempdir()?;
        let evidence = cache(root.path(), "crates/dna-sanger=abc;")?;
        evidence.document("aa", || Ok::<_, CacheError>(document("aa")))?;
        let key = FormsKey::with_code("t", "r", [("w", 1, 2)], "forms.rs=abc;")
            .ok_or("a code identity gives a key")?;
        let forms = FormsCache::open(&root.path().join("cache"), key.clone())?;
        forms.store("table\n");
        let legacy = root.path().join("cache").join("legacy");
        fs::create_dir(&legacy)?;
        write_entry(&legacy.join(format!("bb{SUFFIX}")), &document("bb"))?;
        let listing = list(&root.path().join("cache"))?;
        let shown = listing
            .generations
            .iter()
            .map(|generation| (generation.kind, generation.entries))
            .collect::<Vec<_>>();
        assert_eq!(
            shown,
            [(Some(Kind::Evidence), 1), (Some(Kind::Forms), 1), (None, 1)]
        );
        assert_eq!(
            listing.generations[1].key,
            Some(serde_json::to_value(&key)?)
        );
        assert!(listing.generations.iter().all(|generation| {
            generation
                .last_used
                .is_some_and(|time| time > UNIX_EPOCH + Duration::from_secs(1_000_000_000))
        }));
        let table = listing_table(&listing);
        assert!(table.starts_with("generation\tkind\tlast_used\tentries\tbytes\tkey\n"));
        assert!(table.contains("\tunlabelled\t"), "{table}");
        assert!(table.ends_with("# 0 orphaned temporary file(s), 0 bytes\n"));
        Ok(())
    }

    /// Pruning keeps the most recently used generations of each kind, and
    /// removes the temporary files of runs no longer alive; a dry run
    /// removes nothing.
    #[test]
    fn prunes_per_kind_and_sweeps_orphaned_temporary_files()
    -> Result<(), Box<dyn std::error::Error>> {
        let root = tempfile::tempdir()?;
        let directory = root.path().join("cache");
        let old = cache(root.path(), "crates/dna-sanger=old;")?;
        let new = cache(root.path(), "crates/dna-sanger=new;")?;
        let placed = EvidenceCache::with_code(
            &directory,
            &root.path().join("dna.toml"),
            (Kind::Placement, "placement;"),
        )?;
        fs::write(old.generation.join(LAST_USED), "100\n")?;
        fs::write(new.generation.join(LAST_USED), "200\n")?;
        fs::write(placed.generation.join(LAST_USED), "50\n")?;
        // A writer that no longer runs, and this one.
        let orphan = new.generation.join(".aa.4294967295.0.tmp");
        let live = new
            .generation
            .join(format!(".bb.{}.0.tmp", std::process::id()));
        fs::write(&orphan, b"half")?;
        fs::write(&live, b"half")?;
        let planned = prune(&directory, 1, true)?;
        assert_eq!(
            planned
                .generations
                .iter()
                .map(|(path, kind, _)| (path.clone(), *kind))
                .collect::<Vec<_>>(),
            [(old.generation.clone(), Some(Kind::Evidence))]
        );
        assert_eq!(planned.orphans, [(orphan.clone(), 4)]);
        assert!(old.generation.exists() && orphan.exists());
        assert_eq!(
            pruned_summary(&planned, true),
            "dna-eval: would remove 1 cache generation(s): evidence 1 (0 bytes); 1 orphaned temporary file(s), 4 bytes"
        );
        prune(&directory, 1, false)?;
        assert!(!old.generation.exists() && !orphan.exists());
        assert!(new.generation.exists() && placed.generation.exists() && live.exists());
        Ok(())
    }

    #[test]
    fn writes_times_in_utc() {
        assert_eq!(utc(UNIX_EPOCH), "1970-01-01T00:00:00Z");
        assert_eq!(
            utc(UNIX_EPOCH + Duration::from_secs(951_825_599)),
            "2000-02-29T11:59:59Z"
        );
        assert_eq!(
            utc(UNIX_EPOCH + Duration::from_secs(1_792_627_200)),
            "2026-10-22T00:00:00Z"
        );
    }
}
