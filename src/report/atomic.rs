//! Atomic no-overwrite publication of one completed JSON file.

use std::fs::{self, File, OpenOptions};
use std::io::Write;
use std::path::{Path, PathBuf};

use crate::error::{Error, Result};

struct TemporaryFile {
    path: PathBuf,
    published: bool,
}

impl Drop for TemporaryFile {
    fn drop(&mut self) {
        if !self.published {
            let _ = fs::remove_file(&self.path);
        }
    }
}

/// Writes bytes to a sibling and atomically links them to a new target.
pub(crate) fn publish(path: &Path, bytes: &[u8]) -> Result<()> {
    if path.exists() {
        return Err(Error::Path {
            kind: "output",
            path: path.to_path_buf(),
            reason: "target already exists".into(),
        });
    }
    let parent = path
        .parent()
        .filter(|value| !value.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    fs::create_dir_all(parent).map_err(|source| Error::Output {
        path: parent.to_path_buf(),
        source,
    })?;
    let (temporary_path, mut file) = create_temporary(path)?;
    let mut temporary = TemporaryFile {
        path: temporary_path.clone(),
        published: false,
    };
    file.write_all(bytes).map_err(|source| Error::Output {
        path: temporary_path.clone(),
        source,
    })?;
    file.sync_all().map_err(|source| Error::Output {
        path: temporary_path.clone(),
        source,
    })?;
    drop(file);
    fs::hard_link(&temporary_path, path).map_err(|source| Error::Output {
        path: path.to_path_buf(),
        source,
    })?;
    if let Err(source) = fs::remove_file(&temporary_path) {
        let _ = fs::remove_file(path);
        return Err(Error::Output {
            path: temporary_path,
            source,
        });
    }
    temporary.published = true;
    if let Err(source) = File::open(parent).and_then(|directory| directory.sync_all()) {
        let _ = fs::remove_file(path);
        return Err(Error::Output {
            path: parent.to_path_buf(),
            source,
        });
    }
    Ok(())
}

fn create_temporary(path: &Path) -> Result<(PathBuf, File)> {
    for attempt in 0_u16..1024 {
        let mut name = path.as_os_str().to_os_string();
        name.push(format!(".tmp.{}.{attempt}", std::process::id()));
        let temporary_path = PathBuf::from(name);
        match OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&temporary_path)
        {
            Ok(file) => return Ok((temporary_path, file)),
            Err(source) if source.kind() == std::io::ErrorKind::AlreadyExists => continue,
            Err(source) => {
                return Err(Error::Output {
                    path: temporary_path,
                    source,
                });
            }
        }
    }
    Err(Error::Path {
        kind: "temporary output",
        path: path.to_path_buf(),
        reason: "could not reserve a sibling temporary file after 1024 attempts".into(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn skips_a_stale_temporary_name() -> std::result::Result<(), Box<dyn std::error::Error>> {
        let directory = tempfile::tempdir()?;
        let target = directory.path().join("result.json");
        let stale = PathBuf::from(format!("{}.tmp.{}.0", target.display(), std::process::id()));
        fs::write(&stale, b"stale")?;

        publish(&target, b"complete")?;

        assert_eq!(fs::read(&target)?, b"complete");
        assert_eq!(fs::read(&stale)?, b"stale");
        Ok(())
    }
}
