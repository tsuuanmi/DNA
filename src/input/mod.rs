//! Input adapters normalize external sequencing sources into validated DNA models.

use std::path::Path;

use crate::config::{self, Config};
use crate::error::{Error, Result};
use crate::model::reference::Reference;
use crate::profile::Profile;
use crate::reference;

pub(crate) mod sanger;
pub(crate) mod sequence;

fn load_config(path: &Path) -> Result<Config> {
    let config = config::load_path(path)?;
    require_regular_file(&config.source_path, "configuration")?;
    Ok(config)
}

fn load_profile(config: &Config) -> Result<Profile> {
    require_regular_file(&config.profile_path, "profile")?;
    Profile::load(&config.profile_path)
}

/// Loads the reference with the profile's topology and fails closed when it is
/// not the sequence the profile is validated against.
fn load_reference(path: &Path, profile: &Profile) -> Result<Reference> {
    let reference = reference::load(path, profile.topology)?;
    profile.require_reference(&reference)?;
    Ok(reference)
}

fn require_regular_file(path: &Path, kind: &'static str) -> Result<()> {
    let metadata = path.metadata().map_err(|source| Error::Read {
        kind,
        path: path.to_path_buf(),
        source,
    })?;
    if !metadata.is_file() || metadata.len() == 0 {
        return Err(Error::Path {
            kind,
            path: path.to_path_buf(),
            reason: "path must be a non-empty regular file",
        });
    }
    Ok(())
}
