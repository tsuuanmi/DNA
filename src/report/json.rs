//! What every document shares: the provenance projections of the workflow's
//! plugins and the profile identity, and deterministic JSON serialization.

use crate::error::{KernelError, Result};
use crate::model::result::{PluginResult, ProfileResult};
use dna_kernel::plugin::PluginDescriptor;
use dna_kernel::profile::ProfileIdentity;

/// Projects the workflow's plugin identities recorded in result provenance.
pub(crate) fn project_plugins(plugins: &[&PluginDescriptor]) -> Vec<PluginResult> {
    plugins
        .iter()
        .map(|plugin| PluginResult {
            id: plugin.id,
            family: plugin.family.label(),
            version: plugin.version,
        })
        .collect()
}

/// Projects the profile identity recorded in result provenance.
pub(crate) fn project_profile(profile: ProfileIdentity) -> ProfileResult {
    ProfileResult {
        id: profile.id,
        sha256: profile.sha256,
    }
}

/// Serializes any typed result with a trailing newline for stable text files.
pub(crate) fn serialize<T: serde::Serialize>(result: &T) -> Result<Vec<u8>> {
    let mut bytes = serde_json::to_vec_pretty(result)
        .map_err(|error| KernelError::Serialize(Box::new(error)))?;
    bytes.push(b'\n');
    Ok(bytes)
}
