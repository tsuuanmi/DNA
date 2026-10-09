//! Target profiles: knowledge about one sequencing target, kept out of code.
//!
//! A profile names the reference a target is validated against, its topology,
//! the reportable regions, and the optional representation chain (indel
//! placement, nomenclature windows, notation style). Method parameters stay in
//! the scientific configuration that references the profile. Profiles are
//! strict versioned TOML, validated completely at load, and identified in
//! result provenance by their declared id and file SHA-256.

mod raw;
mod window;

use std::fs::File;
use std::io::Read;
use std::path::Path;

use serde::Deserialize;

use crate::checksum::hex_sha256;
use crate::error::{Error, ProfileError, Result};
use crate::model::reference::{Reference, ReferenceTopology};

pub use window::{Anchor, NomenclatureWindow, WindowRule};

/// Largest accepted profile file.
const MAX_PROFILE_BYTES: usize = 1024 * 1024;

/// A validated target profile.
///
/// Load one with [`Profile::load`]; pass it to profile-driven capabilities such
/// as `variant_nomenclature::apply`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Profile {
    pub(crate) identity: ProfileIdentity,
    /// SHA-256 of the normalized reference sequence the profile is validated against.
    pub(crate) reference_sha256: Option<String>,
    /// Topology the reference is loaded with.
    pub topology: ReferenceTopology,
    /// Inclusive 1-based reportable regions.
    pub regions: Vec<[usize; 2]>,
    /// Non-overlapping nomenclature windows in reference order.
    pub windows: Vec<NomenclatureWindow>,
    /// The sample notation chain, when the profile declares one.
    pub notation: Option<Notation>,
    /// Notation conventions checked on represented calls, when declared.
    pub conformance: Option<Conformance>,
}

/// Notation conventions that represented calls are checked against; findings
/// are reported, never applied.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Conformance {
    /// Rules in reporting order.
    pub rules: Vec<ConformanceRule>,
    /// Shortest run of one reference base that counts as a homopolymer.
    pub minimum_run_length: usize,
}

/// One notation convention.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ConformanceRule {
    /// An insertion that lengthens a homopolymer is placed after its 3' base.
    InsertionAtRunEnd,
    /// A base inserted inside a homopolymer is the homopolymer's base; any
    /// other difference there is written as a substitution.
    InsertionMatchesRun,
}

impl ConformanceRule {
    /// Stable machine-readable label recorded in results.
    #[must_use]
    pub fn label(self) -> &'static str {
        match self {
            Self::InsertionAtRunEnd => "insertion_at_run_end",
            Self::InsertionMatchesRun => "insertion_matches_run",
        }
    }
}

/// How each read's calls are represented and rendered for sample notation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Notation {
    /// Indel placement applied before the nomenclature windows.
    pub indel_placement: IndelPlacement,
    /// How represented variants are rendered.
    pub style: NotationStyle,
}

/// Identity of the target profile used to produce a result.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProfileIdentity {
    /// Identifier declared by the profile, for example `"human-mtdna-rcrs"`.
    pub id: String,
    /// SHA-256 of the profile file bytes.
    pub sha256: String,
}

/// Where sequence-equivalent indels are placed before nomenclature.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum IndelPlacement {
    /// 3'/right-most placement, never across the FASTA coordinate seam.
    Right,
}

/// How represented variants are serialized for reviewers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NotationStyle {
    /// One call per changed base: `73G`, `249DEL`, `309.1C`.
    PerBaseDecimal,
}

impl NotationStyle {
    /// Stable machine-readable label recorded in results.
    #[must_use]
    pub fn label(self) -> &'static str {
        match self {
            Self::PerBaseDecimal => "per_base_decimal",
        }
    }
}

impl Profile {
    /// Loads and validates one target profile file.
    ///
    /// # Errors
    ///
    /// Returns [`Error`] when the file cannot be read, is too large, is not
    /// UTF-8 TOML, or violates the profile contract.
    pub fn load(path: &Path) -> Result<Self> {
        let read_error = |source| Error::Read {
            kind: "profile",
            path: path.to_path_buf(),
            source,
        };
        // Reading at most one byte past the cap bounds memory even for
        // endless sources such as FIFOs or character devices.
        let mut bytes = Vec::new();
        File::open(path)
            .map_err(read_error)?
            .take(MAX_PROFILE_BYTES as u64 + 1)
            .read_to_end(&mut bytes)
            .map_err(read_error)?;
        if bytes.len() > MAX_PROFILE_BYTES {
            return Err(ProfileError::TooLarge {
                maximum: MAX_PROFILE_BYTES,
            }
            .into());
        }
        let text = std::str::from_utf8(&bytes).map_err(ProfileError::NotUtf8)?;
        let raw = raw::parse(text).map_err(|source| Error::ProfileParse {
            path: path.to_path_buf(),
            source: Box::new(source),
        })?;
        raw.validate(hex_sha256(&bytes))
    }

    /// Identity recorded in result provenance.
    #[must_use]
    pub fn identity(&self) -> &ProfileIdentity {
        &self.identity
    }

    /// Fails unless `reference` is the sequence the profile is validated
    /// against and carries every nomenclature window sequence at its position.
    ///
    /// # Errors
    ///
    /// Returns `ProfileError` when the reference digest differs from the pinned
    /// one or a window sequence is absent.
    pub fn require_reference(&self, reference: &Reference) -> Result<()> {
        if self
            .reference_sha256
            .as_ref()
            .is_some_and(|expected| *expected != reference.sequence_sha256)
        {
            return Err(ProfileError::ReferenceMismatch {
                profile: self.identity.id.clone(),
            }
            .into());
        }
        if let Some(window) = self.windows.iter().find(|window| {
            reference.sequence.get(window.start..window.end()) != Some(window.sequence.as_str())
        }) {
            return Err(ProfileError::WindowNotInReference {
                window: window.name.clone(),
            }
            .into());
        }
        Ok(())
    }
}

/// Shipped profiles for tests of the crates that consume profiles.
#[cfg(any(test, feature = "test-support"))]
pub mod fixtures {
    use std::path::Path;

    use super::Profile;
    use crate::error::Result;

    /// The shipped human-mtDNA profile, from the workspace configuration.
    ///
    /// # Errors
    ///
    /// Returns an error when the shipped profile cannot be loaded.
    pub fn human_mtdna() -> Result<Profile> {
        Profile::load(
            &Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("../../config/profiles/human-mtdna-rcrs.toml"),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::fixtures::human_mtdna;
    use super::*;

    #[test]
    fn loads_the_shipped_human_mtdna_profile() -> Result<()> {
        let profile = human_mtdna()?;
        assert_eq!(profile.identity().id, "human-mtdna-rcrs");
        assert_eq!(profile.identity().sha256.len(), 64);
        assert_eq!(profile.topology, ReferenceTopology::Circular);
        assert_eq!(profile.regions, vec![[16024, 16365], [73, 340], [438, 576]]);
        assert_eq!(
            profile.notation,
            Some(Notation {
                indel_placement: IndelPlacement::Right,
                style: NotationStyle::PerBaseDecimal,
            })
        );
        let windows: Vec<_> = profile
            .windows
            .iter()
            .map(|window| (window.name.as_str(), window.start, window.rules.len()))
            .collect();
        assert_eq!(
            windows,
            [("HVS-II", 302, 4), ("HVS-III", 512, 1), ("HVS-I", 16180, 2)]
        );
        Ok(())
    }

    #[test]
    fn rejects_a_reference_other_than_the_declared_sequence() -> Result<()> {
        let profile = human_mtdna()?;
        let reference = Reference {
            name: "other".into(),
            sequence: "ACGT".into(),
            topology: ReferenceTopology::Circular,
            sequence_sha256: hex_sha256(b"ACGT"),
        };
        assert!(matches!(
            profile.require_reference(&reference),
            Err(Error::Profile(ProfileError::ReferenceMismatch { profile }))
                if profile == "human-mtdna-rcrs"
        ));
        Ok(())
    }

    #[test]
    fn rejects_a_reference_without_a_window_sequence() -> Result<()> {
        let mut profile = human_mtdna()?;
        profile.reference_sha256 = None;
        let reference = |sequence: String| Reference {
            name: "other".into(),
            sequence_sha256: hex_sha256(sequence.as_bytes()),
            sequence,
            topology: ReferenceTopology::Circular,
        };
        // Too short to hold HVS-III or HVS-I.
        let short = reference(format!("{}CCCCCCCTCCCCC{}", "A".repeat(302), "G".repeat(5)));
        assert!(matches!(
            profile.require_reference(&short),
            Err(Error::Profile(ProfileError::WindowNotInReference { window }))
                if window == "HVS-III"
        ));
        // Long enough, but without the HVS-II sequence.
        let mismatched = reference("A".repeat(16_569));
        assert!(matches!(
            profile.require_reference(&mismatched),
            Err(Error::Profile(ProfileError::WindowNotInReference { window }))
                if window == "HVS-II"
        ));
        Ok(())
    }

    #[test]
    fn rejects_oversized_and_unreadable_profiles()
    -> std::result::Result<(), Box<dyn std::error::Error>> {
        let directory = tempfile::tempdir()?;
        let large = directory.path().join("large.toml");
        std::fs::write(&large, vec![b'#'; MAX_PROFILE_BYTES + 1])?;
        assert!(matches!(
            Profile::load(&large),
            Err(Error::Profile(ProfileError::TooLarge { .. }))
        ));
        assert!(matches!(
            Profile::load(&directory.path().join("missing.toml")),
            Err(Error::Read {
                kind: "profile",
                ..
            })
        ));
        let binary = directory.path().join("binary.toml");
        std::fs::write(&binary, [0xff, 0xfe])?;
        assert!(matches!(
            Profile::load(&binary),
            Err(Error::Profile(ProfileError::NotUtf8(_)))
        ));
        let invalid = directory.path().join("invalid.toml");
        std::fs::write(&invalid, "schema_version = [")?;
        assert!(matches!(
            Profile::load(&invalid),
            Err(Error::ProfileParse { .. })
        ));
        Ok(())
    }
}
