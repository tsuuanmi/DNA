//! Variants documents (`dna.variants/v1`) as input to post-calling plugins
//! (ADR-0069).
//!
//! Only the fields post-calling needs are read: the sample, the reference
//! identity, and each read's eligible variants. The document's own provenance
//! stays with the document, which the consumer identifies by its SHA-256.

use std::fs;
use std::path::{Path, PathBuf};

use serde::Deserialize;

use crate::checksum::hex_sha256;
use crate::config::{Config, MAX_VARIANTS_BYTES};
use crate::error::{Error, Result, VariantsError};
use crate::model::nucleotide::is_canonical;
use crate::model::reference::Reference;
use crate::profile::Profile;
use crate::variant::{Variant, VariantKind};

use super::{load_config, load_profile, load_reference, require_regular_file};

/// The contract this reader accepts.
const SCHEMA_VERSION: &str = "dna.variants/v1";

/// One read of a variants document with its eligible variants.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct DocumentRead {
    pub(crate) name: String,
    pub(crate) sha256: String,
    /// Eligible variants, in document order.
    pub(crate) variants: Vec<Variant>,
}

/// A validated variants document.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct VariantsDocument {
    /// SHA-256 of the document bytes.
    pub(crate) sha256: String,
    pub(crate) reads: Vec<DocumentRead>,
}

/// Scientific inputs for deriving notation from a variants document.
pub(crate) struct NotationInputs {
    pub(crate) config: Config,
    pub(crate) profile: Profile,
    pub(crate) reference: Reference,
    pub(crate) document: VariantsDocument,
}

#[derive(Deserialize)]
struct RawDocument {
    schema_version: String,
    sample_id: String,
    provenance: RawProvenance,
    reads: Vec<RawRead>,
}

#[derive(Deserialize)]
struct RawProvenance {
    reference: RawReference,
}

#[derive(Deserialize)]
struct RawReference {
    sha256: String,
}

#[derive(Deserialize)]
struct RawRead {
    name: String,
    sha256: String,
    variants: Vec<RawVariant>,
}

#[derive(Deserialize)]
struct RawVariant {
    position: usize,
    reference: String,
    alternate: String,
    kind: String,
    eligible: bool,
}

/// Validates paths, loads configuration, profile, and reference, and reads
/// the variants document of `sample_id` against that reference.
pub(crate) fn load_notation(
    sample_id: &str,
    document_path: &Path,
    reference_path: &Path,
    config_path: &Path,
) -> Result<NotationInputs> {
    require_regular_file(document_path, "variants")?;
    require_regular_file(reference_path, "reference")?;
    let config = load_config(config_path)?;
    let profile = load_profile(&config)?;
    if profile.notation.is_none() {
        return Err(VariantsError::NoNotation.into());
    }
    let reference = load_reference(reference_path, &profile)?;
    let document = load(document_path, sample_id, &reference)?;
    Ok(NotationInputs {
        config,
        profile,
        reference,
        document,
    })
}

fn load(path: &Path, sample_id: &str, reference: &Reference) -> Result<VariantsDocument> {
    let read_error = |source| Error::Read {
        kind: "variants",
        path: path.to_path_buf(),
        source,
    };
    let length = fs::metadata(path).map_err(read_error)?.len();
    if length == 0 || length > MAX_VARIANTS_BYTES as u64 {
        return Err(VariantsError::FileSize {
            bytes: length,
            maximum: MAX_VARIANTS_BYTES,
        }
        .into());
    }
    let bytes = fs::read(path).map_err(read_error)?;
    parse(&bytes, path.to_path_buf(), sample_id, reference)
}

fn parse(
    bytes: &[u8],
    path: PathBuf,
    sample_id: &str,
    reference: &Reference,
) -> Result<VariantsDocument> {
    let raw: RawDocument = serde_json::from_slice(bytes).map_err(|error| Error::VariantsParse {
        path,
        source: Box::new(error),
    })?;
    if raw.schema_version != SCHEMA_VERSION {
        return Err(VariantsError::UnsupportedSchema {
            found: raw.schema_version,
        }
        .into());
    }
    if raw.sample_id != sample_id {
        return Err(VariantsError::SampleMismatch {
            expected: sample_id.to_owned(),
            found: raw.sample_id,
        }
        .into());
    }
    if raw.provenance.reference.sha256 != reference.sequence_sha256 {
        return Err(VariantsError::ReferenceMismatch.into());
    }
    let reads = raw
        .reads
        .into_iter()
        .map(|read| {
            let variants = read
                .variants
                .into_iter()
                .filter(|variant| variant.eligible)
                .map(|variant| to_variant(&read.name, &reference.name, variant))
                .collect::<Result<_>>()?;
            Ok(DocumentRead {
                name: read.name,
                sha256: read.sha256,
                variants,
            })
        })
        .collect::<Result<_>>()?;
    Ok(VariantsDocument {
        sha256: hex_sha256(bytes),
        reads,
    })
}

fn to_variant(read: &str, contig: &str, raw: RawVariant) -> Result<Variant> {
    let invalid = || VariantsError::InvalidVariant {
        read: read.to_owned(),
        position: raw.position,
    };
    let kind = match raw.kind.as_str() {
        "SNV" => VariantKind::Snv,
        "INS" => VariantKind::Ins,
        "DEL" => VariantKind::Del,
        _ => return Err(invalid().into()),
    };
    let canonical = |allele: &str| !allele.is_empty() && allele.chars().all(is_canonical);
    if raw.position == 0 || !canonical(&raw.reference) || !canonical(&raw.alternate) {
        return Err(invalid().into());
    }
    Ok(Variant {
        contig: contig.to_owned(),
        position_1based: raw.position,
        reference: raw.reference,
        alternate: raw.alternate,
        kind,
    })
}

#[cfg(test)]
mod tests {
    use crate::model::reference::ReferenceTopology;

    use super::*;

    fn reference() -> Reference {
        Reference {
            name: "rCRS".into(),
            sequence: "ACGT".into(),
            topology: ReferenceTopology::Circular,
            sequence_sha256: hex_sha256(b"ACGT"),
        }
    }

    fn document(schema: &str, sample: &str, kind: &str) -> String {
        format!(
            r#"{{"schema_version":"{schema}","sample_id":"{sample}","provenance":{{"reference":{{"name":"rCRS","topology":"circular","sha256":"{}"}}}},"reads":[{{"name":"HV1","sha256":"00","variants":[{{"position":2,"reference":"C","alternate":"T","kind":"{kind}","eligible":true,"exclusion_reasons":[]}},{{"position":3,"reference":"G","alternate":"A","kind":"SNV","eligible":false,"exclusion_reasons":["read_end"]}}]}}]}}"#,
            hex_sha256(b"ACGT")
        )
    }

    fn read(text: &str, sample: &str) -> Result<VariantsDocument> {
        parse(
            text.as_bytes(),
            "s.variants.json".into(),
            sample,
            &reference(),
        )
    }

    #[test]
    fn keeps_eligible_variants_on_the_reference_contig() -> Result<()> {
        let text = document("dna.variants/v1", "s", "SNV");
        let parsed = read(&text, "s")?;
        assert_eq!(parsed.sha256, hex_sha256(text.as_bytes()));
        assert_eq!(parsed.reads[0].name, "HV1");
        assert_eq!(
            parsed.reads[0].variants,
            [Variant {
                contig: "rCRS".into(),
                position_1based: 2,
                reference: "C".into(),
                alternate: "T".into(),
                kind: VariantKind::Snv,
            }]
        );
        Ok(())
    }

    #[test]
    fn rejects_other_contracts_samples_references_and_variants() {
        assert!(read(&document("dna.sample_evidence/v10", "s", "SNV"), "s").is_err());
        assert!(read(&document("dna.variants/v1", "other", "SNV"), "s").is_err());
        assert!(read(&document("dna.variants/v1", "s", "MNV"), "s").is_err());
        assert!(read("{}", "s").is_err());
        let foreign = document("dna.variants/v1", "s", "SNV").replace(&hex_sha256(b"ACGT"), "ff");
        assert!(read(&foreign, "s").is_err());
    }
}
