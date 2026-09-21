//! Original trace-call mappings derived from selected alignment columns.

use crate::error::{Error, Result};
use crate::model::alignment::AlignmentColumn;
use crate::model::variant::{VariantCallMapping, VariantCallRole};

/// Maps a difference-bearing query column to its original trace call.
pub(crate) fn supporting(column: &AlignmentColumn) -> Result<VariantCallMapping> {
    let call_index_0based = column.original_call_index_0based.ok_or_else(|| {
        Error::Variant("supporting query column lacks an original call index".into())
    })?;
    Ok(VariantCallMapping {
        role: VariantCallRole::Supporting,
        call_index_0based,
        reference_position_0based: column.reference_index_0based,
    })
}

/// Maps a reference-aligned query column used to bound an indel.
pub(crate) fn flanking(column: &AlignmentColumn) -> Result<VariantCallMapping> {
    let call_index_0based = column
        .original_call_index_0based
        .ok_or_else(|| Error::Variant("flanking column lacks an original call index".into()))?;
    let reference_position_0based = column
        .reference_index_0based
        .ok_or_else(|| Error::Variant("flanking call lacks a reference position".into()))?;
    Ok(VariantCallMapping {
        role: VariantCallRole::Flanking,
        call_index_0based,
        reference_position_0based: Some(reference_position_0based),
    })
}

/// Sorts and deduplicates stable call mappings.
pub(crate) fn sort_dedup(mut calls: Vec<VariantCallMapping>) -> Vec<VariantCallMapping> {
    calls.sort_unstable();
    calls.dedup();
    calls
}

/// Validates that SNV calls map to the substituted reference base.
pub(crate) fn validate_snv(calls: &[VariantCallMapping], position: usize) -> Result<()> {
    if calls.is_empty()
        || calls.iter().any(|call| {
            call.role != VariantCallRole::Supporting
                || call.reference_position_0based != Some(position)
        })
    {
        return Err(Error::Variant(
            "SNV calls must support the substituted reference position".into(),
        ));
    }
    Ok(())
}

/// Validates inserted calls and any aligned flanks.
pub(crate) fn validate_insertion(calls: &[VariantCallMapping]) -> Result<()> {
    let mut supporting = 0;
    for call in calls {
        match call.role {
            VariantCallRole::Supporting => {
                supporting += 1;
                if call.reference_position_0based.is_some() {
                    return Err(Error::Variant(
                        "inserted calls must not have a reference position".into(),
                    ));
                }
            }
            VariantCallRole::Flanking if call.reference_position_0based.is_some() => {}
            VariantCallRole::Flanking => {
                return Err(Error::Variant(
                    "insertion flanks must have a reference position".into(),
                ));
            }
        }
    }
    if supporting == 0 {
        return Err(Error::Variant(
            "insertion lacks a supporting trace call".into(),
        ));
    }
    Ok(())
}

/// Validates that deletion evidence consists only of aligned flanks.
pub(crate) fn validate_deletion(calls: &[VariantCallMapping]) -> Result<()> {
    if calls.is_empty()
        || calls.iter().any(|call| {
            call.role != VariantCallRole::Flanking || call.reference_position_0based.is_none()
        })
    {
        return Err(Error::Variant(
            "deletions must have only reference-aligned flanking calls".into(),
        ));
    }
    Ok(())
}
