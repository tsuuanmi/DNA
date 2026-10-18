//! Per-read profile representation for the notation view of sample, call,
//! and notation documents.
//!
//! ADR-0060 composes representation per read before sample reconciliation:
//! each read's eligible calls are normalized and then given the profile's
//! window representation, so sequence-equivalent descriptions from different
//! reads converge on one notation.

use crate::report::{ReadRepresentation, SampleNotation};
use dna_call::model::called_read::CalledRead;
use std::ops::Range;

use crate::error::Result;
use dna_kernel::model::reference::Reference;
use dna_kernel::profile::{Notation, Profile};
use dna_kernel::variant::{CalledVariantSet, Variant};
use dna_normalize::variant_nomenclature::{self, from_normalization};
use dna_normalize::variant_normalization::{NormalizationPolicy, TradeBounds, normalize_read};

/// Represents every called read's eligible calls, or `None` when the profile
/// declares no notation. The reference has already been checked against the
/// profile.
pub(crate) fn represent(
    reads: &[CalledRead],
    reference: &Reference,
    profile: &Profile,
) -> Result<Option<SampleNotation>> {
    represent_variants(
        reads.iter().map(|read| {
            (
                read.input_sha256.clone(),
                read.variants.reported.iter().map(Variant::from).collect(),
                read.alignment
                    .callable_segments
                    .iter()
                    .map(|segment| segment.start_0based..segment.end_0based_exclusive)
                    .collect(),
            )
        }),
        reference,
        profile,
    )
}

/// Represents each read's eligible variants, given with its content identity
/// and its callable reference segments (0-based, half-open).
pub(crate) fn represent_variants(
    reads: impl IntoIterator<Item = (String, Vec<Variant>, Vec<Range<usize>>)>,
    reference: &Reference,
    profile: &Profile,
) -> Result<Option<SampleNotation>> {
    let Some(Notation {
        indel_placement,
        style,
    }) = profile.notation
    else {
        return Ok(None);
    };
    let policy = NormalizationPolicy::for_placement(indel_placement);
    let reads = reads
        .into_iter()
        .map(|(input_sha256, variants, callable)| {
            Ok(ReadRepresentation {
                input_sha256,
                variants: represent_read(reference, profile, policy, variants, &callable)?,
            })
        })
        .collect::<Result<_>>()?;
    Ok(Some(SampleNotation {
        style,
        conventions: profile.conventions.clone(),
        reads,
    }))
}

/// Normalizes one read's calls, then applies the profile windows. A window an
/// edit straddles keeps its normalized form (SRS-NOM-007), and so does one
/// whose representation would place a call where the read reads nothing
/// (SRS-NOM-010): a window rule rewrites the whole window's sequence, which over
/// an unread base is the reference's, not the read's. For the same reason a
/// deletion trade crosses only read bases, and it crosses no window boundary,
/// leaving a window's representation to its rule (SRS-VN-008).
fn represent_read(
    reference: &Reference,
    profile: &Profile,
    policy: NormalizationPolicy,
    variants: Vec<Variant>,
    callable: &[Range<usize>],
) -> Result<Vec<Variant>> {
    let called = CalledVariantSet {
        reference: reference.identity(),
        variants,
    };
    let windows = profile
        .windows
        .iter()
        .map(|window| window.start..window.end())
        .collect::<Vec<_>>();
    let normalized = normalize_read(
        reference,
        &called,
        policy,
        TradeBounds {
            windows: &windows,
            readable: Some(callable),
        },
    )?;
    Ok(variant_nomenclature::apply_read(
        reference,
        profile,
        from_normalization(&normalized),
        callable,
    )?
    .represented_variants)
}

#[cfg(test)]
mod tests {
    use std::cmp::Ordering;
    use std::path::Path;

    use dna_call::model::variant::{Variant as CalledVariant, VariantKind};
    use dna_kernel::model::reference::ReferenceTopology;
    use dna_kernel::profile::fixtures::human_mtdna;
    use dna_kernel::reference;
    use dna_kernel::variant;

    use super::*;

    /// A read of the whole reference, and one of HVS-I.
    const WHOLE: Range<usize> = 0..16_569;
    const READ: Range<usize> = 16_043..16_379;

    type TestResult = std::result::Result<(), Box<dyn std::error::Error>>;

    fn rcrs() -> Result<Reference> {
        Ok(reference::load(
            &Path::new(env!("CARGO_MANIFEST_DIR")).join("references/rCRS.fasta"),
            ReferenceTopology::Circular,
        )?)
    }

    fn called(
        position: usize,
        reference: &str,
        alternate: &str,
        kind: VariantKind,
    ) -> CalledVariant {
        CalledVariant {
            contig: "rCRS".into(),
            position_1based: position,
            reference: reference.into(),
            alternate: alternate.into(),
            kind,
            calls: Vec::new(),
        }
    }

    fn public(position: usize, reference: &str, alternate: &str) -> Variant {
        Variant {
            contig: "rCRS".into(),
            position_1based: position,
            reference: reference.into(),
            alternate: alternate.into(),
            kind: match reference.len().cmp(&alternate.len()) {
                Ordering::Less => variant::VariantKind::Ins,
                Ordering::Greater => variant::VariantKind::Del,
                Ordering::Equal => variant::VariantKind::Snv,
            },
        }
    }

    #[test]
    fn converges_a_left_run_c_insertion_on_309() -> TestResult {
        let reference = rcrs()?;
        let represented = represent_read(
            &reference,
            &human_mtdna()?,
            NormalizationPolicy::RightAligned,
            vec![Variant::from(&called(303, "C", "CC", VariantKind::Ins))],
            &[WHOLE],
        )?;
        assert_eq!(represented, [public(309, "C", "CC")]);
        Ok(())
    }

    #[test]
    fn an_edit_crossing_the_hv2_window_keeps_its_form_and_hv3_is_represented() -> TestResult {
        let reference = rcrs()?;
        // A deletion of 302-303 crosses HVS-II and stays as normalized; the
        // AC-unit deletion at 515-516 still takes HVS-III's 523DEL 524DEL.
        let crossing = public(301, "AAC", "A");
        let represented = represent_read(
            &reference,
            &human_mtdna()?,
            NormalizationPolicy::RightAligned,
            vec![
                Variant::from(&called(301, "AAC", "A", VariantKind::Del)),
                Variant::from(&called(514, "CAC", "C", VariantKind::Del)),
            ],
            &[WHOLE],
        )?;
        assert_eq!(represented, [crossing, public(522, "CAC", "C")]);
        Ok(())
    }

    #[test]
    fn a_window_whose_form_would_call_an_unread_base_keeps_its_form() -> TestResult {
        let reference = rcrs()?;
        // 16183C read before an unread stretch 16184-16191: the HVS-I rule
        // would write it as 16183DEL 16188.1C, a call where the read shows
        // nothing, so the window keeps 16183C.
        let call = called(16183, "A", "C", VariantKind::Snv);
        let partial = represent_read(
            &reference,
            &human_mtdna()?,
            NormalizationPolicy::RightAligned,
            vec![Variant::from(&call)],
            &[16_043..16_183, 16_191..16_379],
        )?;
        assert_eq!(partial, [public(16183, "A", "C")]);
        let whole = represent_read(
            &reference,
            &human_mtdna()?,
            NormalizationPolicy::RightAligned,
            vec![Variant::from(&call)],
            &[READ],
        )?;
        assert_eq!(whole, [public(16182, "AA", "A"), public(16188, "C", "CC")]);
        Ok(())
    }

    #[test]
    fn a_window_whose_form_rests_on_an_unread_base_keeps_its_form() -> TestResult {
        let reference = rcrs()?;
        // 16182C 16183C read, but 16189-16192 unread: the HVS-I rule would
        // take T16189 from the reference and write 16182DEL 16183DEL
        // 16188.1C 16188.2C, a form that holds only if 16189 is T.
        let calls = vec![
            Variant::from(&called(16182, "A", "C", VariantKind::Snv)),
            Variant::from(&called(16183, "A", "C", VariantKind::Snv)),
        ];
        let partial = represent_read(
            &reference,
            &human_mtdna()?,
            NormalizationPolicy::RightAligned,
            calls.clone(),
            &[16_077..16_188, 16_192..16_379],
        )?;
        assert_eq!(partial, [public(16182, "A", "C"), public(16183, "A", "C")]);
        let whole = represent_read(
            &reference,
            &human_mtdna()?,
            NormalizationPolicy::RightAligned,
            calls,
            &[READ],
        )?;
        assert_eq!(
            whole,
            [public(16181, "AAA", "A"), public(16188, "C", "CCC")]
        );
        Ok(())
    }

    #[test]
    fn a_deletion_at_the_hv1_window_end_does_not_trade_out_of_it() -> TestResult {
        let reference = rcrs()?;
        // 16193DEL 16194T is the same sequence as 16193T 16194DEL, but the
        // trade would carry the deletion past HVS-I's last base, 16193.
        let calls = vec![
            Variant::from(&called(16192, "CC", "C", VariantKind::Del)),
            Variant::from(&called(16194, "A", "T", VariantKind::Snv)),
        ];
        let represented = represent_read(
            &reference,
            &human_mtdna()?,
            NormalizationPolicy::RightAligned,
            calls,
            &[READ],
        )?;
        assert_eq!(
            represented,
            [public(16192, "CC", "C"), public(16194, "A", "T")]
        );
        Ok(())
    }

    #[test]
    fn produces_notation_only_when_the_profile_declares_it() -> TestResult {
        let reference = rcrs()?;
        let mut profile = human_mtdna()?;
        assert!(represent(&[], &reference, &profile)?.is_some());
        profile.notation = None;
        assert!(represent(&[], &reference, &profile)?.is_none());
        Ok(())
    }
}
