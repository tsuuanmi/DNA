//! The registered plugins and the plugin composition of every workflow,
//! validated at compile time (ADR-0069).

use crate::input::sequence::PLUGIN as SEQUENCE;
use crate::report::REVIEW;
use dna_core::consensus::PLUGIN as CONSENSUS_PLUGIN;
use dna_core::consensus::decided::PLUGIN as DECIDED;
use dna_core::read_call::PLUGIN as CORE;
use dna_empop::conformance::PLUGIN as CONFORMANCE;
use dna_empop::export::PLUGIN as EMPOP;
use dna_kernel::plugin::{Contract, PluginDescriptor, composition, validate_registry};
use dna_normalize::variant_nomenclature::PLUGIN as NOMENCLATURE;
use dna_normalize::variant_normalization::PLUGIN as NORMALIZATION;
use dna_sanger::read_processing::PLUGIN as SANGER;

/// Every plugin this build contains.
const REGISTRY: &[&PluginDescriptor] = &[
    &SANGER,
    &SEQUENCE,
    &CORE,
    &CONSENSUS_PLUGIN,
    &DECIDED,
    &NORMALIZATION,
    &NOMENCLATURE,
    &CONFORMANCE,
    &EMPOP,
    &REVIEW,
];

/// Plugins of the `basecall` workflow, in execution order.
pub(crate) const BASECALL: &[&PluginDescriptor] = composition(REGISTRY, &[], &[&SANGER]);
/// Plugins of the `analyze` workflow, in execution order.
pub(crate) const ANALYZE: &[&PluginDescriptor] = composition(REGISTRY, &[], &[&SANGER, &CORE]);
/// Plugins of the `sample` workflow without notation, in execution order.
pub(crate) const SAMPLE: &[&PluginDescriptor] = composition(REGISTRY, &[], &[&SANGER, &CORE]);
/// Plugins of the `consensus` workflow, in execution order.
pub(crate) const CONSENSUS: &[&PluginDescriptor] =
    composition(REGISTRY, &[], &[&SANGER, &CORE, &CONSENSUS_PLUGIN]);
/// Plugins of the `call` workflow without notation, in execution order.
pub(crate) const CALL: &[&PluginDescriptor] = composition(REGISTRY, &[], &[&SEQUENCE, &CORE]);
/// Plugins of the `call` workflow with notation, in execution order.
pub(crate) const CALL_WITH_NOTATION: &[&PluginDescriptor] = composition(
    REGISTRY,
    &[],
    &[&SEQUENCE, &CORE, &NORMALIZATION, &NOMENCLATURE],
);
/// Plugins of the `call` workflow over a consensus document, in execution order.
pub(crate) const CALL_CONSENSUS: &[&PluginDescriptor] =
    composition(REGISTRY, &[Contract::Consensus], &[&DECIDED]);
/// Plugins of the `call` workflow over a consensus document with notation, in
/// execution order.
pub(crate) const CALL_CONSENSUS_WITH_NOTATION: &[&PluginDescriptor] = composition(
    REGISTRY,
    &[Contract::Consensus],
    &[&DECIDED, &NORMALIZATION, &NOMENCLATURE],
);
/// Plugins of the `notation` workflow over a variants document, in execution order.
const NOTATION: &[&PluginDescriptor] = composition(
    REGISTRY,
    &[Contract::CalledVariants],
    &[&NORMALIZATION, &NOMENCLATURE],
);
/// Plugins of the `notation` workflow with conformance checks, in execution order.
const NOTATION_WITH_CONFORMANCE: &[&PluginDescriptor] = composition(
    REGISTRY,
    &[Contract::CalledVariants],
    &[&NORMALIZATION, &NOMENCLATURE, &CONFORMANCE],
);
/// Plugins of the `notation` workflow with the EMPOP export, in execution order.
const NOTATION_WITH_EMPOP: &[&PluginDescriptor] = composition(
    REGISTRY,
    &[Contract::CalledVariants],
    &[&NORMALIZATION, &NOMENCLATURE, &EMPOP],
);
/// Plugins of the `notation` workflow with conformance checks and the EMPOP
/// export, in execution order.
const NOTATION_WITH_CONFORMANCE_AND_EMPOP: &[&PluginDescriptor] = composition(
    REGISTRY,
    &[Contract::CalledVariants],
    &[&NORMALIZATION, &NOMENCLATURE, &CONFORMANCE, &EMPOP],
);

/// Plugins of the `notation` workflow, by the post-calling plugins the
/// profile declares.
pub(crate) const fn notation(
    conformance: bool,
    empop: bool,
) -> &'static [&'static PluginDescriptor] {
    match (conformance, empop) {
        (false, false) => NOTATION,
        (true, false) => NOTATION_WITH_CONFORMANCE,
        (false, true) => NOTATION_WITH_EMPOP,
        (true, true) => NOTATION_WITH_CONFORMANCE_AND_EMPOP,
    }
}
/// Plugins of the `sample` workflow with notation, in execution order.
pub(crate) const SAMPLE_WITH_NOTATION: &[&PluginDescriptor] = composition(
    REGISTRY,
    &[],
    &[&SANGER, &CORE, &NORMALIZATION, &NOMENCLATURE],
);

/// Plugins of the `review` workflow without notation, in execution order: the
/// consensus workflow, the call of its decided edits, then the review list.
const REVIEW_LIST: &[&PluginDescriptor] = composition(
    REGISTRY,
    &[],
    &[&SANGER, &CORE, &CONSENSUS_PLUGIN, &DECIDED, &REVIEW],
);
/// Plugins of the `review` workflow with notation, in execution order.
const REVIEW_WITH_NOTATION: &[&PluginDescriptor] = composition(
    REGISTRY,
    &[],
    &[
        &SANGER,
        &CORE,
        &CONSENSUS_PLUGIN,
        &DECIDED,
        &NORMALIZATION,
        &NOMENCLATURE,
        &REVIEW,
    ],
);
/// Plugins of the `review` workflow with notation and the EMPOP query, in
/// execution order.
const REVIEW_WITH_EMPOP: &[&PluginDescriptor] = composition(
    REGISTRY,
    &[],
    &[
        &SANGER,
        &CORE,
        &CONSENSUS_PLUGIN,
        &DECIDED,
        &NORMALIZATION,
        &NOMENCLATURE,
        &EMPOP,
        &REVIEW,
    ],
);

/// Plugins of the `review` workflow, by the post-calling plugins the profile
/// declares, as `SampleAnalysis::call` runs them: the EMPOP query needs the
/// notation.
pub(crate) const fn review(notation: bool, empop: bool) -> &'static [&'static PluginDescriptor] {
    match (notation, empop) {
        (false, _) => REVIEW_LIST,
        (true, false) => REVIEW_WITH_NOTATION,
        (true, true) => REVIEW_WITH_EMPOP,
    }
}

const _: () = validate_registry(REGISTRY);

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    use dna_kernel::plugin::provided_by;

    use super::*;

    #[test]
    fn registry_sections_are_exactly_the_shipped_configuration_sections() {
        let shipped: toml::Table = toml::from_str(include_str!("../../config/dna.toml"))
            .unwrap_or_else(|error| panic!("shipped configuration parses: {error}"));
        let sections: BTreeSet<&str> = shipped
            .iter()
            .filter(|(_, value)| value.is_table())
            .map(|(key, _)| key.as_str())
            .collect();
        let owned: BTreeSet<&str> = REGISTRY
            .iter()
            .flat_map(|plugin| plugin.config_sections.iter().copied())
            .collect();
        assert_eq!(sections, owned);
    }

    #[test]
    fn checks_requirements_against_earlier_plugins() {
        assert!(provided_by(CORE.requires, &[&SANGER]));
        assert!(!provided_by(CORE.requires, &[&NORMALIZATION]));
        assert!(!provided_by(NOMENCLATURE.requires, &[&SANGER, &CORE]));
        assert!(provided_by(NOMENCLATURE.requires, SAMPLE_WITH_NOTATION));
    }
}
