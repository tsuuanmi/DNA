//! The registered plugins and the plugin composition of every workflow,
//! validated at compile time (ADR-0069).

use crate::conformance::PLUGIN as CONFORMANCE;
use crate::input::sequence::PLUGIN as SEQUENCE;
use crate::plugin::{Contract, PluginDescriptor, composition, validate_registry};
use crate::read_call::PLUGIN as CORE;
use crate::read_processing::PLUGIN as SANGER;
use crate::variant_nomenclature::PLUGIN as NOMENCLATURE;
use crate::variant_normalization::PLUGIN as NORMALIZATION;

/// Every plugin this build contains.
const REGISTRY: &[&PluginDescriptor] = &[
    &SANGER,
    &SEQUENCE,
    &CORE,
    &NORMALIZATION,
    &NOMENCLATURE,
    &CONFORMANCE,
];

/// Plugins of the `basecall` workflow, in execution order.
pub(crate) const BASECALL: &[&PluginDescriptor] = composition(REGISTRY, &[], &[&SANGER]);
/// Plugins of the `analyze` workflow, in execution order.
pub(crate) const ANALYZE: &[&PluginDescriptor] = composition(REGISTRY, &[], &[&SANGER, &CORE]);
/// Plugins of the `sample` workflow without notation, in execution order.
pub(crate) const SAMPLE: &[&PluginDescriptor] = composition(REGISTRY, &[], &[&SANGER, &CORE]);
/// Plugins of the `call` workflow without notation, in execution order.
pub(crate) const CALL: &[&PluginDescriptor] = composition(REGISTRY, &[], &[&SEQUENCE, &CORE]);
/// Plugins of the `call` workflow with notation, in execution order.
pub(crate) const CALL_WITH_NOTATION: &[&PluginDescriptor] = composition(
    REGISTRY,
    &[],
    &[&SEQUENCE, &CORE, &NORMALIZATION, &NOMENCLATURE],
);
/// Plugins of the `notation` workflow over a variants document, in execution order.
pub(crate) const NOTATION: &[&PluginDescriptor] = composition(
    REGISTRY,
    &[Contract::CalledVariants],
    &[&NORMALIZATION, &NOMENCLATURE],
);
/// Plugins of the `notation` workflow with conformance checks, in execution order.
pub(crate) const NOTATION_WITH_CONFORMANCE: &[&PluginDescriptor] = composition(
    REGISTRY,
    &[Contract::CalledVariants],
    &[&NORMALIZATION, &NOMENCLATURE, &CONFORMANCE],
);
/// Plugins of the `sample` workflow with notation, in execution order.
pub(crate) const SAMPLE_WITH_NOTATION: &[&PluginDescriptor] = composition(
    REGISTRY,
    &[],
    &[&SANGER, &CORE, &NORMALIZATION, &NOMENCLATURE],
);

const _: () = validate_registry(REGISTRY);

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    use crate::plugin::provided_by;

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
