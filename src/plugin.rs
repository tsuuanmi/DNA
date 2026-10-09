//! Static plugin registry of the composition kernel (ADR-0069).
//!
//! DNA composes three coarse plugin families at compile time: modality
//! plugins that produce per-read evidence, the core caller, and post-calling
//! plugins. Each descriptor names the plugin, its family, its method version,
//! the contracts it provides and requires, and the configuration sections it
//! owns. The registry and every workflow composition are validated at compile
//! time: identities and configuration sections are unique, and each plugin's
//! required contracts are provided by a plugin that runs before it.

/// Coarse plugin family.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum PluginFamily {
    /// Turns raw data of one sequencing modality into per-read evidence.
    Modality,
    /// Calls variants from per-read evidence and aggregates samples.
    Core,
    /// Consumes called variants.
    PostCalling,
}

impl PluginFamily {
    /// Public label used in result provenance.
    pub(crate) const fn label(self) -> &'static str {
        match self {
            Self::Modality => "modality",
            Self::Core => "core",
            Self::PostCalling => "post_calling",
        }
    }
}

/// Data contract exchanged between plugins.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Contract {
    /// Per-read modality evidence (`read_evidence::ReadEvidence`).
    ReadEvidence,
    /// Called variants (`variant::CalledVariantSet`).
    CalledVariants,
    /// Haplotype-preserving normalized variants.
    NormalizedVariants,
    /// Target nomenclature representation.
    Nomenclature,
    /// Notation-convention findings.
    Conformance,
}

/// Identity, interface, and configuration ownership of one plugin.
#[derive(Debug)]
pub(crate) struct PluginDescriptor {
    /// Stable plugin identity.
    pub(crate) id: &'static str,
    pub(crate) family: PluginFamily,
    /// Method version: increases whenever the plugin's output can change for
    /// the same inputs and configuration.
    pub(crate) version: u32,
    provides: &'static [Contract],
    requires: &'static [Contract],
    /// Top-level configuration sections the plugin owns.
    config_sections: &'static [&'static str],
}

/// Sanger ABIF modality: basecalling, signal processing, callability, quality
/// control, and the evidence adapter.
pub(crate) const SANGER: PluginDescriptor = PluginDescriptor {
    id: "sanger",
    family: PluginFamily::Modality,
    version: 1,
    provides: &[Contract::ReadEvidence],
    requires: &[],
    config_sections: &[
        "basecalling",
        "signal_processing",
        "callability",
        "quality_control",
        "sanger_evidence",
    ],
};

/// Reviewed consensus sequences in FASTA.
pub(crate) const SEQUENCE: PluginDescriptor = PluginDescriptor {
    id: "sequence",
    family: PluginFamily::Modality,
    version: 1,
    provides: &[Contract::ReadEvidence],
    requires: &[],
    config_sections: &[],
};

/// Core caller: evidence-profile alignment, per-read variant calling, and
/// sample aggregation.
pub(crate) const CORE: PluginDescriptor = PluginDescriptor {
    id: "core",
    family: PluginFamily::Core,
    version: 1,
    provides: &[Contract::CalledVariants],
    requires: &[Contract::ReadEvidence],
    config_sections: &["alignment", "sample_reconciliation", "variant_calling"],
};

/// Haplotype-preserving variant normalization.
pub(crate) const NORMALIZATION: PluginDescriptor = PluginDescriptor {
    id: "normalization",
    family: PluginFamily::PostCalling,
    version: 1,
    provides: &[Contract::NormalizedVariants],
    requires: &[Contract::CalledVariants],
    config_sections: &[],
};

/// Profile-driven target nomenclature.
pub(crate) const NOMENCLATURE: PluginDescriptor = PluginDescriptor {
    id: "nomenclature",
    family: PluginFamily::PostCalling,
    version: 1,
    provides: &[Contract::Nomenclature],
    requires: &[Contract::NormalizedVariants],
    config_sections: &[],
};

/// Notation-convention checks declared by the target profile.
pub(crate) const CONFORMANCE: PluginDescriptor = PluginDescriptor {
    id: "conformance",
    family: PluginFamily::PostCalling,
    version: 1,
    provides: &[Contract::Conformance],
    requires: &[Contract::Nomenclature],
    config_sections: &[],
};

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
pub(crate) const BASECALL: &[&PluginDescriptor] = composition(&[], &[&SANGER]);
/// Plugins of the `analyze` workflow, in execution order.
pub(crate) const ANALYZE: &[&PluginDescriptor] = composition(&[], &[&SANGER, &CORE]);
/// Plugins of the `sample` workflow without notation, in execution order.
pub(crate) const SAMPLE: &[&PluginDescriptor] = composition(&[], &[&SANGER, &CORE]);
/// Plugins of the `call` workflow without notation, in execution order.
pub(crate) const CALL: &[&PluginDescriptor] = composition(&[], &[&SEQUENCE, &CORE]);
/// Plugins of the `call` workflow with notation, in execution order.
pub(crate) const CALL_WITH_NOTATION: &[&PluginDescriptor] =
    composition(&[], &[&SEQUENCE, &CORE, &NORMALIZATION, &NOMENCLATURE]);
/// Plugins of the `notation` workflow over a variants document, in execution order.
pub(crate) const NOTATION: &[&PluginDescriptor] = composition(
    &[Contract::CalledVariants],
    &[&NORMALIZATION, &NOMENCLATURE],
);
/// Plugins of the `notation` workflow with conformance checks, in execution order.
pub(crate) const NOTATION_WITH_CONFORMANCE: &[&PluginDescriptor] = composition(
    &[Contract::CalledVariants],
    &[&NORMALIZATION, &NOMENCLATURE, &CONFORMANCE],
);
/// Plugins of the `sample` workflow with notation, in execution order.
pub(crate) const SAMPLE_WITH_NOTATION: &[&PluginDescriptor] =
    composition(&[], &[&SANGER, &CORE, &NORMALIZATION, &NOMENCLATURE]);

const _: () = validate_registry(REGISTRY);

/// Fails compilation unless identities and configuration sections are unique
/// and every required contract is provided by a registered plugin.
const fn validate_registry(registry: &[&PluginDescriptor]) {
    let mut index = 0;
    while index < registry.len() {
        let plugin = registry[index];
        let mut other = index + 1;
        while other < registry.len() {
            assert!(
                !same(plugin.id, registry[other].id),
                "plugin identities must be unique"
            );
            let mut section = 0;
            while section < plugin.config_sections.len() {
                assert!(
                    !contains_section(
                        registry[other].config_sections,
                        plugin.config_sections[section]
                    ),
                    "a configuration section must have one owning plugin"
                );
                section += 1;
            }
            other += 1;
        }
        assert!(
            provided_by(plugin.requires, registry),
            "every required contract must be provided by a registered plugin"
        );
        index += 1;
    }
}

/// Returns a workflow composition after checking at compile time that its
/// plugins are registered and that each plugin's required contracts are read
/// from the workflow's input documents or provided by a plugin running before it.
const fn composition(
    inputs: &[Contract],
    plugins: &'static [&'static PluginDescriptor],
) -> &'static [&'static PluginDescriptor] {
    let mut index = 0;
    while index < plugins.len() {
        let plugin = plugins[index];
        assert!(
            registered(plugin.id),
            "a composition may use only registered plugins"
        );
        let (earlier, _) = plugins.split_at(index);
        let mut required = 0;
        while required < plugin.requires.len() {
            let contract = plugin.requires[required];
            assert!(
                listed(inputs, contract) || provided_by(&[contract], earlier),
                "a plugin's required contracts must be provided before it runs"
            );
            required += 1;
        }
        index += 1;
    }
    plugins
}

const fn listed(contracts: &[Contract], contract: Contract) -> bool {
    let mut index = 0;
    while index < contracts.len() {
        if contracts[index] as u8 == contract as u8 {
            return true;
        }
        index += 1;
    }
    false
}

const fn registered(id: &str) -> bool {
    let mut index = 0;
    while index < REGISTRY.len() {
        if same(REGISTRY[index].id, id) {
            return true;
        }
        index += 1;
    }
    false
}

const fn provided_by(required: &[Contract], plugins: &[&PluginDescriptor]) -> bool {
    let mut index = 0;
    while index < required.len() {
        let mut found = false;
        let mut plugin = 0;
        while plugin < plugins.len() {
            let mut offered = 0;
            while offered < plugins[plugin].provides.len() {
                if plugins[plugin].provides[offered] as u8 == required[index] as u8 {
                    found = true;
                }
                offered += 1;
            }
            plugin += 1;
        }
        if !found {
            return false;
        }
        index += 1;
    }
    true
}

const fn contains_section(sections: &[&str], section: &str) -> bool {
    let mut index = 0;
    while index < sections.len() {
        if same(sections[index], section) {
            return true;
        }
        index += 1;
    }
    false
}

const fn same(left: &str, right: &str) -> bool {
    let (left, right) = (left.as_bytes(), right.as_bytes());
    if left.len() != right.len() {
        return false;
    }
    let mut index = 0;
    while index < left.len() {
        if left[index] != right[index] {
            return false;
        }
        index += 1;
    }
    true
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;

    use super::*;

    #[test]
    fn registry_sections_are_exactly_the_shipped_configuration_sections() {
        let shipped: toml::Table = toml::from_str(include_str!("../config/dna.toml"))
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

    #[test]
    fn compares_identities_and_sections_bytewise() {
        assert!(same("core", "core"));
        assert!(!same("core", "cores"));
        assert!(contains_section(SANGER.config_sections, "callability"));
        assert!(!contains_section(CORE.config_sections, "callability"));
        assert!(registered("nomenclature"));
        assert!(!registered("ngs"));
    }
}
