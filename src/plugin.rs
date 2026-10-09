//! Plugin descriptors and their compile-time validation (ADR-0069).
//!
//! DNA composes three coarse plugin families at compile time: modality
//! plugins that produce per-read evidence, the core caller, and post-calling
//! plugins. Each plugin declares one descriptor: its identity, family, method
//! version, the contracts it provides and requires, and the configuration
//! sections it owns. The composer lists the registered descriptors and the
//! workflow compositions and validates both in constant evaluation:
//! identities and configuration sections are unique, and each plugin's
//! required contracts are read from the workflow's input documents or provided
//! by a plugin that runs before it.

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
pub(crate) enum Contract {
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
    /// Contracts the plugin produces.
    pub(crate) provides: &'static [Contract],
    /// Contracts the plugin needs before it runs.
    pub(crate) requires: &'static [Contract],
    /// Top-level configuration sections the plugin owns.
    pub(crate) config_sections: &'static [&'static str],
}

/// Fails compilation unless identities and configuration sections are unique
/// and every required contract is provided by a registered plugin.
pub(crate) const fn validate_registry(registry: &[&PluginDescriptor]) {
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
pub(crate) const fn composition(
    registry: &[&PluginDescriptor],
    inputs: &[Contract],
    plugins: &'static [&'static PluginDescriptor],
) -> &'static [&'static PluginDescriptor] {
    let mut index = 0;
    while index < plugins.len() {
        let plugin = plugins[index];
        assert!(
            registered(registry, plugin.id),
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

const fn registered(registry: &[&PluginDescriptor], id: &str) -> bool {
    let mut index = 0;
    while index < registry.len() {
        if same(registry[index].id, id) {
            return true;
        }
        index += 1;
    }
    false
}

pub(crate) const fn provided_by(required: &[Contract], plugins: &[&PluginDescriptor]) -> bool {
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
    use super::*;

    const MODALITY: PluginDescriptor = PluginDescriptor {
        id: "modality",
        family: PluginFamily::Modality,
        version: 1,
        provides: &[Contract::ReadEvidence],
        requires: &[],
        config_sections: &["modality"],
    };
    const CALLER: PluginDescriptor = PluginDescriptor {
        id: "caller",
        family: PluginFamily::Core,
        version: 1,
        provides: &[Contract::CalledVariants],
        requires: &[Contract::ReadEvidence],
        config_sections: &["caller"],
    };

    #[test]
    fn compares_identities_and_sections_bytewise() {
        assert!(same("core", "core"));
        assert!(!same("core", "cores"));
        assert!(contains_section(MODALITY.config_sections, "modality"));
        assert!(!contains_section(CALLER.config_sections, "modality"));
        assert!(registered(&[&MODALITY, &CALLER], "caller"));
        assert!(!registered(&[&MODALITY], "caller"));
    }

    #[test]
    fn accepts_requirements_from_inputs_or_earlier_plugins() {
        assert!(provided_by(CALLER.requires, &[&MODALITY]));
        assert!(!provided_by(CALLER.requires, &[]));
        assert!(listed(&[Contract::ReadEvidence], Contract::ReadEvidence));
        assert_eq!(
            composition(&[&MODALITY, &CALLER], &[Contract::ReadEvidence], &[&CALLER]).len(),
            1
        );
    }
}
