# Target Profiles

Owns loading and validation of target profiles: strict versioned TOML that
holds what is known about one sequencing target — the reference it is validated
against and its topology, the reportable regions, and the optional
representation chain (indel placement, nomenclature windows, notation style).

Key children: `raw.rs` parses the TOML records and validates every rule;
`window.rs` holds the validated nomenclature windows. Each window rule carries
the structure it needs, so an invalid rule/structure pairing cannot be
constructed. The public `Profile` type is opaque and built only by
`Profile::load`; `ProfileIdentity` (`id` plus file SHA-256) is what results
record.

Consumers: `input::sanger` loads the profile named by the configuration and
checks the reference against it; `variant_calling` receives its regions;
`variant_nomenclature` runs its windows; `pipeline::sample_notation` follows its
normalization and notation choices; `report` records its identity.

This module is in the target-data layer (ADR-0064): it depends only on `model`,
`error`, `checksum`, and the reference-length cap in `config`. It names its own
vocabulary (`IndelPlacement`, `NotationStyle`); `variant_normalization` maps
the placement to its `NormalizationPolicy`. It does not run any rule, load
references, or read the scientific configuration.

See [profile requirements](../../docs/requirements/profiles.md),
[profile reference](../../docs/reference/profiles.md), and
[ADR-0063](../../docs/decisions/adr/0063-target-profiles.md).
