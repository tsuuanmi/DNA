# ADR-0016: Defer ML Feature Boundary to a Separate Versioned Training Contract

- **Status:** Accepted
- **Date:** 2026-08-26

## Context

The compact `dna.analysis/v5` and `dna.basecalls/v1` contracts are deliberately small, deterministic, and privacy-constrained (ADR-0014, ADR-0015). A future machine-learning path needs richer internal evidence—per-call peaks, rolling windows, spacing, and alignment context—than those production results expose. Adding bulk fields to the compact contracts would violate their one-result, no-compatibility-output boundary and increase identifying payload.

The repository has no approved truth-labeled corpus, no feature registry, no chosen first prediction target, and no trained model (ADR-0013, [`data.md`](../data.md)). AGENTS.md prohibits speculative features, single-use abstractions, and unused configuration. Unused Rust scaffolding would fail the clippy dead-code lint and the `docs/src` one-to-one mirror gate.

## Options

1. Add ML feature fields directly to `dna.analysis/v5` or `dna.basecalls/v1`.
2. Scaffold a `src/features/` module with typed extraction and engineering structs now, without a schema, CLI command, caller, or approved task.
3. Record the feature-boundary direction as a roadmap section and this ADR, keep production contracts unchanged, and defer all implementation until a task, corpus, feature registry, label registry, and closed schema are approved.

## Decision

Choose option 3. The ML-ready JSON direction in [`docs/roadmap.md`](../roadmap.md) documents a future opt-in, independently versioned training envelope such as `dna.training-example/v1`. One explicit feature subsystem sits after the existing scientific pipeline and before a `report::training` JSON projection:

- `features::extract` receives the complete immutable output bundle for one successful command and produces a private `ExtractedEvidence` value;
- `features::engineer` consumes only `ExtractedEvidence` plus a versioned feature definition and produces one canonical `FeatureSet`;
- `report::training` serializes the training example without recomputing or filtering.

`ExtractedEvidence` and `FeatureSet` are feature-subsystem types owned by `src/features/`, not `src/model/`. The `model/` module is the shared scientific domain vocabulary and serializable production result records, not transient post-pipeline algorithmic state. The top-level `crate::features` module is namespace-isolated from the existing private `crate::signal_processing::features` submodule (rolling-SNR calculation); `features` is the chosen name for the post-pipeline subsystem, and the conceptual collision is resolved by module path.

Labels remain independently governed data joined by an opaque example identity. The training example embeds a minimal, privacy-reviewed provenance and scientific projection, not the unchanged full command result, because the basecall result contains complete identifying sequences that contradict [`data.md`](../data.md).

No `features` module, training schema, CLI command, or training output is implemented until the roadmap delivery phases approve a task, corpus, feature registry, label registry, and closed schema.

## Consequences

The production contracts remain the single authority for `signal analyze` and `signal basecall`. No speculative code, unused module, or compatibility output is introduced. The roadmap records the intended boundary, feature catalog, leakage safeguards, privacy rules, and phased gates so a future implementation can proceed from an approved contract rather than ad hoc scaffolding.

A training exporter conflicts with the current one-file CLI contract: `SRS-IN-001`, `SRS-IN-010`, and `SRS-OUT-001` fix each command to one input and exactly one command-specific JSON output with no duplicate compatibility output. An opt-in exporter therefore requires a future ADR amendment or successor and an SRS update establishing its invocation contract (separate subcommand, output path, overwrite semantics) before any `report::training` or pipeline integration. The compact v5 and basecalls-v1 schemas are already closed (`additionalProperties: false`), so a training-keyed production document is already rejected at the schema-field level by the existing validator; no additional schema-field guard test is needed. CLI emission-boundary tests remain a future-implementation concern.

Implementation requires a later ADR or amendment, a closed Draft 2020-12 schema, synthetic example, `docs/src` mirror for any new Rust source, focused tests, and biological validation before any model alters scientific results. All new Rust sources, their mirrors, CLI wiring, and tests must land in one atomic commit so the docs-mirror and clippy dead-code gates stay green.

## Supersession

None. This ADR records direction only and defers to a future implementation ADR.
