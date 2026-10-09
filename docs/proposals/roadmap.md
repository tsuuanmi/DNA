# DNA Roadmap

The roadmap is non-normative. It records current validation priorities and intentionally deferred product directions. Current behavior is defined by the [SRS](../requirements/SRS.md), [methods](../design/README.md), [contracts](../reference/README.md), and accepted [ADRs](../decisions/adr/README.md).

## Current objective

Strengthen evidence-backed confidence in the implemented Sanger analysis path before expanding biological interpretation.

The core confidence floor remains the read-level path defined by [ADR-0021](../decisions/adr/0021-scientific-core-confidence-floor.md), but the current supported baseline already includes intentional capabilities beyond that floor. The roadmap does not redefine those capabilities.

## Architectural evolution

[PROP-0001](0001-modular-dna-analysis-platform.md), accepted by [ADR-0058](../decisions/adr/0058-canonical-contracts-and-modular-analysis-composition.md), defines the long-term direction for DNA as a modular analysis platform with a public Rust API, canonical typed contracts between capabilities, multiple input modalities such as Sanger and NGS, independently replaceable alignment and variant-calling implementations, replaceable scientific data providers, and downstream modules such as haplogroup and nomenclature.

The proposal is not current production architecture. Its implementation is intentionally incremental: preserve the validated Sanger path, introduce stable public and canonical boundaries first, add alternative implementations only where independent variation is real, and extract crates only when dependency or lifecycle boundaries justify them.

[PROP-0002](0002-plugin-first-architecture.md), accepted by [ADR-0069](../decisions/adr/0069-plugin-first-modality-core-post-calling.md), sequences that direction as three plugin families (modality, core caller, post-calling) with phased exit criteria. All phases are implemented: the modality → core evidence contract, the plugin registry with per-plugin configuration, a modality-neutral sample aggregation, a core-only `call` over reviewed consensus sequences, a post-calling `notation` with conformance findings, and the workspace split into plugin-family crates (ADR-0070).

## Validation priorities

### Approved real-trace baseline

Establish and maintain approved real-AB1 evidence for the core path. For each approved trace, retain the provenance required by [data governance](../governance/data.md), including source context, AB1/reference/configuration identity, expected sequence or independently established truth, expected straightforward variants, and explained disagreements.

### Stage confidence

Prioritize evidence that demonstrates:

- ABIF/channel/PLOC decoding matches the intended container and scientific-tag semantics;
- signal-derived re-calling is deterministic and preserves unresolved evidence;
- trimming removes justified tails without silently rewriting internal evidence, and read callability masks internal dephased or mixed stretches without rewriting evidence;
- forward/reverse placement is correct and ambiguous placement fails explicitly;
- reported SNVs and supported indels map back to the observed trace evidence and reference strand;
- sample evidence preserves independent read observations, coverage, overlap, differential loci, and normalized variant support without prematurely turning them into consensus or genotype claims.

Detailed acceptance rules belong to the SRS and method/contract documents rather than being copied here.

## Current production behavior

The roadmap does not maintain a second list of current capabilities. Current behavior is persistent production documentation owned by the [SRS](../requirements/SRS.md), [methods](../design/README.md), [contracts](../reference/README.md), and [../architecture/invariants](../architecture/README.md).

## Production-readiness track

Engineering maturity is developed in parallel with scientific confidence. The repository targets a production-grade delivery baseline with:

- pinned Rust release/MSRV toolchains and an explicit Ubuntu 24.04 runner baseline;
- warnings-clean all-feature Rust CI, strict dependency/source hygiene, and locked builds;
- Rust-only production source boundaries with Python isolated to research/validation/test tooling;
- cargo-deny, RustSec audit, and pull-request dependency review;
- immutable GitHub Action pins maintained by Dependabot;
- actionlint/zizmor workflow analysis, CodeQL, and OpenSSF Scorecard;
- ABIF fuzz smoke plus scheduled adversarial campaigns;
- protected-branch/tag governance and immutable-release policy;
- self-contained auditable release bundles with SBOMs, checksums, and cryptographic attestations.

These controls make the software build and supply chain defensible. They do not replace scientific release evidence. The production-ready milestone is reached only when an exact tagged revision also has approved real-AB1 ground-truth regression, difficult-locus/disagreement analysis, and runtime/peak-memory evidence recorded with the release template.

## Deferred product directions

The following remain outside the current production interpretation boundary until separately specified, decided, implemented, and validated:

- new or substantially more complex indel models;
- phase recovery or shadow-ladder deconvolution behind a detected shift, and calibrated weighting from callability features (detection and state are production behaviour under ADR-0067; research under `docs/research/phase-recovery/`);
- signal denoising and baseline correction (research spike under `docs/research/denoising/`);
- sample-level consensus and adjudicated sample variants;
- quantitative heteroplasmy;
- mixed-template or length-mixture decomposition;
- haplogroup-based QC or inference;
- calibrated quality/error probabilities;
- ML-based calling, correction, or training export;
- SCF, VCF/BCF, multi-contig references, and genome-scale indexing.

## Promotion rule

Exploratory work belongs under `docs/research/<topic>/`. A roadmap or research item becomes production behavior only through the documentation-governance promotion path: requirements/decision/method/contract changes as applicable, implementation and tests, validation evidence, and release evidence when user-visible.
