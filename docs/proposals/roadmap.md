# DNA Roadmap

The roadmap is non-normative. It records current validation priorities and intentionally deferred product directions. Current behavior is defined by the [SRS](../requirements/SRS.md), [methods](../design/README.md), [contracts](../reference/README.md), and accepted [ADRs](../decisions/adr/README.md).

## Current objective

Strengthen evidence-backed confidence in the implemented Sanger analysis path before expanding biological interpretation.

The core confidence floor remains the read-level path defined by [ADR-0021](../decisions/adr/0021-scientific-core-confidence-floor.md), but the current supported baseline already includes intentional capabilities beyond that floor. The roadmap does not redefine those capabilities.

## Validation priorities

### Approved real-trace baseline

Establish and maintain approved real-AB1 evidence for the core path. For each approved trace, retain the provenance required by [data governance](../governance/data.md), including source context, AB1/reference/configuration identity, expected sequence or independently established truth, expected straightforward variants, and explained disagreements.

### Stage confidence

Prioritize evidence that demonstrates:

- ABIF/channel/PLOC decoding matches the intended container and scientific-tag semantics;
- signal-derived re-calling is deterministic and preserves unresolved evidence;
- trimming removes justified tails without silently rewriting internal evidence;
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
- production mtDNA poly-C detector/state/recovery/weighting;
- sample-level consensus and adjudicated sample variants;
- quantitative heteroplasmy;
- mixed-template or length-mixture decomposition;
- haplogroup-based QC or inference;
- calibrated quality/error probabilities;
- ML-based calling, correction, or training export;
- SCF, VCF/BCF, multi-contig references, and genome-scale indexing.

## Promotion rule

Exploratory work belongs under `docs/research/<topic>/`. A roadmap or research item becomes production behavior only through the documentation-governance promotion path: requirements/decision/method/contract changes as applicable, implementation and tests, validation evidence, and release evidence when user-visible.
