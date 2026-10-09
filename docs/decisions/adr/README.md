# Architecture Decision Records (ADR)

ADRs are historical decision records. Current-state truth lives in requirements, architecture, design, reference, source, and validation documentation.

## ADR lifecycle and deduplication

ADRs record durable architectural or scientific decisions, not every implementation increment, schema revision, or research step. The repository preserves historical ADRs, but each current decision scope must have one clear authority.

Before creating a new ADR:

1. Search the SRS, architecture/invariants, methods/contracts, and this index for an existing decision with the same scope.
2. If the work only implements, documents, validates, or clarifies an accepted decision, update the owning SRS/method/contract/tests instead of creating another ADR.
3. If the same decision is being refined without replacing its core choice, amend the existing ADR with a concise dated revision or follow-up section.
4. If a material choice changes, create a new ADR and explicitly mark the previous ADR as `Superseded` or `Superseded in part`, with links in both directions.
5. Create a new ADR only when the decision has an independently reviewable boundary, alternatives, tradeoffs, and consequences.
6. Do not leave two `Accepted` ADRs claiming authority over the same decision scope. Split the scopes explicitly or resolve the relationship through supersession.

ADR numbers are historical identifiers and are never reused. Gaps in numbering are expected.

## Current decision families

Use this map to find current rationale without treating the chronological ADR list as a second specification.

| Concern | Current decision chain | Current production authority |
|---|---|---|
| biological claim boundary | ADR-0009, ADR-0019 | SRS + scientific-state invariants |
| circular/reference placement | ADR-0010, ADR-0029, ADR-0047 | alignment SRS/method + alignment invariants |
| variant calling / canonicalization / nomenclature | ADR-0047 → ADR-0057 → ADR-0060 | current variant SRS/method; platform direction in ADR-0060 |
| signal/locus evidence | ADR-0013, ADR-0028, ADR-0031, ADR-0046 | signal-processing SRS/method + evidence invariants |
| signal-derived read callability | ADR-0067 | SRS-CALL-* + callability method + evidence/pipeline invariants |
| sample boundary and sparse evidence | ADR-0023, ADR-0025 | sample SRS/method + sample contract |
| overlap/coverage/support topology | ADR-0030, ADR-0032, ADR-0033, ADR-0035, ADR-0068 | sample SRS/method + sample contract |
| sample call-signal projection | ADR-0037 (partially supersedes ADR-0034 and ADR-0036) | sample method + sample contract |
| nucleotide contribution/profile geometry | ADR-0038 through ADR-0043 | sample SRS/method; public subset in sample contract |
| reviewer-facing public signal evidence | ADR-0026, ADR-0053 | analysis/sample contracts |
| release/readiness | ADR-0018, ADR-0021 | release operations + roadmap/validation evidence |
| documentation governance | ADR-0006 → ADR-0022 | documentation governance + source-local README policy |
| modular analysis/public API boundaries | ADR-0058 → ADR-0069 | proposals + architecture/design/reference as implemented |
| plugin families and the modality → core evidence contract | ADR-0069 (supersedes ADR-0058 §4 and ADR-0060 §1/§6 in part) | PROP-0002 + interface architecture + variant-calling method |
| ecosystem reuse / dependency implementation policy | ADR-0059 | SRS-NFR + dependency policy + owning design |
| operational logging mechanism | ADR-0007 → ADR-0061 | SRS-OUT-008 + interface architecture |
| Sanger variant eligibility | ADR-0027, ADR-0062 → ADR-0067 | SRS-VAR-012/013 + variant-calling method |

The production authority column describes **current truth**. ADRs explain why that truth exists; they should not be copied into new production docs verbatim.

| ADR | Decision | Status |
|---|---|---|
| [0001](0001-mvp-vertical-slice.md) | End-to-end MVP first | Accepted |
| [0002](0002-single-crate-layering.md) | Single-crate layering | Accepted |
| [0003](0003-behavioral-compatibility.md) | Apollo behavioral evidence | Superseded in part by ADR-0009 |
| [0004](0004-rcrs-direct-alignment.md) | Direct rCRS alignment | Superseded by ADR-0010 |
| [0005](0005-versioned-output-contracts.md) | JSON plus VCF | Superseded by ADR-0008 |
| [0006](0006-source-documentation-mirroring.md) | Source/manual mirroring | Superseded by ADR-0022 |
| [0007](0007-configuration-and-environment.md) | Strict TOML and environment path | Accepted |
| [0008](0008-json-only-auditable-analysis.md) | JSON-only auditable output | Superseded in part by ADR-0011 |
| [0009](0009-biological-semantics.md) | Biologically explicit semantics | Accepted |
| [0010](0010-circular-rcrs-alignment.md) | Circular rCRS direct alignment | Accepted |
| [0011](0011-compact-variant-focused-json.md) | Compact variant-focused JSON | Superseded in part by ADR-0012 |
| [0012](0012-concise-mapped-variant-calls.md) | Concise mapped variant calls | Superseded in part by ADR-0013 |
| [0013](0013-observational-signal-quality.md) | Observational rolling signal quality | Accepted |
| [0014](0014-compact-result-summary.md) | Compact v5 result summary | Superseded by ADR-0026 |
| [0015](0015-reference-free-basecalling.md) | Reference-free basecall JSON | Accepted |
| [0016](0016-defer-ml-feature-boundary.md) | Defer ML feature boundary to separate training contract | Accepted |
| [0017](0017-primary-sample-peak-colocalization.md) | Gate secondary calls at the primary peak sample | Accepted |
| [0018](0018-production-readiness-release-contract.md) | Production readiness is an explicit release contract | Accepted |
| [0019](0019-scientific-evidence-hierarchy.md) | Separate signal evidence, read interpretation, and biological claims | Accepted |
| [0020](0020-rust-as-correctness-architecture.md) | Use Rust as correctness architecture, not only as an implementation language | Accepted |
| [0021](0021-scientific-core-confidence-floor.md) | Define a scientific core confidence floor | Accepted |
| [0022](0022-documentation-knowledge-system.md) | Govern documentation as an executable knowledge system | Accepted |
| [0023](0023-evidence-derived-read-placement-and-sample-boundary.md) | Derive read placement from evidence and reconcile samples from read observations | Accepted |
| [0024](0024-reference-coordinate-sample-evidence.md) | Aggregate independently placed reads in reference-coordinate/variant space | Superseded in part by ADR-0025 |
| [0025](0025-compact-sample-evidence.md) | Factor sample evidence into a read registry and sparse differences | Accepted |
| [0026](0026-reviewer-facing-signal-evidence.md) | Prefer reviewer-facing signal evidence over implementation call coordinates | Accepted |
| [0027](0027-mixed-supporting-signal-snv-eligibility.md) | Treat mixed supporting signal as observed evidence, not a clean SNV | Accepted |
| [0028](0028-basecall-independent-locus-evidence.md) | Derive locus evidence and profiles independently of basecall verdicts | Accepted |
| [0029](0029-profile-aware-gotoh.md) | Align reference placement from basecall-independent evidence profiles | Accepted |
| [0030](0030-tracy-derived-sample-overlap-admission.md) | Add Tracy-derived pairwise overlap admission before sample consensus | Accepted |
| [0031](0031-trace-integrity-evidence.md) | Preserve PLOC and signal-integrity evidence without artifact reclassification | Accepted |
| [0032](0032-sample-coverage-topology.md) | Expose sample coverage topology before consensus | Accepted |
| [0033](0033-variant-support-topology.md) | Factor normalized-variant support by eligibility and orientation | Accepted |
| [0034](0034-sample-evidence-profiles.md) | Preserve reference-oriented evidence profiles in sample evidence | Superseded in part by ADR-0037 |
| [0035](0035-locus-support-topology.md) | Factor differential-locus support topology before consensus | Accepted |
| [0036](0036-sample-local-noise-context.md) | Preserve local noisy-region context in sample evidence | Superseded in part by ADR-0037 |
| [0037](0037-sample-call-signal-evidence.md) | Unify reference-oriented call signal evidence at sample scope | Accepted |
| [0038](0038-locus-profile-availability.md) | Expose differential-locus profile availability before consensus | Accepted |
| [0039](0039-nucleotide-contribution-eligibility.md) | Define structural nucleotide contribution eligibility | Accepted |
| [0040](0040-unweighted-nucleotide-profile-support.md) | Accumulate eligible profiles with unit read mass | Accepted |
| [0041](0041-mean-nucleotide-evidence-profiles.md) | Derive mean nucleotide evidence profiles | Accepted |
| [0042](0042-profile-heterogeneity-decomposition.md) | Decompose nucleotide profile heterogeneity | Accepted |
| [0043](0043-directional-profile-distance.md) | Measure directional nucleotide-profile distance | Accepted |
| [0046](0046-nearest-locus-event.md) | Anchor nucleotide evidence to the nearest locus event | Accepted |
| [0047](0047-canonical-right-aligned-mtdna-gaps.md) | Canonical right-aligned mtDNA gap placement | Superseded in part by ADR-0060 |
| [0053](0053-public-differential-locus-signal-evidence.md) | Expose concise differential-locus signal evidence in sample output | Accepted |
| [0057](0057-haplotype-correctness-and-variant-nomenclature.md) | Separate haplotype correctness from variant nomenclature | Accepted |
| [0058](0058-canonical-contracts-and-modular-analysis-composition.md) | Canonical contracts and modular analysis composition | Superseded in part by ADR-0069 |
| [0059](0059-reuse-ecosystem-machinery-behind-dna-contracts.md) | Reuse ecosystem machinery behind DNA-owned contracts | Accepted |
| [0060](0060-separate-variant-canonicalization-nomenclature.md) | Separate variant calling, canonicalization, and nomenclature | Superseded in part by ADR-0069 |
| [0061](0061-tracing-for-operational-logging.md) | Use `tracing` for operational logging | Accepted |
| [0062](0062-read-callability.md) | Read callability from the read's own calls | Superseded in part by ADR-0067 |
| [0063](0063-target-profiles.md) | Target knowledge in versioned profiles | Accepted |
| [0064](0064-crate-ready-module-layering.md) | Crate-ready module layering | Accepted |
| [0065](0065-result-comparison-downstream.md) | Result comparison belongs to downstream pipelines | Accepted |
| [0066](0066-python-limited-to-repository-tooling.md) | Python in DNA is limited to repository tooling | Accepted |
| [0067](0067-signal-derived-read-callability.md) | Signal-derived read callability | Accepted |
| [0068](0068-variant-opposition-evidence.md) | Variant opposition evidence | Accepted |
| [0069](0069-plugin-first-modality-core-post-calling.md) | Plugin-first composition of modality, core, and post-calling plugins | Accepted |
