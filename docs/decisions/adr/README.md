---
family-labels:
  - claim-boundary=biological claim boundary
  - placement=circular/reference placement
  - variant-nomenclature=variant calling / canonicalization / nomenclature
  - signal-evidence=signal/locus evidence
  - callability=signal-derived read callability
  - sample-boundary=sample boundary and sparse evidence
  - consensus=sample consensus (adjudication above sample evidence)
  - sample-topology=overlap/coverage/support topology
  - call-signal=sample call-signal projection
  - profile-geometry=nucleotide contribution/profile geometry
  - public-signal=reviewer-facing public signal evidence
  - release=release/readiness
  - documentation=documentation governance
  - modular-api=modular analysis/public API boundaries
  - crates=crate structure
  - comparison=comparison with reviewer truth and other tools
  - plugins=plugin families and the modality → core evidence contract
  - ecosystem-reuse=ecosystem reuse / dependency implementation policy
  - logging=operational logging mechanism
  - sanger-eligibility=Sanger variant eligibility
---

# Architecture Decision Records (ADR)

ADRs are historical decision records. Current-state truth lives in requirements, architecture, design, reference, source, and validation documentation.

## ADR lifecycle and deduplication

ADRs record durable architectural or scientific decisions, not every implementation increment, schema revision, or research step. The repository preserves historical ADRs, but each current decision scope must have one clear authority.

Before creating a new ADR:

1. Search the SRS, architecture/invariants, methods/contracts, and this index for an existing decision with the same scope.
2. If the work only implements, documents, validates, or clarifies an accepted decision, update the owning SRS/method/contract/tests instead of creating another ADR.
3. If the same decision is being refined without replacing its core choice, amend the existing ADR with a concise dated revision or follow-up section.
4. If a material choice changes, create a new ADR whose front matter `supersedes` names the previous record (with a section qualifier for a partial supersession, for example `ADR-0058 §4`) and set the previous record's `status` to `superseded` or `superseded-in-part`. `docs_graph.py sync` writes the headers and indexes.
5. Create a new ADR only when the decision has an independently reviewable boundary, alternatives, tradeoffs, and consequences.
6. Do not leave two `Accepted` ADRs claiming authority over the same decision scope. Split the scopes explicitly or resolve the relationship through supersession.

ADR numbers are historical identifiers and are never reused. A fully superseded record may be retired to the [ledger](../../governance/retired.md); the ledger also records numbers that were never committed.

## Current decision families

The map below is generated from each record's `family` and from the current
documents' `decided-by`; it shows where a concern's rationale lives and which
current documents it governs. It is not a second specification.

<!-- docs:begin families:decision -->
| Concern | Decisions († superseded) | Governs |
|---|---|---|
| biological claim boundary | [ADR-0009](0009-biological-semantics.md), [ADR-0019](0019-scientific-evidence-hierarchy.md) | [architecture/invariants/scientific-state](../../architecture/invariants/scientific-state.md), [proposals/roadmap](../../proposals/roadmap.md), [requirements/SRS](../../requirements/SRS.md) |
| circular/reference placement | [ADR-0010](0010-circular-rcrs-alignment.md), [ADR-0029](0029-profile-aware-gotoh.md), [ADR-0047](0047-canonical-right-aligned-mtdna-gaps.md)†, [ADR-0072](0072-run-structure-re-expression-of-masked-repeat-stretches.md), [ADR-0086](0086-unplaced-unreliable-tie.md), [ADR-0091](0091-anchor-kept-in-re-expressed-repeats.md) | [architecture/invariants/alignment](../../architecture/invariants/alignment.md), [architecture/variant-lifecycle](../../architecture/variant-lifecycle.md), [design/alignment](../../design/alignment.md), [design/signal-processing](../../design/signal-processing/README.md), [design/signal-processing/locus-evidence](../../design/signal-processing/locus-evidence.md), [design/variant-calling](../../design/variant-calling.md) +10 more |
| variant calling / canonicalization / nomenclature | [ADR-0047](0047-canonical-right-aligned-mtdna-gaps.md)†, [ADR-0057](0057-haplotype-correctness-and-variant-nomenclature.md), [ADR-0060](0060-separate-variant-canonicalization-nomenclature.md)† | [architecture/data-flow](../../architecture/data-flow.md), [architecture/invariants/alignment](../../architecture/invariants/alignment.md), [architecture/variant-lifecycle](../../architecture/variant-lifecycle.md), [design/alignment](../../design/alignment.md), [design/variant-calling](../../design/variant-calling.md), [design/variant-nomenclature](../../design/variant-nomenclature.md) +6 more |
| signal/locus evidence | [ADR-0013](0013-observational-signal-quality.md), [ADR-0028](0028-basecall-independent-locus-evidence.md), [ADR-0031](0031-trace-integrity-evidence.md), [ADR-0046](0046-nearest-locus-event.md), [ADR-0088](0088-electrophoretic-spikes.md) | [architecture/invariants/evidence](../../architecture/invariants/evidence.md), [design/abif-decoding](../../design/abif-decoding.md), [design/callability](../../design/callability.md), [design/signal-processing](../../design/signal-processing/README.md), [design/signal-processing/spikes](../../design/signal-processing/spikes.md), [design/variant-calling](../../design/variant-calling.md) +2 more |
| signal-derived read callability | [ADR-0067](0067-signal-derived-read-callability.md), [ADR-0076](0076-bounded-run-read-end-exemption.md), [ADR-0080](0080-reference-only-calls.md), [ADR-0082](0082-phase-mixtures.md), [ADR-0084](0084-trace-crate.md) | [architecture/invariants/evidence](../../architecture/invariants/evidence.md), [architecture/invariants/pipeline](../../architecture/invariants/pipeline.md), [architecture/overview](../../architecture/overview.md), [design/callability](../../design/callability.md), [design/consensus](../../design/consensus.md), [design/length-mixture](../../design/length-mixture.md) +26 more |
| sample boundary and sparse evidence | [ADR-0023](0023-evidence-derived-read-placement-and-sample-boundary.md), [ADR-0025](0025-compact-sample-evidence.md) | [design/sample-evidence](../../design/sample-evidence/README.md), [reference/sample-evidence](../../reference/sample-evidence/README.md), [requirements/sample-evidence](../../requirements/sample-evidence/README.md), [requirements/sample-evidence/callability](../../requirements/sample-evidence/callability.md), [requirements/sample-evidence/core](../../requirements/sample-evidence/core.md), [requirements/sample-evidence/overlap-topology](../../requirements/sample-evidence/overlap-topology.md) +2 more |
| sample consensus (adjudication above sample evidence) | [ADR-0073](0073-sample-consensus.md), [ADR-0074](0074-run-structure-consensus.md), [ADR-0083](0083-set-aside-dissent.md), [ADR-0085](0085-not-analysed-region-edge.md), [ADR-0089](0089-insertions-need-both-strands.md), [ADR-0090](0090-consensus-over-n-reads.md) | [architecture/overview](../../architecture/overview.md), [design/consensus](../../design/consensus.md), [design/length-mixture](../../design/length-mixture.md), [design/sample-evidence](../../design/sample-evidence/README.md), [design/sample-evidence/read-registry](../../design/sample-evidence/read-registry.md), [engineering/improvement-method](../../engineering/improvement-method.md) +11 more |
| overlap/coverage/support topology | [ADR-0030](0030-tracy-derived-sample-overlap-admission.md), [ADR-0032](0032-sample-coverage-topology.md), [ADR-0033](0033-variant-support-topology.md), [ADR-0035](0035-locus-support-topology.md), [ADR-0068](0068-variant-opposition-evidence.md) | [design/sample-evidence](../../design/sample-evidence/README.md), [design/sample-evidence/coverage-overlap](../../design/sample-evidence/coverage-overlap.md), [design/sample-evidence/variants](../../design/sample-evidence/variants.md), [reference/sample-evidence](../../reference/sample-evidence/README.md), [reference/sample-evidence/variants-signal](../../reference/sample-evidence/variants-signal.md), [requirements/sample-evidence/callability](../../requirements/sample-evidence/callability.md) +4 more |
| sample call-signal projection | [ADR-0037](0037-sample-call-signal-evidence.md) | [design/sample-evidence/call-signal](../../design/sample-evidence/call-signal.md), [reference/sample-evidence/variants-signal](../../reference/sample-evidence/variants-signal.md) |
| nucleotide contribution/profile geometry | [ADR-0038](0038-locus-profile-availability.md), [ADR-0039](0039-nucleotide-contribution-eligibility.md), [ADR-0040](0040-unweighted-nucleotide-profile-support.md), [ADR-0041](0041-mean-nucleotide-evidence-profiles.md), [ADR-0042](0042-profile-heterogeneity-decomposition.md), [ADR-0043](0043-directional-profile-distance.md) | [design/sample-evidence/locus-profiles](../../design/sample-evidence/locus-profiles.md), [requirements/sample-evidence/profile-aggregation](../../requirements/sample-evidence/profile-aggregation.md) |
| reviewer-facing public signal evidence | [ADR-0026](0026-reviewer-facing-signal-evidence.md), [ADR-0053](0053-public-differential-locus-signal-evidence.md) | [design/sample-evidence](../../design/sample-evidence/README.md), [governance/retired](../../governance/retired.md), [reference/analysis](../../reference/analysis/README.md), [reference/sample-evidence](../../reference/sample-evidence/README.md), [requirements/sample-evidence/callability](../../requirements/sample-evidence/callability.md), [requirements/sample-evidence/core](../../requirements/sample-evidence/core.md) +3 more |
| release/readiness | [ADR-0018](0018-production-readiness-release-contract.md), [ADR-0021](0021-scientific-core-confidence-floor.md) | [engineering/ci-cd](../../engineering/ci-cd.md), [engineering/release](../../engineering/release.md), [governance/versioning](../../governance/versioning.md), [operations/production-readiness](../../operations/production-readiness.md), [operations/release-evidence-template](../../operations/release-evidence-template.md), [proposals/roadmap](../../proposals/roadmap.md) +1 more |
| documentation governance | [ADR-0022](0022-documentation-knowledge-system.md)†, [ADR-0081](0081-typed-documentation-graph-and-retention.md) | [governance/documentation](../../governance/documentation.md), [governance/retired](../../governance/retired.md) |
| modular analysis/public API boundaries | [ADR-0058](0058-canonical-contracts-and-modular-analysis-composition.md)†, [ADR-0069](0069-plugin-first-modality-core-post-calling.md), [ADR-0070](0070-workspace-split-by-plugin-family.md) | [architecture/data-flow](../../architecture/data-flow.md), [architecture/interfaces](../../architecture/interfaces.md), [architecture/overview](../../architecture/overview.md), [design/alignment](../../design/alignment.md), [design/callability](../../design/callability.md), [design/pipeline](../../design/pipeline.md) +20 more |
| crate structure | [ADR-0064](0064-crate-ready-module-layering.md)†, [ADR-0070](0070-workspace-split-by-plugin-family.md), [ADR-0075](0075-evaluation-harness-in-workspace.md) | [architecture/context](../../architecture/context.md), [architecture/interfaces](../../architecture/interfaces.md), [architecture/overview](../../architecture/overview.md), [engineering/ci-cd](../../engineering/ci-cd.md), [engineering/dependencies](../../engineering/dependencies.md), [engineering/development](../../engineering/development.md) +14 more |
| comparison with reviewer truth and other tools | [ADR-0065](0065-result-comparison-downstream.md)†, [ADR-0075](0075-evaluation-harness-in-workspace.md) | [architecture/context](../../architecture/context.md), [architecture/interfaces](../../architecture/interfaces.md), [architecture/overview](../../architecture/overview.md), [design/callability](../../design/callability.md), [design/variant-nomenclature](../../design/variant-nomenclature.md), [engineering/ci-cd](../../engineering/ci-cd.md) +12 more |
| plugin families and the modality → core evidence contract | [ADR-0069](0069-plugin-first-modality-core-post-calling.md) | [architecture/data-flow](../../architecture/data-flow.md), [architecture/interfaces](../../architecture/interfaces.md), [architecture/overview](../../architecture/overview.md), [design/alignment](../../design/alignment.md), [design/callability](../../design/callability.md), [design/pipeline](../../design/pipeline.md) +15 more |
| ecosystem reuse / dependency implementation policy | [ADR-0059](0059-reuse-ecosystem-machinery-behind-dna-contracts.md) | [architecture/interfaces](../../architecture/interfaces.md), [architecture/overview](../../architecture/overview.md), [design/abif-decoding](../../design/abif-decoding.md), [design/alignment](../../design/alignment.md), [engineering/dependencies](../../engineering/dependencies.md), [requirements/quality-attributes](../../requirements/quality-attributes.md) |
| operational logging mechanism | [ADR-0007](0007-configuration-and-environment.md), [ADR-0061](0061-tracing-for-operational-logging.md) | [architecture/interfaces](../../architecture/interfaces.md), [engineering/dependencies](../../engineering/dependencies.md), [requirements/output](../../requirements/output.md) |
| Sanger variant eligibility | [ADR-0027](0027-mixed-supporting-signal-snv-eligibility.md), [ADR-0062](0062-read-callability.md)†, [ADR-0067](0067-signal-derived-read-callability.md), [ADR-0071](0071-run-length-edits-need-resolved-run-ends.md), [ADR-0080](0080-reference-only-calls.md), [ADR-0088](0088-electrophoretic-spikes.md) | [architecture/invariants/evidence](../../architecture/invariants/evidence.md), [architecture/invariants/pipeline](../../architecture/invariants/pipeline.md), [design/callability](../../design/callability.md), [design/consensus](../../design/consensus.md), [design/second-peak-cause](../../design/second-peak-cause.md), [design/signal-processing](../../design/signal-processing/README.md) +18 more |
<!-- docs:end families:decision -->

## All decisions

<!-- docs:begin index:decision -->
| ADR | Decision | Status |
|---|---|---|
| [ADR-0001](../../governance/retired.md) | Deliver an End-to-End MVP First | Retired |
| [ADR-0002](../../governance/retired.md) | Use a Single Layered Crate for the MVP | Retired, superseded by ADR-0070 |
| [ADR-0003](0003-behavioral-compatibility.md) | Target Behavioral Compatibility with Apollo C++ | Superseded in part by [ADR-0009](0009-biological-semantics.md) |
| [ADR-0004](../../governance/retired.md) | Align Directly to Bundled rCRS in the MVP | Retired, superseded by ADR-0010 |
| [ADR-0005](../../governance/retired.md) | Use Versioned JSON and Standards-Compliant VCF | Retired, superseded by ADR-0008 |
| [ADR-0006](../../governance/retired.md) | Mirror Rust Source Documentation One-to-One | Retired, superseded by ADR-0022 |
| [ADR-0007](0007-configuration-and-environment.md) | Use Strict TOML Configuration and a Safe Environment Policy | Accepted |
| [ADR-0008](0008-json-only-auditable-analysis.md) | Emit a Single JSON-Only, Auditable Analysis Document | Superseded in part by [ADR-0011](0011-compact-variant-focused-json.md) |
| [ADR-0009](0009-biological-semantics.md) | Report Primary-Sequence Differences, Not Genotypes or Heteroplasmy | Accepted |
| [ADR-0010](0010-circular-rcrs-alignment.md) | Align to the Circular rCRS Reference with Origin Wrapping | Accepted |
| [ADR-0011](0011-compact-variant-focused-json.md) | Emit Compact Variant-Focused JSON | Superseded in part by [ADR-0012](0012-concise-mapped-variant-calls.md) |
| [ADR-0012](0012-concise-mapped-variant-calls.md) | Use Concise Coordinates and Direct Mapped Variant Calls | Superseded in part by [ADR-0013](0013-observational-signal-quality.md) |
| [ADR-0013](0013-observational-signal-quality.md) | Add Observational DNA-Quality Windows Before DNA Cleaning | Accepted |
| [ADR-0014](../../governance/retired.md) | Project a Compact Analysis Result Summary | Retired, superseded by ADR-0026 |
| [ADR-0015](0015-reference-free-basecalling.md) | Add One Reference-Free Basecall JSON Contract | Accepted |
| [ADR-0016](0016-defer-ml-feature-boundary.md) | Defer ML Feature Boundary to a Separate Versioned Training Contract | Accepted |
| [ADR-0017](0017-primary-sample-peak-colocalization.md) | Gate Secondary Calls at the Primary Peak Sample | Accepted |
| [ADR-0018](0018-production-readiness-release-contract.md) | Production readiness is an explicit release contract | Accepted |
| [ADR-0019](0019-scientific-evidence-hierarchy.md) | Separate signal evidence, read interpretation, and biological claims | Accepted |
| [ADR-0020](0020-rust-as-correctness-architecture.md) | Use Rust as correctness architecture, not only as an implementation language | Accepted |
| [ADR-0021](0021-scientific-core-confidence-floor.md) | Define a scientific core confidence floor | Accepted |
| [ADR-0022](0022-documentation-knowledge-system.md) | Govern documentation as an executable knowledge system | Superseded in part by [ADR-0081](0081-typed-documentation-graph-and-retention.md) (§4/§6/§8) |
| [ADR-0023](0023-evidence-derived-read-placement-and-sample-boundary.md) | Evidence-derived read placement and sample boundary | Accepted |
| [ADR-0024](0024-reference-coordinate-sample-evidence.md) | Reference-coordinate sample evidence contract | Superseded in part by [ADR-0025](0025-compact-sample-evidence.md) |
| [ADR-0025](0025-compact-sample-evidence.md) | Factor sample evidence into a read registry and sparse differences | Accepted |
| [ADR-0026](0026-reviewer-facing-signal-evidence.md) | Prefer reviewer-facing signal evidence over implementation call coordinates | Accepted |
| [ADR-0027](0027-mixed-supporting-signal-snv-eligibility.md) | Mixed supporting signal is not a clean SNV | Accepted |
| [ADR-0028](0028-basecall-independent-locus-evidence.md) | Locus evidence and evidence profiles are basecall-independent | Accepted |
| [ADR-0029](0029-profile-aware-gotoh.md) | Align reference placement from basecall-independent evidence profiles | Accepted |
| [ADR-0030](0030-tracy-derived-sample-overlap-admission.md) | Add Tracy-derived pairwise overlap admission before sample consensus | Accepted |
| [ADR-0031](0031-trace-integrity-evidence.md) | Preserve PLOC and signal-integrity evidence without artifact reclassification | Accepted |
| [ADR-0032](0032-sample-coverage-topology.md) | Expose sample coverage topology before consensus | Accepted |
| [ADR-0033](0033-variant-support-topology.md) | Factor normalized-variant support by eligibility and orientation | Accepted |
| [ADR-0034](0034-sample-evidence-profiles.md) | Preserve reference-oriented evidence profiles in sample evidence | Superseded in part by [ADR-0037](0037-sample-call-signal-evidence.md) |
| [ADR-0035](0035-locus-support-topology.md) | Factor differential-locus support topology before consensus | Accepted |
| [ADR-0036](0036-sample-local-noise-context.md) | Preserve local noisy-region context in sample evidence | Superseded in part by [ADR-0037](0037-sample-call-signal-evidence.md) |
| [ADR-0037](0037-sample-call-signal-evidence.md) | Unify reference-oriented call signal evidence at sample scope | Accepted |
| [ADR-0038](0038-locus-profile-availability.md) | Expose differential-locus profile availability before consensus | Accepted |
| [ADR-0039](0039-nucleotide-contribution-eligibility.md) | Define structural nucleotide contribution eligibility | Accepted |
| [ADR-0040](0040-unweighted-nucleotide-profile-support.md) | Accumulate eligible profiles with unit read mass | Accepted |
| [ADR-0041](0041-mean-nucleotide-evidence-profiles.md) | Derive mean nucleotide evidence profiles | Accepted |
| [ADR-0042](0042-profile-heterogeneity-decomposition.md) | Decompose nucleotide profile heterogeneity | Accepted |
| [ADR-0043](0043-directional-profile-distance.md) | Measure directional nucleotide-profile distance | Accepted |
| [ADR-0044](../../governance/retired.md) | (number skipped) | Retired |
| [ADR-0045](../../governance/retired.md) | (number skipped) | Retired |
| [ADR-0046](0046-nearest-locus-event.md) | Anchor nucleotide evidence to the nearest locus event | Accepted |
| [ADR-0047](0047-canonical-right-aligned-mtdna-gaps.md) | Canonical right-aligned mtDNA gap placement | Superseded in part by [ADR-0060](0060-separate-variant-canonicalization-nomenclature.md) |
| [ADR-0048](../../governance/retired.md) | (number skipped) | Retired |
| [ADR-0049](../../governance/retired.md) | (number skipped) | Retired |
| [ADR-0050](../../governance/retired.md) | (number skipped) | Retired |
| [ADR-0051](../../governance/retired.md) | (number skipped) | Retired |
| [ADR-0052](../../governance/retired.md) | (number skipped) | Retired |
| [ADR-0053](0053-public-differential-locus-signal-evidence.md) | Expose concise differential-locus signal evidence in sample output | Accepted |
| [ADR-0054](../../governance/retired.md) | (number skipped) | Retired |
| [ADR-0055](../../governance/retired.md) | (number skipped) | Retired |
| [ADR-0056](../../governance/retired.md) | (number skipped) | Retired |
| [ADR-0057](0057-haplotype-correctness-and-variant-nomenclature.md) | Separate haplotype correctness from variant nomenclature | Accepted |
| [ADR-0058](0058-canonical-contracts-and-modular-analysis-composition.md) | Canonical contracts and modular analysis composition | Superseded in part by [ADR-0069](0069-plugin-first-modality-core-post-calling.md) (§4), [ADR-0070](0070-workspace-split-by-plugin-family.md) (§9) |
| [ADR-0059](0059-reuse-ecosystem-machinery-behind-dna-contracts.md) | Reuse ecosystem machinery behind DNA-owned contracts | Accepted |
| [ADR-0060](0060-separate-variant-canonicalization-nomenclature.md) | Separate variant calling, canonicalization, and nomenclature | Superseded in part by [ADR-0069](0069-plugin-first-modality-core-post-calling.md) (§1, §6) |
| [ADR-0061](0061-tracing-for-operational-logging.md) | Use `tracing` for operational logging | Accepted |
| [ADR-0062](0062-read-callability.md) | Read callability from the read's own calls | Superseded in part by [ADR-0067](0067-signal-derived-read-callability.md) |
| [ADR-0063](0063-target-profiles.md) | Target knowledge in versioned profiles | Accepted |
| [ADR-0064](0064-crate-ready-module-layering.md) | Crate-ready module layering | Superseded in part by [ADR-0070](0070-workspace-split-by-plugin-family.md) |
| [ADR-0065](0065-result-comparison-downstream.md) | Result comparison belongs to downstream pipelines | Superseded in part by [ADR-0075](0075-evaluation-harness-in-workspace.md) |
| [ADR-0066](0066-python-limited-to-repository-tooling.md) | Python in DNA is limited to repository tooling | Accepted |
| [ADR-0067](0067-signal-derived-read-callability.md) | Signal-derived read callability | Accepted |
| [ADR-0068](0068-variant-opposition-evidence.md) | Variant opposition evidence | Accepted |
| [ADR-0069](0069-plugin-first-modality-core-post-calling.md) | Plugin-first composition of modality, core, and post-calling plugins | Accepted |
| [ADR-0070](0070-workspace-split-by-plugin-family.md) | Workspace split by plugin family | Accepted |
| [ADR-0071](0071-run-length-edits-need-resolved-run-ends.md) | Run-length edits need resolved run ends | Accepted |
| [ADR-0072](0072-run-structure-re-expression-of-masked-repeat-stretches.md) | Run-structure re-expression of masked repeat stretches | Accepted |
| [ADR-0073](0073-sample-consensus.md) | Sample consensus as the adjudication layer | Accepted |
| [ADR-0074](0074-run-structure-consensus.md) | Run-structure consensus and how a run length is known | Accepted |
| [ADR-0075](0075-evaluation-harness-in-workspace.md) | Evaluation harness as a workspace crate | Accepted |
| [ADR-0076](0076-bounded-run-read-end-exemption.md) | Bounded runs before a post-repeat phase loss are exempt from `read_end` | Accepted |
| [ADR-0077](0077-normalize-and-empop-crates.md) | Post-calling split into normalize and EMPOP crates | Accepted |
| [ADR-0078](0078-call-the-consensus-decisions.md) | Call the consensus's decided edits without placing it again | Accepted |
| [ADR-0079](0079-core-crates-by-step.md) | One crate per core step | Accepted |
| [ADR-0080](0080-reference-only-calls.md) | A modality may assert a call confirms only the reference | Accepted |
| [ADR-0081](0081-typed-documentation-graph-and-retention.md) | Document as a typed graph with generated indexes and ledgered retention | Accepted |
| [ADR-0082](0082-phase-mixtures.md) | A modality may assert a phase mixture, and the core compares its components | Accepted |
| [ADR-0083](0083-set-aside-dissent.md) | Set aside dissent from a reference decision that the other strand contradicts | Accepted |
| [ADR-0084](0084-trace-crate.md) | Trace signal engineering is a library crate the modality calls | Proposed |
| [ADR-0085](0085-not-analysed-region-edge.md) | The core decides the analysed range at a region edge from the reads' unreliable spans | Accepted |
| [ADR-0086](0086-unplaced-unreliable-tie.md) | A placement tie confined to a read's unreliable end leaves those calls unplaced | Accepted |
| [ADR-0087](0087-placed-read-contract.md) | Placed reads as a step contract that dna-align writes and reads | Accepted |
| [ADR-0088](0088-electrophoretic-spikes.md) | A call inside an electrophoretic spike shows no base | Accepted |
| [ADR-0089](0089-insertions-need-both-strands.md) | Insertions need both strands | Accepted |
| [ADR-0090](0090-consensus-over-n-reads.md) | Consensus over N reads per site and strand | Accepted |
| [ADR-0091](0091-anchor-kept-in-re-expressed-repeats.md) | A re-expressed repeat stretch keeps the read's own anchor base | Accepted |
<!-- docs:end index:decision -->
