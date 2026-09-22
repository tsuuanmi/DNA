# Documentation and Implementation Traceability

This map helps developers and agents move from intent to implementation without duplicating the specification.

| Requirement family | Current method / architecture | Owning implementation | Primary tests / evidence | Public contract |
|---|---|---|---|---|
| [`SRS-IN-*`](../requirements/input.md) | [`design/pipeline.md`](../design/pipeline.md), [`architecture/system.md`](../architecture/overview.md) | `src/trace/`, `src/reference/`, `src/cli/` | parser/reference/CLI tests, synthetic ABIF | CLI/config contracts |
| [`SRS-CFG-*`](../requirements/configuration.md) | [`reference/configuration.md`](../reference/configuration.md) | `src/config/` | strict TOML/config tests | `config/dna.toml` semantics |
| [`SRS-BC-*`](../requirements/basecalling.md) | [`design/basecalling.md`](../design/basecalling.md) | `src/basecalling/`, `src/model/basecalls.rs` | basecalling unit tests + approved trace comparison | `dna.basecalls/v2`, mapped calls in analysis |
| [`SRS-SIG-*`](../requirements/signal-processing.md) | [`design/signal-processing/`](../design/signal-processing/README.md), ADR-0028 + ADR-0046 | `src/locus.rs`, `src/model/locus_evidence.rs`, `src/signal_processing/` | locus geometry/profile tests + nearest-event regression + signal feature boundary tests | internal evidence profile + public integrity/noisy-region projections |
| [`SRS-QC-*`](../requirements/quality-control.md) | [`design/quality-control.md`](../design/quality-control.md) | `src/quality_control/` | quality/trim unit + real-read review | trim + reviewer-facing `quality` |
| [`SRS-ALN-*`](../requirements/alignment.md) | [`design/alignment.md`](../design/alignment.md), ADR-0029 + ADR-0047 | `src/alignment/canonical.rs`, `src/alignment/`, `src/model/locus_evidence.rs`, `src/model/alignment.rs`, `src/model/read_observation.rs` | fixed-point profile scoring, one-hot compatibility, unresolved-profile, orientation/traceback/circular tests; score-verified homopolymer/tandem-repeat right-gap and origin-seam fixtures | analysis alignment summary |
| [`SRS-SAMPLE-*`](../requirements/sample-evidence/README.md) | [`design/sample-evidence/`](../design/sample-evidence/README.md), ADR-0023, ADR-0030, ADR-0053 | `src/sample/`, `src/model/sample_evidence.rs`, `src/model/sample_result.rs`, `src/report/sample.rs`, `src/pipeline/sample.rs` | deterministic coverage topology, overlap/admission, differential-locus support topology, public reference-oriented profile/noisy-context projection, structural nucleotide-contribution eligibility, unweighted eligible-profile support aggregation, arithmetic mean nucleotide-profile derivation, threshold-free profile heterogeneity decomposition, directional Total Variation geometry, normalized-variant support-topology unit tests, sparse sample evidence tests, end-to-end CLI/schema validation, local non-committed multi-read review | `dna.sample_evidence/v8` |
| [`SRS-VAR-*`](../requirements/variants.md) | [`design/variant-calling.md`](../design/variant-calling.md) | `src/variant_calling/`, `src/model/variant.rs` | SNV/indel/normalization mapping tests + mixed-supporting-signal eligibility + real truth | analysis variants |
| [`SRS-OUT-*`](../requirements/output.md) | [contracts](../reference/README.md), architecture | `src/report/`, `src/pipeline/` | schema/example + publication tests | JSON schemas |
| [`SRS-BAT-*`](../requirements/batch.md) | [batch operations](../operations/runbooks/batch-analysis.md) | `scripts/analyze_samples.py` | Python batch tests | external orchestration behavior |
| [`SRS-COMPAT-*`](../requirements/constraints.md) | [`governance/reference-validation.md`](../governance/reference-validation.md) | cross-cutting | approved differential evidence | reference-validation policy |
| [`SRS-NFR-*`](../requirements/quality-attributes.md) | [`architecture/invariants/`](../architecture/invariants/README.md), validation | cross-cutting | CI, fuzz/property/release evidence as adopted | release evidence |
| [`SRS-VAL-*`](acceptance-criteria.md) | validation and operations docs | cross-cutting | required repository and release gates | release evidence |

## Navigation rule

For a behavior change:

1. locate the requirement family;
2. read the linked design/architecture document;
3. read relevant accepted ADRs;
4. read the nearest owning source-directory `README.md` and source;
5. inspect the linked tests/reference contracts;
6. update every affected layer in the same change.

Research documents are intentionally absent from this table because they are not current production authority.
