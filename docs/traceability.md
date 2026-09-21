# Documentation and Implementation Traceability

This map helps developers and agents move from intent to implementation without duplicating the specification.

| Requirement family | Current method / architecture | Owning implementation | Primary tests / evidence | Public contract |
|---|---|---|---|---|
| `SRS-IN-*` | `pipeline.md`, `architecture.md` | `src/trace/`, `src/reference/`, `src/cli/` | parser/reference/CLI tests, synthetic ABIF | CLI/config contracts |
| `SRS-CFG-*` | `configuration.md` | `src/config/` | strict TOML/config tests | `config/dna.toml` semantics |
| `SRS-BC-*` | `pipeline.md` Stage 2 | `src/basecalling/`, `src/model/basecalls.rs` | basecalling unit tests + approved trace comparison | `dna.basecalls/v2`, mapped calls in analysis |
| `SRS-SIG-*` | `signal-processing.md`, ADR-0028 + ADR-0046 | `src/locus.rs`, `src/model/locus_evidence.rs`, `src/signal_processing/` | locus geometry/profile tests + nearest-event regression + signal feature boundary tests | internal evidence profile + public integrity/noisy-region projections |
| `SRS-QC-*` | `pipeline.md` Stage 4 | `src/quality_control/` | quality/trim unit + real-read review | trim + reviewer-facing `quality` |
| `SRS-ALN-*` | `pipeline.md` Stage 5, ADR-0029 + ADR-0047 + read-observation boundary | `src/alignment/canonical.rs`, `src/alignment/`, `src/model/locus_evidence.rs`, `src/model/alignment.rs`, `src/model/read_observation.rs` | fixed-point profile scoring, one-hot compatibility, unresolved-profile, orientation/traceback/circular tests; score-verified homopolymer/tandem-repeat right-gap and origin-seam fixtures | analysis alignment summary |
| `SRS-SAMPLE-*` | ADR-0023, ADR-0030, ADR-0053, architecture invariants | `src/sample/`, `src/model/sample_evidence.rs`, `src/model/sample_result.rs`, `src/report/sample.rs`, `src/pipeline/sample.rs` | deterministic coverage topology, overlap/admission, differential-locus support topology, public reference-oriented profile/noisy-context projection, structural nucleotide-contribution eligibility, unweighted eligible-profile support aggregation, arithmetic mean nucleotide-profile derivation, threshold-free profile heterogeneity decomposition, directional Total Variation geometry, normalized-variant support-topology unit tests, sparse sample evidence tests, end-to-end CLI/schema validation, local non-committed multi-read review | `dna.sample_evidence/v8` |
| `SRS-VAR-*` | `pipeline.md` Stage 6 | `src/variant_calling/`, `src/model/variant.rs` | SNV/indel/normalization mapping tests + mixed-supporting-signal eligibility + real truth | analysis variants |
| `SRS-OUT-*` | output docs, architecture | `src/report/`, `src/pipeline/` | schema/example + publication tests | JSON schemas |
| `SRS-BAT-*` | data/batch docs | `scripts/analyze_samples.py` | Python batch tests | external orchestration behavior |
| `SRS-NFR-*` | architecture/invariants, validation | cross-cutting | CI, fuzz/property/release evidence as adopted | release evidence |

## Navigation rule

For a behavior change:

1. locate the requirement family;
2. read the linked method/architecture document;
3. read relevant accepted ADRs;
4. read the owning `docs/src` manual and source;
5. inspect the linked tests/contracts;
6. update every affected layer in the same change.

Research documents are intentionally absent from this table because they are not current production authority.
