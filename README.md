# DNA

DNA is a deterministic Rust platform for DNA analysis. The current production implementation analyzes Sanger sequencing traces stored as ABIF, but the platform architecture is not restricted to Sanger or mitochondrial DNA.

The current Sanger capability reads analyzed A/C/G/T chromatogram channels, projects ABIF locus metadata into canonical Sanger loci, re-calls bases from signal, performs read-quality handling, aligns the retained read to a reference, and reports auditable primary-sequence differences through typed Rust results and versioned JSON contracts.

DNA is designed as scientific software rather than as a generic sequence-conversion utility. The goal is not maximum feature count. The goal is to make signal processing, biological interpretation, and software correctness explicit enough to inspect, test, validate, and evolve safely.

## Design principles

DNA is developed through three complementary lenses:

- **signal processing** — preserve what the chromatogram actually measured and make derived transformations explicit;
- **biology** — report only claims supported by the available evidence and preserve unresolved states instead of guessing;
- **engineering** — encode important invariants in types, module boundaries, schemas, tests, and release gates wherever practical.

Rust is an architectural choice, not a branding choice. DNA uses Rust to move correctness from developer discipline into the programming model: validated states, explicit errors, immutable evidence, ownership boundaries, exhaustive state handling, and machine-checked invariants.

The compiler cannot prove biological correctness. Real scientific claims still require independent data and validation.

DNA's platform scope treats biological target, sequencing modality, file format, and analysis capability as separate dimensions. Mitochondrial DNA is one biological target; nuclear/genomic DNA and targeted loci are future targets. SNP analysis/genotyping is an analysis use case or variant class, not a sequencing modality.

## Current status

The JSON-based pipeline is implemented with production-oriented engineering controls. A release is not called production-ready until the exact revision also satisfies the scientific validation and release-evidence contract in ADR-0018.

Current supported behavior includes:

- strict bounded Sanger ABIF decoding;
- canonical analyzed A/C/G/T channels using ABIF channel-order metadata;
- projection of `PLOC.2` into canonical Sanger locus positions at the input boundary;
- signal-derived re-calling at validated canonical Sanger loci;
- explicit primary and ambiguity states;
- observational Sanger-integrity and rolling signal-to-noise annotations;
- signal-derived read callability: phase-state segments, a typed mask, and the callable span, published as observation (ADR-0067);
- deterministic read-quality scoring and end trimming;
- forward/reverse profile-aware semi-global alignment to one short reference;
- linear and circular reference handling;
- primary-sequence SNVs and supported small insertions/deletions;
- reviewer-facing reference-oriented A/C/G/T peak and quality evidence for reported variants;
- run-length total/forward/reverse coverage topology, Tracy-derived pairwise overlap/admission evidence, and factorized normalized-variant support topology across independently placed sample reads;
- closed versioned JSON schemas;
- atomic no-overwrite result publication;
- typed failures and bounded resource use.

The core confidence floor is deliberately simpler than the full current implementation:

```text
Sanger ABIF
 ↓
validated Sanger evidence decode
 ↓
signal-derived base re-calling
 ↓
basic QC / trimming
 ↓
forward-or-reverse evidence-profile alignment
 ↓
primary-sequence variant calling
 ↓
versioned JSON
```

The confidence floor defines what must be understood and validated first. It is not a reason to remove known-good capabilities that already exceed it.

## Scientific interpretation

DNA distinguishes source evidence from interpretation.

```text
analyzed chromatogram channels
        ↓
per-locus observations
        ↓
read interpretation
        ↓
reference differences
```

A single chromatogram does **not** establish genotype, quantitative heteroplasmy, phase, contamination, pathogenicity, or clinical significance.

Important boundaries:

- `PBAS.2` / `PCON.2` are optional vendor evidence; they do not determine DNA's final call;
- `PLOC.2` is interpreted only by the ABIF format layer and projected to canonical Sanger locus positions before downstream scientific stages;
- rolling SNR and relative quality are not Phred-calibrated error probabilities;
- secondary or mixed signal is an observation, not automatically heteroplasmy;
- unresolved evidence remains unresolved;
- normalized variant representation must not erase the trace evidence from which it was observed.

See [ADR-0019](docs/decisions/adr/0019-scientific-evidence-hierarchy.md) and the [system invariants](docs/architecture/invariants/README.md).

## Quick start

Reference-free base calling:

```bash
cargo run --release -- basecall sample.ab1
```

Reference-guided analysis:

```bash
cargo run --release -- analyze sample.ab1 \
  --reference references/rCRS.fasta
```

Multi-read sample evidence:

```bash
cargo run --release -- sample AB0442 read1.ab1 read2.ab1 \
  --reference references/rCRS.fasta
```

DNA reads `DNA_CONFIG` or `config/dna.toml` and, for reference-guided commands, the [target profile](docs/reference/profiles.md) it names (`config/profiles/human-mtdna-rcrs.toml` by default).

Successful core commands publish exactly one command-specific JSON result without overwriting an existing result:

```text
basecall -> results/<trace-stem>.basecalls.json
analyze  -> results/<trace-stem>.json
sample   -> results/<sample-id>.sample.json
```

Operational logs are separate append-only sidecars under `logs/` by default. Standalone `basecall`/`analyze` operations use `<trace-stem>.log`; `sample` uses one `<sample-id>.log` containing the nested processing events for all traces in that sample.

Running DNA over many samples, converting its results, and comparing them with other sources is done by downstream pipelines that drive the CLI or library (ADR-0065, ADR-0066).

## Output contracts

Current public result contracts are:

- `dna.basecalls/v3` — reference-free primary/ambiguity/retained read result with the read callability view;
- `dna.analysis/v9` — compact reference-guided analysis result with reviewer-facing four-channel peak evidence and the read callability view;
- `dna.sample_evidence/v10` — compact multi-read coverage and overlap evidence plus sparse differential loci that preserve factorized support topology, per-read callability, A/C/G/T evidence profiles/noisy context, normalized-variant evidence, and explicit eligibility reasons.

The schemas, examples, coordinate conventions, and human-readable semantics live under [docs/reference](docs/reference/README.md).

Public schemas are versioned contracts. Incompatible output changes require a new schema version rather than silent mutation of an existing version.

## Documentation

Start with [docs/README.md](docs/README.md).

The documentation system is organized as a knowledge lifecycle:

```text
requirements -> research/proposal -> decision
            -> architecture/design/reference
            -> source + tests -> validation
            -> engineering/release -> operations
```

Current truth lives in requirements/architecture/design/reference and executable source. Decisions preserve why; proposals explore change; research provides evidence; validation proves behavior; operations keep the released system supportable.

Key entry points:

- [requirements / SRS](docs/requirements/SRS.md)
- [architecture](docs/architecture/README.md)
- [system invariants](docs/architecture/invariants/README.md)
- [decisions / ADRs](docs/decisions/README.md)
- [design](docs/design/README.md)
- [reference / contracts](docs/reference/README.md)
- [source modules](src/README.md)
- [engineering](docs/engineering/README.md)
- [traceability](docs/validation/traceability.md)
- [research](docs/research/README.md)
- [roadmap](docs/proposals/roadmap.md)
- [documentation architecture standard](docs/governance/documentation-architecture.md)

## Development

The repository root is a Rust project. Executable source under `src/` is Rust; source-local `README.md` files document module ownership and boundaries. Python is isolated under `tools/python/` and is used only for repository checks, tests, measurement, and research (ADR-0066); it never produces DNA results.

The release Rust toolchain is pinned by `rust-toolchain.toml`; `Cargo.toml` separately declares the minimum supported Rust version (MSRV). GitHub-hosted Linux verification and delivery jobs pin Ubuntu 24.04 rather than following the moving `ubuntu-latest` label.

Create the locked Python tooling environment only when those companion tools are needed:

```bash
uv sync --project tools/python --frozen
```

Required repository checks:

```bash
cargo fmt --all --check
cargo shear --deny-warnings
cargo check --locked --all-targets --all-features
cargo clippy --locked --all-targets --all-features -- -D warnings
cargo test --locked --all-targets --all-features
RUSTDOCFLAGS="-D warnings" cargo doc --locked --no-deps --all-features
cargo build --locked --release

cd tools/python
uv run ruff format --check scripts tests
uv run ruff check scripts tests
uv run basedpyright scripts tests
uv run python -m unittest discover -s tests -p 'test_*.py'
uv run python scripts/validate_result_schemas.py
uv run python scripts/validate_rust_source_policy.py
uv run python scripts/validate_docs_structure.py
```

CI additionally verifies GitHub Actions syntax/security, the declared MSRV, Rust-only production source, dependency/source hygiene, dependency policy/review, RustSec, CodeQL, schemas/reference data, an ABIF fuzz smoke campaign, and a release-package smoke. Mandatory CI jobs feed an aggregate `CI success` check for branch protection. Third-party Actions are pinned to immutable commits and Dependabot maintains those pins.

Tagged `v*` releases rerun required Rust/security gates, require the tagged commit to belong to `main`, build the explicit `x86_64-unknown-linux-gnu` target as an auditable Rust binary, preserve and verify embedded dependency metadata after stripping, bundle the authoritative config, target profiles, and rCRS reference with checksums, generate an SPDX SBOM and SHA-256 checksums, attest the verified artifacts, then publish the supported Linux artifact.

Longer scientific validation—approved real-AB1 ground-truth comparison (run by downstream pipelines on DNA's published results), extended fuzzing, and runtime/resource evidence—remains release evidence rather than being conflated with ordinary software CI.

See [CI/CD](docs/engineering/ci-cd.md), [repository governance](docs/governance/repository.md), [release evidence](docs/operations/release-evidence-template.md), [production readiness](docs/operations/production-readiness.md), and [ADR-0018](docs/decisions/adr/0018-production-readiness-release-contract.md).

## Agent development

Coding agents should start with [AGENTS.md](AGENTS.md). It is intentionally limited to routing and repository-wide invariants; project knowledge remains in canonical docs, source-directory README files, source, and tests.

The repository intentionally treats documentation as part of the correctness system, not as an after-the-fact description of the code.
