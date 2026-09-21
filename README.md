# DNA

DNA is a deterministic Rust system for DNA ABIF/AB1 analysis.

It reads analyzed A/C/G/T chromatogram channels, re-calls bases at validated ABIF PLOC loci, performs read-quality handling, aligns the retained read to a short reference in either orientation, and reports auditable primary-sequence differences through versioned JSON contracts.

DNA is designed as scientific software rather than as a generic sequence-conversion utility. The goal is not maximum feature count. The goal is to make signal processing, biological interpretation, and software correctness explicit enough to inspect, test, validate, and evolve safely.

## Design principles

DNA is developed through three complementary lenses:

- **signal processing** — preserve what the chromatogram actually measured and make derived transformations explicit;
- **biology** — report only claims supported by the available evidence and preserve unresolved states instead of guessing;
- **engineering** — encode important invariants in types, module boundaries, schemas, tests, and release gates wherever practical.

Rust is an architectural choice, not a branding choice. DNA uses Rust to move correctness from developer discipline into the programming model: validated states, explicit errors, immutable evidence, ownership boundaries, exhaustive state handling, and machine-checked invariants.

The compiler cannot prove biological correctness. Real scientific claims still require independent data and validation.

## Current status

The JSON-based single-trace pipeline is implemented and is being hardened toward an evidence-backed production release profile.

Current supported behavior includes:

- strict bounded ABIF/AB1 decoding;
- canonical analyzed A/C/G/T channels using ABIF channel-order metadata;
- signal-derived re-calling at validated `PLOC.2` loci;
- explicit primary and ambiguity states;
- observational trace-integrity and rolling signal-to-noise annotations;
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
AB1
 ↓
validated chromatogram decode
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
- `PLOC.2` is currently the locus authority for the re-calling method;
- rolling SNR and relative quality are not Phred-calibrated error probabilities;
- secondary or mixed signal is an observation, not automatically heteroplasmy;
- unresolved evidence remains unresolved;
- normalized variant representation must not erase the trace evidence from which it was observed.

See [ADR-0019](docs/adr/0019-scientific-evidence-hierarchy.md) and the [system invariants](docs/architecture/invariants.md).

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

DNA reads `DNA_CONFIG` or `config/dna.toml`.

Successful core commands publish exactly one command-specific JSON result without overwriting an existing result:

```text
basecall -> results/<trace-stem>.basecalls.json
analyze  -> results/<trace-stem>.json
sample   -> results/<sample-id>.sample.json
```

Operational logs are separate append-only sidecars under `logs/` by default. Standalone `basecall`/`analyze` operations use `<trace-stem>.log`; `sample` uses one `<sample-id>.log` containing the nested processing events for all traces in that sample. The batch runner persists only the sample log while keeping per-trace JSON results.

The external batch runner `scripts/analyze_samples.py` keeps per-trace results and
the aggregate together:

```text
results/<sample-id>/
├── <trace-stem>.json
├── ...
└── <sample-id>.json
```

The final `<sample-id>.json` is generated only when every selected trace for that
sample succeeds.

## Output contracts

Current public result contracts are:

- `dna.basecalls/v2` — reference-free primary/ambiguity/retained read result;
- `dna.analysis/v7` — compact reference-guided analysis result with reviewer-facing four-channel peak evidence;
- `dna.sample_evidence/v8` — compact multi-read coverage and overlap evidence plus sparse differential loci that preserve factorized support topology, per-read A/C/G/T evidence profiles/noisy context, normalized-variant evidence, and explicit eligibility reasons.

The schemas, examples, coordinate conventions, and human-readable semantics live under [docs/contracts](docs/contracts/README.md).

Public schemas are versioned contracts. Incompatible output changes require a new schema version rather than silent mutation of an existing version.

## Documentation

Start with [docs/README.md](docs/README.md).

The documentation system is organized by authority:

```text
SRS
 ↓
architecture + invariants
 ↓
ADRs
 ↓
current methods + public contracts
 ↓
docs/src implementation mirror
 ↓
source + tests
 ↓
validation + release evidence
```

Exploratory work lives under `docs/research/<topic>/` and is non-normative until promoted into the root production SRS/ADR/contract system.

Key entry points:

- [SRS](docs/srs/README.md)
- [architecture](docs/architecture/README.md)
- [system invariants](docs/architecture/invariants.md)
- [ADR index](docs/adr/README.md)
- [current methods](docs/methods/README.md)
- [contracts](docs/contracts/README.md)
- [source mirror](docs/src/README.md)
- [development and release operations](docs/operations/README.md)
- [traceability](docs/traceability.md)
- [research](docs/research/README.md)
- [roadmap](docs/roadmap.md)

## Development

Create the locked development environment:

```bash
uv sync --locked
```

Required repository checks:

```bash
uv run ruff format --check scripts/
uv run ruff check scripts/
uv run basedpyright scripts/
uv run python scripts/validate_result_schemas.py
uv run python scripts/validate_rust_source_policy.py

cargo fmt --all --check
cargo check --all-targets
cargo clippy --all-targets -- -D warnings
cargo test --all-targets
cargo doc --no-deps
```

Longer-running or release-oriented validation such as fuzzing, dependency audit, mutation testing, performance measurement, and approved real-AB1 regression belongs to the extended validation/release lanes rather than being added mechanically to every pull request.

See [CI and verification lanes](docs/operations/ci.md) and [production readiness](docs/adr/0018-production-readiness-release-contract.md).

## Agent development

Coding agents should start with [AGENTS.md](AGENTS.md). It defines a reusable discover → understand → plan → implement → verify → reconcile → review → report workflow, then expects the agent to discover this repository's own requirements, architecture, contracts, tests, and validation sources rather than relying on hard-coded file paths.

The repository intentionally treats documentation as part of the correctness system, not as an after-the-fact description of the code.
