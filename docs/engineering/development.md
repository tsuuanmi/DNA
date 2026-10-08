# Development

## Environment

Create the locked Python tooling environment:

```bash
uv sync --project tools/python --locked
```

Rust dependencies are locked by `Cargo.lock`. Use the repository toolchain and configuration rather than ad-hoc local substitutions.

## Fast local Rust iteration

Use the dedicated `local` Cargo profile when repeatedly exercising real Sanger
workloads during development:

```bash
cargo build --profile local
cargo run --profile local -- analyze sample.ab1 --reference references/rCRS.fasta
cargo test --profile local <focused-test-filter>
```

The `local` profile inherits development settings, keeps incremental
compilation enabled, and applies moderate optimization so CPU-heavy alignment is
much faster than an unoptimized dev binary without paying the full rebuild cost
of the release profile. Release validation still uses the repository's required
release build and CI gates.

For repeated local experiments, build once with `cargo build --profile local`
and run `target/local/dna` directly. Running DNA over a corpus is orchestrated by
downstream pipelines, not by tooling in this repository (ADR-0066).

## Required checks

Before finalizing a repository change, run the required gates defined in [CI/CD](ci-cd.md). During iteration, run the narrowest checks that exercise the changed failure modes first.

## Documentation with code

When behavior, ownership, public contracts, architecture, or validation expectations change, update the affected canonical documentation in the same change. Source-directory responsibility changes also update the nearest `src/**/README.md`.

See [documentation governance](../governance/documentation.md).
