# Development

## Environment

Create the locked Python tooling environment:

```bash
uv sync --locked
```

Rust dependencies are locked by `Cargo.lock`. Use the repository toolchain and configuration rather than ad-hoc local substitutions.

## Required checks

Before finalizing a repository change, run the required gates defined in [CI/CD](ci-cd.md). During iteration, run the narrowest checks that exercise the changed failure modes first.

## Documentation with code

When behavior, ownership, public contracts, architecture, or validation expectations change, update the affected canonical documentation in the same change. Source-directory responsibility changes also update the nearest `src/**/README.md`.

See [documentation governance](../governance/documentation.md).
