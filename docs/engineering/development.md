---
decided-by: [ADR-0064, ADR-0066, ADR-0070, ADR-0075]
---

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
cargo test --profile local --workspace <focused-test-filter>
```

The `local` profile inherits development settings, keeps incremental
compilation enabled, and applies moderate optimization so CPU-heavy alignment is
much faster than an unoptimized dev binary without paying the full rebuild cost
of the release profile. Release validation still uses the repository's required
release build and CI gates.

For repeated local experiments, build once with `cargo build --profile local`
and run `target/local/dna` directly. Running DNA over a corpus is orchestrated by
downstream pipelines, not by tooling in this repository (ADR-0066).

## Evaluation probes

A probe runs `dna-eval` over a panel to measure a candidate change. Build it with
the `probe` profile, which inherits `release` without link-time optimisation and
with incremental compilation:

```bash
cargo build --profile probe -p dna-eval
target/probe/dna-eval run --manifest … --output …
```

Its outputs equal the release build's: on development panel v5 every table
`dna-eval run` writes was byte-identical, and the run took the same time with
about a quarter more CPU. After a first build, rebuilding for an edit to one
crate takes seconds instead of minutes (3 s against 124 s for an edit to
`dna-consensus`), so keep one probe worktree and its target directory and switch
branches in it rather than building a fresh checkout per probe. Approvals,
readiness uses, and releases use `--release`.

A probe that changes nothing of the read stage (the Sanger plugin, `dna-trace`,
the kernel, or the facade's evidence writer) can skip it with
`dna-eval run --evidence-cache <dir>`. Every read is then analysed from its
`dna.read_evidence/v1` document, read from the cache or made and stored there,
so the outputs equal a run from the traces. On a 410-sample subset of panel v5
every table was byte-identical from the traces, a cold cache, and a warm one,
and a warm run took a fifth of the CPU. The cache needs a harness built from
committed read-stage sources, since a plugin version does not move on a
work-in-progress commit; `run.json` records its key and hits. Compare a probe with its base run with
`dna-eval diff --base <run> --run <run>`: it refuses runs that measure
different things, and exits 2 with `GATE FAIL` when more passed samples or
zones carry a non-het unreviewed block, or when the two runs failed on
different samples (rerun a side that failed transiently). Keep one cache
directory per machine, and remove old generations with
`dna-eval cache prune --dir <dir> --keep <n>` (about 3 GB per generation for
panel v5).

## Required checks

Before finalizing a repository change, run the required gates defined in [[engineering/ci-cd|CI/CD]]. During iteration, run the narrowest checks that exercise the changed failure modes first.

## Documentation with code

When behavior, ownership, public contracts, architecture, or validation expectations change, update the affected canonical documentation in the same change. Source-directory responsibility changes also update the nearest `src/**/README.md`.

See [[governance/documentation|documentation governance]].
