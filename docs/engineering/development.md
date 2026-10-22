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
branches in it rather than building a fresh checkout per probe. A commit
alone (the clean build a cached probe needs) only relinks the thin
`dna-eval` binary, which keeps the build stamps, and not the harness
library `dna-eval-lib`: 0.9 s instead of 7.0 s in the probe profile, 41 s
instead of 56 s in the release one, whose link-time optimisation dominates.
Approvals, readiness uses, and releases use `--release`.

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
different things, and exits 2 with `GATE FAIL` when, in the current period
(`--current-since`), more samples, passed samples, or passed zones carry a
non-het unreviewed block, or when the two runs failed on different current
samples (rerun a side that failed transiently); het blocks and the earlier
period are reported, not gated. Keep one cache
directory per machine; `dna-eval cache list --dir <dir>` shows its
generations by kind and last use, and
`dna-eval cache prune --dir <dir> --keep <n> [--dry-run]` keeps the `n` most
recently used of each kind. An evidence generation holds
every read run through it: about 3 GB for panel v5 alone, about 7 GB when
readiness-v2 runs through the same one. The placement cache shares the
directory, under 1 GB a generation.
Before any panel run, run the artefact bench on the probe
(`dna-eval bench --cases <cases.tsv> --manifest <panel.tsv>…`, with the same
cache): it runs DNA on the adjudicated artefact cases alone, in seconds, and
exits 2 with `BENCH FAIL` when a probe breaks one. Every signal change keeps
it passing; add a case for each newly adjudicated artefact, and promote a
known failure to a guard once it passes.
For a first pass, run the base and the probe on a quick subset
(`dna-eval subset`, [[validation/evaluation-harness]]); a subset verdict is a
pre-filter only, and a candidate must still pass `diff` on the full panel.

## Keeping a branch current

Merge `main` into a branch rather than rebasing a branch others may have
fetched. Two kinds of file used to conflict on almost every merge:

- `CHANGELOG.md` merges with git's built-in `union` driver
  (`.gitattributes`): every side's new lines are kept, so check the merged
  entries for order and duplicates before committing.
- Generated blocks (`<!-- docs:begin … -->`, such as the traceability
  coverage and the ADR index) are never resolved by hand. Take either side,
  then regenerate them; also run it after a clean merge, since a block
  generated on each side alone goes stale once both sides are combined:

  ```bash
  git merge origin/main
  (cd tools/python && uv run python scripts/docs_graph.py sync)
  ```

Wait for a pull request's checks with `gh pr checks <pr> --watch --fail-fast`.
It exits non-zero when a check fails. `CI success` is the aggregate of every
mandatory job ([[engineering/ci-cd#Aggregate required check]]). Once branch
rules mark it required ([[governance/repository]]), add `--required` to watch
only the required checks; `main` has no such rule yet, and with none
`--required` reports nothing.

## Required checks

Before finalizing a repository change, run the required gates defined in [[engineering/ci-cd|CI/CD]]. During iteration, run the narrowest checks that exercise the changed failure modes first.

## Documentation with code

When behavior, ownership, public contracts, architecture, or validation expectations change, update the affected canonical documentation in the same change. Source-directory responsibility changes also update the nearest `src/**/README.md`.

See [[governance/documentation|documentation governance]].
