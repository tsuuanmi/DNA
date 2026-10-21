---
decided-by: [ADR-0075]
covers: [crates/dna-eval/src/bench.rs]
---

# Artefact bench

Part of the [[validation/evaluation-harness|evaluation harness]]: a
regression bench of adjudicated cases, run before any panel run.

```bash
cargo run --release -p dna-eval -- bench \
  --cases <cases.tsv> --manifest <panel.tsv> [--manifest <panel.tsv>...] \
  --config config/dna.toml --reference references/rCRS.fasta \
  --threads <1-256> [--evidence-cache <dir>] [--output <new-directory>]
```

## What it checks

A case is a Sanger artefact or a truth error whose expected outcome was
settled from the raw signal: an extra peak, an electrophoretic spike, a dye
blob, a truth range cut short. `bench` runs DNA once on each sample the
cases name, from the traces its panel manifest lists (or their cached
evidence with `--evidence-cache`), and checks each case's expectation. It
reads no truth: the expectation is the adjudicated outcome. A signal change
that breaks known-good behaviour on a case therefore shows in seconds,
before a panel run, and two rules that interact badly show on the case they
break. It does not replace `diff`, which still gates a candidate on the
full development panel.

## Case list

The case list is a TSV with these columns:

| Column | Content |
|---|---|
| `case` | a unique name |
| `sample_id` | the sample |
| `panel` | `<panel>-<version>` of a manifest header, for example `dev-v5`; every panel a case names needs its `--manifest` |
| `reads` | the reads that show the artefact; descriptive, since DNA runs on every trace of the sample |
| `class` | the artefact or truth-error class, which the report counts by |
| `first`, `last` | the span, 1-based and inclusive, at most 1,000 positions |
| `expect` | the outcome DNA must keep over the span |
| `status` | `guard`, or `known` for a known failure |
| `source` | where the case was adjudicated |

`expect` holds alternatives separated by `|`, each a list of terms
separated by `;` that must all hold:

| Term | Holds when |
|---|---|
| `alleles=<a>,<b>` | the alleles anchored in the span are exactly these, in reviewer notation (`573.1C`, `249DEL`); `alleles=` means none |
| `absent=<a>,<b>` | none of these alleles is called anywhere in the sample |
| `no_indel` | no insertion or deletion is anchored in the span |
| `check` / `check=<label>` | a check item of `dna review` overlaps the span (carrying the label) |
| `called` / `uncalled` | every / no position of the span lies in DNA's callable intervals |

For example, an extra peak expects `no_indel`; a spike, `uncalled |
alleles=;called` (withheld, or the reference called); a dye blob,
`check=reference_second | alleles=;called`; a truth error where DNA is
right, `alleles=573.1C,573.2C`.

## Verdicts

- `PASS`: a guarded case holds.
- `FAIL`: a guarded case fails, a regression.
- `KNOWN` / `FIXED`: a known failure still fails, or now passes; a fixed
  case should become a guard.
- `ERROR`: DNA could not analyse the sample, so its cases cannot be judged.

The report (stdout, and `bench.md` with `--output`) lists every case that
is not a plain pass with what DNA showed over its span (alleles, overlapping
checks, callable positions), then the counts per class; `bench.tsv` holds
every case. The verdict is `BENCH FAIL`, with exit status 2, when a guarded
case fails or any case is `ERROR`; otherwise `BENCH PASS`.

## Data

The case list names samples, so it is data: it stays with the panels under
`results/` and is never committed ([[governance/data]]). The repository holds
only the runner and this format, and its tests use synthetic outcomes. A
committed fixture such as each case's read evidence was not chosen: it is
sample-derived sequence, which the data policy keeps out of the repository.
Add a case for every newly adjudicated artefact, and keep its source with
it.
