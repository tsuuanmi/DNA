---
implements: [SRS-API, SRS-CONS]
---

# DNA Run Context JSON

`dna run-context <run-id> --traces <list>|--evidence <list> --reference
<fasta>` writes one deterministic document,
`results/<run-id>.run_context.json`, identified as `dna.run_context/v1`
(SRS-API-021, SRS-CONS-014, [[PROP-0014]]). It records what a laboratory
batch's own samples say about the analysis of each of them: additive
statistics only, so a consumer takes its own sample's contribution out. The
unit is the batch: a sample's forward and reverse reads of one amplicon often
sit on different instrument runs, and the length statistic needs both
strands.

The authoritative schema is
[`schemas/run-context-v1.schema.json`](schemas/run-context-v1.schema.json);
the synthetic example is
[`examples/run-context-v1.example.json`](examples/run-context-v1.example.json).

## Command

The list holds one `sample_id<TAB>path` line per read of the batch: its
traces with `--traces`, or their `dna.read_evidence/v1` documents with
`--evidence`. Every sample's reads are processed as `dna consensus` processes
them, up to the two-strand rule's views. A read the sample commands set aside
or cannot place counts in `rejected_reads` and contributes nothing. A listed
file that is missing, not a regular file, empty, or unreadable is left out
and recorded in `provenance.skipped` with its sample, file name, and reason
(`missing`, `not_a_file`, `empty`, `unreadable`); the sample's other reads
still count. Every sample command fails on such a file, so a caller analyses
the sample without it, and all of the reads it analyses are in the
context. The
output does not depend on the order of the list. Build it over the batch's
full folder: a sample analysed later with the context must have all of its
reads in it, or none. The operation log is `logs/<run-id>.log`.

## Fields

- `schema_version`, `run_id`: the caller's batch identifier.
- `provenance`: the reference, the configuration's SHA-256, the profile, the
  plugins, and `skipped`, the listed files left out.
- `reads`: the content SHA-256 of every contributing read, sorted;
  `rejected_reads`; `samples`, the samples that contributed.
- `length` (owner `dna-consensus`, method 1): the profile's recorded anchored
  windows the pairs come from (`windows`, on rCRS `[[303, 315]]`), the τ the
  differences were normalised with (`noise_forward`, `noise_reverse`), the
  histogram's bin `edges` (`consensus.length_context_edges`), and per run of
  a window (run 1 before the anchor, run 2 after it):
  - `pairs`: cross-strand pairs. A sample gives a run one when each strand
    has exactly one view of it that is not in phase, both with the same top
    length, at least one with a next length, and both with a standard error
    (SRS-CONS-013 names the views);
  - `sum_z2`: the sum of `d² / (SE_F² + SE_R² + τ_F² + τ_R²)`, with
    `d = p_F − p_R`;
  - `sum_z4`: the sum of `z⁴`, so the variance of `z²`, and with it how far
    a batch's dispersion varies by sampling alone, is known (for the prior's
    fit and for control limits);
  - `histogram`: the pairs by normalised difference `d / √(…)`: below the
    first edge, between each two edges, and from the last edge on.

## Use

`dna consensus` and `dna review` take `--run-context <file>`
([[reference/consensus|consensus]], [[reference/review|review]]). The
consensus takes the sample's own pairs out when all its reads are among
`reads`, and keeps them when none is. It rescales τ of each run whose
remaining pairs reach `consensus.length_context_minimum_pairs` by
`max(m, √(φ̃ / φ_ref))`, at most ½, where `φ̃` is the run's mean `z²` shrunk
towards `φ_ref` by `consensus.length_context_prior_pairs` pseudo-pairs
(SRS-CONS-014). A run with fewer pairs keeps the configured τ. A context made
against another reference, over other windows, or with another τ or other
edges, or one that holds only some of the sample's reads, fails the command.
The document's provenance then records `run_context`: the run identifier,
the context's SHA-256, `left_out`, and `length_scale` per run (`null` where
the pairs were too few).

In process, `SampleAnalysis::with_run_context` applies a context that
`RunContext::read` loaded ([[reference/rust-api-sample|sample API]]);
`dna-eval run --run-context <dir>` reads `<dir>/<batch>.run_context.json`
for each sample's batch ([[validation/evaluation-harness|evaluation
harness]]).
