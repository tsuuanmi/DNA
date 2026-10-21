---
id: PROP-0014
status: accepted
created: 2026-10-21
summary: "a run-level context stage, `dna.run_context/v1`, computed once per laboratory batch from all of the batch's reads and read by the modality and the consensus: per-batch length noise for the anchored-window runs, then the run signal, read-type context shadows, and blob baselines as further sections. It holds additive statistics only, so a consumer leaves the sample out itself, and every consumer falls back to today's fixed levels without one."
implements: [SRS-API, SRS-CALL, SRS-CONS]
decided-by: [ADR-0065, ADR-0069, ADR-0087]
---

# Proposal: Run context

## Decision

Accepted by the user on 2026-10-21, with increment 1 to implement (contract,
`dna run-context`, the `length` section, and τ calibration). The one-strand
scale gate of the research note is not part of it: PR #204's one-strand
in-phase confirmation was merged unconditionally after its re-measurement,
and dev truth gives the gate no support. The gate stays a deferred
candidate, to revisit only if the readiness-v4 test set shows a one-strand
confirmation failure.

## Problem

Several of DNA's levels are fixed, although the data say they depend on the
run the sample was sequenced in:

- **Length noise.** Rule A's replicate noise τ (SRS-CONS-013) comes from 743
  re-sequenced samples. Run 1 (HVS-II `303-309`) shows real between-batch
  spread: its dispersion has a between-batch sd of 0.31, and its shrunk scale
  ranges 0.79-1.34 (research note `results/design-notes/run-calibrated-length.md`,
  2026-10-21).
- **Second-channel levels.** `dna.signal/v1` ([[PROP-0009]]) already judges
  a second channel against its own run, but only per 5-mer context and
  reference-free.
- **Context shadows.** Some shadows are tied to a read type and a reference
  position, and their level is batch-dependent. Example: the HV3R
  `phase_shift` at 320.
- **Blob and spike baselines.** These are per instrument chemistry, and are
  measured today only within one read.

Each of these is per-run evidence that a single sample cannot estimate.
Today each one is fixed in configuration, or, for the run signal, built by a
separate command with its own contract and flag. A fourth or fifth such
input would mean a fourth or fifth contract, cache key, and flag.

## Goals

- **One stage.** One run-level stage and one contract, `dna.run_context/v1`,
  computed once per laboratory batch from all of the batch's reads. Every
  consumer reads it through one input.
- **Additive sections.** Each section holds additive sufficient statistics
  only (counts, sums, histograms), so a consumer can take a sample's own
  contribution out (leave-one-out) and apply its own validated calibration.
- **Fail closed.** Without a context, or below a section's minimum
  evidence, every consumer behaves exactly as today.
- **Owner per section.** Each section is computed by the crate that owns its
  method, with pure functions. The facade only assembles the sections.
- **Testable and cheap.** Byte-identical behaviour without a context,
  deterministic output independent of read order, and a cache that keeps
  probes fast.

## Non-goals

- New scientific rules in this proposal. Each consumer's use of a section is
  its own increment, measured with `dna-eval diff` and decided on its own.
- An instrument-run unit. The research found the lab batch is the right unit
  (below).
- Changes to `mtdna_raw` itself. A handoff note describes its step.

## Evidence and research

From `results/design-notes/run-calibrated-length.md` (prototype
`feat/run-calibrated-length` b63f05cb, not merged):

- **Unit.** The lab batch, not the instrument run. A sample's HV2F and HV3R
  sit on different instrument runs 45 % of the time, and the length
  statistic needs both strands. Grouping by instrument run gives the same
  spread.
- **Length.** Run 1 has between-batch variance (sd 0.31; method-of-moments
  prior n0 ≈ 6 pairs). Run 2's batch differences are sampling noise (n0 ≈ 64
  by HV3R instrument run). Calibrating τ per batch moves it by at most about
  ±35 %, and was inert on dev v5 (`dna-eval diff` byte-identical in checks).
  Gating #204 (one-strand in-phase confirmation) on the batch scale is not
  supported by dev truth.
- **Shadows.** HV3R calls at 316-333 are scored yet show unexplained seconds
  of 22-48 % in some batches (D-increment trace, 2026-10-21). The 320
  `phase_shift` level is batch-dependent.
- **Run signal.** [[PROP-0009]] built the reference-free per-run, per-context
  second-channel levels. It is the seed of the `signal` section.

## Proposed design

### Contract `dna.run_context/v1`

One deterministic document per batch: `results/<run-id>.run_context.json`.

- `schema_version`, `run_id` (the caller's batch identifier).
- `provenance`:
  - the configuration's SHA-256 and the plugins;
  - the reference and profile identities (for placed sections);
  - per section, its `method` version and the parameters its statistics
    were normalised with (e.g. the τ the length pairs used).

  A consumer refuses a section whose method or parameters differ from its
  configuration's, as `dna review --signal` refuses another flank today.
- `reads`: the content SHA-256 of every read evidence document that
  contributed, sorted, and the rejected count. A consumer leaves its own
  reads out by these identities.
- **Sections**, each optional. A run context may carry one or several.
  - `length` (placed; owner `dna-consensus::noise`): the profile's recorded
    anchored windows the pairs come from (`windows`, each `[first, last]`),
    the τ and bin edges the differences were normalised and binned with, and
    one entry per window run index, `{run, pairs, sum_z2, histogram}`, over
    those windows. A consumer refuses a section made over other windows, or
    against another reference.
    - `pairs` counts cross-strand pairs. A sample contributes when each
      strand has exactly one decoded view of the run with the same top
      length and a defined `SE`.
    - `sum_z2` is the sum of `d² / (v + τ_F² + τ_R²)`.
    - `histogram` counts the pairs by their normalised `d` in fixed bins.
      The bin edges come from configuration and are recorded in the
      provenance, so a consumer can use a robust dispersion (median,
      trimmed) without new contexts. Run 2's heavy tail needs that: its
      mean statistic does not separate batches.

    Every field is additive, and a sample's own contribution is
    recomputable, so leave-one-out needs no per-sample rows.
  - `signal` (reference-free; owner `dna-sanger::run_signal`): today's
    `dna.signal/v1` keys (context, base, second, share histogram), moved
    here unchanged.
  - `shadows` (placed; owner `dna-call`, which owns per-read trust at
    reference positions): per `(read type, reference position)` the counts
    of in-phase calls with a second peak, by second-base share bin, and how
    many the modality left unexplained.
    - A read type is the amplicon and strand, from the read's placement,
      never its filename.
    - Shares and counts only, never cause labels: the core stays free of
      modality reasons ([[ADR-0080]]).
    - Interpretation belongs to the consumer: the consensus's confirmation
      guard, or the review's unexplained-second test.
  - `blobs` (reference-free; owner `dna-sanger`): per channel and
    scan-region bin, the counts of dye-blob and spike detections, as
    baselines for the read-level detectors.
- The schema and an example go in `docs/reference/schemas/` with a
  reference page. Each section names its owner and method version.

### Where it is computed

- **Command.** `dna run-context <run-id> <traces...|--evidence <docs>...>
  [--placed <docs>...] --reference <fasta>` processes every read of the batch
  as `dna sample` does, up to each section's input:
  - evidence, for `signal` and `blobs`;
  - placement and the consensus's own views, for `length` and `shadows`.

  Placed-read documents (`dna.placed_read/v1`, [[ADR-0087]]) let the
  command and the per-sample commands share one placement per read.
- **Library.** The facade assembles the sections. Each owner exposes a pure
  function from its inputs to its section (e.g. `noise::cross_strand`,
  `noise::LengthStatistics`), and a consumer function that applies
  calibration with leave-one-out (e.g. `RuleParams::calibrated`). A
  section's types live with its only reader (the `length` section in
  `dna-consensus`); a section the modality and the core both read moves its
  types to `dna-kernel` ([[ADR-0069]]: plugins depend only on the kernel).
- **Consumers.** `dna consensus|review --run-context <file>` and the
  in-process API (`SampleAnalysis::with_run_context`). The
  `--signal` flag and `dna.signal/v1` are removed once `signal` lives here,
  leaving one implementation.

### Calibration (consumers' configuration)

- **Length.** New `[consensus]` keys:
  - `length_context_minimum_pairs` (10);
  - `length_context_prior_pairs` ([4.3, 64]: run 1 refitted by the
    method of moments from the 1,000 dev batch contexts over every read;
    run 2's between-batch variance is not distinguishable from zero, so
    its research value stays);
  - `length_context_reference` (the corpus mean z² over every read of the
    dev batches' folders, [0.4577, 0.0542]; HV2F and HV3R alone, as the
    research measured, give [0.4643, 0.0563]);
  - `length_context_minimum_scale` (1.0: τ never below the replicate
    estimate).

  The calibrated noise is `τ_run = τ · max(minimum, √(φ̃/ref))`, where
  `φ̃ = (n0·ref + n·φ)/(n0 + n)` and `n`, `φ` are computed without the
  sample's own pairs. With fewer than the minimum pairs, τ is the
  configured one. Each key is documented with its measured basis.
- **Other sections.** Calibration for the other sections is defined in
  their own increments. Every one must be inert without a context.

### Cache keys

- **Reference-free sections** (`signal`, `blobs`): keyed by the sorted
  evidence-document SHA-256s, the configuration's SHA-256, and the read-stage
  code identity (as the evidence cache).
- **Placed sections** (`length`, `shadows`): additionally the reference's
  SHA-256 and the placement and consensus code identity.
- **`dna-eval`.** `--run-context <dir>` reads each batch's
  `<batch>.run_context.json`, as `dna run-context` writes them over the
  batches' full folders; the directory is the cache. A panel holds too few
  of a batch's samples (about 16, or 2.5 run-1 pairs) to build a context in
  process. It reports context coverage (`run_context.tsv`): per batch month,
  how many samples had a context and in how many it calibrated each run, so
  the cost of failing closed is visible.

### `mtdna_raw` batch integration (handoff note)

The DNA adapter runs `dna review` per sample in a process pool
(`src/tools/dna/pipeline.py::run_batch`). With a run context it becomes:

1. `dna evidence` per trace of the batch folder (all samples, cached by the
   trace's SHA-256);
2. `dna run-context <batch-id> --evidence …` once per batch, writing
   `<batch-id>.run_context.json` (with `--placed`, after a placement pass, to
   place each read once);
3. `dna review <sample> --evidence … --run-context <file>` per sample.

The batch id is the adapter's `request.batch_id`. Without step 2, step 3
runs as today. DNA does not change `mtdna_raw`; this is a handoff note.

### Increments

1. **Contract, command, and `length`.** The contract, `dna run-context`, the
   `length` section (with its histogram), the consensus's τ calibration with
   leave-one-out, and `dna-eval --run-context` with its coverage report. This
   promotes the research prototype. Acceptance:
   - outputs byte-identical to main without a context;
   - with contexts, the dev v5 τ-only variant of the research note, and the
     coverage table per period;
   - then one logged test-set run (readiness-v3, when it exists).
2. **`signal` moves in.** `dna.signal/v1` and `--signal` are retired
   (`docs_graph.py retire`); `dna review` reads the `signal` section.
   Byte-identical review output.
3. **`shadows`.** Research first: per read type and position, batch levels
   against truth. Then a consumer (review or the modality's mixed verdict)
   only where the data support it.
4. **`blobs`.** Baselines for the blob and spike detectors, likewise
   research-first.
5. **Handoff.** The `mtdna_raw` adapter note, and an end-to-end test on one
   batch folder.

## Alternatives considered

- **One contract per run-level statistic** (`dna.signal/v1`,
  `dna.length_noise/v1`, …). Rejected: each adds a contract, a cache key, a
  flag, and an adapter step, for statistics computed from the same batch.
- **The instrument run as the unit** (ABIF `RunN1`). Rejected: the length
  statistic is cross-strand, and strands sit on different instrument runs;
  the spread is the same by either unit.
- **Corpus-wide fixed levels only** (today). Kept as the fallback, but run 1
  has real between-batch spread.
- **Per-sample estimation.** Rejected: a sample has one or two pairs per
  run.
- **Computing the context inside `dna review` from sibling files.**
  Rejected: every sample would recompute the batch, and the result would
  depend on which files a call happens to see.

## Validation plan

- **Unit tests per section.** Additivity (merge of two halves equals the
  whole), leave-one-out equal to the context without the sample,
  order-independence, refusal of mismatched method or parameters.
- **No context.** Byte-identical outputs to main on dev v5 (`dna-eval diff`)
  for every consumer.
- **With contexts.** Each consumer increment is measured on dev v5 per period
  and source, then one logged readiness-v2 use. The gate (non-het) must
  pass.
- **Performance.** The context command's cost per batch, with and without
  placed-read documents.

## Rollout and rollback

Every consumer is opt-in through `--run-context`. Removing the flag (or the
file) restores today's behaviour exactly. Section method versions let a
consumer refuse a stale context.

## Observability and operations

- `dna run-context` logs reads, rejected reads, and per section its counts
  and minimum-evidence status.
- Consumers log which context and sections they used, whether leave-one-out
  applied, and any refusal.
- `dna-eval` records contexts in `run.json`.

## Security and privacy

A run context holds aggregate statistics and read-document hashes, no
sample identity beyond the batch's. It is validated like any input.

## Compatibility and migration

- Adds a contract and a command.
- Increment 2 removes `dna.signal/v1` and `--signal` (pre-release; no
  legacy path kept).
- Consumers' configuration gains keys (a schema bump per increment).

## Decision

Proposed; awaiting the user's review.

## Implementation

None yet. Research prototype for the `length` section:
`feat/run-calibrated-length` b63f05cb (Researcher).
