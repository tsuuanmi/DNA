---
id: PROP-0011
status: implementing
created: 2026-10-19
summary: "every crate's output is a versioned, schema-validated contract producible on its own (`dna.read_evidence/v1`, `dna.analysis/v10` with per-call trust and component calls, schemas for the evaluation outputs), with a crate → contract → command map, so a wrong result is traceable to the step that made it."
implements: [SRS-OUT, SRS-API-014]
decided-by: [ADR-0069, ADR-0075, ADR-0079, ADR-0080]
implementation: [PR-160, PR-161]
---

# Proposal: A contract at every step

## Problem

The workspace is split into one crate per step (ADR-0069, ADR-0079). The point of that split
is to find the failing step when a result is wrong. That only works if each step's output
can be inspected on its own. Today several boundaries have no serialized contract, so finding
the faulty step means patching code.

Two examples from one day of work (2026-10-19):

- **Why a read-start call was withheld.** This needed a probe build that printed each call's
  placement, vetoes, trust and withholders, and later one that dumped every aligned call with
  its quality and height. The modality → core boundary, `ReadEvidence`, exists only in memory.
- **The HVS-II "one-call lag" of the length decomposition.** The fault was localised step by
  step:
  - the signal and DNA's primaries were right (296/300 reads, checked by re-reading AB1 files
    in Python, because `dna.basecalls/v3` carries no per-call data);
  - the component sequences were right;
  - the evaluation's measurement was wrong.

  Every step had to be checked by hand-built scripts.

| step (crate) | output | contract today |
|---|---|---|
| basecalling, QC, callability (dna-sanger) | primaries, trim, phase segments | `dna.basecalls/v3`, without per-call signal (heights, relative quality, second peak) |
| per-call evidence (dna-sanger → dna-kernel `ReadEvidence`) | verdicts per scope, reasons, `reference_only`, second base and cause, `second_unexplained`, `PhaseMixture` | **none** |
| placement (dna-align) | orientation, columns | partly in `dna.analysis/v9` (`alignment`) |
| trust and per-read calls (dna-call) | per-call trust, withholders, eligibility, variants, component calls | variants in `dna.analysis/v9`; **trust and component calls: none** |
| sample reconciliation (dna-sample) | sample evidence | `dna.sample_evidence/v10` |
| consensus (dna-consensus) | segments, sites, `reference_only` | `dna.consensus/v1` |
| call of decisions (dna-consensus `decided`) | variants | `dna.variants/v1` |
| normalization, nomenclature, EMPOP | notation | `dna.notation/v1` |
| review (dna) | review list | `dna.review/v1` |
| run signal (dna-sanger) | per-run background | `dna.signal/v1` |
| evaluation (dna-eval) | `run.json`, TSV tables | **none** (convention only) |

## Goals

- Every step's output is a versioned, schema-validated contract that can be produced on its
  own, for one read or one sample, from the same code path as production.
- The data flow is legible: for each crate, its input contract(s), its output contract, and
  the command that writes it, listed in one place.
- A check fails when a crate's output has no declared contract.

## Non-goals

- New scientific behaviour. Contracts expose what each step already computes.
- Stable public APIs for every internal type. A contract is the serialized evidence a step
  hands on, not its Rust types.
- Per-call arrays in operational logs (SRS-CALL-006 stays).

## Evidence and research

The examples above. Both investigations needed scratchpad probe builds (`DNA_PROBE_DUMP`,
review `probe` fields) that were not production code paths. The same information, written by a
contract, would have localised each fault in one command.

## Proposed design

Increments are merged one at a time. Each is display-only: no call changes, and existing
contracts stay byte-identical apart from their own new fields.

### Increment 1: `dna.read_evidence/v1` (modality → core)

- Command `dna evidence <trace.ab1>`, reference-free. It writes the read's `ReadEvidence`
  exactly as the core receives it.
- Content per call:
  - index, primary and ambiguity, scan position;
  - primary peak height and relative quality;
  - second base, share and cause, and `stands_out`;
  - placement, with withheld reasons per scope (identity, indel, bound), edge;
  - `second_unexplained`, `reference_only`.
- Content per read: informative interval, read ends, phase segments, and `mixtures`
  (`PhaseMixture`, ADR-0081).
- Owner: dna-sanger produces the values; the facade serializes. The schema follows the
  kernel contract, so another modality writes the same shape.

### Increment 2: `dna.analysis/v10` (placement and per-read calling)

`dna analyze` gains, per aligned call:
- its reference position and column;
- the core's trust (`trusted`, `exempt`, withholders, scope);
- the read's component calls (span, sequence, edits) for each mixture.

A read then shows the whole chain sanger → align → call in one document.

### Increment 3: evaluation contracts

`run.json` and the harness tables dna-eval writes (`scores.tsv`, `components.tsv`, the
workload tables) get documented, versioned schemas. Evaluation outputs are then checked the
same way as production outputs. The 2026-10-19 lag was a measurement fault.

### Increment 4: contracts the next step reads

A contract that is only written can be checked but not replayed. Increment 4 makes each
step able to start from the previous step's file, so a later step can be rerun, compared,
or fed an edited input without the traces:

- 4a. `dna analyze --evidence <read_evidence.json> --reference <fasta>` places and calls a
  read from its `dna.read_evidence/v1` alone. The contract then carries everything the core
  needs: bases, profiles, verdicts, bounds, flags, informative interval, ends, and mixtures.
- 4b. `dna sample`, `dna consensus` and `dna review` start from a sample's read-evidence
  documents in place of its traces (`--evidence`) and run the core in process, as
  `dna analyze --evidence` does. A read the Sanger plugin sets aside gets a document with
  its rejection, so the sample documents still list it with its cause.
- The acceptance test for every step: running through the files is byte-identical to
  running in process, on the synthetic tests and on real traces. A field missing from a
  contract then shows as a difference, not as a silent default.

### Data-flow map and enforcement

- `docs/architecture/data-flow.md` gains a table: crate → input contract(s) → output
  contract → command.
- Each crate's README names its output contract.
- `validate_docs_structure` (or the module-layer check) fails when a crate README names no
  output contract or names one without a schema.

## Alternatives considered

- **Debug logging behind an environment variable.** A non-contract path that drifts from
  production, and the log rules forbid per-call arrays.
- **Probe builds per investigation.** Today's practice. Slow, and not reproducible by others.
- **One giant debug document.** It would mix steps and defeat the purpose of locating the
  failing one.

## Validation plan

- Schemas, synthetic examples, and validator cases for each new contract and version.
- For each increment, every existing contract is byte-identical to main on panel v5, apart
  from version fields.
- Reproduce the two 2026-10-19 investigations from the new contracts alone, without a probe
  build:
  - the read-start withholding classes, from `dna.read_evidence/v1`;
  - the HVS-II lag localisation, from `dna.analysis/v10` and the eval schema.

## Rollout and rollback

Each increment is a PR. New commands and versions are additive. A version bump of
`dna.analysis` follows the versioning policy; consumers pin the version.

## Observability and operations

This proposal is itself the observability change: per-step outputs replace probe builds.
Operational logs are unchanged.

## Security and privacy

The outputs are per-trace scientific evidence, of the same kind as today's analysis
documents. Nothing new leaves the machine.

## Compatibility and migration

`dna.analysis/v9` → `v10` is a new version. `dna.read_evidence/v1` is new. No existing
contract changes meaning.

## Decision

Accepted by the user, 2026-10-19: "each crate should have an input and output contract, so the
output of that step can be checked and the flow is clear".

## Implementation

- Increment 1, `dna.read_evidence/v1` and `dna evidence`: #160 (SRS-API-017,
  `docs/reference/read-evidence.md`).
- Data-flow map and enforcement: the "Step contracts" table in
  `docs/architecture/data-flow.md`, an "Output contract:" line in every crate
  README, and the check in `validate_docs_structure.py`.
- Increment 2, `dna.analysis/v10`: each aligned call's reference column and
  the core's trust in it; component calls follow PROP-0010.
- Increment 3, evaluation contracts: `dna.eval_run/v1`
  (`docs/reference/schemas/eval-run-v1.schema.json`), one contract for the
  `dna-eval run` directory: which files a run writes, `run.json`, and every
  table's columns, formats, and meaning. `validate_result_schemas.py` checks
  the synthetic example and, with `--eval-run`, a run directory; dna-eval's
  tests check the headers and enumerated cells it writes.
- Increment 4 (contracts the next step reads): accepted by the user on 2026-10-19.
  4a, `dna analyze --evidence` (SRS-API-018), is byte-identical to `dna analyze <trace>`
  on 999 of 1,000 real traces; the remaining read is rejected by the Sanger plugin both
  ways. Its evidence now carries full-precision profile weights: rounding them to three
  decimals changed 2 of 1,000 analyses. 4b, `dna sample|consensus|review --evidence`
  (SRS-API-019), is byte-identical to the trace-based commands on 1,120 real samples
  (30 picked with every rejection kind, the 90 samples of one run reviewed with
  `--signal`, and 1,000 random): sample, consensus document and FASTA, and review. For it
  `dna.read_evidence/v1` changed in place (pre-release): the trace's file name in
  provenance (the read's name in the sample documents), each call's window peaks
  (`signal.channels`, in place of `height`, which the run signal counts), second shares,
  mixture fractions, and residuals at full precision, and a document with its `rejection`
  for a read the Sanger plugin sets aside, which `dna evidence` now writes instead of
  failing.
