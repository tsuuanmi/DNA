---
decided-by: [ADR-0069, ADR-0075, ADR-0080, ADR-0089, ADR-0090]
---

# Improvement method

How the team changes DNA's scientific behaviour: **census first, then fix the
logic**. It applies to every change that can move a call, a label, a check, or
a pass. The tools it names are specified in the
[[validation/evaluation-harness|evaluation harness]], the
[[validation/artefact-bench|artefact bench]], and the
[[validation/simulated-bench|simulated bench]]; this document says how to use
them, not what they compute.

## 1. Census first

Before changing a rule, classify what actually blocks or errs:

- by label and by evidence state (which reads, which strand, what each read
  shows), against the truth;
- on development data, the current period first;
- counting every member of the class, not the examples that prompted it.

Then fix the logic the census points to. A rule written from a few examples
fits those examples; a census shows whether the class is what it seemed, how
large it is, and what the rule would cost elsewhere.

Examples:

- **A hypothesis disproved before any code.** "Downstream `phase_shift`
  seconds are explained by the decoded length mixture" was the proposed fix
  for the window@HV2 checks. Its census on dev v5 showed that second peaks
  did not cause those checks: `phase_shift` there is a context label, and only
  1 of 1,030 window items had an unexplained `possible_heteroplasmy`. The
  real blocker was the window-wide `reported` rule, and the fix became
  per-run reporting (#226).
- **Length notes.** A census of C-tract blocks separated count-only
  differences from structure differences. Forensic practice (ISFG, EMPOP)
  notes a length heteroplasmy whose dominant length is determined rather
  than comparing it, so the review and the harness now report those as
  length notes, and structure differences stay checks (#222).
- **One-strand gates.** A census of decided insertions by the other strand's
  reading ([[ADR-0089]]) found the class where one strand alone errs.

Signal hypotheses (denoising, decomposition, new models of the trace) come
after the logic is right: a logic error measured as a signal problem sends the
work to the wrong crate.

## 2. Localize by crate

The crate split exists to find the failing step. Before changing a rule,
inspect the stage outputs in pipeline order and stop at the first that is
wrong:

1. the modality, `dna-sanger`, with its `dna-trace` decomposition: the read
   evidence;
2. `dna-align`: the placement;
3. `dna-call`: the read assertion;
4. `dna-consensus`: the decided calls;
5. the review ([[reference/review]]): labels, checks, and zones;
6. `dna-eval`: the comparison and its classes.

Fix the step that first goes wrong, not a later one that could hide it.
Ownership of each step is in the [crate router](../../crates/README.md).

## 3. Precision first

- **A flag beats a wrong pass.** An error the review flags costs a reviewer a
  look; an unflagged error reaches a result.
- **The gate** compares the sets of escapes in the current period, the
  samples whose `batch_date` is at or after `--current-since`
  ([[validation/evaluation-harness#Comparing two runs]]): per zone kind, the
  non-het blocks and the truth alleles the report omits outside a noted
  tract, in zones no check item touches, that a change adds and removes. A
  change must not add more than it removes, and every escape it adds is
  adjudicated before a merge.
- **Heteroplasmy is optional** but always reported: het changes go in every
  result, beside the gate.
- **Measure before adopting.** An earlier blanket rule, "every `phase_shift`
  second is explained", cleared more samples in its probe but hid real
  heteroplasmies and added non-het escapes, so it was rejected.
- **Truth is guidance, not an oracle.** A difference the raw signal shows to
  be a truth error goes to the re-check list or to the `truth_n` and
  `length_note` classes ([[validation/evaluation-workload]]). It is never fitted:
  no rule is tuned to reproduce a reviewer's habit.

## 4. Data-driven levels

No hard-coded positions or thresholds in a rule. A level is either
configuration (validated, documented in [[reference/configuration]], and
recorded in provenance) or derived from the data: from the read itself, or
from its batch through the run context ([[PROP-0014]]). Format limits and
numeric tolerances are not levels.

## 5. N reads, modality vocabulary

- A rule handles any number of reads per site and per strand; none assumes
  one forward and one reverse read ([[ADR-0090]]).
- Modality vocabulary stays in the modality. The core reads neutral
  assertions; a Sanger cause such as a second peak's shape is the Sanger
  plugin's to name and act on ([[ADR-0069]], [[ADR-0080]]).

## 6. Panels and their uses

| Panel | Use |
|---|---|
| development: dev v5 and readiness-v2 | tuning, censuses, per-sample inspection, `dna-eval diff` |
| the untouched test set: readiness-v6 | one logged use per decision, aggregates only, each `source` reported separately; never inspected per sample while tuning |
| the reserved hold-out remainder: the approval hold-out less the batches readiness-v6 took | none: it stays reserved. `dna-eval` refuses a hold-out panel; the approval path that will read it is designed later, with endpoints per zone kind among auto-passed zones on recent batches (evaluation validity review, recommendation 5) |

- Test sets are built from recent batches only.
- A test-set result that drives a rule change makes that panel development
  data; the next test set needs batches it does not share.
- Every test-set use is recorded in the readiness log (`readiness-log.md`
  under `results/`), which is not in git.

## 7. The loop

1. **Probe** from a clean committed tree, and record its revision; a dirty
   build is not evidence ([[engineering/development#Evaluation probes]]).
2. **Artefact bench** (`dna-eval bench`): the adjudicated cases first, in
   seconds.
3. **Simulated bench** (`dna-eval sim-bench`) for any tract-decoder change
   ([[PROP-0015]]).
4. **`dna-eval diff`** on the development panels, base against candidate. A
   subset run is a pre-filter only.
5. **Equal failures.** Check that the samples that failed to score are the
   same on both sides; otherwise the two runs measure different samples.
6. **Combined stacks.** Measure the set of changes that will land together,
   not only each pull request on its own. Rules interact: an end-norm change
   combined with a downstream-confirmation change produced a false deletion
   that neither made alone. Measuring the stack caught it, and the artefact
   bench now guards it.
7. **Merge** only after CI has run on the pull request combined with the
   current `main`.
8. **Bump** the version of every plugin whose output changes
   ([[governance/versioning#Plugins]]).

## 8. What to report

Current period first, then the earlier period and each source. For each
change:

- clear samples;
- checks by zone;
- non-het changes (the gate) and het changes;
- costs: samples no longer clear, new checks, and run time when it moves.
