# ADR-0071: Run-length edits need resolved run ends

- **Status:** Accepted
- **Date:** 2026-10-10
- **Related decisions:** [ADR-0062](0062-read-callability.md),
  [ADR-0065](0065-result-comparison-downstream.md),
  [ADR-0067](0067-signal-derived-read-callability.md),
  [ADR-0069](0069-plugin-first-modality-core-post-calling.md)

## Context

The core reproduces reviewer calls from reviewed consensus sequences, so on
the same samples it is a ceiling for the Sanger path. A downstream comparison
of the two paths on 160 held-out samples (ADR-0065) put most of the remaining
Sanger differences, 19 of 26, in the poly-C stretches of HVS-I and HVS-II.

In those reads, a homopolymer run ends where the read stops being in phase.
Length heteroplasmy and polymerase slippage in the run make the signal after
it dephased, and callability masks it as `dephased`, an anchoring mask
(ADR-0067). The last in-phase calls of the run then over- or under-count it
by about one base. The resulting edits are an extra run base inserted at the
run end, or the base after the run read as the run base. `read_end`
(ADR-0062) does not catch them, because anchoring calls count as informative.
The read shows a run but not where the run ends.

## Decision

A variant that changes the length of a run of one base is ineligible with the
core reason `run_boundary` when the read does not resolve the end of the run
that the edit touches.

- **The run** is the read's maximal stretch of the edited base around the
  edit:
  - the substituted base;
  - the inserted bases, when they are all one base;
  - or, for a deletion of one repeated base, what the deletion leaves of the
    run, read by a flank.
- **Not run-length edits:**
  - a lone substituted or inserted base, with no neighbouring call of the same
    base;
  - an edit of more than one distinct base;
  - the deletion of a whole run.
- **A resolved run end** is the call just beyond the run. It lies inside the
  read's informative interval and reads another resolved base, either unmasked
  or anchoring-masked. A dephased call that still reads the next base resolves
  the end. A masked call that reads the run base or `N`, an unresolved mask,
  and an unvouched physical read end do not. A vouched read end does
  (ADR-0069).
- **Insertions and deletions** need both run ends resolved and every call of
  the run unmasked.
- **Substitutions** need only the ends they lie on.
- **Placement in the core:** the rule uses only `ReadEvidence` and the read's
  placement, so it lives in the core and applies to every modality. It has no
  parameter. A sweep of a minimum run length from one to eight calls changed
  no Sanger result.
- **Reason order:** the region reason, support vetoes, `read_end`,
  `run_boundary`, then mask reasons.

## Alternatives

- **Vetoing every call next to a dephased segment.** On the same comparison
  this would discard about 28 correct supports to remove 2 false ones.
- **Requiring an unmasked call beyond the run.** On the 320-sample held-out set
  this removed 14 false calls but lost 5 true HVS-II length insertions
  (`309.1C`, `309.2C`), whose run end is a dephased call that still reads `T`.
- **Excluding all of a read's edits in a nomenclature window when any is
  ineligible.** On the same set this lost only true calls.
- **Narrowing `read_end_margin`.** Supports that `read_end` alone excludes have
  equally clean signal whether they are true or false. A narrower margin
  recovers fewer true calls than it admits false ones.

## Consequences

- **320-sample held-out set:** 9 false length calls are removed (`16193.1C`,
  `16194C`, and a `573.1C`–`573.5C` run) and no true call is lost. Precision
  rises from 0.9870 to 0.9894; recall stays at 0.9800.
- **Consensus calls** from `dna call` are unchanged.
- **Remaining HVS-II differences:** when every read meets the 303–315 run at its
  phase boundary, the run's length is still unread. Edits split between
  eligible and masked calls in that window can still be named as another
  haplotype. Repeat-length evidence remains open work
  ([known limitations](../../validation/known-limitations.md)).
- **Contracts:** the sample-evidence and variants contracts gain the core
  reason label `run_boundary`.
- **Sample consensus:** the same test keeps an unbounded insertion or deletion
  from deciding a run's length in `dna consensus` (ADR-0073). The consensus
  applies it to runs rather than edits ([ADR-0074](0074-run-structure-consensus.md)).
  It counts a run end as resolved only when the call beyond the run reads
  another base, and it labels how each length is known.
