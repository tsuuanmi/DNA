# ADR-0078: Call the consensus's decided edits without placing it again

- **Status:** Accepted
- **Date:** 2026-10-17
- **Related decisions:** [ADR-0047](0047-canonical-right-aligned-mtdna-gaps.md),
  [ADR-0071](0071-run-length-edits-need-resolved-run-ends.md),
  [ADR-0073](0073-sample-consensus.md),
  [ADR-0074](0074-run-structure-consensus.md)
- **Proposal:** [PROP-0006](../../proposals/0006-call-consensus-decisions.md)

## Context

The recommended workflow is `consensus` then `call`. Until now `call`, and
`SampleAnalysis::call` in process, read the consensus FASTA, placed each
record on the reference again, and extracted variants from that new
alignment. The consensus had already placed every base and decided every
edit. The FASTA cannot carry where a base sits or why a base is `N`, so around
undecided bases the re-alignment was free to choose differently from the
consensus. Each consequence was guarded separately:
- a record placed away from its interval (SRS-VAR-017);
- an AC deletion beside `N`s (SRS-ALN-017, SRS-VAR-018).

## Decision

1. `call` given a consensus document (`dna.consensus/v1`), and
   `SampleAnalysis::call`, call each consensus segment from the consensus's
   decisions, without placing it on the reference again (SRS-VAR-021).
   - Positions outside every site read the reference.
   - A site's call is aligned against its own reference interval only. An
     `N` stands for the position it replaces and is not callable.
   - Every insertion and deletion then takes its 3'-most sequence-equivalent
     placement over the segment, the move of ADR-0047.
2. The decided variants are judged by the same record rules as a placed
   record (SRS-VAR-013 to SRS-VAR-020). Waiving `run_boundary` for lengths
   the consensus decided by run structure was measured and rejected: it
   reported lengths whose bound the consensus itself could not read.
3. One rule is the segment's own. A length change of the reference run at a
   segment's first or last position is excluded as `run_boundary`, since no
   read shows that run's outer end (ADR-0071).
4. A FASTA sequence keeps the placement path, with its guards, because
   nothing else knows its placement.
5. SRS-API-012 changes from "exactly as `call` reads the consensus FASTA" to
   "exactly the consensus's decided edits", the result of `call` given the
   consensus document.

## Consequences

- One path from consensus decisions to called variants serves the command
  and the library. The plugin is `decided`, version 1, with
  `dna.consensus/v1` as its input contract.
- The `dna.variants/v1` reads of a consensus input are the consensus
  segments. Their alignment summary describes the segment's own columns.
  This is revised in place under the unreleased rule.
- Development panel v5 against main (`haplotype_lod`):
  - 14 samples become concordant and none discordant;
  - 8 samples that failed to call now call;
  - no exported EMPOP range covers a position the consensus left undecided
    (16 positions before, 0 now).
- A length estimate (`phase_loss`, `reference_frame`) is still reported
  where a record of the same sequence would report it. Making it uncallable
  is PROP-0008 step 3.
