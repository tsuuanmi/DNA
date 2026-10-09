# ADR-0065: Result comparison belongs to downstream pipelines

- **Status:** Accepted
- **Date:** 2026-10-08
- **Related decisions:** [ADR-0005](0005-versioned-output-contracts.md),
  [ADR-0019](0019-scientific-evidence-hierarchy.md),
  [ADR-0058](0058-canonical-contracts-and-modular-analysis-composition.md)

## Context

Validating DNA means comparing its calls with other sources: reviewer-checked
Sequencher calls, legacy pipelines, or other tools. Each source has its own
format, notation, analyzed ranges, and sample identities. Converting them, and
deciding what counts as agreement, is project- and lab-specific work.

Downstream pipelines such as `mtdna_raw` already own that work: their tools read
each source, convert it into a canonical sample representation, and compare
those representations.

## Decision

DNA produces evidence and publishes versioned results (`dna.analysis`,
`dna.basecalls`, `dna.sample_evidence`). It contains no logic that compares its
results with external call sets, reviewer truth, or other tools, and no parsers
for their formats.

Downstream pipelines consume DNA's published results, convert them with the
other sources into their canonical sample representation, and compare there.
Their comparison outputs are validation and release evidence about a DNA
revision; they are not DNA code, fixtures, or contracts.

## Consequences

- DNA's public contracts must carry everything a downstream comparison needs,
  such as provenance, profile identity, eligibility, and notation, so that no
  comparison has to reach into DNA internals.
- Metrics quoted in DNA decisions and design (for example the Sequencher
  concordance in ADR-0062 and ADR-0063) are evidence measured outside the
  repository against a named DNA revision.
- Tests in DNA assert DNA's own behaviour against synthetic expectations; they
  do not implement a comparison method.
