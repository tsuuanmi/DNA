# ADR-0009: Report Primary-Sequence Differences, Not Genotypes or Heteroplasmy

- **Status:** Accepted
- **Date:** 2026-08-22

## Context

The MVP analyzes Sanger chromatograms of mitochondrial DNA. A naive reading of
the output could imply genotype, allele-fraction, or heteroplasmy calls. The
implemented pipeline derives variants from a single primary basecall sequence
and performs no two-allele decomposition, no allele-frequency estimation, and no
heteroplasmy inference. The output must not be mistaken for such a call.

## Options

1. Emit variants with implied genotype or heteroplasmy semantics.
2. Emit only normalized primary-sequence differences and document the
   limitation explicitly.

## Decision

Choose option 2. Variants are classified `primary_sequence_difference` and are
derived solely from the primary sequence alignment. The pipeline does not
estimate allele fractions, genotype likelihoods, or heteroplasmy levels, and it
does not decompose mixed or two-allele signals. The documentation
(`docs/pipeline.md`) states these limitations explicitly.

## Consequences

The output is honest about what the MVP computes and does not overstate
biological meaning. Consumers who need genotype or heteroplasmy calls must use a
different tool. The variant record carries supporting and flanking call indices
so the evidence for each primary-sequence difference is auditable.

## Supersession

A future ADR may add two-allele decomposition or heteroplasmy estimation as a
distinct, clearly labeled output; it must not be conflated with the
primary-sequence difference report.
