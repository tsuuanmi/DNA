# Kernel model

Owns the vocabulary every plugin family shares. `nucleotide` is the single
authority for canonical-base parsing and A/C/G/T channel order. `reference` is
the validated reference model and the one place that derives its
`ReferenceIdentity` (`Reference::identity`). `placement` is the one
primitive of 3'/right-most indel placement (ADR-0047): the sequence-equivalent
one-base move of a gap that alignment canonicalization in the core and
variant normalization both use. `runs` is the one scan for maximal runs of
equal items (`maximal_runs`); each caller keeps its own rule for runs of
unresolved items.

See [dna-kernel](../README.md).
