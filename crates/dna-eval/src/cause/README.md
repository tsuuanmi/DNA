# cause

The cause of a sample's discrepancy
([discrepancy causes](../../../../docs/validation/evaluation-harness.md#discrepancy-causes)):
whether DNA, the lab (its assay, runs, or reviewers' truth), or the
sequencing technology accounts for each difference, read from the evidence,
so the difference reaches the one who can act on it.

- `mod.rs` — the rules per discrepancy class, the findings (cause, scope,
  evidence, position or window), and the sample's cause and scope.
- `loci.rs` — per region position, what the admitted reads and the consensus show: covered,
  observed cleanly on each strand, something other than the reference,
  masked for phase loss or for weak signal, unresolved, a tie the consensus
  left `N`, a run stretch it left undecided, and where a clean but
  untrusted call stands, read from the core's own verdict
  (`SampleAnalysis::trust`): beside phase loss, at the read's own end, with a
  weak or mixed support veto, or inside; and the faint second
  bases (3 % or more). The runner takes it from the sample analysis.
- `forms.rs` — how the reviewers write each sequence over the whole truth:
  per repeat window and window sequence, the samples per written form, and
  outside the windows the samples per set of calls (a proxy); written as
  `forms.tsv` and read back for the same truth database.

It depends on `classes`, `compare`, and `lod` for the sample's class,
comparisons, and limit-of-detection signals, on `truth` for the written
forms, and on DNA's sample analysis and the `dna-core` and `dna-sanger`
evidence types for the reads; `report`, `runner`, and `main` (the forms) use
it. It reads no configuration: the core's verdict carries the margins and
thresholds DNA applied. It never rewrites the truth.
