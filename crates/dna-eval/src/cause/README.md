# cause

The cause of a sample's discrepancy
([discrepancy causes](../../../../docs/validation/evaluation-harness.md#discrepancy-causes)):
whether DNA, the lab (its assay, runs, or reviewers' truth), or the
sequencing technology accounts for each difference, read from the evidence,
so the difference reaches the one who can act on it.

- `mod.rs` — the rules per discrepancy class, the findings (cause, evidence,
  position or window), and the sample's cause.
- `loci.rs` — per region position, what the admitted reads show: covered,
  observed cleanly on each strand, something other than the reference,
  masked for phase loss or for weak signal, unresolved; and the faint second
  bases (3 % or more). The runner takes it from the sample analysis.
- `forms.rs` — how the reviewers write each sequence over the whole truth:
  per repeat window and window sequence, the samples per written form, and
  outside the windows the calls each sample carries (a proxy).

It depends on `classes`, `compare`, and `lod` for the sample's class,
comparisons, and limit-of-detection signals, and on DNA's sample analysis
for the reads; `report` and `runner` use it. It never rewrites the truth.
