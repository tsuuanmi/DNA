# approve

The single pre-registered approval run of triage pass rules on a hold-out
panel, as the [approval protocol](../../../../docs/proposals/0005-evaluation-feedback-loop.md)
defines it, and the power projection that sizes it beforehand. Only
`dna-eval approve` may read a hold-out panel; `run` refuses one.

- `prereg.rs` — parses the pre-registration (TOML): the frozen inputs and
  their checksums, the bootstrap draws and seed, the utility floor, the
  gates (measure and bound), and the fixed sequence of hypotheses (policy,
  test condition `first` or `previous_pass`, optional marginal endpoint with
  its recorded projection; the sequence stops at the first failure). A value
  still written as `"pending"` is an error, every policy must keep the `n`
  cause, a marginal bound that 0 leaks among its projected samples cannot
  meet is an error, and each policy is pinned by the SHA-256 of its full
  definition.
- `stats.rs` — the bounds: Clopper–Pearson upper and lower bounds for
  fractional counts through the beta quantile, the Korn–Graubard upper and
  lower bounds of a stratified rate, the binomial distribution function, and the batch-cluster bootstrap with a local `SplitMix64` seeded by
  draw index. Pure and deterministic.
- `mod.rs` — `approve`: checks the frozen inputs before DNA runs (committed
  pre-registration, build revision without local changes, manifest,
  configuration, truth copy, and panel truth digest), runs DNA once, decides
  and scores the hypotheses in order, and writes `approval.json`,
  `approval.md` (aggregates only), `leaks.tsv` (for adjudication), and a copy
  of the pre-registration under a new directory in `results/`.
- `power.rs` — `power`: projects the passes of an approval run per stratum
  from a development run's `triage.tsv` and the panel's stratum sizes, at the
  stratified development rate and at its joint Korn–Graubard lower bound,
  with the attainable bound at 0–10 leaks. With a base and a candidate sweep
  policy it projects, from `triage_policies.tsv`, the samples the candidate
  newly passes, which a marginal endpoint counts, with the attainable bound
  at 0–3 leaks and, for a true rate, the chance that each given bound is met.
  It reads no truth of the panel.

It depends on `../triage` for the policies and decisions, `../report.rs` for
the compared rows, `../estimate.rs` for the stratified rate, and
`../runner.rs` for the DNA pass; nothing depends on it.
