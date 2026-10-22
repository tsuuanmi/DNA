# regress

`dna-eval regress`: the versioned regression set
([regression set](../../../../docs/validation/regression-set.md)). `check`
runs DNA on a `panel=regress` manifest through the harness's shared run
(`Loaded`, `Ran`, and the report of `main.rs` and `report.rs`), writes the
golden output, and judges it against the head of the version ledger;
`record` appends a checked run of `main` to the ledger. It reads the facade's
plugin registry (`dna::plugins()`) and the contracts `build.rs` records for
the version key. Nothing else in the harness depends on it.

- `mod.rs` — the command line, the run header `run.json` (`Header`, with the
  `VersionKey`), the comparability of a check with the head, the intent
  file's git checks (`Git`: committed and unmodified at the harness's
  revision; a recorded revision that is `origin/main`'s current commit), the
  judgement and the `run.json` fields it derives (`Judged`, shared by `check`
  and `record`'s re-judgement), `regress.md`, and `record`'s refusals and
  staged writes.
- `golden.rs` — `golden.tsv`: the golden kinds, built from DNA's outcomes and
  the run's report, their canonical order, the strict reader, the file,
  kind, and stratum digests, and the ledger's counts.
- `intents.rs` — the intent file: its header, rows, `where` forms (positions,
  zone, region), and the glob of `match`.
- `judge.rs` — the changed lines of a candidate against the head, the rules
  in order (the forbidden rules, an escape by a rise of a sample's escape
  count in a zone, the version bump, the intents' bounds, unexplained lines,
  a rebase's truth-free changes), the escape lines changed without a rise,
  the verdict, and `changes.tsv`.
- `ledger.rs` — the local `index.tsv`, the committed baseline files' rows and
  their strict reader, entry ids, and the UTC date of a record.
- `tests.rs` — the run header, comparability, and `record`'s refusals on
  synthetic headers and ledgers, and `check` and `record` end to end on a
  synthetic set: a fixture truth database, a derived `panel=regress`
  manifest, and traces DNA cannot read.

Tests use synthetic outcomes, golden lines, and ledgers only; the manifest,
golden outputs, and the local ledger name samples and stay under `results/`.
