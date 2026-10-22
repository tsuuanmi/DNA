# dna-eval

The `dna-eval` binary of the evaluation harness
([ADR-0075](../../../docs/decisions/adr/0075-evaluation-harness-in-workspace.md)).
It holds no harness logic: `main.rs` passes what the binary was built from
to the harness library ([dna-eval-lib](../../dna-eval-lib/src/README.md))
and runs it. It depends on `dna-eval-lib` alone, and nothing depends on it.

The build stamps change with every commit and every `git add`, so they live
here: a new commit recompiles this file and relinks the binary, and the
harness is not recompiled.

Output contract: none of its own; the harness writes
[`dna.eval_run/v1`](../../../docs/reference/schemas/eval-run-v1.schema.json)
([dna-eval-lib](../../dna-eval-lib/src/README.md)).

- `main.rs` — the command: the package name and version and the build
  script's stamps, as `dna_eval_lib::Identity`.
- `../build.rs` — the stamps: the git revision and whether tracked files
  differ from it (`run.json`'s `source`), and each cached stage's code
  identity, the git object ids at the commit of the sources `../stages.txt`
  lists for the stage, or nothing when one of them differs from the commit.
- `../stages.txt` — the sources of each cached stage: the read stage (a
  trace to its read evidence document), the placement stage (a read
  evidence document to its placed-read document), and the index of the
  truth's written forms. `tools/python/scripts/validate_module_layers.py`
  checks that each list is closed under the module graph, so a source a
  stage newly uses cannot be left out of its cache key.

The caches and their keys are specified in the
[evaluation harness](../../../docs/validation/evaluation-harness.md) document.
