---
implements: [SRS-NOM]
decided-by: [ADR-0069, ADR-0075, ADR-0090]
covers: [crates/dna-eval-lib/src/regress/, crates/dna-eval-lib/baselines/regress-*]
---

# Regression set

Part of the [[validation/evaluation-harness|evaluation harness]]: a frozen
set of samples whose DNA output is recorded per version in a ledger, so that
every output change on known samples is declared, recorded, and tied to a
version. `dna-eval regress check` compares a build's output with the ledger
head; `dna-eval regress record` appends a checked run of `main`.

| Instrument | Unit | Expected outcome | Role |
|---|---|---|---|
| [[validation/artefact-bench|artefact bench]] | adjudicated cases | the adjudicated outcome | correctness of known artefacts |
| [[validation/simulated-bench|simulated bench]] | simulated reads | the known mixture | the tract decoder |
| **regression set** | **fixed samples** | **the previous version's golden output** | **a version ledger** |
| `dna-eval diff` on the development panel | panel samples | base against candidate | the gate |

## The set

The set is development data: its samples may be inspected one by one, and it
uses no test budget. It is enriched for hard samples by strata taken from one
build's labels, so it reports counts and changes only, never rates; every
weight is 1. It is never the gate: in the loop
([[engineering/improvement-method#7. The loop]]) it runs after the benches
and before `diff`.

- **Membership is frozen.** A manifest version is never edited; a withdrawn
  sample is retired by `dna-eval panel`. A later version adds rows and is a
  new ledger series (`regress-v2`), seeded by its own first entry.
- **regress-v1** has 589 samples: 562 from recent batches (2026 onwards) of
  the development panels, selected by strata, and 27 adjudicated artefact
  bench samples from before 2026, a labelled exception (stratum
  `bench_pre2026`). Its manifest, checksums, and selection scripts are under
  `results/test-panel/`; it is derived against the truth pin
  `statistics-2026-10-21.sqlite` and checked with the frozen re-check list of
  its selection.
- **Sample identifiers stay out of git** ([[governance/data]]). The
  repository holds the command, the aggregate ledger rows, the per-stratum
  digests, and the intent files, none of which names a sample. The manifest,
  golden outputs, change lists, and the local ledger stay under `results/`.

## Commands

```bash
cargo run --release -p dna-eval -- regress check \
  --manifest <regress-v1.tsv> --truth <statistics.sqlite> \
  --recheck <frozen re-check list> \
  --config config/dna.toml --reference references/rCRS.fasta \
  --ledger results/regress/ledger --threads <1-256> \
  [--evidence-cache <dir>] [--run-context <dir>] \
  [--review-window <bp>] [--zone-after <bp>] [--zone-edge <bp>] [--zone-overrun <bp>] \
  [--intent crates/dna-eval-lib/baselines/regress-intents/<entry>.tsv] [--rebase] \
  --output <new-directory>

cargo run --release -p dna-eval -- regress record \
  --run <check output directory> --ledger results/regress/ledger --entry <rNNNN> \
  --baseline crates/dna-eval-lib/baselines/regress-v1-ledger.tsv \
  --strata crates/dna-eval-lib/baselines/regress-v1-strata.tsv
```

- **`check`** reads only a manifest whose header says `panel=regress`; any
  other is an input error. It runs DNA on every row exactly as `run` does:
  the same runner, caches, run contexts, review options, comparison, and pass
  classes, with the options `run` takes under the same names. It writes the
  golden output and compares it with the ledger head, the last entry of
  `<ledger>/index.tsv`; with no head, it writes the golden output and the
  verdict `REGRESS INITIAL`. It writes nothing to the ledger.
- **`record`** never runs DNA. It refuses a run whose build was dirty or
  records no revision, whose revision is not the current commit of
  `origin/main` (fetch first; a run of an older `main` is checked again),
  whose verdict is `REGRESS FAIL`, whose head is no longer the ledger head, or
  whose `--entry` already exists (in the index, the ledger directory, or a
  baseline file); it also refuses when either committed baseline does not end
  at the local head, or when the intent file changed since the check.
- **`record` re-judges** rather than trusting `run.json`: it reads the
  manifest and the configuration the check names (refusing either if it
  changed since), places the zones from the recorded widths, checks
  comparability with the head again, and judges the run's `golden.tsv`
  against the head's with the intent file. A verdict, rule, digest, count, or
  change count that the re-judgement does not reproduce is refused.
- **`record` writes atomically:** once everything is checked, it stages the
  entry directory (the check, with the intent file as `intent.tsv`) and the
  new `index.tsv` and baseline files as hidden `.<name>.partial` siblings,
  then renames each into place, `index.tsv` last, so an interrupted record
  leaves no half-written file and the head unchanged. A leftover staged file
  is refused until removed. The pull request then commits the ledger
  baseline's new row and the strata baseline's row per stratum.
- **Exit status:** 0 for `SAME`, `EXPLAINED`, `INITIAL`, `REBASE`, and a
  recorded entry; 2 for `FAIL`; 1 for an input error.

## Comparability with the head

A series starts on the truth its manifest pins: the first check refuses a
truth database whose panel digest or per-sample digests differ from the
manifest's. Before recording `r0001`, verify the frozen files against their
checksums, from `results/test-panel/`:

```bash
grep ' (file)$' SHA256SUMS-regress-v1 | sed 's/ (file)$//' | sha256sum -c
```

Every later check measures against the head's truth instead, which only a
rebase changes.

`check` compares, against the head entry's `run.json`: the manifest body
SHA-256 and the series; the panel truth digest; the re-check list SHA-256;
the profile and reference SHA-256; the review window and zone widths; and the
run-context digest. Any difference is an input error, refused before DNA
runs, with two exceptions.

- **`--rebase`**: the truth digest or the re-check list may differ, and one
  of them must; nothing else, the configuration included. The truth-free
  lines must not change (`rebase_dna_change`); the evaluation lines change
  freely, and the verdict is `REGRESS REBASE` (entry kind `rebase`).
- **A new configuration** is a DNA change, not a rebase: its SHA-256 is part
  of the version key, the check accepts it, and the change goes through
  intents like a code change.

## Golden output

`golden.tsv` has the columns `sample_id`, `kind`, `key`, `value`. Samples
come in manifest order; within a sample, kinds in the order of the table;
within a kind, keys ascending by the number they start with (keys without
one last), then by text, then by value. A key a sample repeats is compared
occurrence by occurrence. Cells are UTF-8 with no tab or line end and no
trailing space; an empty cell is `-`. The file is byte-identical for any
thread count, cache state, and host.

| Kind | Group | Key | Value |
|---|---|---|---|
| `failure` | truth-free | `-` | the failure kind (`trace`, `rejected`, `dna`, `panic`, `notation`, `reference`); a sample DNA could not analyse has only this line |
| `allele` | truth-free | the allele in reviewer notation (`309.1C`, `249DEL`) | `recorded`, or `not_recorded` (SRS-NOM-017) |
| `callable` | truth-free | `start-end` (1-based, closed) | `-` |
| `check` | truth-free | `first-last` of a check item | its label names, sorted, space-separated |
| `note` | truth-free | `<tract>:run<run>` | `first-last status called reading convention`, `.` for none, as in `notes.tsv` |
| `review` | truth-free | `clear`, `zones`, `positions` | `Y`/`N`; the checked zones, sorted, comma-joined; the positions to review |
| `empop` | truth-free | `ranges`, `mitotype` | as in `scores.tsv`, when the profile declares an EMPOP query |
| `rejected` | truth-free | `-` | the traces set aside |
| `block` | evaluation | `first-last:<text>` of an unreviewed block or omitted truth allele | `class zone zone_checked`, one line per zone, as in `unreviewed.tsv` |
| `status` | evaluation | each measure; `class`; `cause` | the status; the class and its regions; the cause, scope, and evidence, as in `scores.tsv` |

The review lines read the check items alone, so they are truth-free and
written for every analysed sample. The `golden_sha256` is the SHA-256 of the
whole file; a kind digest is the SHA-256 of that kind's lines in file order;
a stratum digest is the SHA-256 of the lines of the samples whose manifest
`stratum` is that stratum, in file order, and names no sample.

## Changes and intents

`changes.tsv` lists every changed line, in golden order, with the columns
`sample_id`, `stratum`, `kind`, `key`, `change` (`added`: only in the
candidate; `removed`: only in the head; `changed`: another value), `before`
and `after` (`-` when absent), and `intent`: the id of the first intent row
that matches, `rebase` for an evaluation line of a rebase, `forbidden:<rule>`,
or `unexplained`. It has its header only when nothing changed.

A pull request that changes output declares the changes in an **intent
file**, committed as `crates/dna-eval-lib/baselines/regress-intents/<entry>.tsv`.
It names no sample: `check` refuses a file that holds a sample identifier of
the manifest.

```text
# regress-intent series=regress-v1 base=<head entry> pr=<#n>
id	kind	change	where	match	max_samples	posthoc	reason
note-hv2	check	removed	window@HV2	tract_length	-	0	#222 count-only C-tract differences become notes
note-hv2b	note	added	303-316	*	-	0	#222 the same samples gain a note
```

| Column | Content | A changed line matches when |
|---|---|---|
| `id` | unique in the file, `[a-z0-9_-]+` | – |
| `kind` | a golden kind, or `*` | equal, or `*` |
| `change` | `added`, `removed`, `changed`, or `any` | equal, or `any` |
| `where` | `a-b` (positions), a zone such as `window@HV2`, a region of the profile such as `HV1`, or `*` | positions: the line's positions overlap `a-b` (an allele's anchor, a callable interval, a check item, a block, or a note's span); zone: a block's zone, or a zone a check item touches; region: the line's positions overlap the region. A line without positions (`review`, `status`, `empop`, `rejected`, `failure`) matches a zone or region only when its sample's changed blocks or check items lie there |
| `match` | a glob (`*`, `?`) | the key, the value before, or the value after matches; for a check item, also any one of its labels |
| `max_samples` | a positive integer, or `-` | – (a bound: `intent_exceeded`) |
| `posthoc` | `0` when declared before the first check of the change; `1` when added after a check showed it | – |
| `reason` | a pull request, ADR, or SRS reference and a short phrase | – |

The header's `series` and `base` must be the check's series and head entry.
The file must be committed and unmodified at the harness's revision: `check`
compares its git blob with the blob at the revision the harness was built
from and records the blob id. Pre-declaration rests on git history: the
commit that adds a `posthoc=0` row must not be later than the first regress
output posted for that change, which reviewers check; the ledger counts the
`posthoc` rows per entry. Without an intent file, every change is
unexplained.

## Verdict

The rules are evaluated in this order, and `REGRESS FAIL` lists every rule
that fired. A line a forbidden rule names carries it as its label.

| Rule | Fires when |
|---|---|
| `forbidden:escape` | a `block` line of class `non_het` with `zone_checked` `no` is added, or a block changes into that state; whatever the intents. A rebase only reports such blocks, since it changes the truth they are measured against |
| `forbidden:flag_lost` | a sample whose head `haplotype_lod` status is `discordant` goes from `review clear N` to `Y` while still `discordant`; whatever the intents |
| `forbidden:failures` | a `failure` line is added or removed: the set of failed samples differs |
| `forbidden:no_bump` | a truth-free line changed while the candidate's version key equals the head's |
| `intent_exceeded` | an intent row matched lines of more distinct samples than its `max_samples` |
| `unexplained` | a changed line matches no intent row |
| `rebase_dna_change` | with `--rebase`, a truth-free line changed |

Otherwise the verdict is `REGRESS SAME` when no line changed, `REGRESS
EXPLAINED` when every changed line is explained, and `REGRESS REBASE` for a
rebase. Reported, never the verdict: a version key that moved without a
changed line (the change may lie outside the set), intent rows that matched
nothing (`unused intent`), `posthoc` rows, and the version-key fields that
moved.

The **version key** is the workspace version, every plugin `id=version` of
the static registry ([[ADR-0069]]; [[governance/versioning#Plugins]]), the
result contracts DNA emits (every `dna.<name>/v<n>` the published schemas
declare, but the harness's own), and the configuration SHA-256. A truth-free
change therefore needs a plugin version bump, a new contract version, or a
configuration change.

## Outputs

`check` writes `dna.eval_regress/v1`
([schema](../reference/schemas/eval-regress-v1.schema.json),
[example](../reference/examples/eval-regress-v1.example.json)):

| File | Content |
|---|---|
| `golden.tsv` | the golden output |
| `changes.tsv` | the changed lines |
| `regress.md` | the verdict and the rules that fired, warnings, the head and candidate revisions and version keys, the changes by kind and change, by stratum (changed samples of its members), and by intent row (lines, samples), the first 40 unexplained or forbidden lines, and the counts; the listed lines name samples, so the file stays local |
| `run.json` | the harness and source revision and dirty flag, the manifest fields, the truth and re-check identities, the configuration, profile, and reference identities, the review options, the run contexts, the version key, cache usage, the head entry, the intent file (blob id, rows, `posthoc` rows), the verdict and rules, the warnings, the golden, kind, and stratum digests, the counts, and the timings |

## Ledger

- **Local**, under `results/regress/ledger/`: `index.tsv` (`entry`, `kind`,
  `base`, `revision`, `golden_sha256`, `verdict`, `recorded`) and one
  directory per entry with a copy of its check and its intent file.
- **Committed, append-only**: `crates/dna-eval-lib/baselines/regress-v1-ledger.tsv`,
  one row per entry: `entry` (`r0001`, … per series), `date`, `kind`
  (`initial`, `change`, or `rebase`), `base`, `revision`, the version key
  (`dna_version`, `plugins`, `contracts`, `configuration_sha256`),
  `profile_sha256`, `reference_sha256`, `manifest_body_sha256`,
  `truth_digest`, `recheck_sha256`, `golden_sha256`, `kinds_sha256`
  (`kind=<sha256>` in kind order, `;`-joined), the counts (`samples`,
  `failed`, `clear`, `with_check`, `checks`, `notes`, `alleles`,
  `blocks_non_het`, `escapes`, `blocks_het`, `blocks_length_note`,
  `blocks_truth_n`, `blocks_recheck`), the changes against the base
  (`changed_samples`, `changed_lines`, `explained_lines`, `posthoc_rows`),
  `intent_blob`, and `verdict`. `crates/dna-eval-lib/baselines/regress-v1-strata.tsv`
  has one row per entry and stratum: `entry`, `stratum`, `samples`,
  `golden_sha256`, `changed_samples`. Each file opens with the comment line
  `# regress-ledger series=regress-v1` or `# regress-strata
  series=regress-v1`.

`clear` and `with_check` count samples; `checks`, `notes`, and `alleles`
count golden lines; the `blocks_*` count block lines by class (the
omission classes have no count); `escapes` counts the non-het blocks in
unchecked zones.

**For each output-changing pull request:**

1. Commit the intent file before the first regress check of the change.
2. Run `check`; resolve or declare (`posthoc=1`) every change; post the
   verdict of `regress.md` (its aggregate sections) in the pull request.
3. After the merge, run `check` and then `record` on a clean build of `main`.
4. Commit the two baseline rows; the intent file is already on `main` from
   step 1.
5. Cite the entry id in the CHANGELOG.

The first entry `r0001` is an `INITIAL` on the `main` revision that
introduces the command.

## Limits

- **Not yet shared with `diff`** (specification section 3.1, item 3): `run`,
  `bench`, and `regress check` share the options that run DNA, and `run` and
  `regress check` one block classifier, but `diff` still classifies passes
  with its own severities; unifying them waits for the set-based rework of
  the `diff` gate, which rewrites that classifier. `--current-since`
  stays `run`'s alone: it orders `run`'s summary and the gate of `diff`, and
  neither `regress check` nor `bench` reports a current period, so it would
  be an unused option there.

- The set guards escapes mostly through its flagged discordant samples and
  `forbidden:flag_lost`: few recent samples carry a non-het escape.
- The strata are a coverage device from one build's labels; membership stays
  frozen whatever the labels become.
- A change outside the set's samples shows only as a moved version key with
  no changed line.
