# Changelog

All notable user-visible and production-readiness changes are recorded here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
DNA uses semantic versions for tagged releases, while scientific JSON contracts
remain independently versioned and are never silently changed in place.

## [Unreleased]

Plugin versions rise with every change to a plugin's output, also before a
release (`docs/governance/versioning.md`). Two changes are recorded without
a bump of their own (2026-10-17, user: no retroactive bump): the
`nomenclature` window form resting only on a wholly read window (#77) shares
the `nomenclature` 2 of #76, and the bounded deletion trade (#79) shares the
`normalization` 2 of #78.

### Added

- `dna.placed_read/v1` (ADR-0087, **proposed**): one read's placement as
  `dna-align` selected it, which `dna-align` alone writes and reads back,
  re-checking it against the read's evidence and the reference (PROP-0011).
  `dna_call::read_call::call_read_placed` starts the core from it;
  `dna::recorded::placed_read_document` and `analyze_placed_outcome` write
  and read it in memory. `dna-eval --evidence-cache` also keeps a placement
  cache keyed by the placement stage's code identity and the reference
  (`run.json` `placement_cache`). On development panel v5 every table is
  byte-identical with and without it, and a consensus-only probe takes
  29-36 % less CPU.
- `dna-eval subset --manifest <panel> --rows <n> --minimum <n> --salt <s>
  --output <file>` writes a deterministic stratified quick subset of a
  panel, every stratum kept, as a manifest of its own, for a first-pass
  `dna-eval diff` before the full panel. Its verdict is a pre-filter, never
  the gate. Development panel v5's 2,111-row subset (`--rows 2000 --minimum
  30`) runs in about 29 s on a warm cache against about 75 s for the panel
  on the same loaded machine.
- `dna-eval diff --base <run> --run <run>` compares a run with a base run of
  the same panel from their directories alone, and gates it.
  - **Comparability:** it refuses runs with another manifest, truth, re-check
    list, current month, review window, or zones, or a dirty build (`--force`
    warns instead), and runs whose files disagree with their `run.json`.
  - **Report:** per group (all, period, source), samples, clear samples,
    unreviewed samples, passed samples by their most severe block, and
    positions to review; per changed zone, units needing a check and passed
    units with a non-het block; coverage changes by allele; status changes
    per measure; and the changed samples, in `samples.tsv` and
    `statuses.tsv`.
  - **Verdict:** it exits 2 with `GATE FAIL` when non-het escapes or
    unreviewed samples rise in the current period (`batch_date` in or after
    the runs' current-data month; all samples when the runs record none), or
    when the samples the harness could not score there differ between the
    runs: the runs then compare different samples (2026-10-21: transient
    trace-read failures under load cost one run 47 samples and went
    unnoticed). The earlier period and the sources, which mix periods, are
    reported in full and never gate (user decision 2026-10-21: the LAB keeps
    improving, so recent data decides).
  - **New table:** `dna-eval run` writes `unreviewed.tsv`, one row per
    unreviewed block with its class, zone, and whether a check touches the
    zone (`dna.eval_run/v1`, in place, pre-release).

- A reportable region's poor start is not analysed (SRS-CONS-012,
  ADR-0085, user decision 2026-10-20). `dna.consensus/v1` publishes
  `unreliable_only`, the positions only reads' unreliable spans observe. The
  call leaves out, at each end of a region, those up to the last one the
  consensus does not settle as the reference, and excludes changes there as
  `poor_start_signal`. `dna.variants/v1` lists the intervals as
  `not_analysed`, and the EMPOP query cuts them as `not_analysed`.
  `dna.review/v1` lists them as `not_analysed` with no item or zone
  (`summary.positions.not_analysed`), makes a check of any difference within
  `review.edge_zone` of such an end (`beside_not_analysed`), and its
  `poor_signal` label reads the consensus's `unreliable_only`.
  `review.edge_zone` defaults to 40 (was 20). On development panel v5's
  current slice, clear samples rise from 73.28 to 75.99 %, and non-het passed
  zones are unchanged in every zone and period. 73G is no longer analysed in
  117 current samples, where it becomes a coverage item. Plugins `consensus`
  17, `decided` 4, `review` 20, `empop` 2.
- `dna-eval run --evidence-cache <dir>` analyses every read from its
  `dna.read_evidence/v1` document instead of its trace, read from the cache
  or made and stored, and `dna-eval cache prune --keep <n>` removes old
  generations (PROP-0011).
  - **Key:** the trace's SHA-256, the configuration's SHA-256, the `sanger`
    version, and the read stage's code identity, which `dna-eval`'s build
    script records as git object ids. A build from modified read-stage
    sources refuses the cache.
  - **Storage:** entries are gzip-compressed and written atomically; a
    corrupt entry is remade.
  - **Record:** `run.json` gains `evidence_cache` (`dna.eval_run/v1`, in
    place, pre-release).
  - **Facade:** gains `read_evidence_document` and
    `analyze_evidence_outcome`. A sample's rejected-read count now comes from
    its sample evidence, so it does not depend on the path.
  - **Measured:** on 410 samples of panel v5 every table is byte-identical
    from the traces, a cold cache, and a warm one; a warm run takes a fifth
    of the CPU. The read stage is about 90 % of a run's CPU since the tract
    decoder.

- The Sanger plugin decodes a read's homopolymer tracts from its own signal
  (SRS-CALL-025, PROP-0012 increment 2, ADR-0084), as `ReadEvidence::tracts`,
  evidence only. Before each call where callability sees the read leave
  phase, `dna-trace` locates the tract from the read's primaries and chooses
  its structure: a single run, or two runs around a mark of another base.
  It decodes the tract as a mixture of molecules that share one downstream
  sequence, from the analysed channels at each call's locus sample. It
  reports the first run's and the total's length marginals, each with
  standard errors, an autocorrelation-deflated margin, and its state
  (`clear`, `near_equal`, or `undefined`, which is never clear), and refits them with the read's own base signatures and leakage.
  - The decode matches the research reference on 13 vectors.
  - It is published as `read.tracts` in `dna.read_evidence/v1` (changed in
    place, pre-release) and as `tracts.tsv` in `dna.eval_run/v1`.
  - The facade gains `SampleAnalysis::tracts`.
  - The constants are the new `[trace]` section: configuration schema 19,
    plugin `sanger` 17.
  - On development panel v5, every scores, checks, windows, estimates,
    anchors, components, and forms table is byte-identical: no decision
    changes. 67,567 tracts are decoded. At 303-315, HV2F reads choose the
    `T` mark on 95 % of their tracts and HV3R reads the `A` mark on 95 %.
  - A panel run takes about 4.6 times as long (50.6 s to 234 s on the same
    loaded machine).

- `dna sample`, `dna consensus`, and `dna review` take `--evidence
  <read_evidence.json>...` in place of a sample's traces (SRS-API-019,
  PROP-0011 increment 4b): the core places every read from its
  `dna.read_evidence/v1` document in process, and the commands write their
  documents byte-identical to the trace-based commands for the same traces
  and configuration, `dna review --signal` included. Traces and `--evidence`
  exclude each other. On 1,120 real samples (30 picked with every rejection
  kind, the 90 samples of one run reviewed with `--signal`, and 1,000
  random) every sample, consensus, FASTA, and review document is
  byte-identical. `dna.read_evidence/v1` changes in place (pre-release) to
  carry everything the sample documents read: the trace's file name in
  `provenance.input.name`, from which the documents name the read; each
  call's selected window peaks `signal.channels` in place of `height`, which
  the run signal counts; second shares, mixture fractions, and residuals at
  full precision; and the second-peak `cause` read back with the call's
  second base. `dna evidence` now writes a document for a read the sample
  commands set aside before the core, instead of failing: its `rejection`
  (`too_short`, or `callable_calls_below_minimum` with its integrity and
  callability) as the sample documents report it; `dna analyze --evidence`
  fails on it as on its trace. The read evidence reader accepts at most
  16 MiB, as documented (it accepted 64 MiB). The `sample` log's channel
  counts now count positive profile weights (`*_positive_profile_channels`)
  in place of positive corrected amplitudes and SNRs, which the documents do
  not carry. No plugin version changes: no plugin's output changes.
- The evaluation harness's output contract `dna.eval_run/v1`
  (`docs/reference/schemas/eval-run-v1.schema.json`, PROP-0011 increment 3):
  which files `dna-eval run` writes, `run.json`, and every table's columns,
  formats, and meaning, with a synthetic example.
  `validate_result_schemas.py --eval-run <directory>` checks a run, and
  dna-eval's tests check the headers it writes. Contracts only: every output
  is byte-identical.
- `dna-eval run` reports what passes without review, on the current data and
  on every sample its review workload counts: the auto-pass samples (no check
  item) and those with an unreviewed difference; the regions (sample ×
  profile region) and zones (sample × zone) needing a check, overall and per
  region, per zone, and per zone kind; and the passed regions and zones
  holding an unreviewed difference by its first position. Unreviewed
  samples, regions, and zones are split into non-het and
  het only (the truth writes an IUPAC code of two or three bases in every
  such block). `summary.md` gains a "Passes" table in both review workload
  sections, and `run.json` a `workload` object (`all`, and `current` with
  `--current-since`); `dna.eval_run/v1` changes in place (pre-release), and
  `validate_result_schemas.py --eval-run` checks `workload` against
  `scores.tsv`. Every table is byte-identical. On development panel v5,
  current data (10,695 samples): 7,786 auto-pass (72.80 %), 281 of them
  with an unreviewed difference (46 non-het); regions needing a check
  3,496 of 32,085 (10.90 %), 360 passed regions with an unreviewed
  difference (65 non-het); zones needing a check 4,455 of 160,425 (2.78 %),
  and 423 passed zones with an unreviewed difference (74 non-het).
- `dna-eval run --recheck <file.tsv>` reads a list of difference blocks
  whose truth is under re-check (`sample_id`, `first`, `last`, `reason`,
  `source`; a malformed list stops the run), so that the pass measures no
  longer count them as DNA errors. Each unreviewed block now has one class:
  re-check when it overlaps a listed span of its sample, else truth N when
  the truth writes `N` in it, else het, else non-het; a passed sample,
  region, or zone counts under its most severe block (non-het, truth N,
  het, re-check). `run.json` `workload` splits every count into `non_het`,
  `truth_n_only`, `het_only`, and `recheck_only`, and records the list's
  path, SHA-256, and rows (`recheck`); the "Passes" tables gain the same
  columns, the unreviewed blocks on the list, and the list's rows that
  overlap no difference block. `dna.eval_run/v1` changes in place
  (pre-release) and `validate_result_schemas.py --eval-run` checks the split.
  Every table is byte-identical. Without `--recheck`, the non-het counts
  fall: a block with a truth `N` used to count as non-het and is now its
  own class. On development panel v5, current data, the passed zones with
  an unreviewed non-het difference fall from 74 to 43 (31 truth N only), and
  with the 186-row re-check list of 2026-10-20 to 20 (5 truth N only, 49
  re-check only); `window@HV2` keeps 17 of its 19.
- The Sanger plugin asserts a read's poor start as an unreliable span
  (SRS-CALL-024, PROP-0012 increment 5): from its first call to its first
  run of `start_identity_margin` consecutive sound calls. `ReadEvidence`
  gains `unreliable`, published in `dna.read_evidence/v1` (`read.unreliable`)
  and in `dna review` (`reads[].unreliable`). The review labels an item
  covered by such a span with no other read sound there `poor_signal`, and a
  difference there is a check. Plugins `sanger` 16 and `review` 18. On
  development panel v5 against `main` 69e1471, passed zones with a non-het
  unreviewed difference fall 74 -> 72 on current data and 177 -> 175 in all
  (`end@HV3`; no zone rises), and 4 current samples (5 in all) lose their
  auto-pass (73.16 -> 73.12 %).
- `dna analyze --evidence <read_evidence.json> --reference <fasta>` places
  and calls a read from its `dna.read_evidence/v1` document alone, without
  the trace (SRS-API-018, PROP-0011 increment 4a), and writes the analysis
  document byte-identical to `dna analyze <trace>` for the same trace and
  configuration, under the trace's `results/<trace-stem>.json`. The trace and
  `--evidence` exclude each other. The reader refuses, typed and without
  output, another schema version, configuration, or evidence plugin, a call
  out of order, a reason label the Sanger plugin does not raise, and Sanger
  evidence inconsistent with the calls; it rebuilds `ReadEvidence` through
  the kernel's constructors. JSON numbers are now parsed correctly rounded
  (`serde_json` `float_roundtrip`), so written values read back bit for bit.
  On 1,000 real traces, 999 analyses are byte-identical and one read the
  Sanger plugin rejects fails both ways alike.
- The consensus sets aside dissent from a reference decision that a clean
  read of the other strand contradicts (SRS-CONS-011, ADR-0083): opposing
  reads whose second peak stands out and is explained by its modality, or
  whose only difference is an insertion or deletion within
  `consensus.dissent_edge` (40) calls of their own ends. The decided value
  never changes. On panel v5, `strands_disagree` items on the current slice
  fall from 271 to 170 and 110 samples become clear, with unreviewed and
  non-het unchanged in every period. New key `dissent_edge` in
  `[consensus]`; configuration schema 18; plugin `consensus` 13.

- `dna analyze` writes, per call of the selected alignment, where the
  alignment places it and what the core decided about it (PROP-0011,
  increment 2; SRS-OUT-003): a new `calls` array, in call order, with each
  call's `index`, its `column` (1-based reference position and insertion
  ordinal), and `trust` (`placement`, `trusted`, `exempt`, `withheld_by`,
  `scope`), the same trust `dna review` reads. One document then shows the
  chain sanger → align → call for its read; it joins `dna.read_evidence/v1`
  by call index. Every other field is unchanged; component calls of a phase
  mixture (PROP-0010) are pending. **Breaking:** `dna.analysis/v9` → `v10`.
  No plugin version changes: no plugin's output changes.
- `dna evidence <trace>` writes `results/<trace-stem>.evidence.json`, `dna.read_evidence/v1` (SRS-API-017, PROP-0011): a read's per-call evidence exactly as the core receives it from the Sanger plugin, with the signal each verdict rests on, so the modality → core boundary can be checked on its own.
- `docs/architecture/data-flow.md` maps each step to its input and output contracts and the command that writes them (PROP-0011); every crate's source README names its output contract, and `validate_docs_structure.py` fails when one is missing or names a schema that does not exist.
- The Sanger plugin asserts a phase mixture for every dephased segment that
  a repeat run opens (SRS-CALL-023, ADR-0082, PROP-0010 increment 1): the
  stretch fitted by non-negative least squares as the read's own calls plus
  copies offset by whole calls, with each component's offset and fraction,
  the residual, the unexplained calls, and whether the components explain
  the stretch, all against the read's own in-phase noise. The kernel carries
  it as `ReadEvidence::mixtures`. `dna-call`'s `component_calls` compares
  each component of an explained mixture with the reference on its own
  (SRS-VAR-022), exposed as `SampleAnalysis::components`, and `dna-eval run`
  writes `components.tsv`. New key `decompose_unexplained_quantile` (0.99)
  in `[signal_processing]`; configuration schema 17; plugin `sanger` 13.
  No call, consensus, or review output changes. On panel v5's current
  slice, an explained component's entered run in the HVS-II C-tract is the
  truth's in 95.6 % of HV2F and 99.3 % of HV3R reads.

- Each Sanger call with a second peak above `signal_processing.second_floor`,
  or masked as mixed, is labelled with the peak's cause from its shape
  against the call's neighbours (SRS-CALL-021): `noise`, `neighbour_tail`,
  `off_centre`, `dye_blob`, `phase_shift`, `second_template`, or
  `mixed_stretch`. With a run signal, a `second_template` within its key's
  run level is `context`. The labels sit beside the core's records on the
  read's Sanger attachment and decide nothing. New keys
  `second_cause_noise_ratio` (3.0), `second_cause_tail_offset` (0.4),
  `second_cause_blob_width` (1.8), and `second_cause_width_calls` (6) in
  `[signal_processing]`; configuration schema 16. No output changes.

- `dna-eval run` counts the review workload in zones, the functional parts
  of a region a reviewer opens as a whole: repeat window, after a window,
  amplicon start and end, and interior, from the configuration's target
  profile and `--zone-after`, `--zone-edge`, and `--zone-overrun` (recorded
  in `run.json`). `scores.tsv` gains `review_zones` and
  `review_zone_count`; the summary reports zones per sample, the samples
  with none, the samples with a check in each zone, and the unreviewed
  blocks by zone kind.

- `dna review --signal <run.signal.json>` (PROP-0009, increment ii; SRS-API-014,
  SRS-CALL-020): `possible_heteroplasmy` is judged against the sample's run
  signal, with the sample's own reads taken out. A run signal holding none of
  the sample's traces fails the command. Provenance records the run signal's
  `run_id`, SHA-256, and selection; each qualifying read carries its
  `threshold`, and `calls` when its key's run level was usable.
- `dna-eval run --run-signal batch|month [--run-signal-batch-reads <n>]`
  builds run signals in process (per manifest batch with enough reads, else
  per batch month) and reviews every sample against its group's. Every run
  now also writes `checks.tsv` (each check item with its labels), and with
  `--run-signal`, `signals.tsv` (per group: reads, calls, keys, keys with at
  least 20, 40, and 100 calls).

- The run signal (PROP-0009, increment i; SRS-CALL-019, SRS-API-016):
  `dna signal <run-id> <trace>...` writes `results/<run-id>.signal.json`, a
  new contract `dna.signal/v1`. Per read-strand 5-mer, call base, and other
  channel, it records how often each share of that channel occurs among the
  run's selected calls (every call counts under its three other channels), with the run's trace hashes and the selection version. The
  selection is one function (`dna_sanger::run_signal::contributions`) for
  building and for taking a read out. No command reads a run signal yet, so
  no other output changes.

- `dna-eval run --current-since <YYYY-MM>` reports the current data first:
  the review workload of the samples whose batch date is in or after that
  month, and with a baseline their changes of status class per measure.

- The review summary counts positions (`review` plugin 4, rules version 4;
  `dna.review/v1` revised in place with one summary object):
  `summary.positions.reportable`, the profile's reportable positions, and
  `summary.positions.check`, those a `check` item holds. A review tool can
  show a sample as the few positions to look at and leave the rest as
  decided. On development panel v5 a sample's median is 2 positions to check,
  of 749.
- The evaluation harness measures the review workload: per sample its
  check items, positions to review in windows of `--review-window`
  positions (default 10), whether it is clear, and the `haplotype_lod`
  difference blocks no check item overlaps (`unreviewed`); `windows.tsv`
  ranks the windows by unreviewed differences, and `--forced-windows` reports
  the workload again with a list of windows always reviewed. The review
  items are `dna review`'s own, through the new
  `SampleAnalysis::review_checks` (SRS-API-015). Existing columns and
  outputs are unchanged; the new columns follow `failure`.

- Two review labels (`review` plugin 2, rules version 2; `dna.review/v1`
  revised in place: two label kinds and two summary counts, and `rules` is
  any version from 1). `read_disagrees` makes an item a check when a read
  the consensus counted for its decision, or set aside, still shows the whole
  item as another sequence with every call free of the modality's reasons.
  `not_called` makes an item a check when the consensus decided a difference
  that DNA does not call within one position (an excluded variant). On
  development panel v5 (16,413 samples): `clear` 6,460 → 6,424, leaks among
  clear 238 → 234; `not_called` on 287 samples, mostly `315.1C` excluded at
  311–316.
- The review list (PROP-0004, SRS-API-014): `dna review <sample-id>
  <trace>... --reference <fasta>` runs the sample analysis, the consensus,
  and its call once and writes `results/<sample-id>.review.json`, a new
  contract `dna.review/v1`. It lists the positions a reviewer should look at,
  in priority order, each with what DNA decided, its calls, evidence labels
  (`strands_disagree`, `possible_heteroplasmy`, `estimated_length`,
  `undecided`, `one_read`, `one_strand`, `noisy_signal`, `not_covered`,
  `phase_shift`, `read_edge`), and every covering read's bases, call indices,
  trace scan positions, relative qualities, and second-base shares; a summary
  counts the labels and says whether the sample is `clear`. Every sample gets
  a document: one whose reads are all rejected, or whose call fails, carries
  that status. New post-calling plugin `review` 1 (labelling rules version 1)
  and kernel contract `ReviewList`. Watched positions come later.

- The evaluation harness reports DNA's callable yield: per sample and HV
  region, the truth's analysed positions that DNA's callable intervals cover
  (`callable_yield` in `scores.tsv`), and per region in the summary the
  positions, their share, and the samples short of full yield. It feeds no
  status (PROP-0008 step 0), so a change that makes calls uncallable shows
  what it costs.

- The EMPOP query (SRS-NOM-018, PROP-0007 increment 2): when the profile
  declares `[empop]`, as the human-mtDNA profile does, `dna notation` adds an
  offline `empop` query to `dna.notation/v1`, and `SampleAnalysis::call`
  gives it as `ConsensusCall::empop`. It is sent nowhere. The query has
  interpretation ranges over the positions exactly one record calls (no `N`,
  no excluded variant), the represented calls inside them with an insertion
  written by its anchor, as EMPOP numbers it, and the record positions cut,
  with their reason. New post-calling plugin `empop` 1 in `dna-empop`.
  **Breaking:** `dna.notation/v1` adds the optional `empop`, revised in place
  (unreleased); the profile gains `[empop]`. The evaluation harness writes
  each query to `scores.tsv` (`empop_ranges`, `empop_mitotype`), with the
  positions inside its ranges that the consensus left undecided or decided by
  a tie (`empop_undecided`, `empop_ties`), and summarises them.

- Every `technology` finding of the evaluation harness's discrepancy causes
  is deferred, to signal processing in `dna-sanger`: the core is judged on
  good input, or input the lab made poor. A defect heads a sample before a
  deferred item, so a sample with both a `lab` and a `technology` finding is
  the lab's. The summary's cause table lists defects first, then deferred
  items, each in cause order (`technology (deferred)` replaces
  `technology`). Statuses and classes are unchanged.

- Report conventions (SRS-NOM-017, PROP-0007 increment 1): a profile's
  `[[nomenclature.conventions]]` labels the rendered calls the reviewers
  write differently, and keeps them. The human-mtDNA profile declares
  `hvs1_c_length` (`not_recorded`): the length calls of the uninterrupted
  HVS-I C run that `16189C` creates (C10 and longer). Every notation view
  carries the label as `convention` on a labelled call, and
  `SampleAnalysis::call` carries it on `NotationCall`. The `nomenclature`
  plugin is version 3. **Breaking:** `dna.variants/v1`,
  `dna.sample_evidence/v10`, and `dna.notation/v1` add the optional
  `notation.calls[].convention`, revised in place (unreleased). The
  evaluation harness gives a `sequence` difference made only of such calls,
  with no truth-only call, the cause `lab`, deferred
  (`convention_not_recorded`): 31 samples on development panel v5, all
  previously `undetermined`. The summary's cause table adds a
  `lab (deferred)` column.

- The evaluation harness's haplotype measures compare only inside the
  truth's analysed intervals: a DNA call where the reviewers did not analyse
  no longer counts against DNA there (the primary, strict, and
  IUPAC-compatible measures, the lab's rules, are unchanged). A sample whose
  only primary differences are such calls takes the new class `unanalysed`,
  cause `undetermined` (`outside_truth_intervals`)
  ([measures](docs/validation/evaluation-harness.md#measures)).
- `local_shadow` minor kind and minor rule sets in the draft triage
  (evaluation harness only, #70): a read's neighbour-compatible minor channel
  counts as a shadow only where at least 0.05 of the read's calls within 10
  show the same shadow, and a policy's `minor` is a set of rules any of which
  marks (`upstream@100+local_shadow@80`); the pre-registration parses
  `+`-joined rules.
- The `haplotype_lod` measure (evaluation harness only, #71): the haplotype
  measure, except that a truth two-base IUPAC substitution the reads do not
  observe at the 10 % limit of detection compares IUPAC-compatibly.
  `scores.tsv` adds its status, code, and blocks and `below_lod`; it is the
  approval safety endpoint
  ([measures](docs/validation/evaluation-harness.md#measures)).
- Readiness panels (evaluation harness only, #72): `dna-eval run` reads a
  `readiness` panel as well as `dev`, and `power` reads a development or
  readiness run; any other (hold-out) panel is still refused.

- Discrepancy causes in the evaluation harness: `dna-eval run` gives each
  sample with a discrepancy a cause, `dna`, `lab`, `technology`, or
  `undetermined`, a scope (`defect`, or `deferred` for IUPAC codes DNA does
  not write by design), and its findings (`cause`, `cause_scope`,
  `cause_evidence` in `scores.tsv`), read from the reads' evidence at each position, the
  limit-of-detection signals, and how the reviewers write each window
  sequence over the whole truth. The summary counts causes per class and
  samples per cause and evidence
  ([discrepancy causes](docs/validation/evaluation-harness.md#discrepancy-causes)).
  A sample whose every trace DNA rejects now fails as `rejected`. A
  coverage gap read cleanly yet untrusted is placed by where the call
  stands, with the margins and thresholds read from DNA's configuration:
  beside phase loss or a read end DNA attributes to a repeat run
  (technology), at the read's own end or with a weak support veto (lab), a
  weak veto on a call that reads the reference (DNA, `veto_out_of_scope`),
  mixed channels (undetermined), or otherwise DNA's callability policy; a
  truth IUPAC code a read shows is judged as under heteroplasmy. Each run
  writes the written-form index as `forms.tsv`, keyed by the truth
  database's SHA-256, and `--forms` reads an earlier run's instead of the
  whole database.
- `dna-eval convention` (evaluation harness only): runs every IUPAC-free
  truth sample's calls through DNA's normalization and nomenclature windows
  and counts the samples that come back unchanged, with the change patterns
  of the rest ([notation convention](docs/validation/evaluation-harness.md#notation-convention)).
- Discrepancy classes in the evaluation harness: `dna-eval run` gives every
  sample one class from its statuses (`failed`, `sequence`, `heteroplasmy`,
  `representation`, `coverage`, `not_comparable`, or `none`), naming the part of DNA that owns
  its difference from the truth, with the windows involved (`class`,
  `class_regions` in `scores.tsv`). The summary counts the classes per window
  and, against a baseline, the samples whose class changed
  ([discrepancy classes](docs/validation/evaluation-harness.md#discrepancy-classes)).
- One-strand exception policies in the draft triage (evaluation harness
  only): `region.one_strand` can apply everywhere, nowhere, or only inside
  the repeat windows and near reference homopolymers (R6m, R6m-h4), and a
  single-strand minor-signal rule applies where one strand may pass. The
  approval pre-registration accepts them, and the marginal endpoint reports
  apart the leaks inside the exception.
- `dna-eval approve` and `dna-eval power` (evaluation harness only). `approve`
  runs the pre-registered approval of triage pass rules on a hold-out panel
  once: it checks every frozen input before DNA runs, decides the fixed
  sequence of hypotheses, and bounds each gate's leak rate three ways
  (count-based Clopper–Pearson, Korn–Graubard, batch-cluster bootstrap).
  `power` projects the passes and attainable bounds from a development run
  without the hold-out truth. `dna-eval run` now refuses a hold-out panel
  ([approval](docs/validation/evaluation-harness.md#approval)).
- Haplotype measures in the evaluation harness (PROP-0005).
  - **Harness:** `dna-eval run` also scores every sample under `haplotype`
    and `haplotype_iupac_compatible`: the primary and IUPAC-compatible rules
    applied only to alleles that do not describe the same sequence, so a
    haplotype written or placed differently (`16190.1C` against `16193.1C`,
    `248DEL` against `249DEL`) is not a difference. Scores list the blocks
    where the sequences differ; estimates, anchors, the baseline diff, the
    summary, and the draft triage report both measures
    ([measures](docs/validation/evaluation-harness.md#measures)). The primary,
    strict, and IUPAC-compatible measures are unchanged. The harness now
    loads the reference it passes to DNA, and a sample DNA analyses on another
    reference sequence fails (`reference`).
- Draft triage measurement for PROP-0004 phase 0 (evaluation harness only).
  - **Harness:** `dna-eval run --triage` takes a compact triage evidence
    record per sample from the same analysis, runs draft stage-1 detectors
    (`trace.too_short`, `trace.weak`, `trace.unplaced`, `trace.short`,
    `region.unreliable`, `region.one_strand`, `region.gap`), decides pass,
    review, or rerun per region, trace, and sample, and scores A/B/C, leaks
    with Clopper–Pearson bounds, catches per detector, and a policy sweep
    against the truth (`triage*` outputs,
    [draft triage](docs/validation/evaluation-harness.md#draft-triage)).
    Without `--triage` a run is unchanged. The detectors are research drafts,
    not DNA rules; their results are in
    [PROP-0004](docs/proposals/0004-review-triage.md#phase-0-results).
  - **Rust API:** `SampleAnalysis::regions()` returns the reportable regions
    of the target profile the analysis ran with (SRS-API-010).
  - **Harness:** draft `region.minor_signal` detectors for point
    heteroplasmy over each clean read's minor channel at every covered
    position: cross-strand rules with and without a shadow test, a noise
    ratio, single-read and one-strand variants, swept at several thresholds
    and measured against the same policy without them (`triage_minor.tsv`).
    The draft applies `upstream@100`; results are in
    [PROP-0004](docs/proposals/0004-review-triage.md#minor-signal-results).
  - **Rust API:** `SampleAnalysis::observations(first_1based, last_1based)`
    gives every admitted read's observation of each covered position of an
    interval, as the sample evidence records it at a differential locus,
    whether or not a read differs there, without support aggregates
    (SRS-API-013). A zero or reversed interval is a typed
    `SampleError::InvalidInterval`; a span across the circular origin is two
    calls. `Orientation::reference_base` (`dna-core`) is public.
  - **Rust API:** `sample_analysis::analyze_sanger_outcome` runs the same pass
    as `analyze_sanger` but returns `SampleOutcome::Rejected` with every
    rejected read and its cause when none is admitted, instead of the
    `NoAdmittedReads` error, so a caller can name the traces to rerun
    (SRS-API-009).
- Sample analysis in one placement pass (PROP-0005, SRS-API-009/010/011).
  - **Rust API:** `dna::sample_analysis::analyze_sanger(traces, reference,
    config, workers)` places each trace once on up to `workers` threads and returns a typed `SampleAnalysis`: the
    core sample evidence, the consensus, each admitted read's Sanger evidence,
    and each rejected read with its cause. A caller that needs both no longer
    places every trace twice.
  - `dna sample` and `dna consensus` run the same pass; their documents and
    FASTA are unchanged.

- Evaluation harness `dna-eval` (PROP-0005 harness v0, ADR-0075,
  SRS-COMPAT-003).
  - **Crate:** `crates/dna-eval`, a workspace binary outside the `dna` binary
    and the release. It is the only place for truth formats and comparison;
    no crate depends on it.
  - **Run:** it verifies a panel manifest, runs `consensus` then `call` for
    every sample in process on a thread pool, reads the reviewer truth from
    the lab's statistics database (read-only, `rusqlite`), and classifies each
    sample `Y`, `N`, or `N/A` with the lab's comparison rules, plus a strict
    measure that counts DNA's coverage gaps and an IUPAC-compatible measure.
  - **Panel versions:** `dna-eval panel` derives a new panel version against
    the database, keeping every row, retiring samples without truth, and
    recording per-sample and panel truth digests; a run reports the samples
    whose truth changed since the panel was built.
  - **Outputs:** per-sample comparisons, stratified weighted estimates with
    95 % intervals, anchors, the diff against a baseline run, a Markdown
    summary, and a run header
    ([evaluation harness](docs/validation/evaluation-harness.md)).
  - **Rust API:** `SampleAnalysis::call(sample_id)` calls a sample
    analysis's consensus as `dna call` reads the consensus FASTA, completing
    the recommended `consensus` then `call` workflow in process. It returns a
    `ConsensusCall`: the callable segments and structured `NotationCall`s
    (SRS-API-011, SRS-API-012).
- Sample consensus (PROP-0003 increment 1, ADR-0073, SRS-CONS).
  - **Command:** `dna consensus <sample-id> <trace>... --reference
    <reference.fasta>` adjudicates between a sample's placed reads.
  - **Outputs:** `results/<sample-id>.consensus.json` (`dna.consensus/v1`:
    segments, notable sites with their supporting, opposing and
    uninformative reads, and decision counts) and
    `results/<sample-id>.consensus.fasta`. The FASTA is ready for `dna call`.
  - **Rule:** only clean observations decide. Differences are decided as
    stretches of whole reference runs (as a whole, run by run, or position by
    position). One read against one read takes the reference and is marked
    `contested`.
  - **Measured on two disjoint 320-sample held-out subsets** (fair scale):
    exact samples rise from 288 to 290 and from 290 to 292 against the union
    of reads, and precision from 0.9924 to 0.9967 and from 0.9934 to 0.9964.
  - `dna sample` and `dna call` are unchanged.
  - **Rust API:** `dna_core::consensus` and `dna_core::model::consensus` are
    new; `dna_kernel::plugin::Contract` gains `Consensus`.

- Post-calling notation of a variants document (ADR-0069 phase 4, PROP-0002,
  SRS-IN-014, SRS-OUT-011, SRS-NOM-016, SRS-PRF-008). The command is
  `dna notation <sample-id> <variants.json> --reference <reference.fasta>`.
  - It reads a `dna.variants/v1` document and runs normalization and
    nomenclature over each read's eligible variants.
  - It writes `results/<sample-id>.notation.json` (`dna.notation/v1`): the
    per-read notation (equal to the in-process notation), the source
    document's identity, and the findings of the new `conformance` plugin.
  - That plugin reports represented calls that break notation conventions
    declared in a new profile `[conformance]` section
    (`insertion_at_run_end`, `insertion_matches_run`); it never changes a
    call.
  - The shipped human-mtDNA profile declares both rules, so its file SHA-256
    recorded in result provenance changes.
  - **Rust API:** `Error` gains `Variants(VariantsError)` and
    `VariantsParse`.
- Core-only calls from reviewed consensus sequences (ADR-0069 phase 3,
  PROP-0002, SRS-IN-013, SRS-OUT-010, SRS-VAR-015). The command is
  `dna call <sample-id> <sequences.fasta>... --reference <reference.fasta>`.
  - Each FASTA record becomes modality-neutral `ReadEvidence` through the new
    `sequence` modality plugin: one-hot profiles, shared-weight profiles for
    IUPAC codes, no masks or vetoes, and vouched read ends.
  - The core caller alone places each record and calls its variants into
    `results/<sample-id>.variants.json` (`dna.variants/v1`), with per-read
    alignment summaries, observed variants with eligibility, and optional
    notation.
  - Sanger results are unchanged.
  - **Rust API:** `Error` gains `Sequence(SequenceError)`, and
    `analyze_sanger` opens a `read_evidence` stage span, for the Sanger
    evidence adapter, between `quality_control` and `alignment`.
- Signal-derived read callability (ADR-0067, SRS-CALL-001 to SRS-CALL-008): a
  new `callability` stage derives, from each read's own signal in trace order,
  phase-state segments (`in_phase`, `dephased`, `mixed`, `weak`, `irregular`),
  a per-position mask, and the callable span, and publishes them as
  observation. Trim bounds, alignment, warnings, and variant eligibility are
  unchanged in this increment. **Breaking:** `dna.basecalls/v2` → `v3`,
  `dna.analysis/v8` → `v9`, and `dna.sample_evidence/v9` → `v10` gain the
  required `read.callability` / `reads[].callability` view; configuration
  schema 6 → 7 adds `[callability]` (`window_calls`, `onset_defect_fraction`,
  `exit_defect_fraction`, `minimum_main_share`, `maximum_far_share`,
  `minimum_shadow_share`, `weak_amplitude_fraction`) and
  renames `quality_control.trim_window_size` to `penalty_window_size`.
  **Breaking (Rust API):** `Error` gains `Callability(CallabilityError)`, and
  `analyze_sanger` opens a `callability` stage span between
  `signal_processing` and `quality_control`.
- Target profiles (ADR-0063, SRS-PRF-001 to SRS-PRF-007): target knowledge —
  reference identity and topology, reportable regions, indel placement,
  nomenclature windows and their rules, notation style — moves out of code and
  configuration into strict versioned TOML profiles. The shipped profile is
  `config/profiles/human-mtdna-rcrs.toml`; a profile that pins a reference fails
  closed on another sequence. Scientific output for human mtDNA is unchanged.
  **Breaking:** configuration schema 6 gains a root `profile` key and loses
  `[reference]` and `variant_calling.regions`; `dna.analysis/v8` and
  `dna.sample_evidence/v9` record `provenance.profile` (`id`, file `sha256`);
  sample `notation.policy` becomes `notation.style = "per_base_decimal"`;
  the variant exclusion reason `outside_configured_region` becomes
  `outside_target_region`.
  **Breaking (Rust API):** new `dna::profile::{Profile, ProfileIdentity}`;
  `variant_nomenclature::mtdna::apply_control_region` is replaced by
  `variant_nomenclature::apply(reference, &profile, input)`;
  `NormalizationPolicy::MtDnaRightAligned` is renamed `RightAligned`;
  `VariantAnalysisResult` gains `profile`; `Error` gains `Profile` and
  `ProfileParse`; `ConfigError::RegionOutOfBounds` moves to `ProfileError`;
  `NomenclatureError` window names are `String`s.
- Read callability (ADR-0062, SRS-VAR-013): variants with a call near either
  read end or right after a long homopolymer are ineligible with reasons
  `read_end` / `post_homopolymer`. **Breaking:** configuration schema 6 adds
  `variant_calling.read_end_margin`, `homopolymer_min_length`, and
  `post_homopolymer_window`.
- `dna.sample_evidence/v9`: when the target profile declares notation, an
  optional `notation` view publishes each read's eligible calls after the
  profile's right alignment and nomenclature windows (HVS-II, HVS-III, HVS-I for
  human mtDNA), rendered per base (`73G`, `249DEL`, `309.1C`) with supporting
  reads (SRS-NOM-010 to SRS-NOM-012). All other v8 fields are unchanged.
- Production-oriented CI/CD and supply-chain verification, including strict cargo-shear dependency/source hygiene and an explicit Ubuntu 24.04 runner baseline.
- Declare the crate proprietary (`license = "LicenseRef-Proprietary"`) so SBOMs
  and audits record its ownership.
- Crate-wide Rust lint policy in `Cargo.toml` (`missing_docs`, `unreachable_pub`,
  Clippy `pedantic`) and a tuned release profile (thin LTO, one codegen unit).
- Rust/Python source-boundary enforcement.
- Dependency policy, dependency review, CodeQL, OpenSSF Scorecard, fuzzing,
  self-contained explicit-target auditable release bundles with post-strip metadata verification, release SBOMs, and artifact attestations.

### Changed

- The review shows no length of DNA's for a run of an anchored repeat
  window it does not report (a window position not callable): such a run's
  `length_unconfirmed` (or `length_not_recorded`) `called` is null, the
  estimate is not judged against the reads (a run both strands read alike is
  `undecided_call`, no longer confirmed by the estimate), and it stays only in
  `estimated_length` with its `length_evidence`. The run carries `floor`
  instead, the reads' floor reading: the least length any read gives as its
  top or next, a suggestion beside the near-equal pair's `convention`
  (SRS-API-020, user decision 2026-10-21). On development panel v5 the
  floor matched the truth in 354 of 406 current unreported HVS-II run 1s
  whose reads top at different lengths, the phase-loss estimate in 41.
  `dna.review/v1` gains the field in place (pre-release); plugin `review` 27.
  Against `main` 60252579 (before #204) the check items keep their positions (25 HVS-II
  `estimated_length` items also carry `length_unconfirmed`), and the scores,
  unreviewed blocks, and clear samples are byte-identical (diff gate pass).

- The review confirms an anchored window run that only one strand gives a
  length when that strand gives one clear length, DNA's call, and at least
  one of its reads shows the run in phase (SRS-API-020, ADR-0074): an
  in-phase read shows the length directly, and the other strand having no
  view is absence of evidence, not a contradiction. A run that only decodes
  give still needs both strands, and a read whose structure differs from
  DNA's still never confirms. Plugin `review` 26 (rules version 26).
  - **Measured.** Against main 3cad890 (`dna-eval diff`, frozen re-check
    list). Development panel v5, current data: window@HV2 checks
    1,974 → 1,664, clear samples 7,728 → 7,824 of 10,695, no zone or
    sample with a non-het block newly passed; GATE PASS. Readiness-v2 (now
    development data), current data: window@HV2 checks 2,900 → 2,483, clear
    13,383 → 13,511 of 17,307, passed non-het samples 11 → 11 and passed
    window@HV2 zones with a non-het block 4 → 4. The six current samples
    that newly pass with an unreviewed block that is not het hold only
    `573.1C 573.2C` past the truth's HVS-III interval (truth N, #211); their
    window@HV2 runs, now confirmed, match the truth.
  - **Earlier data.** Readiness-v2 gains two window@HV2 zones with a
    non-het block, both in one 2024-11 batch: DNA's `309.1C`, which the
    truth lacks, is shown by eight clean, evenly spaced C peaks on HVS-II
    forward (and on HVS-III reverse in one), so the truth is the likely
    error; the diff gate fails on them.

- `dna-eval` classes an unreviewed block whose every difference is an
  insertion past the truth's analysed interval (anchored at the last base of
  the truth's own interval, e.g. DNA's `573.1C 573.2C` where the truth's
  HVS-III interval ends at 573) as `truth_n`, not `non_het`: the truth records
  nothing there, and the allele measure already leaves it uncompared (user
  decision 2026-10-21). A block that ends only where DNA's callable interval
  ends stays `non_het`. On development panel v5 against main 60252579, 27
  blocks move (25 at 573, 2 at 309); non-het unreviewed blocks fall 71 → 44
  (current 34 → 20), passed samples with a non-het block 24 → 12 (current
  8 → 2), and `end@HV3` passed non-het zones 24 → 1. DNA's output is
  unchanged.
- A Sanger call before a run that closes a stop of the read's in-phase
  signal (the informative end, or a phase loss after a repeat run inside the
  interval) is judged against the medians of the last `end_norm_calls` (30)
  scored single-channel calls before that stop, with a peak of at least
  `end_height_fraction` (0.3) of theirs, when SRS-CALL-017 decides whether
  its own signal takes its base off the signal edge. A read's signal decays
  over its last calls and its peaks vary with their dye context, so the
  read's own median (`edge_height_fraction` 1.0, kept at the informative
  ends) held back clean calls. Such a call reads its outer neighbour's base
  only inside a run that a vouched call of another base closes. The reverse
  read that ends a few calls past the HVS-I C-tract, or regains phase after
  it, now reads 16194-16205. A call further in than the run's bound (the
  call right before the run) qualifies only behind a vouched bound: an
  ambiguous bound leaves the read's phase at the run in doubt, and without
  it the consensus completed `566-573` from the reverse HVS-III read's start
  (`567C` in five samples). On development panel v5 against `main` 85cab63
  the diff gate passes: unreviewed samples (618) and passed samples with a
  non-het block unchanged in every period and source, current `after@HV1`
  checks 254 → 88, 131 more samples clear (current 71.49 → 72.13 %), and
  positions to review 164,199 → 158,668. Under the haplotype
  IUPAC-compatible measure 3 samples turn discordant and 2 concordant; the
  three new calls (`573.1C`-`573.3C` twice, where both strands show the
  longer run past the truth's HVS-III interval end at `573`, and `309DEL`)
  lie inside check items. **Breaking:** configuration schema 23 requires
  `sanger_evidence.end_norm_calls` and `end_height_fraction`;
  plugin `sanger` 22.
- A junction whose spacing marks an extra peak, a peak the basecaller called
  twice, now withholds the two calls beside it like a missing peak's
  junction does (SRS-CALL-026, plugin `sanger` 21): they place no indel and
  show no run's end. It is the mirror of the missing-peak test on the same
  measure and level: narrower than the local median by
  `callability.spacing_defect` of it, and still so together with either
  neighbouring junction, so a peak merely placed early is no extra peak. No
  new configuration. Such a doubled call (4-6 scans after its neighbour
  against a spacing of about 12) gave one-read insertions near read starts
  that passed unflagged: `16049.1G`, `16025.1T`, `16042.1G`, `16060.1G`,
  `16360.1C`, `560.1C`, `87.1A`. On development panel v5 against the
  previous commit the diff gate passes: samples passed with a non-het
  difference 23 → 19, passed non-het zones `start@HV1` 6 → 1, `end@HV1`
  2 → 1, `end@HV3` 12 → 11 (one new: a one-read `563.1A` at an HVS-III
  reverse read's irregular start); zones needing a check fall or hold
  everywhere but `end@HV1` (+1); 9 more samples clear; primary concordance
  14 samples better, 1 worse.
- The review labels `reference_second` a position where a difference rests
  on one read whose call there reads another base than the reference while
  showing the reference base as its second peak above the call's
  second-channel threshold (SRS-CALL-020): a check, decision unchanged
  (rules version 25, plugin `review` 25). No new configuration. On
  development panel v5 such bases were mostly dye blobs over the
  reference's peak at reverse read starts (`527A`, `530G`, `16326G`,
  `16328G`). Against 85cab63 the diff gate passes: samples passed with a
  non-het difference 18 → 15 (`after@HV3` 2 → 0, `end@HV1` 2 → 1, all
  earlier data), het-only 340 → 338, unreviewed samples 618 → 581, no status
  or coverage changes; 25 more zones need a check (13 current) and 6 samples
  lose their clear state (2 current), one of them for a correct call.
- The review's `length_unconfirmed` runs carry `convention` (null in
  `length_not_recorded`): when both strands read a run as a near-equal adjacent pair
  `n`, `n + 1`, the lab's reading `n` (the shorter member, the run ending at
  the reverse read's first mixed call), as a suggestion beside DNA's call,
  never applied (SRS-API-020, user decision 2026-10-20). No decision, label,
  or kind changes. `dna.review/v1` gains the field in place (pre-release);
  plugin `review` 24.
- The consensus decides a run it would leave as an estimate
  (`anchored_end`, `phase_loss`, or `reference_frame`) when the reads' tract
  decodes on both strands put one length on top beyond the strands' spread
  and replicate noise (SRS-CONS-013, ADR-0074 amendment, PROP-0012
  increment 4). This applies to either of two runs of one base around a
  single-base anchor, in a stretch whose call keeps that reference
  composition. The window comes from the stretch's own structure, not from
  the profile.
  - **Evidence.** Such a run carries `length_evidence` `two_strand`
    (`dna.consensus/v1`, in place, pre-release). It is observed, so its 3'
    end stays callable. An `in_phase` or `confirmed_end` length is never
    replaced.
  - **Review.** A `two_strand` run counts as confirmed in
    `length_unconfirmed`, unless a read's structure differs from DNA's.
  - **Configuration.** The constants are new `[consensus]` keys
    `length_confidence` (2.326), `length_noise_forward` ([0.065, 0.111]),
    and `length_noise_reverse` ([0.111, 0.128]): replicate noise per strand
    and run, from 743 re-sequenced samples
    (`docs/design/length-mixture.md`). Configuration schema 22 (21 is taken by `signal_spike`, merged first).
  - **Shared code.** The non-negative least squares solver moves from
    `dna-trace` to `dna_kernel::nnls`, unchanged.
  - **Measured.** On development panel v5 against increment 3
    (`dna-eval diff`, GATE PASS), the effect is near zero:
    - one call changes (a C8 at 303-309 where truth has `N`), and no error
      is introduced;
    - window@HV2 checks go from 3,073 to 3,072;
    - 31 more run 1s become `two_strand` at the length DNA already had, but
      none of their checks goes, because run 2 stays one-strand there.
    Most flagged run 1s are a real C8/C9 mix (61 of 69 near-equal on both
    strands), which the rule refuses.
  - Plugins `consensus` 18 and `review` 23 (rules 23).

- A Sanger call inside an electrophoretic spike shows no base (SRS-SIG-012,
  SRS-CALL-027, ADR-0088, plugin `sanger` 20). A spike is a narrow event in
  all four analysed channels at once, each far above the read's signal and
  staggered by the dyes' mobility offsets; the call whose primary peak lies
  in it reads the spike's channel. Signal processing finds the spikes from
  the read's own channels and calls; the Sanger evidence withholds their
  calls for the base and for indels (`signal_spike`), with no run's end and
  no reference confirmation. New keys `spike_height` (3.0), `spike_span`
  (1.5), and `spike_level_calls` (6) in `[signal_processing]`.
  **Breaking:** configuration schema 21; `signal_spike` joins the
  `dna.sample_evidence/v10` exclusion reasons (in place, pre-release). On
  development panel v5 against `main` 85cab63 the diff gate passes: the one
  unflagged non-het error it targets (`186T` in 2538054, `interior@HV2`)
  becomes an undecided check, samples with an unreviewed block 618 → 617,
  passed non-het units unchanged in every other zone and period, five more
  samples clear (current 7,646 → 7,651), positions to review 164,199 →
  164,139, no coverage change, and one sample concordant under every
  measure.
- A read whose two best placements tie only on calls at one end that the
  read itself asserts unreliable (inside an unreliable span, or masked) is
  no longer rejected as `ambiguous_placement`: those calls are left unplaced
  and the rest of the read is aligned again (SRS-ALN-006, ADR-0086, plugin
  `core` 12). A tie that disputes any other call is still rejected. On
  development panel v5 all 31 such reads are placed (23 current: HV1R and
  HV3R poor starts, two dephased C-run ends). Against `main` 85cab63 the
  diff gate passes: non-het passed zones and unreviewed samples unchanged in
  every period, 11 more samples clear (current 7,646 → 7,654), positions to
  review 164,199 → 162,335, and 20 truth alleles newly covered, 19 called
  as the truth holds; `16192Y` in one sample is called `16192T`, inside a
  check.
- `dna-eval diff` gates on non-het blocks only: more samples with a
  non-het unreviewed block fail it, while samples whose only unreviewed
  blocks are het are reported in their own column and never gate
  (heteroplasmy is optional, user decision). The `Unreviewed` column splits
  into `Unreviewed non-het` and `Unreviewed het only`.
- A mixed Sanger call whose second peak is a neighbour's tail or a dye blob
  may confirm the reference beside a neighbour that is itself mixed, when
  that neighbour's second peak is also a tail or a blob: a blob is wider than
  one call (SRS-CALL-022, ADR-0080 amended, plugin `sanger` 19). Every core
  guard of a confirmation is unchanged. On development panel v5 against the
  previous commit the diff gate passes: unreviewed and non-het passed zones
  unchanged, no status changes, 84 more samples clear (current 76.19 → 76.51
  %, earlier 64.78 → 65.63 %), and positions to review 141,588 → 136,579.
- The Sanger evidence vouches for the two calls beside a junction whose
  spacing marks a missing peak by their bases only (SRS-CALL-026, plugin
  `sanger` 18): they place no indel and show no run's end. Callability's
  spacing measure finds the junctions: wider than the local median by
  `callability.spacing_defect` of it, and still so together with either
  neighbouring junction, so a peak merely placed late is no gap. No new
  configuration. A one-read deletion at such a junction, opposed by a clean
  read of the other strand, had been a check; on development panel v5 every
  one was truth reference. Against the previous commit the diff gate
  passes: unreviewed 668 → 667, non-het passed zones unchanged, clear
  72.08 → 72.21 % (current 76.06 → 76.19 %), current `interior@HV1` checks
  311 → 303, and four samples become concordant under every measure.
- The review's `edge_deletion` label reads its own zone,
  `review.deletion_zone` (default 20, its width before ADR-0085), instead of
  `review.edge_zone`, which keeps 40 for `beside_not_analysed` (ADR-0085,
  amended). The wider zone had flagged the common true `106-111DEL`, read by
  one HVS-II forward read. On development panel v5 against `feat/eval-diff`
  ae98063 the diff gate passes: unreviewed and non-het passed zones
  unchanged everywhere, current `start@HV2` checks 285 → 245 and
  `interior@HV2` 257 → 221, clear 75.99 → 76.06 % on current data and 12
  samples overall. Configuration schema 20; plugin `review` 22.
- The review flags the runs of an anchored repeat window that the reads'
  tract decodes do not confirm (SRS-API-020, PROP-0012 increment 3). This is
  flag-only: no decision changes.
  - **Credit.** Each read's decode is credited to the run its first call is
    aligned to. A read gives the length of the run it enters first; with
    the anchor as its mark and clear margins, it also gives the next run's
    length, from its total. A read the consensus reads in phase through a
    run gives that length too (ADR-0074, `dna_consensus::read_runs_over`).
  - **Confirmation.** A run is confirmed when both strands give it one
    clear, equal length that is DNA's call.
  - **`length_unconfirmed`** (check, rank 2): every other run that a decode
    reaches. It lists the runs with each read's length, source, share,
    margin, and the agreement.
  - **`length_not_recorded`** (context, rank 4): in HVS-I, whose run
    lengths the profile's convention leaves unrecorded.
  - **Scope.** HVS-III's 568-573 run has no profile window and is not
    flagged.
  - **Measured.** On development panel v5 against the same build without
    it (`dna-eval diff`, GATE PASS):
    - passed window@HV2 zones with a non-het unreviewed difference fall from
      50 to 8 (19 → 1 on current data);
    - passed samples with a non-het difference fall from 69 to 30;
    - window@HV2 checks rise from 1,664 to 3,073;
    - clear samples on current data fall from 76.0 % to 71.1 %;
    - 1,363 of the 1,409 newly checked window@HV2 units hold no truth
      difference.
  - Review rules 21.

- The review's `possible_heteroplasmy` is `explained` context, not a check,
  when every qualifying read of one strand shows a second peak the Sanger
  evidence explains by its shape (a neighbour's tail or a dye blob, not
  `second_unexplained`, SRS-CALL-022); a position outside every site that
  holds only such a label is a `note` item (`review` plugin 19, rules
  version 19). On the development panel none of the 23 such flags was a
  heteroplasmy or hid a difference: 14 checks go and 4 samples clear, with
  no unreviewed difference added.
- A Sanger call before a run that the read loses phase after
  (`post_homopolymer`) is no longer at a signal edge when its own signal is
  unambiguous (SRS-CALL-017): before the run every molecule shows the same
  sequence. The reverse read now reads 16194-16205 before the HVS-I C-tract;
  on panel v5, 33 more samples are clear, with unreviewed and non-het
  unchanged. Plugin `sanger` 14.

- More Sanger calls may confirm the reference (ADR-0080, SRS-CALL-022;
  plugin `sanger` 15, after the previous entry's 14). A single-channel call withheld only for its quality,
  or at a benign signal edge, may show a second base standing out when its
  cause is a neighbour's tail or a dye blob. An isolated mixed tail or blob
  call keeps its confirmation at a benign signal edge (the same test), may
  also be withheld for the peak floor or its relative quality, and may have a
  neighbour masked as weak before the trim interval. An off-centre,
  second-template, or phase-shift second never qualifies, a read's physical
  end is never benign, and every core guard is unchanged. On panel v5
  against `main` 20fc7ca, samples with an unreviewed or non-het unreviewed
  difference are unchanged in every period and no status changes; clear
  samples rise 71.00 → 72.05 % on current data (+113), 48.59 → 49.29 % in
  2025-08..11 (+19), and 69.41 → 70.84 % older (+43); the current
  `start@HV1` zone falls from 7.09 to 4.76 % of samples.
- An unexplained second peak on one strand alone no longer withdraws a
  reference confirmation made by a read of the other strand (SRS-CONS-010,
  ADR-0080): it still does on the confirming strand, where the read's second
  base does not stand out, or where a read of the confirming strand shows
  that second base standing out too. On panel v5, unreviewed and non-het are
  unchanged in every period, no status changes, and 156 samples become clear
  (current 71.00 -> 71.39 %, 2025-08..11 48.59 -> 52.15 %, older 69.41 ->
  70.01 %). Plugin `consensus` 14.
- A call within `variant_calling.read_end_margin` of an unvouched read end
  confirms the reference when its modality vouches for its base (scored,
  off its edges, withheld for nothing) and it shows no second base standing
  out and no unexplained second peak (SRS-CONS-010, ADR-0080). Every guard of
  a confirmation applies; such a call still bounds no run and supports no
  variant. On panel v5 against the stack of #169-#173, unreviewed and
  non-het are unchanged in every period, no status changes, and 48 samples
  held only by HVS-I read-start checks become clear (current 72.80 ->
  73.16 %, older 71.63 -> 71.96 %). Plugin `consensus` 15.
- An opposing read whose dissent from a called variant is only a run's
  length is set aside under the explained-second rule (SRS-CONS-011,
  ADR-0083): every call it shows there is resolved, it shows the decided
  runs in the same order, it shows a second base standing out there that its
  modality explains, and a clean read of the other strand without a second
  peak supports the decision. The indel-at-a-read-end rule stays
  reference-only, and a substitution against a called variant is never set
  aside. This lifts the HV1F overcall of the G run 66-71 beside `73G`. On
  panel v5 against `main` 69e1471, unreviewed and non-het are unchanged in
  every period, no status changes, and 19 samples become clear (current
  73.16 -> 73.31 %, 2025-08..11 53.38 -> 53.41 %, older 71.96 -> 72.00 %).
  Plugin `consensus` 16.
- Phase-mixture decomposition moves from `dna-sanger` into the new library
  crate `dna-trace` (PROP-0012 increment 1, ADR-0084), which reads plain
  per-call data and depends on `dna-kernel` only; `dna-sanger` calls it. No
  output changes, so no plugin version changes.
- `dna.read_evidence/v1` carries everything the analysis document reports
  beside the core's products, in place (pre-release): `read.trim`,
  `signal_quality`, and `warnings` as in `dna.basecalls/v3`, and each call's
  `signal.peaks` (A/C/G/T heights at its primary event, read strand). Profile
  weights are written at full precision: rounded to three decimals they
  changed the quantized alignment weight of 27 % of calls and the analysis
  of 2 reads in 1,000. No plugin version changes: no plugin's output changes.
- The review shows each admitted read's phase mixtures (ADR-0082): for
  every stretch after a run that the Sanger plugin decomposes into the read's
  own sequence at whole-call offsets, `reads[].mixtures` lists the dominant
  and the other sequences with their fractions, each compared with the
  reference on its own and written in the profile notation, and the calls
  they leave unexplained (`review` plugin 17, rules version 17). Display
  only: no item, label, or `clear` changes.
- The review judges a site the reportable regions clip over its reportable
  columns only (`review` plugin 16, rules version 16): the item's `sequence`
  is the decision there, and an opposing read that shows it there, its
  difference outside the regions, no longer makes `strands_disagree`. On the
  current development data 32 such checks go (66 at position 73, where a
  wrapped HV1F read differs at 66-72), clearing 30 samples, and 1,789 clipped
  items now show their clipped sequence.
- Read-start calls may confirm the reference too (ADR-0080, SRS-CALL-022,
  SRS-CONS-010). The Sanger plugin marks as reference-only a single-channel
  call without a second base standing out or unexplained, withheld only for
  its peak floor or relative quality, or lying at a signal edge whose every
  cause is the untrimmed read start or a weak mask clear of repeat runs. The
  core lets such a call confirm beyond the read's unvouched ends (new
  `Eligibility::clear_of_ends`).
  - Every confirmation is now also withdrawn:
    - at or beside any read's observed difference, clean or not;
    - in a stretch;
    - in a reference run of at least `consensus.seeded_run` across a
      junction no read observes cleanly.
  - Plugins `sanger` 12 and `consensus` 12.
  - Development panel v5 against `main` b49ff87:
    - samples with an unreviewed difference are unchanged in every period;
    - clear samples go from 68.7 to 70.2 % on current data;
    - no strict or `haplotype_lod` changes.
- The review lists where DNA decided by reference confirmations alone:
  `reference_only` holds the reportable positions the consensus published as
  `reference_only` (ADR-0080), and `summary.positions.reference_only` counts
  them (`review` plugin 15, rules version 15). Such a position has no item,
  so this is the only place a reviewer can find it.
- A mixed Sanger call whose second peak is a neighbour's tail or a dye blob,
  isolated between in-phase single-channel calls and withheld for nothing
  else, may now confirm the reference at its position (ADR-0080,
  SRS-CALL-022, SRS-CONS-010). It keeps its `mixed_supporting_dna` reason
  and its run role. The kernel contract gains three modality assertions:
  `CallEvidence.reference_only`, `CallEvidence.second_unexplained`, and
  `SecondBase.stands_out`.
  - The consensus withdraws a confirmation where any read shows an
    unexplained second peak, where the other strand shows the same second
    base standing out, and within or beside a site that differs from the
    reference.
  - It publishes the positions only confirmations decided as `reference_only`
    in `dna.consensus/v1` (revised in place). The call step lets them bound
    no neighbouring run.
  - Plugins `sanger` 11, `consensus` 11, and `decided` 3.
  - Development panel v5 against `main` 65c221e:
    - samples with an unreviewed difference are unchanged on current data
      and in 2025-08..11; older data gains one, a truth `N` at 538;
    - clear samples go from 66.1 to 68.7 % on current data, 42.3 to 47.4 %
      in 2025-08..11, and 60.0 to 67.2 % older;
    - no strict status changes.
- The review lists the sample's zones: `zones` tiles the reportable regions
  into maximal stretches of one kind, from the profile's repeat windows and
  the sample's own reads (`repeat`, `not_callable`, `dephased`, `masked`,
  `clean`), each with the `check` items touching it, and `summary.zones`
  counts them (`review` plugin 14, rules version 14). No width or threshold
  places a zone. On the current development data a called sample has 0.93
  zones to check on average (p90 3), and the `dephased` zones hold almost no
  checks.
- The review marks a long deletion one read decides inside a run: a new
  check label `run_deletion` for a site resting on a single read whose call
  is at least `review.run_growth` bases shorter than its reference, where the
  reference holds a run of at least `review.phase_run` bases (`review` plugin
  13, rules version 13). On the current development data it turns the 4
  auto-passed sites whose deletion was wrong (2550268's HVS-I run, three
  compressed HVS-II runs) into checks, with 2 more checks on true deletions;
  older slices add one check.
- The review shows why a call has a second peak: an evidence position's
  `second` carries `cause` (SRS-CALL-021) when its share passes
  `signal_processing.second_floor` or the call is masked as mixed, and with
  `--signal` an unexplained second the run's reads share there reads
  `context` (`review` plugin 12, rules version 12). Items, calls, and checks
  are unchanged.
- A read's `masked` stretches no longer list a leading call the modality
  vouches for at a read start (SRS-CALL-014): such a call is scored and the
  read shows its base (`review` plugin 11, rules version 11). On current
  development data this removes 5,167 of 486,695 masked positions, all at
  the HV1F and HV2F read starts; nothing else changes.
- A leading call before the callable span is no longer vouched for by its base
  (SRS-CALL-014) when its read is masked after a repeat run
  (`post_homopolymer`): that stretch is the read's dephased tail, not the
  start noise the rule is for. With it, no base, insertion, or deletion comes
  from a masked call; a masked call shows at most a run's end (`anchored_end`,
  ADR-0071), which the review labels `dephased_end` (#145). Plugin `sanger`
  10. Development panel v5: no output changes.
- The run signal's context width and share resolution are configuration:
  `signal_processing.run_signal_flank` (2, a 5-mer key) and
  `signal_processing.run_signal_resolution` (1000, thousandths), configuration
  schema 15. `dna.signal/v1` records both in `provenance` (`flank`,
  `resolution`), revised in place before its first release, and `dna review
  --signal` refuses a run signal whose values differ from the
  configuration's. Plugin `sanger` 9. At the default values run signals,
  reviews, and calls are unchanged.
- The review records where each read does not show its bases, and which
  run lengths lean on a dephased tail (`review` plugin 10, rules version 10;
  `dna.review/v1` revised in place):
  - every admitted read lists `masked`, its stretches of reference positions
    inside the reportable regions covered while its signal is masked, with the
    callability phase state and whether the stretch follows a repeat run,
    taken from the read's own signal, never from fixed positions;
  - a context label `dephased_end` marks a site whose run length rests on an
    anchored end, the first masked call after the run. Reading the dephased
    tail itself is later work; the 2,119 such auto-passed lengths on current
    development data were 98.4 % right, so they stay decided.
- A start call may read the base of the call before it when that call is
  vouched for its base (SRS-CALL-018), judged in read order, so the runs of a
  read's start (`TT` at 16024-16025, `TTT` at 16027-16029, `GGG` at
  16034-16036 in HVS-I) chain from their first call; the start's own
  second-channel share is `sanger_evidence.start_second_share` (0.15), set
  apart from the read end's `edge_second_share` (0.10); configuration schema
  14. Plugin `sanger` 8. Development panel v5 after #141, current data
  (batch date 2026-01 on): samples with a check in the HVS-I amplicon start
  16.9 to 13.2 %, clear samples 63.88 to 65.69 %, zones to check per sample
  0.62 to 0.58, samples with an unreviewed difference unchanged in every
  period (4.14 % current); two more correct `73G` calls. The share was chosen
  on the development half (0.20 adds a false `77C`). The one new false call,
  a single-read `16038DEL` in a 2025-10 batch, is an `edge_deletion` check.

- The Sanger peak floor no longer withholds a call whose relative quality
  reaches `sanger_evidence.peak_rescue_quality` (50) and whose second base
  holds less than `sanger_evidence.peak_rescue_second_share` (0.20) of the two
  strongest peaks (`sanger` plugin 7, configuration schema 13); a low peak
  beside a strong second base is the split signal of two templates. A low absolute peak at a weak read end or in a
  weak dye context (the reverse read's C at 317 and 320 after the HVS-II
  run) is not a doubtful call when the call stands clear of its
  neighbourhood. On the current development data the floor-withheld calls at
  that quality matched the truth at every one of 11,756 positions outside
  heteroplasmies and indels. It decides about 8,300 undecided checks: clear
  samples rise from about 56 % to about 64 %, review zones per sample fall
  from 0.84 to 0.62, 133 calls are gained that the truth has; one two-base
  heteroplasmy (195Y) is no longer inside a check, against eight without the
  second-share bound.
- `estimated_length` no longer lists a run whose length a profile report
  convention leaves unrecorded (`review` plugin 9, rules version 9): on the
  human profile, the uninterrupted HVS-I C run of ten or more bases that
  16189C creates, whose length the lab's reviewers do not record. The run
  must reach the bound without the bases the site's call adds, so an
  estimate that lengthens a shorter run past it stays listed. Such a site's
  other labels still apply. On the current development data it turns
  1,810 checks into variants and clears about 750 more samples (7 points);
  the checks it removes held no non-heteroplasmy difference of their own,
  but about 40 two-base heteroplasmies next to the run (16192Y, 16183M) are
  no longer inside a check.
- `read_disagrees` also fires on one end of an item that a read does not
  show whole (`review` plugin 8, rules version 8; configuration schema 12
  with `review.phase_run`). At a decided item without a length change, the
  end a read's strand reaches first is compared with the decision, base by
  base, up to where the read has read `review.phase_run` (4) equal bases in a
  row; past that the run's far end is a length question dephasing leaves
  open. A clean, full-trust base there that differs makes the item a check.
  It catches a run exit that one strand's dephased tail decided while the
  other strand reads it cleanly in phase (`16194C` beside a clean reverse
  `A`). On the current development data it adds 17 checks over 10,697
  samples, 7 of them on an unreviewed non-heteroplasmy difference.
- The review list marks two kinds of decided indels that one strand cannot
  assure (`review` plugin 7, rules version 7; configuration schema 11 with a
  new `[review]` section). Where both strands cover a site, both must agree
  before DNA decides it alone:
  - `long_run`: a run at least `review.run_growth` (3) bases longer than its
    reference run whose length rests on an anchored end or on a single read;
  - `edge_deletion`: a deletion resting on a single read within
    `review.edge_zone` (20) positions of a reportable region's edge.
  On the current development data such decisions were 89 % and 88 % wrong.
  As checks they take the samples with an unreviewed non-heteroplasmy
  difference from 1.49 % to about 0.77 %, adding about one check per 100
  samples.
- The core and post-calling method thresholds are configuration too
  (SRS-CFG-009), with defaults equal to the former constants, so outputs
  change only in the configuration checksum and the profile identity:
  - new `[consensus]` section (configuration schema 10): `seeded_run` (4,
    was a method constant of SRS-CONS-004), `maximum_undecided_gap` and
    `minimum_resolved_bases` (20 each, set apart from
    `alignment.minimum_callable_bases`), and the decided site alignment's
    `site_mismatch_cost`, `site_gap_open_cost`, and
    `site_gap_extension_cost` (2, 3, 2, with their ordering checked at load);
  - `alignment.max_stretch` (64, SRS-ALN-016);
  - target profile schema 2: `normalization.max_trade` (20, SRS-VN-008),
    carried by `NormalizationPolicy::RightAligned { max_trade }`, and a
    tandem-repeat window's `read_copies` (2), the motif copies `motif_shift`
    asks a read to show.
  `dna-eval run` records its measurement constants in `run.json`
  (`constants`).

- `possible_heteroplasmy` fires only when both strands read the position
  cleanly and confirm the same second base (`review` plugin 6, rules
  version 6). A site that only one strand reads gives no label. On current
  development data (floor 0.20), two-strand flags were 90 % precise (128 of
  142) and one-strand flags 3 % (76 of 2,272), so one-strand flags cost
  reviewers far more than they found. This follows the user's rule that a
  heteroplasmy flag must not add false positives, and the forensic practice
  of confirming point heteroplasmy on both strands.
- Every Sanger method threshold is configuration (SRS-CFG-009; schema 9),
  with defaults equal to the former constants, so outputs change only in
  the configuration checksum:
  - `[sanger_evidence]`: the read-edge rules' `edge_second_share` (0.10),
    `edge_width_factor` (1.5), `edge_spacing_minimum`/`maximum` (0.5, 1.5),
    `edge_height_fraction` and `leading_height_fraction` (1.0),
    `start_norm_calls` (30), `start_height_fraction` (0.5), and
    `post_repeat_minimum_calls` (2);
  - `[callability]`: `spacing_defect` (0.5) and `spacing_half_window` (8);
  - `[signal_processing]`: `minimum_noise_sigma` (1.0).
  Keys that drove several rules are split, so one tuning moves one rule:
  - `variant_calling.read_end_margin` is now the core's alone; the Sanger
    plugin has `signal_edge_margin`, `post_repeat_margin`,
    `edge_identity_margin`, `start_identity_margin`, and
    `quality_control.trim_context_margin` (all 12);
  - `callability.double_peak_ratio` (0.33) is set apart from
    `basecalling.secondary_peak_ratio`;
  - `callability.repeat_window_calls` (16) is set apart from `window_calls`.
  The shadow model's offset range, minimum fit size, and solver tolerances
  stay method constants of its definition.

- `signal_processing.second_floor` defaults to 0.20 (was 0.16), the floor
  that meets the heteroplasmy policy (2026-10-18, user): a
  `possible_heteroplasmy` flag is worth raising only without more false
  flags for reviewers. Chosen on the development half as the lowest floor
  whose false flags (no two-base IUPAC code in the truth at the item) are
  at most today's and whose precision is at least today's, with and
  without a run signal. Development panel v5, current data, test half,
  without a run signal: 1,172 flags, 99 true, 1,073 false (precision
  8.4 %), against 1,610, 123, and 1,487 (7.6 %) before PROP-0009 (ii); with
  batch run signals 1,022 flags, 99 true, 923 false (9.7 %). Truth
  heteroplasmy flagged falls from 25.5 % to 20.5 %; samples with an
  unreviewed difference rise from 4.03 % to 4.48 % (4.50 % with batch run
  signals), and clear samples from 47.8 % to 48.8 %.

- One second-channel scale and threshold (PROP-0009, increment ii;
  SRS-CALL-020). `dna_sanger::run_signal::threshold` is the configured
  `signal_processing.second_floor`, raised to the call's key run level when a
  run signal holds the key with enough reads and calls. The review's
  `possible_heteroplasmy` now judges a clean call's second peak share
  (SRS-CALL-015) against it in place of a fixed `0.15` profile weight, with
  the neighbour-shadow and strand rules unchanged, and the review evidence
  `second` moves to the peak scale. Configuration schema 8 adds
  `second_floor` (0.16) and the `run_signal_*` keys; plugin `review` 5.
  Development panel v5, current data (batch date 2026-01 on), test half, at
  that first floor of 0.16: without a run signal 3.96 % of samples with an
  unreviewed difference (4.03 % before), 18.9 positions to review per
  sample (18.3), 47.7 % clear (47.8 %); with batch run signals 4.01 %,
  18.4, and 48.1 %. The default floor is now 0.20 (above).

- Start-of-read calls vouched for by their own signal (SRS-CALL-018). At a
  read's start, up to the edge margin past the callable span's start, a
  scored edge call whose signal is unambiguous against the start's own norms
  (the first 30 in-phase calls of the callable span; a primary peak at least
  half their median) has its verdict on its base taken off the edge, even
  behind masked start noise. Plugin `sanger` 6. Development panel v5: clear
  samples 40.56 to 43.79 % (45.11 to 48.22 % on batches from 2026),
  positions to review 3.20 to 3.07 %, samples with an unreviewed difference
  unchanged, `haplotype_lod` unchanged.

- Edge calls vouched for by their own signal (SRS-CALL-017). At an end of a
  read's informative interval, a scored call at a signal edge whose signal is
  unambiguous (one qualifying channel, a primary peak at least the read's
  median, a second share below 0.10, a peak no broader than 1.5 times the
  read's median width, regular spacing, no phase loss beyond the end, and not
  the base of its neighbour towards the end) has its verdict on its base taken
  off the edge; its indel verdict keeps the edge. Signal processing now
  measures each channel's peak width. Plugin `sanger` 5. Development panel
  v5: positions to review 3.28 to 3.20 %, clear samples 39.36 to 40.56 %,
  samples with an unreviewed difference unchanged; `315.1C` called in 126
  more samples, `73G` in 56; `strict` 33 better and 1 worse.

- `Alignment` is read-only outside `dna-align` (ADR-0079 follow-up): its
  orientation, segments, origin flag, and metrics are getters, and the cuts
  of its callable segments (an excluded variant's positions, SRS-VAR-019; a
  placement conflict, SRS-VAR-017) are its own methods, so `dna-align` is the
  only writer of an alignment. Outputs are unchanged.

- The core is one crate per step (PROP-0008 step 5, ADR-0079): `dna-core`
  becomes `dna-align`, `dna-call`, `dna-sample`, and `dna-consensus`, a chain
  on the kernel that `validate_module_layers.py` enforces. Only what the next
  step reads crosses a boundary: a read-only view of the alignment columns,
  a validating `Alignment::placed` for the consensus's decided reads,
  `Eligibility` as the one trust rule, `CallConfig` (alignment and variant
  calling) instead of `CoreConfig`, and a kernel contract `CalledReads` that
  the consensus requires. Plugin identities and every output are unchanged.
  Operation-log records of the core now carry the target
  `dna_call::read_call` instead of `dna_core::read_call`; a log filter that
  matches on the old target must change.

- Sanger failures belong to the Sanger plugin (PROP-0008 step 4, part 2c).
  `dna-sanger` owns `SangerError` with the ABIF, re-calling, signal,
  callability, and quality-control failures, and `ProductsError` for one
  read's Sanger products disagreeing while its evidence is assembled; the
  kernel keeps only modality-neutral failures. Outputs, command-line
  messages, and log text are unchanged. **Breaking (Rust API):**
  `dna::error::Error` is now `Kernel(KernelError)`, `Sanger(SangerError)`, or
  `OperationAndLog`; the kernel's former `Error` is re-exported as
  `dna::error::KernelError`, and `Profile::load`,
  `variant_normalization::normalize`, and `variant_nomenclature::apply`
  return it; `VariantError` loses the unused `MissingQuality`.

- Verdict scope (PROP-0008 step 4, part 2b). The consensus reads the
  modality's reasons by scope, as variant calling does (SRS-CONS-003); an
  insertion or deletion on a call vouched for by its base only is excluded
  with the new core reason `base_only` instead of `read_end` (SRS-VAR-013);
  overlap agreement compares only bases the core trusts (SRS-SAMPLE-008); a
  verdict refuses a reason it cannot hold instead of dropping it. A
  candidate excluded before it is observed (an indel no flank reads, an
  unresolved allele) leaves its positions undecided (SRS-VAR-019), and a
  consensus site that changes the length stays undecided while a position
  beside it is undecided, unless reads show its runs' lengths
  (SRS-CONS-009). Plugins `core` 11 and `consensus` 7.
  `dna.sample_evidence/v10` lists `base_only` among the exclusion reasons
  (revised in place, unreleased). Development panel v5: `haplotype_lod` 3
  better and none worse.

- `ReadEvidence` v2 (PROP-0008 step 4, part 2a). Each call carries its
  placement (scored, or masked as anchoring or unresolved with its reason)
  and the modality's verdicts on its base, on an indel, and on a run's end,
  each with its reasons and signal edge; the read lists the modality's
  reasons in reporting order. The Sanger plugin marks the bounded runs of
  ADR-0076 itself (SRS-CALL-016). Outputs are unchanged. **Breaking (Rust
  API):** the kernel drops `CallMask`, `MaskOnset`, `MaskedAlignment`,
  `VetoSet`, `SupportVeto`, and `VetoScope`; `CallEvidence::new` becomes
  `CallEvidence::scored` and `masked`; `ReadEvidence::new` takes the reason
  order; `CallTrust.mask` becomes `placement`, and `Withholder::Masked`
  carries the reason.

- Signal-edge verdicts and the read start (PROP-0008 step 4, part 1). The
  core's `read_end` margin around trim ends and masked calls is now the
  modality's verdict per call (`CallEvidence::edge`); the core keeps only the
  margin at a read end the modality does not vouch for (SRS-VAR-013,
  SRS-CALL-013). Before the callable span, a masked call with one qualifying
  channel and a primary peak at least the read's median is vouched for by its
  base only (SRS-CALL-014): it supports a substitution, never an insertion,
  deletion, or run end, and not where it may sit one position off
  (`unresolved_flank`, SRS-VAR-018). Every call carries its second base and
  share (SRS-CALL-015). Plugins `core` 9, `consensus` 6, `sanger` 4.
  **Breaking (Rust API):** `CallEvidence` gains `edge`, `identity_only`, and
  `second`, with `CallEvidence::new`; `CallTrust` gains `scope`
  (`TrustScope`) and `second` and is no longer `Eq`. Development panel v5:
  truth `73G` uncallable in 192 samples instead of 745; `haplotype_lod` 2
  better and 3 worse.

- Consensus lengths (PROP-0008 step 3). A reference homopolymer of at least
  four bases that no read shows cleanly with both neighbours, while some read
  differs on it, is decided by its run structure, so its length is labelled
  (SRS-CONS-004). A run of length evidence `phase_loss` or `reference_frame`
  keeps its site even at the reference length (SRS-CONS-007). `call` of the
  consensus document leaves such a run's 3' end, where a change of its length
  is written, not callable and excludes a change there as `run_boundary`; the
  run's composition stays callable (SRS-VAR-021). Plugins `consensus` 5 and
  `decided` 2. Development panel v5: `haplotype_lod` 677 samples better and
  5 worse; strict 549 worse, all truth insertions at `309` after a run whose
  length is an estimate.

- `call` given a sample's consensus document (`dna call <id>
  <id>.consensus.json`), and `SampleAnalysis::call`, report the edits the
  consensus decided without placing its segments on the reference again
  (PROP-0006, ADR-0078, SRS-VAR-021). Each site's call is aligned against its
  own interval only, gaps take their 3'-most placement, the record rules
  apply, and a length change of the run at a segment's end is excluded as
  `run_boundary`. New plugin `decided` 1; a FASTA sequence keeps the placement
  path. **Breaking:** SRS-API-012 now means the consensus's decided edits, and
  the `dna.variants/v1` reads of a consensus input are its segments (revised
  in place, unreleased). Development panel v5: `haplotype_lod` 14 samples
  better and none worse, 8 samples that failed now call, and no EMPOP range
  covers a position the consensus left undecided.

- The evaluation harness's own measures (strict, IUPAC-compatible, and the
  haplotype measures, the approval endpoint `haplotype_lod` included) judge
  insertions and deletions at 309, 455, 463, and 573 like any other: DNA is
  target-agnostic, and only the primary measure, the lab's rules, keeps the
  lab's exception for these positions. The haplotype measures no longer
  build each interval with the other side's special-position indels.
  Development panel v5: the primary measure is unchanged; 738 samples become
  discordant on `haplotype_lod` (695 at 309, 44 at 573, 1 at 455/463), its
  corpus not-concordant rate goes from 0.0650 to 0.1112, and the `sequence`
  class from 346 to 1,125 samples, mostly the HVS-II 303-309 run DNA writes
  at the reference length without a label.

- The evaluation harness's discrepancy causes read the core's own trust
  verdict on each call (`SampleAnalysis::trust`, `CallTrust`: mask, vetoes,
  `read_end` trust, ADR-0076 exemption, and what withholds trust) instead of
  re-deriving it from four configuration keys. A read's end counts as phase
  loss only where DNA attributes it so (a phase-loss mask label, or a segment
  after a repeat run); development panel v5 moves `beside_phase_loss` from
  141 to 43 findings, mostly to `read_edge` (45 to 112), and statuses and
  classes are unchanged.
- Right-aligned normalization takes an insertion across a base that a
  substitution of the same haplotype makes equal to its first inserted base,
  until the insertion follows a reference run of that base: with `16172C`, an
  inserted C is written `16174.1C` rather than `16171.1C`, as the lab's
  reviewers write it, while `16262.1C 16263C` stays. A deletion that cannot
  shift further moves up to 20 bases 3' across substitutions when the same
  haplotype needs no more edits, every base it crosses is read, and it
  crosses no nomenclature-window boundary (`16293DEL 16298C` becomes
  `16293C 16298DEL`). A candidate shift onto a substituted base is rejected
  rather than failing the read. The `normalization` plugin version rises to 2.
- Window representation now applies per window rather than per read
  (SRS-NOM-007, SRS-NOM-010). An edit that crosses one nomenclature window
  leaves that window in its normalized form, but no longer prevents the
  read's other windows from taking their forensic form: for example a
  deletion of 512-513 beside HVS-III no longer leaves HVS-II's
  `309.1C 315.1C` as `309.1C 309.2T 310C`. **Breaking for library
  consumers:** `dna::variant_nomenclature::apply` no longer fails with
  `NomenclatureError::WindowCrossing`, which is removed. The `nomenclature`
  plugin version rises to 2.
- The nomenclature windows follow the lab's reviewer convention (ADR-0077,
  SRS-NOM-005, SRS-NOM-014, SRS-NOM-015). HVS-II writes any C/T haplotype as
  the left run's length at 309 and the rest left-anchored with C insertions
  after 315 (`312T 315.1C 315.2C 315.3C`); HVS-III now starts at 512 and
  writes whole AC-copy losses and gains at its 3' end (`513A 523DEL 524DEL`,
  `524.1A 524.2C`); HVS-I now starts at 16180 and gains the `anchored_tract`
  rule for the A-tract and C-tract (`16183DEL`, `16188.1C`,
  `16182C 16183C 16189C`, `16189C 16190T 16193.1C 16193.2C`). In `sample`,
  `call`, and `notation`, a window keeps its normalized form when its
  representation would place a call where the read reads nothing, or rests
  on a base the read does not read (for the phylogenetic HVS-I form, any base
  of the window). Every window sequence of at least 5
  reviewed corpus samples now takes the form most reviewers wrote; on
  development panel v5 the primary measure gains 60 samples, strict 56, and
  IUPAC-compatible 80, with no regression on any measure.
- **Breaking for library and log consumers:** the post-calling crate
  `dna-post` is split into `dna-normalize` (representation, normalization,
  nomenclature) and `dna-empop` (conformance) (ADR-0077). Rust paths under
  `dna_post::` move to `dna_normalize::` or `dna_empop::`, and the operation
  log target `dna_post` becomes `dna_normalize` and `dna_empop`. CLI outputs,
  document contracts, the `dna::variant_normalization` and
  `dna::variant_nomenclature` re-exports, and plugin identifiers and versions
  are unchanged.
- Bounded runs before a post-repeat phase loss are exempt from `read_end`
  (ADR-0076, SRS-VAR-013, SRS-CONS-003). An insertion or deletion of a run the
  read bounds with unmasked calls, and the consensus's observations of that
  run's matching and inserted calls and of its reference-matching bounds, are
  no longer distrusted when the margin's edge is the read losing phase right
  after a repeat run. Substitutions, and bounds that read another base, keep
  the margin. `ReadEvidence` masks carry a `MaskOnset`; the `core` plugin
  version rises to 4, the `sanger` plugin version to 3, and the `consensus`
  plugin version to 3.
- Faster per-sample analysis with byte-identical results (SRS-ALN-014).
  - **Seeded pruning** now covers reads with masked or mixed calls (wildcard
    loci) and windows across the circular origin, and an orientation
    evaluated exactly sets the threshold that rejects the other. On 1,197
    real reads, 96 % avoid full-reference Gotoh, up from none. A read whose
    seed probes cover most of the reference goes straight to Gotoh.
  - **Gotoh** fills its matrix in one pass per row with one traceback write
    per cell; **signal statistics** take medians by selection.
  - `dna consensus` drops from about 4 s to 0.13 s per four-trace sample.
- Faster sample commands with byte-identical results (SRS-SAMPLE-001,
  SRS-OUT-008).
  - **Signal statistics** of each sample span are estimated once per read and
    shared by the windowed SNR features and the locus contexts: signal
    processing takes about half the time.
  - **Reads run concurrently** in `dna sample` and `dna consensus`, up to the
    available parallelism. The first failing read in trace order still decides
    the failure, and the operation log still lists each read's records in
    trace order, through that read.
  - On 2,618 panel samples every `consensus` JSON and FASTA and every `sample`
    JSON is byte-identical to the previous release. On a shared server, mean
    `dna consensus` wall time per sample drops from about 134 ms to 72 ms, and
    from 123 ms to 102 ms when pinned to one CPU.
- `consensus` method version 2 decides a stretch by its run structure and
  labels how each run length is known (ADR-0074, SRS-CONS-003/004/008).
  - **Stretches.** A stretch also takes in a neighbouring run that a
    substituted base or a deletion joins, as `16189C` does with the HVS-I
    C-tract.
  - **Composition, then lengths.** The ordered run bases are decided first.
    Run lengths follow, voted as whole sequences when reads show every run's
    length. A length comes only from reads that show both ends of the run
    (ADR-0071); an insertion or deletion of a run the read does not bound is
    never clean.
  - **A run no read shows both ends of** takes one of two lengths:
    - `reference_frame`: the reference positions between the ends that reads
      locate, when every base between is read cleanly;
    - otherwise `phase_loss`: the in-phase length before a read loses phase.
  - **Measured** on the fair scale against method version 1:
    - exact samples rise from 290 to 297 and from 292 to 299 on the two
      320-sample subsets;
    - recall rises from 0.9866 to 0.9899 and from 0.9901 to 0.9928;
    - precision moves from 0.9967 to 0.9970 and from 0.9964 to 0.9967;
    - reviewer-consensus calls identical on 145 of 160 samples (was 142);
    - the local Sequencher set has 28 exact samples (was 26).
    - The false `16182.1C` that a phase-lost read added in HVS-I is gone, and
      HVS-II dominant lengths (`309.1C`, `309.2C`) are kept.
  - **Breaking:** `dna.consensus/v1` is revised in place. Each site gains
    `runs` (`base`, `length`, `length_evidence`), and the summary gains
    `phase_loss_runs` and `reference_frame_runs`.
  - `dna sample` and `dna call` are unchanged apart from the `core` version.

- The `core` plugin method version rises to 2. It records the output changes
  of `run_boundary` (ADR-0071) and run-structure re-expression (ADR-0072),
  which earlier entries in this cycle made without raising it
  ([versioning](docs/governance/versioning.md#plugins)).

- A read the core cannot place no longer fails `dna sample` (SRS-SAMPLE-029,
  PROP-0003 increment 0).
  - The covered cases are too few callable columns, identity below
    `alignment.minimum_identity`, more than one best placement, and tied
    orientations.
  - Such a read is recorded in `rejected_reads[]` with its own reason and
    measured values, and the other reads are still aggregated.
  - On a fresh 320-sample held-out subset, this lets 5 samples that failed
    before produce evidence.
  - **Breaking:** the `rejection` object of `dna.sample_evidence/v10` gains
    four reasons (revised in place), and the error for a sample whose every
    read is rejected now reads "all N reads were rejected and none could be
    analyzed".

- The core re-expresses masked repeat stretches of the selected alignment
  (ADR-0072, SRS-ALN-016).
  - A stretch qualifies when its read calls spell the reference's runs in
    order, with every changed run longer or every changed run shorter, and
    it holds a masked call.
  - It is rewritten as one length edit at the 3' end of each changed run.
  - This recovers HVS-II length insertions that the optimal path split into
    a substitution next to a masked insertion.
  - On the 320-sample held-out set, 6 true `309.1C`/`309.2C` calls are
    gained and one false `310C` is removed; precision rises from 0.9894 to
    0.9897 and recall from 0.9800 to 0.9816.
  - Consensus calls from `dna call` are unchanged.
  - Alignment summaries (`identity`, `gap_opens`) and observed evidence
    change for reads with a rewritten stretch.

- A variant that changes the length of a run of one base is ineligible with
  the new core reason `run_boundary` when the read does not show where that
  run ends on the side the edit touches (ADR-0071, SRS-VAR-016).
  - This targets Sanger reads that dephase right after a poly-C run and
    miscount it by one base. A dephased call that still reads the next base
    resolves the run end.
  - On the 320-sample held-out set, 9 false length calls are removed and
    none is gained; precision rises from 0.9870 to 0.9894, and recall stays
    at 0.9800. Consensus calls from `dna call` are unchanged.
  - **Breaking:** the exclusion-reason vocabulary of
    `dna.sample_evidence/v10` and `dna.variants/v1` gains `run_boundary`
    (revised in place).

- Plugin-first architecture, phase 5b (ADR-0070, PROP-0002): DNA is a Cargo
  workspace.
  - `dna-kernel` (shared contracts), `dna-core` (the core caller),
    `dna-sanger` (the Sanger modality), and `dna-post` (post-calling plugins)
    live under `crates/`; the `dna` facade at the root composes them.
  - Each plugin crate builds and tests alone.
  - The public Rust API (`dna::error`, `dna::variant`, `dna::profile`,
    `dna::variant_normalization`, `dna::variant_nomenclature`,
    `dna::variant_analysis`, `dna::cli`, `dna::run`) and every result
    document are unchanged.
  - The operation log's module path for events of moved code now names its
    crate.
  - Cargo commands take `--workspace`, and builds take `-p dna`.

- Plugin-first architecture, phase 5a (PROP-0002): the sources are untangled
  along the future crate boundaries. Each plugin owns its configuration
  sections and its plugin descriptor; the facade composes the configuration
  envelope, the plugin registry, and the workflow compositions. The module
  validator enforces the crate map (kernel, core, sanger, post, dna) instead of
  layers. Configuration files, result documents, and the public Rust API are
  unchanged.

- Plugin-first architecture, phase 2b (ADR-0069, PROP-0002): sample
  aggregation is modality-neutral. It consumes each read's `CalledRead`
  (identity, evidence, alignment, variants), and the sample report joins the
  read's Sanger attachment by read identity and call index. Every result
  document is unchanged.
- Plugin-first architecture, phase 2 (ADR-0069, PROP-0002): a static plugin
  registry (`sanger` modality, `core`, and the `normalization` and
  `nomenclature` post-calling plugins) and the plugin composition of each
  workflow are validated at compile time. Each configuration section has one
  owning plugin. Scientific results are unchanged.
  **Breaking (unreleased, revised in place):**
  - configuration schema 7 moves `minimum_peak_height` and
    `relative_quality_threshold` from `[variant_calling]` to a new
    Sanger-owned `[sanger_evidence]` section;
  - `dna.analysis/v9`, `dna.basecalls/v3`, and `dna.sample_evidence/v10`
    require `provenance.plugins`: the workflow's plugins in execution order,
    each with `id`, `family`, and method `version`.
- Plugin-first architecture, phase 1 (ADR-0069, PROP-0002, SRS-VAR-014):
  alignment and variant calling consume a modality-neutral per-read evidence
  contract (`ReadEvidence`) instead of Sanger types. The Sanger adapter in
  `read_processing` owns the peak, relative-quality, and mixed-signal vetoes
  and the phase-state mask reasons. The module validator rejects any dependency
  from a neutral module on a Sanger module. Every result document is unchanged.
  **Breaking (Rust API):** `Error` gains `Evidence(EvidenceError)`.

- Sample variants publish opposition evidence (ADR-0068, SRS-SAMPLE-028): the
  admitted reads whose callable reference segments cover the variant's evidence
  span without supporting it, with orientation counts. It is evidence only;
  variant calling and notation are unchanged. **Breaking (unreleased, revised
  in place):** `variants[].opposition` is required in
  `dna.sample_evidence/v10`.
- The read-callability mask acts (ADR-0067 increment 2, SRS-CALL-009 to
  SRS-CALL-012). The trim interval is the callable span plus up to
  `read_end_margin` adjacent dephased calls (`dna.callable_span_trim/v1`
  replaces the best-section walk-out). Alignment keeps dephased calls with
  their profile and aligns other masked calls as unresolved. Masked evidence
  calls exclude a variant with the reason of their segment, replacing the
  sequence-only `post_homopolymer` window, and a masked call is no SNV
  candidate. Masked sample observations are kept as `masked` and never retain a
  locus. A mask that starts in the window after a repeat run is pulled back to
  the run's end. `read_end` is measured from every uninformative call,
  including internal unresolved masked segments, and the shipped
  `weak_amplitude_fraction` rises to 0.2. On the local corpus precision rises
  from 0.949 to 0.979, recall from 0.970 to 0.977, and the reads that failed
  the identity gate now
  place. **Breaking (unreleased, revised in place):** configuration schema 7
  removes `quality_control.best_section_fraction`, `trim_stringency`,
  `minimum_retained_bases`, and `variant_calling.post_homopolymer_window`,
  moves `variant_calling.homopolymer_min_length` to
  `callability.repeat_min_length`, adds `callability.minimum_callable_calls`,
  and ships `read_end_margin = 12`; `dna.analysis/v9` and
  `dna.sample_evidence/v10` alignment summaries gain `masked_bases` and
  `callable_reference_segments`; sample evidence gains the `masked` locus
  state, `support_topology.masked_reads`, the exclusion reasons
  `dephased_signal`, `mixed_signal`, `weak_signal`, and `irregular_spacing`, and
  the required `rejected_reads[]` array.
  A `sample` read with too few callable calls is recorded in the new
  `rejected_reads[]` registry, contributes nothing else, and no longer fails
  the operation. **Breaking (Rust API):** `QualityControlError` replaces
  `TooFewCalls` and `RetainedTooShort` with `InvalidCallableSpan`,
  `CallabilityError` gains `TooFewCallableCalls`, and `SampleError` gains
  `NoAdmittedReads`.
- Read callability separates dephased from mixed double-peak segments with a
  shadow model (ADR-0067 calibration amendment): a non-negative least-squares
  fit of each segment's normalized amplitudes to the read's own primary calls
  shifted by −3..+3 calls replaces the single modal shift offset, so two-sided
  slippage is dephased rather than mixed. Segment boundaries, the mask, and the
  callable span are unchanged. **Breaking (unreleased, revised in place):**
  `callabilitySegment` in `dna.basecalls/v3`, `dna.analysis/v9`, and
  `dna.sample_evidence/v10` gains `shadow_offsets`, present exactly on dephased
  segments; configuration schema 7 replaces `callability.shift_coherence` with
  `minimum_main_share`, `maximum_far_share`, and `minimum_shadow_share`.
- **Removed:** the Python batch runner `tools/python/scripts/analyze_samples.py`,
  its `SRS-BAT-*` requirements, and the batch runbook. Python in DNA is limited
  to repository checks, tests, measurement, and research (ADR-0066); running
  DNA over a corpus is the job of downstream pipelines.
- Result comparison with external call sets, reviewer truth, or other tools is
  out of DNA's scope (ADR-0065, SRS-COMPAT-003): downstream pipelines convert
  DNA's published results into their canonical sample representation and
  compare there.
- **Breaking (Rust API):** the canonical called-variant contracts move from
  `dna::variant_analysis` to the core module `dna::variant` (`Variant`,
  `VariantKind`, `ReferenceIdentity`, `CalledVariantSet`). Modules now form
  five acyclic layers (ADR-0064), enforced in CI by
  `validate_module_layers.py`, so a later crate split needs no contract redesign.
- **Breaking (Rust API):** `variant_nomenclature::mtdna::apply_hv2_polyc` is
  replaced by `apply_control_region`, which adds HVS-II `311T 315.1C`,
  `310C 315DEL` and EMPOP `315.1C`, HVS-III `513A 523DEL 524DEL`, and HVS-I
  `16183C 16184A 16189C` / `16189C 16193DEL` (SRS-NOM-013 to SRS-NOM-015).
  `NomenclatureError` window failures name their window.
- **Breaking (Rust API):** the empty public `dna::config` and `dna::model`
  modules are now private; library-produced results (`VariantAnalysisResult`,
  `VariantNormalizationResult`, `VariantNomenclatureResult`) and the enums
  `VariantKind`, `NormalizationPolicy`, `cli::Command`, and `error::Error` are
  `#[non_exhaustive]`.
- Operational logging uses `tracing`: scientific stages emit structured events
  and per-stage spans, and a DNA layer renders them in the unchanged
  per-operation record format (ADR-0061). Library callers can observe stage
  progress through their own `tracing` subscriber.
- **Breaking (Rust API):** `dna::error::Error` stage variants now wrap typed,
  `#[non_exhaustive]` per-stage failure enums (`AbifError`, `FastaError`,
  `ConfigError`, `AlignmentError`, and so on) instead of `String` messages;
  `DNAProcessing` is renamed `Signal`; `Path::reason` is `&'static str`; and
  `ConfigParse`/`Serialize` no longer expose `toml`/`serde_json` types. CLI and
  log error text is unchanged, except two internal-invariant messages that
  valid input cannot reach: `VariantError::NoSupportingCalls` uses the
  uppercase variant label, and a non-ASCII rendered allele reports
  `RepresentationError::InvalidAllele`.
- CI documentation now matches the workflows: CodeQL, fuzzing, MSRV, release
  build, and RustSec audit run on `main` and on schedules, not on pull requests.
- Corrected the minimum supported Rust version to match language features used
  by the codebase.

### Fixed

- A read's clean start or end that holds every run of a stretch's
  composition no longer stands for the whole stretch when the read reads
  further runs there. In the HVS-I C-tract a read reading `C7 T C T C` had the
  clean end `C T C`, the reference composition, and its inserted `C` gave the
  first C run a length of one: a one-read seven-base deletion that
  auto-passed (`16185T` with `16187-16193DEL` against the truth's
  `16189C 16191Y 16191.1C 16192T`). Plugin `consensus` 10. Development panel
  v5 against `main` 3bbf194: samples with an unreviewed difference 4.13 to
  4.11 % on current data, non-heteroplasmic 0.91 to 0.88 %, unchanged in the
  other periods; four samples lose an unreviewed block and one gains a
  heteroplasmic one (`16192Y`).
- A run that no read shows both ends of, but that two reads count alike from
  opposite ends, now takes their count (`length_evidence` `confirmed_end` in
  `dna.consensus/v1` and `dna.review/v1`; plugin `consensus` 9). Each read
  enters the run in phase, reads all of it unmasked, and loses phase at the
  other end on a call that does not read the run's base. Such runs took the
  reference frame before, mostly the HVS-II run `303-309` with the forward
  read losing phase at `310`: DNA wrote `C7` where both strands counted
  `C8`. Development panel v5 against `main` a19d69d: 133 samples from
  discordant to concordant (strict) and none worse; samples with a clear
  review 65.6 to 66.1 % on current data; samples with an unreviewed
  difference unchanged in every period (4.13 % on current data).
- A run read through its exit by one read no longer swallows the exit base
  that another read shows cleanly. A stretch's composition came only from
  reads that read it whole; a read that reads only an edge of it cleanly was
  left out, so a forward read losing phase in the HVS-I C-tract decided
  `16194C` against a reverse read showing `…CA` in phase. Such a read now
  opposes a composition its clean edge contradicts (plugin `consensus` 8).
  Development panel v5 after #142: `haplotype_lod` 11 better and 3 worse
  (HVS-III AC-repeat stretches, now undecided check items), samples with an
  unreviewed difference 4.14 to 4.13 % on current data (4.25 to 4.23 % on the
  development half) and unchanged or lower in every period.

- `read_disagrees` (`review` plugin 3, rules version 3) read a read's
  sequence from its insertion variants, which leave out masked inserted
  calls, so a phase-shifted read whose masked calls agree with the decision
  looked one base short (`2510457`: the forward read's masked `315.1C`). It
  now reads the read's own alignment columns, and counts a read only when it
  also reads both sides of the item cleanly, so it bounds what it shows. On
  development panel v5 the label falls from 441 samples to 6, `clear` returns
  to 6,460, and leaks among clear to 238: the four leaks rules 2 caught were
  caught by the defect.
- A deletion neither of whose flanking calls reads a base (both `N`) made
  `sample`, `consensus`, and `SampleAnalysis` fail for the whole sample
  ("variant has no resolved call"). Such a deletion observes nothing, since
  its gap could lie anywhere among the unresolved calls; it is now recorded
  only as excluded, with `unresolved_flank` (SRS-VAR-018). Plugin `core` 10.

- Variant normalization: a deletion whose one-base 3' shift would overlap a
  substitution of the next base now tries the deletion trade instead of
  stopping, so `16297DEL 16298C` over `TT` is written `16297C 16298DEL`, as
  any other trade (SRS-VN-008). `normalization` is version 3.
- The kernel refuses `run_boundary` as a modality label, as it refuses the
  core's other exclusion reasons.

- Alignment places every gap at its 3'-most sequence-equivalent placement
  whatever the read's profile scores (ADR-0047 §1 amended, SRS-ALN-012): the
  optimal score chooses the explanation, and a move is kept when the
  primary-base score is unchanged. With profile-weighted scores two reads of
  one sequence placed an AC-repeat insertion at `513.1C 513.2A` and at
  `524.1A 524.2C`, and the consensus lost it (118 contested AC sites on
  development panel v5, now 4). The one-base move is shared with variant
  normalization (`dna_kernel::model::placement`). An insertion into a repeat
  the read bounds is judged for support vetoes on its placements along the
  repeat, in variant calling and the consensus alike (SRS-VAR-020), so a
  weak call that happens to be taken as the inserted one no longer decides.
  The `core` plugin version rises to 8 and `consensus` to 4. Development
  panel v5 against main: 232 more samples concordant on haplotype LoD and
  227 on the primary measure; the one LoD regression is a truth C-shadow
  listed for re-check.
- `dna call` and `SampleAnalysis::call` no longer report a position DNA
  excluded a call at as callable (SRS-VAR-019). A record's callable
  reference segments leave out every position of an excluded variant (a
  substitution's position, a deletion's deleted positions, both positions
  beside an insertion's junction), so an excluded call is undecided rather
  than read as the reference, and a record excluded for `placement_conflict`
  contributes no callable segments at all; per-read sample evidence is
  unchanged. The `core` plugin version rises to 7. Development panel v5:
  102 more samples concordant on haplotype LoD and 97 on the primary
  measure, no regression; 287 samples lose 441 callable positions inside the
  HV regions (at most 7).
- `call` writes an AC-repeat deletion that a consensus decided at `523–524`
  there, rather than as `513DEL 514DEL` beside the unresolved bases the
  consensus left at `513–516`. On a reviewed consensus record, whose ends are
  vouched, a deletion beside an `N` without evidence now moves 3' across it
  and the repeat beyond (SRS-ALN-017): the `N` scores the same against any
  base, so the traceback's placement there was arbitrary and could delete a
  base the consensus kept. A variant with a mapped call read as `N` is no
  longer reported, with the new core reason `unresolved_flank`
  (SRS-VAR-018); the exclusion-reason vocabulary of
  `dna.sample_evidence/v10` and `dna.variants/v1` gains it (revised in
  place), and the `core` plugin version rises to 6. The same reason also
  removes three deletions beside an `N` that no truth has (`16197DEL`,
  `533–535DEL`). Development panel v5: 29 more samples concordant on the
  primary measure and 31 on haplotype IUPAC-compatible, no regression.
- `call` no longer reports phantom variants from a consensus record placed
  away from the interval its name declares (SRS-VAR-017). A mostly
  undecided consensus segment named for HVS-I's 5' flank could align onto the
  real HVS-I segment and call its few resolved bases there; on the
  development panel (v5) 24 samples carried 115 such calls, which triage
  could not see. Every variant of such a record is now excluded with the core
  reason `placement_conflict`; the `core` plugin version rises to 5.
- **Breaking:** callable reference segments (`callable_reference_segments`)
  no longer include positions read by an unresolved call (`N`, or an IUPAC
  code a modality presents as `N`), masked or not (SRS-ALN-015). A reviewed
  or sample consensus that writes `N` where no read observes a position
  cleanly therefore no longer reports those positions as callable, so a
  comparison counts truth alleles there as uncovered instead of missed: a
  reclassification of DNA's abstentions, not a change of any call. A read
  with an unresolved call in a variant's evidence span no longer counts as
  opposing it.
  `dna.analysis/v9`, `dna.sample_evidence/v10`, and `dna.variants/v1` are
  edited in place (unreleased cycle); the `core` plugin version rises to 3
  ([versioning](docs/governance/versioning.md)). Calls are unchanged.
- A trace with fewer calls than one signal window no longer fails a whole
  `dna sample` or `dna consensus` (SRS-SIG-007, SRS-SAMPLE-030).
  - Read processing checks the call count before any stage runs. Such a read
    is recorded in `rejected_reads[]` with the new reason `too_short`
    (`calls`, `minimum_calls`), and the other reads continue.
  - `analyze` and `basecall` still fail typed with "N calls are fewer than
    window_size_bases W". A one-locus trace now fails the same way, instead
    of in locus-window geometry during basecalling.
  - On the 2620-sample dev panel, the 2 samples that failed now produce
    results (3 reads rejected as `too_short`); every other consensus
    document and FASTA differs only in the `sanger` plugin version.
  - The `sanger` plugin method version rises to 2, because a short read
    that produced an error now produces a rejection
    ([versioning](docs/governance/versioning.md#plugins)). Every document
    that lists plugins in its provenance records it.
  - **Breaking:** `dna.sample_evidence/v10` and `dna.consensus/v1` are
    revised in place. The rejection reasons gain `too_short`, and a
    `too_short` record in sample evidence has no `integrity` or
    `callability`.
  - **Rust API:** `dna_sanger::model::signal::ShortRead` (with
    `ShortRead::of`, which reads the short-read error) is new;
    `SignalError::TooFewLoci` is removed, because the one length check runs
    before signal processing;
    `dna_sanger::model::attachment::SangerRejection` no longer carries the
    trace name and SHA-256. In `dna::sample_analysis` (new in this cycle),
    `RejectedRead` carries `input_name` and `input_sha256` beside its
    `cause`, and `RejectionCause` gains `TooShort(ShortRead)`; its
    `Callability` and `Placement` variants carry the read's `sanger`
    evidence and their `rejection`, which `RejectionCause::sanger` returns.
- A trace whose `PLOC.2` repeats one position (seen at the first or last
  locus of production traces) no longer fails the whole `analyze`,
  `basecall`, or `sample` operation (SRS-IN-003, ADR-0031). Decode merges the
  repeated position into its predecessor and drops the vendor entries at its
  index; a decreasing position is still rejected, now with "PLOC.2 positions
  must not decrease". **Breaking:** trace integrity in `dna.analysis/v9`,
  `dna.basecalls/v3`, and `dna.sample_evidence/v10` gains the required
  `duplicate_ploc_positions` count (revised in place).
- An indel flanking call whose primary base is unresolved (`N`) is omitted from
  public calls (SRS-VAR-006). Previously a tied flank aborted `analyze` or a
  whole `sample` operation, and a mixed-signal flank was published with the
  schema-invalid base `N`.
- PLOC loci one sample apart, which SRS-IN-003 accepts, no longer fail
  basecalling with an empty locus window.

## [0.1.0] - Unreleased

Initial pre-production development baseline.
