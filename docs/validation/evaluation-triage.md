---
implements: [SRS-ALN, SRS-API, SRS-CALL, SRS-COMPAT, SRS-NOM, SRS-SAMPLE]
decided-by: [ADR-0065, ADR-0075, ADR-0076, ADR-0077]
---

# Evaluation of review triage

Part of the [[validation/evaluation-harness|evaluation harness]]: the review workload, estimates, anchors, and the draft triage evaluated against reviewer truth.

## Review workload

The review workload says how much of each sample DNA leaves to a reviewer,
and which errors no reviewer would see. The review items are `dna review`'s
own, read through `SampleAnalysis::review_checks` (SRS-API-015), so the
harness applies no review rule of its own.

- The HV regions split into windows of `--review-window` positions (default
  10), aligned to each region's first position; a region's last window may
  be shorter. A reviewer looks at every window a check item touches: those
  windows' positions are the sample's positions to review
  (`review_positions`). With a window of 1 they are the positions inside
  check items.
- A `haplotype_lod` differing block that overlaps no check item is
  unreviewed (`unreviewed`): an error the report would carry without anyone
  looking. The safety measure is per sample, the share of samples with an
  unreviewed block, since one such block makes the report wrong.
- `--forced-windows` lists windows always sent to review, one `start-end`
  (1-based, inclusive) per line; the workload is reported again with them
  reviewed (`unreviewed_forced`). Choosing them is a policy decision:
  `windows.tsv` only ranks windows. A list must be learned on one run and
  judged on a hold-out half or later batches; applied to the run it was
  learned on, it overstates its effect.

Only samples compared under `haplotype_lod` with a concordant or discordant
status count. The summary's "Review workload" section gives the positions to
review, the samples with an unreviewed difference overall, by region, and by
batch month (from the manifest's optional `batch_date` column), the clear samples (no check item), the windows with the most
unreviewed differences, the effect of the forced windows, and with a baseline
that reports the workload, its values on the same samples.

**Zones to check.** A reviewer opens a functional part of a region as a
whole, so the workload is also counted in zones. Every position of a
reportable region of the configuration's target profile has one zone kind,
by precedence:

1. `window`: a repeat window, one of the profile's nomenclature windows
   (polyC of HV1 and HV2, the AC repeat of HV3);
2. `after`: within `--zone-after` positions (default 20) of a window on
   either side, since the harness does not see a read's direction;
3. `start` and `end`: the first and last `--zone-edge` positions (default
   40) of the region, its amplicon's start and end;
4. `interior`: the rest.

A zone is a kind in one region (`kind@region`, such as `after@HV2`),
contiguous or not. A check item counts for every zone it touches. An item
that overlaps a window counts its positions within `--zone-overrun`
positions (default 2) of the window for the window, since the consensus
often decides a repeat as one stretch a position or two past it. The
summary gives:

- zones per sample (mean, median, 90th percentile);
- the samples with none;
- the check items that touch more than one zone;
- the samples with a check in each zone;
- the unreviewed blocks by the zone kind of their first position;
- with a baseline that reports zones, its zones per sample on the same
  samples.

The profile is the one the configuration names, and every analysed sample
must have used it. The widths are in `run.json`.

**Passes.** A reviewer works per region and per zone, so the workload also
counts what passes without review: a sample, a region, or a zone that no
check item touches. Over the same samples, the summary's "Passes" table and
`run.json` `workload` (`all`, and `current` with `--current-since`) give:

- auto-pass samples: those without a check item (`review_clear` `Y`), and
  those of them with an unreviewed block;
- regions: one per sample and reportable region of the profile. A region
  needs a check when a check item touches one of its zones; the summary gives
  the regions needing a check and the pass share, overall and per region, and
  the passed regions that hold an unreviewed block, which counts for the
  region of its first position;
- zones: one per sample and zone (a kind in a region, as listed above), and
  those a check item touches; the summary gives the zones needing a check and
  the pass share, overall, per zone, and per zone kind (summed over the
  regions), and the passed zones that hold an unreviewed block, which counts
  for the zone of its first position. A zone a check item touches is not
  passed, even when an unreviewed block in it lies beside the item.

A passed zone may lie in a region that needs a check, so the passed zones
with an unreviewed block can outnumber the passed regions with one; and they
count zones, where the zone kind table above counts blocks. A block outside
every region counts for its sample only.

Each unreviewed block has one class, the first that applies:

1. `re-check`: the re-check list (`--recheck`) names a span of its sample
   that the block overlaps;
2. `truth N`: the truth writes `N` in it (the reviewer abstained), or every
   difference of the block is an insertion past the truth's analysed
   interval: the block reaches the end of the truth's own interval, holds
   only its last base, and one side reads that base as the other does, then
   inserted bases the other has not (DNA's `573.1C 573.2C` where the truth's
   HVS-III interval ends at 573). The truth records nothing past its
   interval, and the allele measure leaves such an insertion uncompared too.
   A block whose interval ends only where DNA's callable interval ends is not
   this case;
3. `het`: the truth writes an IUPAC code of two or three bases (`R`, `Y`,
   `K`, `M`, `S`, `W`, `B`, `D`, `H`, `V`) in it;
4. `non-het`: any other block, a DNA error that neither a heteroplasmy, a
   truth `N`, nor the re-check list explains.

A passed sample, region, or zone with unreviewed blocks counts once, under
its most severe block, in the order non-het, truth N, het, re-check
(`non_het`, `truth_n_only`, `het_only`, `recheck_only`): one non-het block
makes it non-het, and it is re-check only when every block is on the list.
The non-het count is the gate a change must not raise. Het is reported but
counts as a pass, and the truth N and re-check units await a check of the
truth, not of DNA. The sample measure "Samples with an unreviewed
difference" and the TSV tables (`unreviewed` in `scores.tsv`, `windows.tsv`)
count every block whatever its class.

`--recheck <file.tsv>` names the difference blocks whose truth is under
re-check: raw-trace adjudication found the truth wrong or questionable, or
a blind re-check of the truth is pending. The file is tab-separated with the
header `sample_id`, `first`, `last`, `reason`, `source`: per row a sample,
the first and last reference positions of the block (1-based, inclusive, as
in the block notation of `unreviewed`), why its truth is under re-check, and
where the adjudication is recorded. Every cell must be filled and `first`
must not exceed `last`; another header, a row of another width, a position
that is not a number from 1, or a repeated (`sample_id`, `first`, `last`)
stops the run. `run.json` records the list's path, SHA-256, and rows
(`recheck`); the summary names it and gives, with the passes, the
unreviewed blocks on it and the rows that overlap no difference block of
the samples counted, which flags a stale row. The list holds sample
identifiers and stays outside the repository under the
[[governance/data|data policy]]. Without it, no block is re-check.

`--current-since <YYYY-MM>` names the first batch month of the current data:
the samples whose `batch_date` is in or after it. Recent batches come first,
since old batches may predate the laboratory's current practice. The
summary then reports the review workload on the current data before the
whole panel, and with a baseline, the changes of status class per measure on
the current data after those of the whole panel. Without it, nothing changes.

## Estimates

For each weighting (`corpus` = `weight_corpus`, `pool` = `weight`), each
measure with the metrics `not_concordant` (anything but `Y`, failures
included) and `discordant`, and `failed` (measure `all`), over samples of
positive weight that are not retired:

- rate `p = Σ w·y / Σ w`, per stratum and over all;
- stratum size `N_h = Σ_h w`, `W_h = N_h / N`;
- variance `Σ_h W_h² · p_h(1 − p_h) / n_h · max(0, 1 − n_h / N_h)`;
- 95 % interval `p ± 1.96·√variance`, and the unweighted rate.

## Anchors and baseline

Statuses fall into classes so that runs of other conventions stay comparable:
`exact` and `concordant` are good; `error` and `discordant` are bad;
`not_applicable` and `failed` are their own classes.

- Anchor candidates are the manifest's `anchor_candidate = 1` samples, listed
  under each measure. With a baseline, the anchor set is the candidates whose
  baseline `status` was good.
- Changes against the baseline, per measure the baseline reports:
  `regression` (good before, not now), `fix` (the reverse), `changed` (another
  change of class), `added`, `removed`. The summary also tabulates the
  baseline `status` against the primary status.

## Draft triage

`run --triage` measures whether the triage of
[[PROP-0004]] is worth building (phase 0).
Its detectors are research drafts over the evidence DNA already produces:
nothing in DNA reads them, they decide no call, and none is approved.

**Evidence.** Per sample, from the analysis the run already made
(`analyze_sanger_outcome`), kept even when the call then fails:

- the reportable regions of the target profile (`SampleAnalysis::regions`);
- per region position, its clean and masked reads from one source: at a
  consensus site, the site's supporting and opposing reads observe it cleanly
  and its uninformative reads are masked; elsewhere, a read observes it
  cleanly inside its callable reference segments and masked in the rest of its
  reference segments. Clean reads are counted per strand. Sites name reads by
  trace name; a name that two admitted reads share counts a strand only when
  the site names every read of that name or all of them share the strand,
  otherwise as masked;
- the marks of the consensus sites over each position;
- per region position, each read's *minor channel* where that read calls a
  canonical base without a mask (every covered position, through
  `SampleAnalysis::observations`, read once over the span of the regions):
  the strongest channel of its reference-oriented
  A/C/G/T profile other than its own call, kept when its weight is at least
  0.05; the pair (call, minor base); whether the minor base is the base the
  read calls just before or just after in sequencing order (a neighbour's
  shadow); how systematic each of those two shadows is in the read here,
  the share of its calls within 10 on either side (the call itself left
  out) whose minor channel, at 0.05 or more, is the base of the call before
  (or after) them; and the read's local noise, the median minor weight of
  its calls within 10 calls on either side. Weights and shares are kept in
  thousandths, rounded down for shares;
- per input trace: its name, the primer in the name (a `HV<n><F|R>` token),
  whether it was admitted (with its strand) or rejected (too short, too few
  callable calls, or unplaced), and its callable span in calls.

When every read is rejected, the evidence is the traces alone, with no region
evidence. When DNA cannot give an analysed sample's observations, the
stage-1 evidence is kept without minor channels: the sample's call and scores
stand, a policy with a minor-signal rule sends its would-pass regions to
review because the gate could not run, and `triage.tsv` gives the cause
(`minor_unavailable=…`).

**Detectors.**

| Level | Finding | Fires when |
|---|---|---|
| trace | `trace.too_short`, `trace.weak`, `trace.unplaced` | the read was rejected for that cause (SRS-SAMPLE-030, SRS-SAMPLE-027, SRS-SAMPLE-029) |
| trace | `trace.short` | the callable span has fewer than `short_calls` calls (default 100; `0` turns it off) |
| region | `region.unreliable` | a position is `n` (a site without a call or with `N`), `contested`, `masked` (a site with an uninformative read), or `phase_loss` or `reference_frame` (a run of the site has that length evidence); `masked_read` (a covering read reads the position only through a masked call) is measured but not applied by default |
| region | `region.one_strand` | a position is observed cleanly by one strand only, where the policy requires both strands: everywhere (the draft), nowhere (`off`), or only inside the repeat windows ±3 (`outside_repeats`, R6m) and also in or next to a reference homopolymer of at least `k` bases (`outside_repeats_h<k>`, R6m-h4 with 4); the windows are the measures' windows (`compare::Region::WINDOWS`) and the runs are read from the reference |
| region | `region.minor_signal_single_strand` | a position one strand alone reads cleanly shows a minor signal under the policy's single-strand rule (off by default); it applies only where the policy lets one strand pass |
| region | `region.gap` | a position is observed cleanly by no read |
| region | `region.minor_signal` | a position's reads show a second allele under one of the policy's minor rules (default `upstream@100`; off turns it off) |

**Minor rules.** A rule is a kind and a threshold `t` in thousandths; a read
counts when its minor weight (its second channel's share of the read's
normalised A/C/G/T profile, which sums to 1, rounded to thousandths, so 0.0995
counts at `t` = 100) is at least `t` (`>=`) and it passes the kind's shadow
test. The shadow tests compare the minor base with the read's *adjacent
calls* in reference order, not with the reference bases. The cross-strand
kinds mark a position when a forward and a reverse read that count show the
same unordered pair:

| Kind | A read counts when | Marks when |
|---|---|---|
| `strands` | always | a forward and a reverse read show the same pair |
| `unshadowed` | its minor base is neither neighbour's base | as `strands` |
| `upstream` | its minor base is not the base of its 3′-side adjacent call in reference order (the call after it on a forward read, the call before it in sequencing order on a reverse read) | as `strands` |
| `noise` | its minor weight is at least 3 times its local noise | as `strands` |
| `local_shadow` | its minor base is not a neighbour's base whose shadow at least 0.05 of the read's nearby calls show | as `strands` |
| `any` | its minor base is neither neighbour's base | any read counts |
| `lone` | as `unshadowed` | as `unshadowed`, or any read counts where one strand alone reads the position cleanly |

A minor signal never asks for a rerun.

**Blind spot by construction.** A cross-strand kind marks a position only when
a forward and a reverse read count, so none of them can fire where one strand
alone reads cleanly. `upstream` can never mark a heteroplasmy whose alternate
base equals the 3′ neighbour (next reference base), and `unshadowed` and `any`
also miss one equal to the 5′ neighbour. In the merged reviewer truth of about
103,000 samples, 9,368 reference-plus-one-alternate IUPAC sites split as
alternate = previous reference base only 15.5 %, next only 11.3 %, both
10.4 %, neither 62.8 %. `upstream` is therefore blind to about 21.7 % of
corpus heteroplasmy sites (common examples: 146Y, 16189Y, 16183M) and
`unshadowed` to about 37.2 %. The leak bound of a policy with these rules
holds only for the classes they can see; [[PROP-0004]]
gives the dev-panel split.

A trace with any trace finding is *poor*; it stays a diagnostic unless a
region asks for it. A sample without an admitted read is `rerun` with every
poor trace, the clearest rerun case. A sample DNA cannot analyse (a trace that
cannot be read, an invalid input) or whose consensus cannot be called is its
own class, `failed`; after a call failure its trace findings are still
reported.

**Outcomes.** A region is `rerun` when it has a gap, or a one-strand position
while the policy applies `region.one_strand`, and a poor trace could have read
it; `review` when any applied finding remains; `pass` otherwise. A sample takes
its worst region (rerun > review > pass). Which poor trace could have read a
region is a proxy:

- `primer` (default): the poor trace's primer reads the region on a strand
  the region lacks. A primer's strand is its `F` or `R`; its regions are those
  that at least half of the admitted traces with that primer in the run's
  active samples observe cleanly at 50 positions or more (its footprint). A
  sample's own traces are not left out: a poor trace never contributes to a
  footprint, and a sample's admitted traces are a negligible share of a
  primer's;
- `any`: any poor trace of the sample.

**Measurement.** Retired samples are left out; samples without a comparison
count in A/B/C but not in leaks or catches.

- **A/B/C:** the share of samples in each outcome, unweighted and
  corpus-weighted with the stratified estimator (with its interval); and the
  share of regions of decided samples with region evidence, unweighted and
  weighted by their sample's `weight_corpus` (a point estimate).
- **Leaks:** passed samples that are not concordant under each measure, with
  the one-sided 95 % Clopper–Pearson upper bound on the leak rate among passed
  samples (`1 − 0.05^(1/n)`, about `3/n`, at zero leaks) and the
  corpus-weighted rate. Region leaks are discordant positions inside passed
  regions, each with its kind (`dna_unique`, `truth_unique`, `both`,
  `uncovered`), allele class (`iupac`, `insertion`, `deletion`,
  `substitution`), window, and the detectors that mark the position itself.
- **Catches** per region detector: not-concordant samples it flags (catches),
  concordant samples it flags (burden), and primary discordant positions it
  marks itself (hits). For a detector the policy applies, also the catches no
  other applied detector makes (unique) and the concordant samples only it
  holds back (sole burden: clearing it would pass them). For a detector the
  policy does not apply, instead what applying it would add: the
  not-concordant samples it would newly catch and the concordant samples it
  would newly hold back.
- **Sweep:** the draft and policies that each change one thing
  (`short_calls` 0, 50, 150, 200, 300; the `any` proxy; `region.one_strand`
  off; `masked` off; `masked_read` on; `n` and `contested` only; no unreliable
  cause; one-strand and `masked` off; the gap alone; the minor signal off),
  then the draft with each minor rule but `lone` (which equals `unshadowed`
  while `region.one_strand` applies) at thresholds 60, 80, 100, 120, 150, 200,
  and 250, leaving out the draft's own rule, `lone` at the same
  thresholds with `region.one_strand` off, and the draft's rule joined by
  `local_shadow` at the same thresholds;
  each with A/B/C and its leaks under every measure: a risk–coverage table.
  Every policy runs on the evidence already taken, without DNA.
- **Minor signal:** for the draft and each policy with another minor rule,
  against the same policy without the rule, decided like any policy:
  corpus-weighted A without and
  with the rule and its cost; the compared samples the rule alone moves from
  pass to review and the concordant ones among them (burden); passed samples
  and leaks under every measure without and with the rule; and the recall of
  truth IUPAC alleles (primary discordant positions of class `iupac`) that
  the rule marks itself, among those inside a region that passes without the
  rule and among those inside any region.

Phase-0 results are in
[[PROP-0004#phase-0-results|PROP-0004]].
