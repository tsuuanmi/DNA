---
implements: [SRS-ALN, SRS-API, SRS-CALL, SRS-COMPAT, SRS-CONS, SRS-NOM, SRS-SAMPLE]
decided-by: [ADR-0065, ADR-0075, ADR-0076, ADR-0077]
---

# Evaluation of the review workload

Part of the [[validation/evaluation-harness|evaluation harness]]: the review workload and its passes, the estimates, and the anchors.

## Review workload

The review workload says how much of each sample DNA leaves to a reviewer,
and which errors no reviewer would see. The review items are `dna review`'s
own, read through `SampleAnalysis::review_list` (SRS-API-015), so the
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
unreviewed differences, and the effect of the forced windows.

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
- the unreviewed blocks by the zone kind of their first position.

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

**Unreviewed differences.** A difference block is decided per differing
position: the rule positions of its `haplotype_lod` comparison inside it
(correctness audit F2, approved 2026-10-21). A position is reviewed when a
check item touches its reach: the position, or for an indel the reference
run of its anchor base, where the indel may be written elsewhere. A block is
unreviewed when any position is, and it counts once in the zone of each of
its unreviewed positions, with the most severe class of those positions.
Each unreviewed position takes the first class that applies:

1. `re-check`: the re-check list (`--recheck`) names a span of its sample
   that holds the position;
2. `truth N`: every difference of the block is an insertion past the
   truth's analysed interval: the block reaches the end of the truth's own
   interval, holds only its last base, and one side reads that base as the
   other does, then inserted bases the other has not (DNA's `573.1C
   573.2C` where the truth's HVS-III interval ends at 573). The truth
   records nothing past its interval, and the allele measure leaves such an
   insertion uncompared too. A block whose interval ends only where DNA's
   callable interval ends is not this case;
3. `length note`: the block lies in a profile length tract
   (`length_tracts`), both sides read the same sequence once the tract's
   base is taken out, and the count differs by one (user decisions
   2026-10-21; validity review, recommendation 6). Forensic practice notes
   such a length heteroplasmy and leaves it out of comparisons and searches
   (ISFG 2014 Rec. 10, 16), so it is reported and never gated. A difference
   of two or more (a homoplasmic C7 called C9 is DNA's error), a structure
   difference (`T310C`, a `16189C` one side lacks), and `315.1C`, outside
   HVS-II's tract, stay non-het;
4. `truth N`: the truth writes `N` at the position (the reviewer
   abstained);
5. `het`: the truth writes an IUPAC code of two or three bases (`R`, `Y`,
   `K`, `M`, `S`, `W`, `B`, `D`, `H`, `V`) at the position;
6. `non-het`: any other, a DNA error that neither a heteroplasmy, a truth
   `N`, a tract length, nor the re-check list explains.

**Omitted truth alleles.** A truth allele inside the truth's analysed
interval and outside DNA's callable intervals (the strict measure's
`uncovered`) makes no block: the report omits it (correctness audit F1).
When no check item touches its reach, it is an unreviewed difference of its
own, in the zone of its anchor:

- `omitted (not analysed)`: in a stretch at a region's end that DNA did not
  analyse (SRS-CONS-012), mostly a poor read start (`73G`);
- `omitted tract`: in a profile length tract under a `tract_length` note
  whose reading (or, without one, the convention's reading) is the truth's
  count of the tract's base there;
- `note misread`: there, when the note's reading is one unit off the
  truth's count, or the note gives none: a length note (user decision
  2026-10-22); two or more units off is `non-het`;
- `omitted (other)`: anywhere else.

The re-check list applies to an omission as to a position.

A passed sample, region, or zone with unreviewed differences counts once,
under its most severe one, in the order non-het, omitted (other), omitted
(not analysed), truth N, het, re-check, note misread, omitted tract, length
note (`non_het`, `omitted_other`, `omitted_not_analysed`, `truth_n_only`,
`het_only`, `recheck_only`, `note_misread_only`, `omitted_tract_only`,
`length_note_only`). `unreviewed.tsv` lists each difference with its class
and its class without the re-check list (`class_unlisted`), so `diff` can
report its gate both ways. The summary's "Length notes" table gives, per
tract, the items the review notes and the unreviewed length-note blocks,
split by the truth's length: the note's convention reading, DNA's length
less one, or another (`notes.tsv` lists every noted run). The gated classes
are `diff`'s option ([[validation/evaluation-harness#Comparing two runs]]):
by default non-het and the omissions outside a noted tract. Het is reported
but counts as a pass, and the truth N and re-check units await a check of
the truth, not of DNA. The sample measure "Samples with an unreviewed
difference" and the TSV tables (`unreviewed` in `scores.tsv`, `windows.tsv`)
count every unreviewed block whatever its class, and no omission.

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
whole panel, and `dna-eval diff` gates on them alone. Without it, nothing
changes.

## Estimates

For each weighting (`corpus` = `weight_corpus`, `pool` = `weight`), each
measure with the metrics `not_concordant` (anything but `Y`, failures
included) and `discordant`, and `failed` (measure `all`), over samples of
positive weight that are not retired:

- rate `p = Σ w·y / Σ w`, per stratum and over all;
- stratum size `N_h = Σ_h w`, `W_h = N_h / N`;
- variance `Σ_h W_h² · p_h(1 − p_h) / n_h · max(0, 1 − n_h / N_h)`;
- 95 % interval `p ± 1.96·√variance`, and the unweighted rate.

## Anchors

Anchor candidates are the manifest's `anchor_candidate = 1` samples:
`anchors.tsv` lists each under every measure, and the summary gives how many
are concordant under each. Two runs are compared with `dna-eval diff`
([[validation/evaluation-harness#Comparing two runs]]), whose status changes
read statuses by class so that runs of other conventions stay comparable:
`exact` and `concordant` are good; `error` and `discordant` are bad;
`not_applicable` and `failed` are their own classes.
