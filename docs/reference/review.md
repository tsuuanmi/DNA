---
implements: [SRS-API, SRS-CALL, SRS-CONS, SRS-VAR]
decided-by: [ADR-0071, ADR-0074, ADR-0076, ADR-0080, ADR-0082, ADR-0084]
covers:
  - crates/dna-kernel/src/variant.rs
  - crates/dna-normalize/src/
  - src/lib.rs
  - src/model/review_result.rs
  - src/pipeline/review.rs
  - src/pipeline/sample_call.rs
  - src/report/review.rs
  - src/sample_analysis/
  - src/variant_analysis/
---

# DNA Review JSON

`dna review <sample-id> <trace>... --reference <reference.fasta> [--signal <run.signal.json>]`, or with `--evidence <read.evidence.json>...` in place of the traces, writes one
deterministic document, `results/<sample-id>.review.json`, identified as
`dna.review/v1` (SRS-API-014, [[PROP-0004]]).
It lists the positions of one sample a reviewer should look at, each with the
evidence labels that make it quick to settle. It never states the reviewer's
conclusion.

The authoritative schema is
[`schemas/review-v1.schema.json`](schemas/review-v1.schema.json); the
synthetic example is
[`examples/review-v1.example.json`](examples/review-v1.example.json).

## Command

The command runs the in-process path of SRS-API-012 once per sample: the
sample placement pass, the consensus, and `SampleAnalysis::call` (the
consensus's decided edits, then the notation and the EMPOP query the profile
declares). Traces, reference, and configuration (`DNA_CONFIG`) are read as by
[[reference/consensus|`dna consensus`]], and reads are rejected the same way. The
operation log is `logs/<sample-id>.log`. With `--evidence <read.evidence.json>...`
in place of the traces, the command starts from the reads'
[[reference/read-evidence|read evidence documents]], which carry every
Sanger value the review reads (each call's scan position, relative quality,
second peak with its cause, and window peaks for the run signal), and writes
the same document, with or without `--signal` (SRS-API-019).

`--signal` names the sample's run signal ([[reference/signal|`dna.signal/v1`]]).
The sample's own admitted reads are taken out of it, so the sample is never
judged against itself, and `possible_heteroplasmy` and `reference_second`
judge each call against its key's run level (SRS-CALL-020). A run signal that holds none of the
sample's admitted traces is another run's and fails the command without
output, as does one made with another `flank` or `resolution` than the
configuration's. Without it, the configured floor alone sets the threshold.

Every sample gets a document. `status` says how far it got:

- `called`: the analysis and the call succeeded;
- `rejected`: every read was rejected; the document lists the reads with their
  causes and has no items;
- `call_failed`: the analysis succeeded but the call failed; `error` carries
  the failure text, and the items are those of the consensus and the reads,
  without calls.

An unreadable input, an invalid configuration, two traces with the same
reviewer-facing name (the same file stem in two directories), or an existing
target still fails the command without output.

## Fields

- `schema_version`, `sample_id`, `status`, and `error` (only for
  `call_failed`).
- `provenance`:
  - the reference `name`, `topology`, and sequence `sha256`;
  - `configuration_sha256` and the target `profile`;
  - `plugins`: `sanger`, `core`, `consensus`, `decided`, then
    `normalization` and `nomenclature` when the profile declares notation,
    `empop` when it also declares `[empop]`, and `review`;
  - `rules`: the version of the labelling rules below, equal to the `review`
    plugin version. `clear` depends on it;
  - `signal` (with `--signal`): the run signal's `run_id`, file `sha256`, and
    `selection` version.
- `reads`: the admitted reads in the sample's read order (content SHA-256),
  each with `name` (the trace stem), `sha256`, `strand` (`forward` or
  `reverse`, its placement), `admitted: true`, and `masked`; then the
  rejected reads with `name`, `sha256`, `admitted: false`, and `rejection`
  as in [[reference/consensus|`dna.consensus/v1`]]. `masked` lists, in reference
  order, the maximal stretches `{first, last, state, after_repeat}` of
  reference positions inside the reportable regions that the read's aligned
  calls cover while the modality masks its signal and the call stays masked
  (rules version 10; since rules version 11 a leading call the modality
  vouches for at a read start, SRS-CALL-014, is scored and not listed): `state`
  is the callability segment's phase state (`dephased`: superimposed ladders
  offset by slippage, as after a repeat run; `mixed`; `weak`; `irregular`),
  and `after_repeat` whether the segment starts in the window that follows a
  long repeat run. The read does not show its bases there: no call, base, or
  indel comes from such a call, except that a masked call reading another
  base than a run's may mark where that run ends (`dephased_end`). The
  stretches come from the read's own signal, never from fixed positions, so
  they locate what each read left unread for a later reading of the
  dephased tail, also where no item lists it.
  Since rules version 18, an admitted read also has `unreliable`: the maximal
  stretches `{first, last, reason}` of reportable positions that its aligned
  calls cover inside a span its modality asserts unreliable
  (`poor_start_signal`, the read's start before its signal settles,
  SRS-CALL-024).
  Since rules version 17, an admitted read also has `mixtures`: the stretches
  after a run whose signal the modality asserts is the read's own sequence
  superposed with copies offset by whole calls (ADR-0082). Each says whether
  it is `explained`: the components fit the stretch at the read's own noise
  level. Each lists its `components`, the dominant sequence (`offset` 0)
  first, each with its `fraction` of the fitted signal; in an explained
  mixture, a component placed on the reference also gives the interval `first`..`last` it
  covers, its `sequence` there on the reference strand (`first` after `last` when the interval crosses the origin), and its eligible
  edits as profile-notation `calls`. Each component is placed and called as
  a read of its own, so the dominant and the minor sequences are each
  compared with the reference independently (for example `309.1C 315.1C`
  against `309.1C 315.1C 315.2C`). `unexplained` lists the reference positions
  of calls the components do not explain. Mixtures are display only: no
  item, label, or `clear` depends on them.
- `calls` (`called`, under a profile that declares notation): the sample's
  calls in the profile notation, distinct and in reference order, as
  `SampleAnalysis::call` gives them: `call` (`309.1C`), `position`,
  insertion `ordinal`, `reference` (`null` for an inserted base), `alternate`
  (`null` for a deleted base), and the `convention` that labels it, when one
  does ([[reference/notation|notation]]).
- `callable` (`called`): the call's callable reference segments, merged where
  they overlap or touch, as 1-based inclusive `start`/`end` ranges.
- `empop` (`called`, when the profile declares `[empop]`): the EMPOP query,
  as in [`dna.notation/v1`](notation.md#empop-query).
- `reference_only` (`called`, rules version 15): the reportable positions
  (1-based, increasing) that the consensus decided as the reference by
  reference confirmations alone ([ADR-0080](../decisions/adr/0080-reference-only-calls.md),
  `reference_only` in [`dna.consensus/v1`](consensus.md)). No item lists
  them; they are DNA's own decision on weaker evidence than a clean call,
  and are listed so a reviewer can find where DNA decided that way.
- `not_analysed` (`called`, rules version 20): the stretches `{first, last,
  reason}` at the reportable regions' ends that the call does not analyse
  (SRS-CONS-012): positions only reads' unreliable spans observe
  (`poor_start_signal`), up to the last one the consensus does not settle
  as the reference. They are poor data, not a decision: no item, zone, or
  check holds them, the call reports nothing there, and the analysed range
  (`callable`, the EMPOP ranges) starts or ends after them.
- `zones` (`called`, rules version 14): the reportable positions as the
  zones a reviewer opens whole, in reference order, tiling every reportable
  region. Each position takes the first kind that holds: `repeat`, inside a
  nomenclature window of the profile (with the window's `name`);
  `not_callable`, outside `callable`; `dephased`, inside a read's `masked`
  stretch with `after_repeat`; `masked`, inside any other `masked` stretch;
  else `clean`, where every read covering it is in phase or vouched for. A
  zone `{kind, first, last, name, checks}` is a maximal stretch of one kind
  (and window) within one region; `checks` counts the `check` items touching
  it. Every boundary comes from the profile and the sample's own reads, none
  from a width, so a sample whose reads stay in phase has no `dephased`
  zone. A zone without checks is DNA's own decision.
- `items`: the review items, in priority order (below).
- `summary`:
  - `items`, `checks`, `variants`: item counts (`items` includes the
    `note` items);
  - `labels`: for every label, the number of items carrying it, zeros
    included, in the label order below;
  - `positions`: `reportable`, the positions of the profile's reportable
    regions; `check`, those a `check` item holds, which a reviewer looks at
    (rules version 5); `reference_only`, the `reference_only` positions
    (rules version 15; `0` without them); `not_analysed`, the positions of
    `not_analysed` (rules version 20). Every other reportable position is DNA's own decision
    or outside every read, as its items say;
  - `zones`: `total`, the number of `zones`; `check`, those with a `check`
    item, which a reviewer opens (rules version 14; both `0` without
    `zones`);
  - `clear`: `status` is `called` and no item is a `check`.

## Items

An item is one consensus site, or one position or run of positions outside
every site, inside the profile's reportable regions (a site is clipped to
them; since rules version 16, a clipped site's `sequence` is its decision over
the clipped columns, what most of its supporting reads show there). It has:

- `first`, `last`: its first and last reference columns, each a `position`
  and an insertion `ordinal` (`309.1` is `309`/`1`). In version 1 both
  ordinals are `0`, and an item holds the insertions after its last position
  (`315.1C` in `303`–`315`).
- `kind`: `check` when a label leaves the item unassured
  (`strands_disagree`, `read_disagrees`, `possible_heteroplasmy` unless
  `explained`, `reference_second`, `estimated_length`, `length_unconfirmed`, `undecided`,
  `not_called`, `long_run`, `edge_deletion`, `run_deletion`,
  `beside_not_analysed`, `not_covered`), whatever its calls; otherwise `variant`: it
  has calls, or, when there are no notation calls (a profile without
  notation, or `call_failed`), the consensus decided a sequence other than the
  reference; otherwise `note` (rules version 19): it holds an `explained`
  `possible_heteroplasmy`, listed for the reviewer's information only. The
  read labels `one_read`, `one_strand`, `noisy_signal`, `read_edge`, and
  `phase_shift` are context and do not make a check; `poor_signal` makes a
  check of an item with a difference only (rules version 18).
- `decision`: `called` (decided, and every position callable),
  `partly_callable`, or `not_callable` (no position callable).
- `sequence`: the decided sequence: a site's call (`N` where undecided, empty
  when every position is deleted), or the reference base of a position the
  consensus decided outside every site; `null` when nothing was decided, and
  for an item that only holds a call the notation moved out of its site.
- `calls`: the texts of the `calls` whose position lies in
  `first.position`..`last.position`.
- `labels`: the evidence labels, ordered by rank, then name, then position or
  read.
- `evidence`: for every admitted read that covers a position of the item, in
  read order: `read`, `sha256`, `strand`, and `positions`, one per covered
  reference position in order:
  - `position`, `ordinal` (`0`);
  - `base`: the read's base on the reference strand, `-` for a deletion, `N`
    when unresolved;
  - `state`: `reference`, `alternate`, `deletion`, `unresolved`, or `masked`,
    as in [sample evidence](sample-evidence/README.md);
  - `call`: the read's 0-based call index (`null` for a deletion);
  - `scope`: how far the core trusts that call (`SampleAnalysis::trust`):
    `full`, or `base` for a leading call the modality rescued, which shows
    its base only and no run length, end, or indel (`null` for a deletion);
  - `scan`: the trace scan position of that call, so a viewer can jump to the
    peak (`null` for a deletion);
  - `quality`: the call's relative quality score (`null` for a deletion);
  - `second`: the call's second base on the reference strand and its peak
    `share` (SRS-CALL-015), rounded to three decimals; `null` without a second
    peak or for an unresolved base. When the share passes
    `signal_processing.second_floor`, or the call is masked as mixed, it also
    carries `cause` (rules version 12): why the call shows a second peak,
    read from the signal (SRS-CALL-021): `noise` (near the noise floor),
    `neighbour_tail` (the tail of a neighbouring peak of that base),
    `off_centre`, `dye_blob` (a broad single-dye peak), `phase_shift` (the
    base of the previous or next call, centred), `second_template` (none of
    these), or `mixed_stretch` (masked mixed with a small own second). With
    `--signal`, a `second_template` whose share its key's run level holds is
    `context`: the run's reads show it there too (SRS-CALL-020).

Items come from three sources:

1. **Sites.** Every consensus site ([[reference/consensus|consensus]]). A site with no
   call and no label that leaves it unassured is DNA's own decision and is not
   an item.
2. **Second bases.** Every callable position outside every site where
   `possible_heteroplasmy` fires: a one-position item (a `note` when the
   label is `explained`).
3. **Gaps.** Every maximal run of positions outside every site that is not
   callable: `not_covered` when no admitted read's reference segments cover
   any of it, otherwise `undecided` plus the labels of its covering reads.
   For `call_failed`, the admitted reads' callable reference segments stand
   in for `callable`.

A call that no such item holds (for example one the notation moved out of its
site) is a one-position item of its own.

### Priority

Each label has a rank; an item's priority is the lowest rank of its labels,
and `5` for a variant without a label. Items are ordered by priority, then
`first` (position, then ordinal), then `last`.

| Rank | Labels |
|---:|---|
| 0 | `strands_disagree`, `read_disagrees` |
| 1 | `possible_heteroplasmy`, `reference_second` |
| 2 | `beside_not_analysed`, `edge_deletion`, `estimated_length`, `length_unconfirmed`, `long_run`, `not_called`, `run_deletion`, `undecided` |
| 3 | `one_read`, `one_strand` |
| 4 | `dephased_end`, `length_not_recorded`, `noisy_signal`, `not_covered`, `phase_shift`, `poor_signal`, `read_edge`, an `explained` `possible_heteroplasmy` |

## Labels (rules version 25)

Every label is computed from DNA's own evidence, never from a truth. A read's
strand is its placement orientation. A read's call is judged with the core's
trust in it (`SampleAnalysis::trust`): its mask, whether it is trusted
(informative and outside the read-end margin), and whether the bounded-run
exemption lifts the margin (ADR-0076).

- `strands_disagree` `{reads: [{read, strand, role}]}`: the site has
  supporting and opposing reads, and a supporting read lies on another strand
  than an opposing one. Lists the supporting, then the opposing reads.
  On a clipped site (rules version 16) an opposing read that shows the
  clipped decision over the reportable columns no longer opposes the item:
  its difference lies outside the regions. The label stays while an opposing
  read that may differ there, or that does not bound those columns, remains
  on another strand than a supporting one.
- `read_disagrees` `{read, strand, sequence}`: at a decided site that is not
  clipped, a supporting or uninformative read shows the whole item as another
  sequence than the decision. The read's sequence comes from its own
  alignment columns over the item, with every inserted call after a position
  of the item, masked ones included. The read shows the item only when it
  aligns every position of it and reads the reference positions on both
  sides of it, so that it bounds what it shows; each of those calls, and
  each call inside the item, must be A/C/G/T without a modality reason
  against its base, and either anchoring or scored and trusted or exempt.
  One label per such read. Opposing reads are `strands_disagree`'s.
  A read also disagrees through one end of an item it does not show whole
  (rules version 8): at a decided item without a length change, from the
  end its strand reaches first (a forward read's first position, a reverse
  read's last), every observation `reference` or `alternate` with a call of
  full trust, up to where the read has read `review.phase_run` equal bases
  in a row, one of those bases differs from the decision at its position.
  The `sequence` is then that end, in reference order. The far end of such a
  run is a length question that dephasing leaves open, so it is not
  compared. This catches a run whose exit base one read's dephased tail
  decided while the other strand reads the exit cleanly in phase (on
  current development data, `16194C` beside a clean reverse `A`).
- `possible_heteroplasmy` `{position, base, reads: [{read, strand, share,
  threshold, calls?}]}`: at one position, a read reads it cleanly when its
  observation is `reference` or `alternate`, its call is unmasked and trusted
  or exempt, and the call has a second base (SRS-CALL-015): the highest other
  peak, its share that peak's against the primary peak and itself. It
  qualifies when the share exceeds the call's second-channel threshold
  (SRS-CALL-020): `signal_processing.second_floor`, raised to the run level
  of the call's key when `--signal` gives a run signal that holds the key
  with enough calls (`calls`, the sample's own taken out), and the second
  base is not the base the same read shows at the previous or next reference
  position (a neighbour shadow). The neighbour test stays until DNA handles
  noisy regions itself. The label fires only when both strands read the
  position cleanly and a qualifying read on each strand shows the same second
  base: a site only one strand reads cleanly gives no label (rules version 6;
  on the current development data one-strand flags were 3 % precise against
  90 % for two-strand ones, so they cost reviewers far more than they found).
  The base most qualifying reads show is reported (then A/C/G/T order), with
  each of them, its share, and its threshold.
  The label is `explained` (rules version 19) when every qualifying read of
  one strand shows a second peak the modality explains by the signal's own
  shape, a neighbour's tail or a dye blob (not `second_unexplained`,
  SRS-CALL-022): that strand shows no possible second sequence. An
  `explained` label is context and makes no check. On the development panel
  none of the 23 such flags was a heteroplasmy in the truth, and none hid
  a difference from the decision, against 92 % heteroplasmies among the
  flags whose qualifying reads all show a `second_template` peak.
- `reference_second` `{position, read, strand, share, threshold, calls?}`
  (rules version 25): at a difference that rests on one read (`one_read`), a
  position where that read reads another base than the reference while its
  call shows the reference base as its second peak, with a share above the
  call's second-channel threshold (SRS-CALL-020, as for
  `possible_heteroplasmy`). One label per such position. Neither a
  neighbour's base nor the peak's cause (SRS-CALL-021) exempts it: a dye
  blob is wider than one call, and its own shape labels the real peak, not
  the blob. One read cannot tell a base from a broad peak of another dye
  over the reference's, nor from a mixture, so the item is a check. On
  development panel v5, 28 such positions outside checks were truth
  reference (most at a reverse read's start, `527A`, `530G`, `16326G`,
  `16328G`, under a dye blob over the reference's peak), 17 a truth
  heteroplasmy, and 6 the truth's base; every truth-base share was below
  0.25 and every truth-reference share but two above 0.33.
- `estimated_length` `{runs: [{base, length, length_evidence}]}`: the site's
  runs whose length is `phase_loss` or `reference_frame`, so their 3' end is
  not callable. A run whose length a profile report convention leaves
  unrecorded is not listed (rules version 9): a `run_length` convention
  labelled `not_recorded` whose window the site overlaps, of the run's base,
  with the run at least its `minimum_run` long without the bases the site's
  call adds to the reference. An estimate lengthens a run only by adding
  bases, so a run that reaches the bound only by them may be one whose
  length the reviewers record, and stays listed. The reviewers do not record
  that length, so its estimate leaves nothing to assure (on the human
  profile, the HVS-I C run of ten or more that 16189C creates).
- `length_unconfirmed` `{runs: [{run, base, called, agreement, reads}]}`
  (rules version 21, SRS-API-020): runs of an anchored repeat window (on the
  human profile, HVS-II `C7 T C5`) that the reads' tract decodes reach
  ([[PROP-0012]] increment 3) but do not confirm. `run` is 1 before the
  anchor and 2 after it (1 alone when DNA calls one run); `called` is DNA's
  length, null when its call there holds an `N`. Each read gives a length
  with its `source`:
  - `decode`, the first-run marginal of the run the read enters first (a
    forward read enters run 1, a reverse read run 2);
  - `derived`, the next run, from the read's total, when its mark is the
    anchor and both its marginals are clear;
  - `in_phase`, a run the consensus reads in phase in that read (ADR-0074).

  Each also carries its `share`, `next` length, and `margin`. `agreement`
  says why the run is not confirmed, the first that applies:
  - `structure_differs`: a read sees two runs where DNA calls one, or the
    reverse;
  - `disagree`: the strands give different clear lengths;
  - `not_clear`: a strand's lengths are not all clear and equal;
  - `one_strand`: one strand gives no length;
  - `undecided_call`;
  - `differs_from_call`: both strands agree on another length than DNA's.

  `convention` (rules version 24, user decision 2026-10-20) is the lab's
  reading of a near-equal pair, a suggestion and never DNA's call:
  `{length, rule: "shorter_of_near_equal_pair"}` when both strands read the
  run as the adjacent pair `n`, `n + 1` and at least one read's top two are
  near-equal (every read with a next length holds the pair as its top two,
  and every other read gives one of the two), with `length` `n`; the lab
  ends the run at the first mixed call of the reverse read. Null otherwise,
  for a run whose structure differs from DNA's, and in `length_not_recorded`,
  where the lab records no length.

  A run both strands give one clear length that is DNA's is confirmed and not
  listed. So is a run the consensus decided by the two-strand rule
  (`length_evidence` `two_strand`, SRS-CONS-013, rules version 23), unless
  its structure differs. The label is a check: the decodes do not assure
  DNA's length.
- `length_not_recorded` `{runs: [...]}` (rules version 21): the same runs, as
  `length_unconfirmed` lists them, confirmed ones included, of an anchored
  window over which the profile's report convention leaves the run lengths
  unrecorded (on the human profile, HVS-I). Context only: the reviewers do
  not record those lengths, and a check against them would follow a
  convention, not an error.
- `undecided`: the site is `unresolved` or `contested` (a tie, even one the
  consensus settled on the reference and whichever strands the tied reads are
  on), or its call holds an `N`; or a gap is covered but not callable.
- `not_called`: with notation calls, the consensus decided another sequence
  than the reference at the site, and DNA calls nothing within one position
  of it: the call step excluded the difference (its positions are then not
  callable, SRS-VAR-019), so the reviewer sees what the reads showed.
- `long_run` `{runs: [{base, length, reference_length, length_evidence}]}`:
  at a decided site whose call holds the same runs as the reference, a run at
  least `review.run_growth` bases longer than its reference run whose length
  rests on an anchored end (`length_evidence` `anchored_end`) or on a single
  supporting read (rules version 7). Both strands showing the length in
  phase assure it, however long. `length_evidence` is `null` when the site
  was decided position by position. On current development data such
  lengths were 89 % wrong (an anchored end or one read counts one base too
  many in a long run), against 2 % for lengths both strands showed in phase.
- `edge_deletion` `{read}`: a decided deletion that rests on a single read
  and lies within `review.deletion_zone` positions of a reportable region's
  first or last position, where a read's own edge loses or doubles a base
  (rules version 7). On current development data such deletions were 88 %
  wrong, against 0.2 % inside the regions and 0.8 % with both strands. Rules
  version 20 widened the zone with `review.edge_zone`; rules version 22 gives
  it its own key at its earlier width, since the wider zone flagged the
  common true `106-111DEL`, read by one HVS-II forward read 40 to 70 calls
  from its end (41 checks on development panel v5, DNA right in all).
- `beside_not_analysed` (rules version 20): an item with a difference lies
  within `review.edge_zone` positions of a reportable region's end that the
  call does not analyse (`not_analysed`, SRS-CONS-012). The reads there sit
  right after their poor start, where an indel is often wrong: on
  development panel v5, three such differences (an insertion at 111 in two
  samples, a shift at 102 in one) were DNA errors that a check at the cut
  edge had covered before.
- `run_deletion` `{read}`: a decided site that rests on a single read whose
  call is at least `review.run_growth` bases shorter than its reference,
  where the reference holds a run of at least `review.phase_run` equal bases
  (rules version 13). A read that loses or compresses part of a run aligns it
  as a long deletion, and one read cannot assure that. On current
  development data it made 6 sites checks: the 4 auto-passed wrong ones (a
  compressed HVS-II run, and a composition taken from one read's end) and 2
  true deletions.
- `one_read` `{read}`: a variant site with exactly one supporting read.
- `one_strand` `{strand, reads}`: a variant site whose supporting reads all
  lie on one strand. It describes the support of the call, so it can come
  with `strands_disagree` when the opposing reads lie on the other strand.
- For every admitted read covering a position of the item whose call there is
  not clean, one label per label and read:
  - `phase_shift` `{read, from}`: the call is masked but still anchors the
    alignment; `from` is the first reference position of the read's masked
    stretch that reaches the item, walking back in the read's own direction;
  - `noisy_signal` `{read, reason}`: the call is masked otherwise; `reason`
    is the modality's mask reason, verbatim;
  - `read_edge` `{read}`: the call is unmasked but neither trusted nor exempt:
    it lies inside the read-end margin of a read end, an untrimmed call, or a
    masked call.
- `poor_signal` `{read}` (rules version 18): a read's call covering a position
  of the item lies inside the read's unreliable span (`poor_start_signal`,
  SRS-CALL-024), and no read shows that position soundly outside its own
  unreliable span (a scored, trusted call withheld for nothing): the
  position is `unreliable_only` in the consensus (SRS-CONS-012). Inside
  a read's poor start, its differences are often wrong (one-position shifts
  at its first calls), so an item with a difference there is a check; an
  item without one keeps it as context.
- `dephased_end` `{runs: [{base, length, length_evidence}]}`: the site's runs
  whose length rests on an anchored end (`length_evidence` `anchored_end`):
  a read loses phase after the run and its first masked call reading another
  base marks where the run ends (ADR-0071). Context only (rules version 10):
  such lengths were 98.4 % right on current development data, so the user
  kept them; the label records that the length leans on the dephased tail,
  which is read later. `long_run` still makes a long one a check.
- `not_covered`: no admitted read covers the gap.

`phase_shift` and `noisy_signal` read the core's view of the modality's masks;
their counts may move when the modality's per-call verdicts carry their own
reasons ([[PROP-0008]]).

## Not in version 1

- Inserted calls, which have no reference position, are not listed in
  `evidence`; their item's `sequence` and `calls` carry them.
- Watched positions, which a profile declares for every sample
  ([[PROP-0004#review-list-canonical-output|PROP-0004]]),
  come in a later revision with the raw signal they need.

- The circular origin: the neighbour test of `possible_heteroplasmy` and
  the walk of `phase_shift.from` do not wrap from position 16569 to 1, and a
  site is clipped to the reportable regions from its first to its last
  position. The shipped profile's regions do not touch the origin.

## Determinism

The same traces, reference, configuration, run signal, and DNA version give a
byte-identical document, whatever the number of workers. Reads follow content
order, items the priority order above, labels their own order, and evidence
reference order. Shares and thresholds are rounded to three decimals; a
share qualifies when it exceeds its unrounded threshold.

## Schema

The schema reuses `dna.notation/v1` (reference, profile, plugins, call text,
convention, EMPOP query) and `dna.consensus/v1` (read rejection, run) by their
`$id`, so a validator must register those schemas too;
`tools/python/scripts/validate_result_schemas.py` does.
