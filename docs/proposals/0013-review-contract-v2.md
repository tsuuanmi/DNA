---
id: PROP-0013
status: draft
created: 2026-10-20
summary: "`dna.review/v2`: labels with namespaced evidence IDs from one registry, a versioned policy that sorts them into doubt, abstain, and context, clear per region, and every reviewed position shown with its own labels and the data-detected annotations it lies in."
implements: [SRS-API-014, SRS-API-015]
---

# Proposal: Review contract v2

## Problem

`dna.review/v1` ([[reference/review]], [[PROP-0004]]) works, and its labels on
decided calls are precise. Its contract is not yet one a laboratory tool can
build on:

1. **Three things share one version.** A label's predicate, whether it makes
   an item a check, and its rank all move with `provenance.rules` (the
   `review` plugin version, 19 on main at 69ee8be). The document does not say which labels make a
   check; only code and prose do. Moving `not_called` to context means a new
   rules version though no predicate changes.
2. **v1 changed in place.** `zones` (rules 14), `reference_only` (15),
   `mixtures` (17), `unreliable` (18), the `note` kind (19), and
   `run_deletion` (13) entered the closed v1 schema
   without a schema version.
3. **"DNA did not decide" and "DNA decided but doubts it" are one kind.** Both
   are `kind: check`, but they ask the reviewer different things and measure
   differently (below).
4. **Names mix evidence and interpretation.** `one_strand` names evidence;
   `possible_heteroplasmy` and `noisy_signal` name a conclusion, against
   PROP-0004's rule that the document never states one. `undecided` carries
   two meanings (a contested site, and a covered gap with no callable
   position). `noisy_signal` is any mask reason other than a phase shift.
5. **Payloads differ per label** (`{read}`, `{reads: [{read, strand, role}]}`,
   `{strand, reads}`, `{runs}`), so every consumer writes one parser per label.
6. **The document says where a position lies only by `zones`, which flatten.**
   A zone gives each position one kind by precedence, merged over every read.
   It cannot say that the forward read is dephased there while the reverse
   read is clean, which is what decides whether a position is assured.
7. **"Zone" means two things.** The document's `zones` are bounded by the
   profile's repeat windows and the sample's own reads; `dna-eval`'s zones
   (`window`, `after`, `start`, `end`, `interior`) are fixed widths from
   command-line options.

### Terms as they stand

| Term | Owner | Bounded by | Unit of |
|---|---|---|---|
| region | target profile | fixed: HVS-I 16024–16365, HVS-II 73–340, HVS-III 438–576 | reporting; the proposed `clear` |
| site | consensus (`dna.consensus/v1` `sites[]`) | the data: a decided interval that differs from the reference, was not called outright, had disagreeing clean observations, or holds an estimated run length | an item |
| zone (document) | `dna review` | the profile's repeat windows and the reads' masks; one kind per position by precedence | display |
| zone (harness) | `dna-eval` | fixed widths `--zone-after`, `--zone-edge`, `--zone-overrun` | workload measure |

## Evidence and research

Development panel v5, HEAD 69e1471 (rules 17), harness run `results/rebase-2026-10-20/69e1471/r`
(16,413 samples). The hit rate is the share of a label's check items that lie within one
position of a block where DNA differs from the truth under `haplotype_lod`.
It is an approximation computed from `checks.tsv` and `scores.tsv`, not a
harness measure.

| Label | Items | Hit rate | Hit rate when the only check label |
|---|---:|---:|---:|
| `undecided` | 16,262 | 1.2 % | 1.2 % |
| `estimated_length` | 1,782 | 3.7 % | 3.4 % |
| `strands_disagree` | 317 | 6.6 % | 9.4 % |
| `not_covered` | 265 | 0.4 % | 0.4 % |
| `possible_heteroplasmy` | 242 | 23.1 % | 61.5 % |
| `long_run` | 154 | 79.2 % | 83.5 % |
| `not_called` | 147 | 0.7 % | 0.8 % |
| `edge_deletion` | 103 | 84.5 % | 90.5 % |
| `read_disagrees` | 39 | 12.8 % | 14.3 % |
| `run_deletion` | 8 | 12.5 % | 0 % |

- The labels where DNA did not decide hit 0.4–3.7 %: their items are gaps a
  reviewer fills, not errors. The labels on a decided call hit 6.6–84.5 %.
- `undecided` is 86 % of all check items, 92 % of them one position; 2,392 of
  the 4,976 samples that are not clear are so for `undecided` alone.
- On current data (10,695 samples) 7,824 samples pass without review; 281 of
  them hold an unreviewed difference, 235 of these heteroplasmy only.
  Regions: 89.26 % pass, 1.26 % of passed regions hold an unreviewed
  difference.
- Rules 6 dropped one-strand heteroplasmy flags at 3 % precision, before the
  second-peak cause (SRS-CALL-021) and the run signal (SRS-CALL-020) existed.
- Since 69e1471, main added `poor_signal` (rules 18) and the `explained`
  `possible_heteroplasmy` as a `note` (rules 19); the hit rates above predate
  both and are re-measured with increment 1.

## Goals

- One machine-readable registry of labels, and a policy that sorts them,
  each versioned on its own and printed into the document.
- Item kinds that say what the reviewer is asked: confirm a doubt, or fill a
  gap.
- `clear` per region, the unit a reviewer signs off.
- At every reviewed position, both what the position itself shows and the
  data-detected annotations it lies in, per read and per strand.
- Label IDs that name evidence, carry one meaning each, and fit standard
  identifier rules (VCF `FILTER`).

## Non-goals

- No change to any label predicate in increment 1, and no change to calling,
  consensus, or Sanger evidence.
- No reviewer decisions in the document; it still never states a conclusion.
- No threshold choice for demoting a label; that belongs to the
  pre-registration ([[validation/acceptance-criteria]]).
- No compatibility reader for v1 in the `dna` binary (AGENTS.md change
  invariants). The v1 → v2 name map is documentation for downstream tools.

## Proposed design

The draft contract is three files in [0013-contract](0013-contract/README.md); they are
non-normative until increment 1 moves them to `docs/reference/`:

- [`review-v2.schema.json`](0013-contract/review-v2.schema.json):
  the JSON Schema of `dna.review/v2` (`$id` `urn:signal:schema:review:v2`);
- [`review-v2.labels.json`](0013-contract/review-v2.labels.json):
  the label registry, with each label's ID, category, scopes, title,
  description, v1 name, and policy 1 class and rank;
- [`review-v2.example.json`](0013-contract/review-v2.example.json):
  the v1 synthetic example rewritten as v2; it validates against the schema.

### Three layers, three versions

| Layer | What it holds | Version | A change bumps |
|---|---|---|---|
| evidence | reads, calls, windows, annotations, per-position observations | `schema_version` `dna.review/v2` | the schema |
| labels | named predicates over the evidence and their payloads | `provenance.rules` | rules; a new or renamed ID also bumps the schema |
| policy | label → class and rank; the kind and `clear` rules | `provenance.policy.version` | policy only |

The `review` plugin owns the policy. `provenance.policy.labels` prints the
whole table, so a document describes itself.

### Naming

- **Label IDs** are `category.name` in lower-case `snake_case`, valid as VCF
  `FILTER` and `INFO` IDs. An ID has one meaning and is never reused after
  retirement. Titles and descriptions live in the registry, so a laboratory
  may translate them without touching the contract.
- **`annotations`** names the data-detected intervals. The alternatives
  collide: `ranges` is already the EMPOP interpretation ranges that this very
  document embeds (`empop.ranges`, [[reference/notation]]); `coverage` is read
  depth in `dna.sample_evidence/v10` and a label category here; `segments` is
  the callability segment and consensus segment of the core; `zones` is the
  harness's workload unit; `features` is the signal-processing feature set.
  An annotation is a typed, labelled interval, as in genome-browser annotation
  tracks; it carries no functional (VEP-style) meaning here.

### Label IDs

| v1 | v2 | Scope | Policy 1 class | Rank |
|---|---|---|---|---:|
| `strands_disagree` | `conflict.strands` | item | doubt | 0 |
| `read_disagrees` | `conflict.read` | item | doubt | 0 |
| `possible_heteroplasmy` (not `explained`) | `mixture.second_base_both_strands` | position | doubt | 1 |
| `possible_heteroplasmy` (`explained`) | `mixture.second_base_explained` | position | context | 4 |
| — | `mixture.second_base_one_strand` | position | context | 4 |
| `long_run` | `length.run_longer_unconfirmed` | item | doubt | 2 |
| `edge_deletion` | `indel.deletion_single_read_at_edge` | item | doubt | 2 |
| `run_deletion` | `indel.deletion_single_read_in_run` | item | doubt | 2 |
| `undecided` (contested or unresolved site, `N` call) | `decision.contested` | item | abstain | 2 |
| `undecided` (covered gap, nothing callable) | `coverage.not_callable` | item, annotation | abstain | 2 |
| `not_called` | `call.excluded` | item | abstain | 2 |
| `estimated_length` | `length.run_end_unread` | item | abstain | 2 |
| `not_covered` | `coverage.no_reads` | item, annotation | abstain | 4 |
| `one_read` | `support.single_read` | item | context | 3 |
| `one_strand` | `support.single_strand` | item | context | 3 |
| `dephased_end` | `length.run_end_from_dephased_tail` | item | context | 4 |
| `read_edge` | `signal.read_end_margin` | position, annotation | context | 4 |
| `phase_shift` | `signal.masked_anchor` | position, annotation | context | 4 |
| `noisy_signal` | `signal.masked` (`reason` verbatim) | position, annotation | context | 4 |
| `poor_signal` (item without a difference); reads' `unreliable` | `signal.poor_start` | position, annotation | context | 4 |
| `poor_signal` (item with a difference) | `support.difference_in_poor_start` | item | doubt | 4 |
| — | `signal.dephased` | annotation | context | 4 |
| — | `mixture.length` | annotation | context | 4 |
| — | `coverage.single_strand` | annotation | context | 4 |

Ranks are v1's, so item order does not move. **A label's class depends on its
ID only:** where v1 lets one label make a check or not by its payload or its
item (`explained`, `poor_signal`), v2 has two IDs, so the policy stays a plain
table. The title of
`mixture.second_base_both_strands` is "Possible point heteroplasmy" (the ISFG
term); the ID names what the reads show.

### Three scopes

- **position**: what one reference position shows, read by read (base,
  quality, second peak with its share, threshold, and cause, mask reason, the
  core's trust); its labels are in the position's `at`.
- **item**: what a consensus site or a gap shows as a whole; its labels are in
  the item's `labels`.
- **annotation**: a maximal run of columns where the reads show one
  characteristic of their own (dephased tail, masked signal with its reason,
  phase mixture, read-end margin, one strand only, a gap). Found in the data,
  never a fixed width, a profile window, or a whole region. Annotations may
  overlap; a position lists those it lies in under `within`.

Every label has one envelope, `{id, reads: [{read, strand, role}], detail}`.
`role` is `supporting`, `opposing`, or `subject` (the read the label is about);
per-label fields sit in `detail`, whose shape the schema fixes per ID. Per-read
shares and thresholds are not repeated in labels: they are in the position's
observations.

### Document

| Field | For status | Content |
|---|---|---|
| `provenance` | all | v1's, plus `policy {version, labels: [{id, class, rank}]}`; `rules` is the predicate version |
| `reads` | all | admitted reads `{name, sha256, strand, admitted}`, rejected reads with `rejection`; v1's `masked`, `mixtures`, and `unreliable` move to `annotations` |
| `calls`, `callable`, `reference_only`, `empop` | called | as v1 |
| `windows` | called | the profile's repeat windows `{name, first, last}`: profile knowledge, kept apart from data |
| `annotations` | called, call_failed | `{id: "A<n>", label, first, last, reads: [{read, strand, first, last}], strands: {forward, reverse}, detail}` in order of first, then label, then read |
| `items` | called, call_failed | `{id, region, first, last, kind, decision, sequence, calls, labels, positions}` in priority order |
| `summary` | all | `items`, `doubts`, `abstains`, `variants`, `contexts`, `labels` (every ID, zeros included), `positions {reportable, review, reference_only}`, `annotations`, `regions`, `clear` |

v1's `zones` tiling is dropped: it follows from `windows` and `annotations`,
and its name collides with the harness's zones, which keep it.

**`strands`** is `affected` (a read on that strand makes the annotation),
`clean` (a read on that strand covers it callable and unmasked), or `absent`
(no read on that strand covers it). It is the field a reviewer decides on: a
one-strand second base inside an annotation whose other strand is `affected`
cannot be confirmed or refuted there; with the other strand `clean`, that
strand reads cleanly and shows no second base, which is evidence against it.

**Item `id`** is `I-` and the first 12 hex digits of SHA-256 over
`sample_id`, `first`, `last`, and `rules`, so a reviewer decision held outside
DNA can point to the same item across re-runs of the same rules.

**`positions[]`** replaces v1's read-major `evidence`: one entry per reference
position of the item, each with `at`, `within`, and `reads` (every covering
read's `base`, `state`, `call`, `scope`, `scan`, `quality`, `second {base,
share, threshold, cause, run_calls}`, and `mask`).

```json
{
  "position": 16093, "ordinal": 0,
  "at": [{"id": "mixture.second_base_one_strand",
          "reads": [{"read": "HV1F", "strand": "forward", "role": "subject"}],
          "detail": {"base": "C"}}],
  "within": ["A4"],
  "reads": [
    {"read": "HV1F", "strand": "forward", "base": "T", "state": "reference", "call": 212,
     "scope": "full", "scan": 5120, "quality": 41, "mask": null,
     "unreliable": null,
     "second": {"base": "C", "share": 0.312, "threshold": 0.2, "cause": "second_template"}},
    {"read": "HV1R", "strand": "reverse", "base": "T", "state": "masked", "call": 388,
     "scope": null, "scan": 7741, "quality": 12, "mask": "weak", "unreliable": null,
     "second": null}
  ]
}
```

with `A4` = `{"label": "signal.masked", "first": 16071, "last": 16110,
"detail": {"reason": "weak"}, "strands": {"forward": "clean", "reverse":
"affected"}}` (columns shortened). The values are synthetic: point
heteroplasmy at 16093 is possible on the forward read and unconfirmable on
the reverse one, so the reviewer checks it or asks for a rerun.

### Kind and clear

- `kind` is `doubt` when a label of class doubt is present (in `labels` or any
  position's `at`); else `abstain` for class abstain; else `variant` when the
  item holds a call; else `context` (v1's `note`). A `context` item holds only context
  labels, for example a one-strand second base outside every site, so it can
  be seen and measured without making work.
- `summary.regions[]`, one per profile region in profile order: `name`,
  `first`, `last`, `clear` (status `called` and no doubt or abstain item in the
  region), `doubts`, `abstains`. Items are clipped to regions.
- `summary.clear` is true when every region is clear. In the example, HVS-I
  holds a doubt and two abstains while HVS-II and HVS-III are clear.

### Gap items

A gap item (`coverage.not_callable`) of one position is 92 % of v1's
`undecided` items. Gaps are merged before items are ordered:

- a gap that touches an item (adjacent column) joins that item, whose
  `decision` becomes `partly_callable`;
- gaps that lie inside the same annotation, with only positions of that
  annotation between them, form one item spanning them, so a reviewer opens
  one item per data-detected interval;
- a gap in no annotation stays its own item.

Merging changes item counts, never `clear`: an item merged from abstain gaps is
abstain, and the positions it holds were abstain before.

### New label

`mixture.second_base_one_strand` fires when one strand alone reads the
position cleanly, that call's second base stands out with a peak the modality
does not explain (`second_unexplained`, SRS-CALL-022) above the run-signal
threshold of its key (with `--signal`), and the neighbour test passes. It
reads that generic assertion only, never the second-peak cause label, which
may decide nothing (SRS-CALL-021). It is context in policy 1;
promoting it to doubt is a policy change made only when the pre-registered
measurement shows it helps.

### Production readiness

1. **One registry.** `review-v2.labels.json` is the source of the IDs; the
   schema's `labelId` enum and the label section of [[reference/review]] are
   generated from it, and a check fails when they differ.
2. **One envelope and closed shapes.** One parser for every label; the schema
   enumerates the IDs and fixes each `detail`.
3. **Fail closed.** A reader that meets an unknown ID, policy version, or
   class refuses the document.
4. **Stable item IDs** for reviewer decisions kept outside DNA (format out of
   scope).
5. **Golden documents.** The example validates in CI
   (`validate_result_schemas.py`), and a byte-identity test runs on recorded
   `dna.read_evidence/v1` inputs (SRS-API-019), so no trace is needed.
6. **Measured labels.** Every rules or policy bump records, per label, the
   yield, the hit rate when it is the only check label, and the per-region
   leak on a named panel ([[PROP-0005]]). A doubt label under the
   pre-registered threshold moves to context, or the record states why it
   stays; abstain labels are exempt.
7. **Size budget.** Document size p50/p99 on panel v5 is recorded with each
   increment; `positions[]` covers item positions only.

## Alternatives considered

- **Policy in the target profile.** Flexible per laboratory, but it turns
  review classes into calibration by configuration and makes leaks
  incomparable between laboratories. Rejected.
- **Keep `kind: check`, add a class per label.** Keeps v1's shape, but an item
  with labels of both classes stays ambiguous. Rejected.
- **Clear per zone.** Lowest workload (97.26 % of harness zones pass), but
  leaks gather in interior zones (278 of 423) and reviewers sign off regions.
  Kept as a harness measure.
- **Keep `zones` and add read detail to them.** A tiling gives one kind per
  position; annotations must overlap to say "forward dephased, reverse
  clean". Rejected.
- **Call them `ranges` or `coverage`.** Both names are taken (see Naming).
  Rejected.
- **Aliases for v1 names inside DNA.** A compatibility path the specification
  does not require. Rejected; the map above serves downstream tools.

## Validation plan

Each increment is measured on development panel v5 against the main revision
before it:

1. **Equivalence.** v1 and v2 documents for the same inputs hold the same
   label occurrences under the name map, and before gap merging the same
   items: `undecided` equals `decision.contested` plus
   `coverage.not_callable`, every v1 check is a v2 doubt or abstain, every v1
   `note` a v2 `context`, and every v1 `clear` sample is v2 `clear`. Gap
   merging then reports its item counts per region; `clear` must not move.
2. **Region clear.** v2 `summary.regions` agrees with the harness's passed
   regions (89.26 % on current data at 69e1471).
3. **Annotations.** Every v1 `masked` stretch and mixture appears as an
   annotation with the same bounds; `strands` agrees with the reads' callable segments.
4. **New label.** Yield and hit rate of `mixture.second_base_one_strand` with
   and without `--signal`; the auto-pass heteroplasmy leak it would catch if
   promoted.
5. **Determinism and size.** Byte-identical documents for 1 and N workers;
   size p50/p99 recorded.

## Rollout and rollback

Four pull requests, each its own plugin and schema version:

1. Registry, label IDs, policy, `kind`, `summary.regions`, `windows`; v1's
   `zones` removed. `dna.review/v1` → `v2`. No predicate changes.
2. `annotations[]` and item `positions[]` with `at` and `within`, replacing
   `evidence` and the reads' `masked`, `mixtures`, and `unreliable`.
3. `mixture.second_base_one_strand`, as context.
4. `dna-eval` reads v2 and records per-label measures.

Rollback is a revert of the pull request. DNA is pre-release, so no stored v1
documents are migrated.

## Observability and operations

The operation log records the rules and policy versions per sample, as it
records plugin versions now. No new metric or playbook.

## Security and privacy

No new input or trust boundary. Documents name reads by trace stem, as v1
does. Stable item IDs derive from the sample identifier; they add no
identifying data that is not already in the document.

## Compatibility and migration

Breaking for every consumer of `dna.review/v1` and of
`SampleAnalysis::review_checks` (SRS-API-015 label names). The name map in
this proposal is the migration note. `dna-eval`'s workload measures read the
check items and change to doubt plus abstain.

## Decision

Draft. Agreed with the user on 2026-10-20: policy owned by the review plugin
with its own version; kinds `doubt | abstain | variant | context`; `clear` per
region; the renaming above; position and annotation labels; the name
`annotations`; and merging one-position gaps (Gap items). Rebased on main at
69ee8be (rules 19): `explained` and `poor_signal` map to their own IDs, and
v2's first rules version is 20. Open:

- `call.excluded` at 315.1C, where the exclusion matches the reviewers'
  convention: a policy exception or a predicate change;
- the doubt threshold, which needs the pre-registration directory;
- watched positions ([[PROP-0004]]), and whether `mixture.length` ever
  decides anything;
- how this contract consumes the components of [[PROP-0012]] (`dna-trace`),
  which today reach the document as `mixture.length` annotations.

## Implementation

None yet.
