---
id: PROP-0012
status: draft
created: 2026-10-20
summary: "`dna.review/v2`: labels with namespaced evidence IDs from one registry, a versioned policy that sorts them into doubt, abstain, and context, clear per region, and every reviewed position shown with its own labels and the data-detected stretches it lies in."
implements: [SRS-API-014, SRS-API-015]
---

# Proposal: Review contract v2

## Problem

`dna.review/v1` ([[reference/review]], [[PROP-0004]]) works, and its labels on
decided calls are precise. Its contract is not yet one a laboratory tool can
build on:

1. **Three things share one version.** A label's predicate, whether it makes
   an item a check, and its rank all move with `provenance.rules` (the
   `review` plugin version, 17). The document does not say which labels make a
   check; only code and prose do. Moving `not_called` to context means a new
   rules version though no predicate changes.
2. **v1 changed in place.** `zones` (rules 14), `reference_only` (15),
   `mixtures` (17), and `run_deletion` (13) entered the closed v1 schema
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

Development panel v5, HEAD 69e1471, harness run `results/rebase-2026-10-20/69e1471/r`
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

## Goals

- One machine-readable registry of labels, and a policy that sorts them,
  each versioned on its own and printed into the document.
- Item kinds that say what the reviewer is asked: confirm a doubt, or fill a
  gap.
- `clear` per region, the unit a reviewer signs off.
- At every reviewed position, both what the position itself shows and the
  data-detected stretches it lies in, per read and per strand.
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

### Three layers, three versions

| Layer | What it holds | Version | A change bumps |
|---|---|---|---|
| evidence | reads, masks, mixtures, per-position read observations, stretches | `schema_version` `dna.review/v2` | the schema |
| labels | named predicates over the evidence, their payloads | `provenance.rules` | rules (a new or renamed label also bumps the schema) |
| policy | label → class and rank; the `clear` rule | `provenance.policy.version` | policy only |

The `review` plugin owns the policy. The document prints the whole table in
`provenance.policy`, so it describes itself.

### Label IDs

IDs are `category.name`, lower-case `snake_case`, with no space, semicolon, or
other character the VCF `FILTER` and `INFO` key rules forbid. An ID has one
meaning and is never reused after retirement. Human text (`title`,
`description`) lives in the registry, not in the ID, so a laboratory may
translate it.

| v1 | v2 | Scope | Policy v1 class | Rank |
|---|---|---|---|---:|
| `strands_disagree` | `conflict.strands` | item | doubt | 0 |
| `read_disagrees` | `conflict.read` | item | doubt | 0 |
| `possible_heteroplasmy` | `mixture.second_base_both_strands` | position | doubt | 1 |
| — | `mixture.second_base_one_strand` | position | context | 4 |
| `long_run` | `length.run_longer_unconfirmed` | item | doubt | 2 |
| `edge_deletion` | `indel.deletion_single_read_at_edge` | item | doubt | 2 |
| `run_deletion` | `indel.deletion_single_read_in_run` | item | doubt | 2 |
| `undecided` (contested or unresolved site, `N` call) | `decision.contested` | item | abstain | 2 |
| `undecided` (covered gap, nothing callable) | `coverage.not_callable` | stretch | abstain | 2 |
| `not_called` | `call.excluded` | item | abstain | 2 |
| `estimated_length` | `length.run_end_unread` | item | abstain | 2 |
| `not_covered` | `coverage.no_reads` | stretch | abstain | 4 |
| `one_read` | `support.single_read` | item | context | 3 |
| `one_strand` | `support.single_strand` | item | context | 3 |
| `dephased_end` | `length.run_end_from_dephased_tail` | item | context | 4 |
| `read_edge` | `signal.read_end_margin` | position, stretch | context | 4 |
| `phase_shift` | `signal.masked_anchor` | position, stretch | context | 4 |
| `noisy_signal` | `signal.masked` (`reason` verbatim) | position, stretch | context | 4 |
| — | `signal.dephased` | stretch | context | 4 |
| — | `mixture.length` | stretch | context | 4 |
| — | `coverage.single_strand` | stretch | context | 4 |

Ranks are v1's, so item order does not move. The title of
`mixture.second_base_both_strands` is "Possible point heteroplasmy" (the ISFG
term); the ID names what the reads show.

### Three scopes

A label belongs to one scope or more:

- **position**: what one reference position shows, read by read: its base,
  quality, second peak and its cause, and the core's trust. Example:
  `mixture.second_base_both_strands`.
- **item**: what a consensus site or a gap shows as a whole. Example:
  `conflict.strands`, `length.run_longer_unconfirmed`.
- **stretch**: a maximal run of positions where the sample's reads show a
  characteristic of their own. A stretch is found in the data and is never a
  fixed width or a whole region: a read's dephased tail after a repeat run,
  its masked stretch with the mask reason, its phase mixture, its read-end
  margin, positions only one strand reads callably, a covered gap with
  nothing callable, an uncovered gap.

### Stretches

`stretches[]` is a sample-level list in reference order. Stretches may overlap.
Each has:

- `id`: `S1`, `S2`, … in order of `first`, then `label`; deterministic for
  the same inputs;
- `label`, `first`, `last` (columns as in items), `detail` (for example
  `{state: "dephased", after_repeat: true}` or the mask `reason`);
- `reads`: every read whose own signal makes the stretch, with `read`,
  `strand`, and its own `first`..`last`;
- `strands`: for `forward` and `reverse`, `affected` (a read on that strand
  makes the stretch), `clean` (a read on that strand covers the stretch and is
  callable and unmasked there), or `absent` (no read on that strand covers
  it).

`strands` is the field a reviewer decides on. A one-strand second base inside
a stretch where the other strand is `affected` means the other strand cannot
confirm or refute it. The same second base where the other strand is `clean`
means the other strand reads cleanly and shows no second base, which is
evidence against it.

The profile's repeat windows stay as `windows[]` (`name`, `first`, `last`). They
are profile knowledge, not data, and are kept apart from stretches. v1's `zones`
tiling is dropped: it can be derived from `windows` and `stretches`, and its name
collides with the harness's zones, which keep the name.

### Positions inside items

Each item's `evidence` (read-major in v1) becomes `positions[]`
(position-major):

```json
{
  "position": 16093, "ordinal": 0,
  "at": ["mixture.second_base_one_strand"],
  "within": ["S4"],
  "reads": [
    {"read": "HV1F", "strand": "forward", "state": "reference", "base": "T",
     "quality": 41, "scope": "full", "scan": 5120,
     "second": {"base": "C", "share": 0.312, "threshold": 0.200, "cause": "second_template"}},
    {"read": "HV1R", "strand": "reverse", "state": "masked", "base": "T",
     "quality": 12, "scope": null, "scan": 7741, "mask": "weak", "second": null}
  ]
}
```

and the sample's stretch:

```json
{"id": "S4", "label": "signal.masked", "first": {"position": 16071, "ordinal": 0},
 "last": {"position": 16110, "ordinal": 0}, "detail": {"reason": "weak"},
 "reads": [{"read": "HV1R", "strand": "reverse", "first": 16071, "last": 16110}],
 "strands": {"forward": "clean", "reverse": "affected"}}
```

The values are synthetic. Read together: at 16093 the forward read shows a
second base C at 0.312, which its cause labels `second_template`, and the
reverse read is masked as weak over 16071–16110. Point heteroplasmy is
possible and unconfirmable on the reverse strand, so the reviewer checks it
or asks for a rerun. `at` lists the position-scope labels; `within` lists the
stretches the position lies in.

### Item kind and clear

- `kind` is `doubt` when an item has a doubt label, else `abstain` when it has
  an abstain label, else `variant`. A doubt includes the measurable errors,
  so it takes precedence; the reviewer still fills an abstain inside it.
- `summary.regions[]`, one per profile region in profile order: `name`,
  `clear` (the status is `called` and no doubt or abstain item touches the
  region), `doubts`, `abstains`. Items are clipped to regions, so each item
  counts for one region.
- `summary.clear` is true when every region is clear.

### New label

`mixture.second_base_one_strand` `{position, base, read, strand, share,
threshold}` fires when one strand alone reads the position cleanly, that
call's second peak has cause `second_template` (not `context`), its share
exceeds the threshold the run signal gives its key (with `--signal`), and the
neighbour test of `mixture.second_base_both_strands` passes. It is context in
policy 1, so no `clear` changes. Promoting it to doubt is a policy change,
made only when the pre-registered measurement shows it helps.

### Production readiness

1. **One registry, generated docs.** `$defs/labelRegistry` in
   `review-v2.schema.json` lists every ID with `category`, `scope`, `title`,
   `description`, `payload`, and `since_rules`. The label section of
   [[reference/review]] is generated from it, and `validate_docs_structure.py`
   fails when they differ. A stale heading such as "Labels (rules version
   13)" cannot recur.
2. **One label envelope.** Every label is `{id, positions, reads: [{read,
   strand, role}], detail}`; per-label fields sit in `detail`, and thresholds
   used are copied there. Consumers need one parser.
3. **Fail closed.** The schema enumerates the IDs. A reader that meets an
   unknown ID, policy version, or class refuses the document instead of
   guessing.
4. **Stable IDs for items and stretches.** An item ID derives from the sample,
   its first and last columns, and the rules version, so a reviewer decision
   held outside DNA can point to it across re-runs of the same rules. The
   decision format itself is out of scope here.
5. **Golden documents.** Synthetic examples for each label and scope, checked
   by `validate_result_schemas.py`, and a byte-identity test on recorded
   `dna.read_evidence/v1` inputs (SRS-API-019) so no trace is needed in CI.
6. **Measured labels.** Every rules or policy bump records, per label, the
   yield, the hit rate when it is the only check label, and the per-region
   leak on a named panel (`dna-eval`, [[PROP-0005]]). A doubt label under the
   pre-registered threshold moves to context, or the record states why it
   stays. Abstain labels are exempt: they mark gaps.
7. **Size budget.** Document size p50/p99 on panel v5 is recorded with each
   increment; `positions[]` covers item positions only, never the whole
   regions.

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
  position; stretches must overlap to say "forward dephased, reverse clean".
  Rejected.
- **Aliases for v1 names inside DNA.** A compatibility path the specification
  does not require. Rejected; the map above serves downstream tools.

## Validation plan

Each increment is measured on development panel v5 against the main revision
before it:

1. **Equivalence.** v1 and v2 documents for the same inputs hold the same
   items, positions, and label occurrences under the name map. `undecided`
   equals `decision.contested` plus `coverage.not_callable`. Every v1 check is
   a v2 doubt or abstain, and every v1 `clear` sample is v2 `clear`.
2. **Region clear.** v2 `summary.regions` agrees with the harness's passed
   regions (89.26 % on current data at 69e1471).
3. **Stretches.** Every v1 `masked` stretch and mixture appears as a stretch
   with the same bounds; `strands` agrees with the reads' callable segments.
4. **New label.** Yield and hit rate of `mixture.second_base_one_strand` with
   and without `--signal`; the auto-pass heteroplasmy leak it would catch if
   promoted.
5. **Determinism and size.** Byte-identical documents for 1 and N workers;
   size p50/p99 recorded.

## Rollout and rollback

Four pull requests, each its own plugin and schema version:

1. Registry, label IDs, policy, `kind`, `summary.regions`, `windows`; v1's
   `zones` and `evidence` removed. `dna.review/v1` → `v2`. No predicate
   changes.
2. `stretches[]` and item `positions[]` with `at` and `within`.
3. `mixture.second_base_one_strand` (rules 18), as context.
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
with its own version, kinds `variant | doubt | abstain`, `clear` per region,
the renaming above, and position plus stretch labels. Open:

- merging one-position `coverage.not_callable` stretches into their
  neighbouring items (it changes item counts; decide before increment 1);
- `call.excluded` at 315.1C, where the exclusion matches the reviewers'
  convention: a policy exception or a predicate change;
- the doubt threshold, which needs the pre-registration directory;
- watched positions ([[PROP-0004]]), and whether `mixture.length` ever
  decides anything.

## Implementation

None yet.
