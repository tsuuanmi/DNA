---
id: PROP-0007
status: implemented
created: 2026-10-17
summary: "an offline EMPOP query (ranges and mitotype) in `notation/v1` whose ranges cover only positions DNA decided, and profile conventions that label calls the reviewers do not record instead of deleting them."
implements: [SRS-NOM-016, SRS-NOM-017]
decided-by: [ADR-0069, ADR-0077]
implementation: [PR-96, PR-99]
---

# Proposal: EMPOP export and labelled report conventions

## Problem

The laboratory searches a sample's haplotype in EMPOP to estimate its
frequency. DNA writes the calls in the lab's notation but no EMPOP query, so
the query is typed or built outside DNA from the calls and a coverage range.
Two things can go wrong there:

- **An unread position inside a declared range.** EMPOP reads "the omission of
  sequence information within a specified range" as rCRS. A position DNA left
  `N`, or a call DNA excluded, inside a declared range therefore becomes a
  reference base in the search. The search can then match haplotypes the
  sample does not share.
- **A convention that changes content.** The reviewers do not record the
  HVS-I C-tract length (`16193.1C`). DNA does record it, as its evidence shows.
  ADR-0077 §5 rules out handling this in normalization, which must preserve
  the haplotype. On development panel v5, 31 sequence-class samples differ
  from the truth only by such a length call, and the evaluation harness cannot
  yet tell them apart from DNA defects.

## Goals

1. An offline EMPOP query in the `notation` output: interpretation ranges and
   an rCRS-coded mitotype, built from the sample's calls.
2. Every position inside an exported range was decided by DNA. No range covers
   a position DNA left undecided (`N`) or excluded. A tie the consensus
   resolved is exported as DNA's decision, as in the haplotype.
3. A labelled report convention for calls that the reviewers' convention does
   not record. Such a call is classified and kept, never deleted, and the
   haplotype stays complete.
4. Clear ownership per ADR-0077 §5:
   - the core decides what is callable;
   - `dna-normalize` decides how a call is written, and its convention label;
   - `dna-empop` only formats.

## Non-goals

- Querying EMPOP, or sending sample data anywhere. The export is a string the
  user may submit.
- Frequency estimates, haplogroup assignment, or a local EMPOP index.
- Rewriting calls into EMPOP's phylogenetic alignment (ADR-0077: an EMPOP form
  is never guessed from the lab form).
- IUPAC heteroplasmy calls. They are deferred. When DNA reports them, the
  export passes the code through.
- A FASTA-like export (see *Alternatives*).

## Evidence and research

- **EMPOP documentation** ([Using EMPOP for mitotype searches](https://empop-use.readthedocs.io/en/latest/using_empop/)):
  - queries need interpretation ranges;
  - omitted sequence inside a range is read as rCRS;
  - rCRS-coded mitotypes use the extended IUPAC code, and in the default
    pattern match a mixture matches its components;
  - standard query settings disregard length differences in the hotspot
    regions, among them the C-stretches at 16180–16193 and 310–316.
- **The laboratory's current exporter** (`mtdna_raw`, `src/modules/empop/profile.py`):
  - ranges come from the sample's per-region intervals, in the order HVS-I,
    HVS-II, HVS-III (`16024-16365 73-340`);
  - tokens are `<pos><base>`, `<pos>DEL`, `<pos>.<k><base>`, distinct and
    sorted by position.
- **Development panel v5:**
  - consensus `N` positions inside truth intervals: 43,182 in 7,423 samples
    (main `711c43e`), before the planned support-veto scope change;
  - 47 sequence-class samples differ only inside HVS-I 16180–16194. In 31
    of them, the only difference is an HVS-I C-tract length call (main
    `4e6f8d1`).

## Proposed design

### 1. EMPOP export (`dna-empop`, format only)

The `notation` workflow runs a new post-calling plugin, `empop`, after
`conformance`, when the target profile declares an `[empop]` section. The
shipped human-mtDNA profile declares it. The plugin reads:

- the represented calls of the nomenclature plugin;
- from the variants document:
  - each record's `callable_reference_segments`, which already exclude
    unresolved and masked positions;
  - the reference positions of every excluded variant
    (`exclusion_reasons`, for example `placement_conflict` or
    `unresolved_flank`).

It decides nothing about the sequence. The rules:

- **Ranges.** The ranges are:
  - the profile's regions intersected with the positions exactly one record
    calls, which are its callable segments;
  - less every position of an excluded variant, and both positions of a
    junction with an excluded insertion. The callable segments already leave
    these out (SRS-VAR-019);
  - merged where they touch.

  They are written in the profile's region order and, within a region, in
  position order (`16024-16365 73-340 438-576`).
- **Mitotype.** The mitotype is built as follows:
  - every eligible represented call inside the ranges, in the profile notation
    (`263G`, `249DEL`, `309.1C`);
  - distinct, sorted by position and then by insertion ordinal;
  - never a reference-equal token.

  An eligible insertion is written when its anchor position lies inside a
  range: EMPOP numbers `X.k` at `X`. A range that includes `X` therefore
  claims the complete list of insertions after `X`. Leaving `309.1C` out of
  `73-309` would state "no insertion after 309", the very misstatement the
  export exists to prevent. The list is complete: an excluded insertion cuts
  `X` itself (SRS-VAR-019), so where `X` is in a range, every insertion DNA
  observed after `X` is eligible and written. An insertion whose anchor is
  cut is not written, even when `X+1` is in a range.
- **IUPAC.** A code passes through unchanged when DNA writes one. None is
  written today.
- **Document.** The plugin adds an `empop` object to `notation/v1`:
  - `ranges`: 1-based inclusive pairs, plus the rendered `ranges_text`;
  - `mitotype`: the space-separated tokens;
  - `cut`: the excluded positions with their reason (`unresolved`,
    `excluded_variant`, `overlap`), so a reviewer sees why a range is split.

  The contract change follows the
  [[governance/versioning#public-json-contracts|versioning policy]]:
  `notation/v1` has not been released, so it is revised in place and
  recorded as **Breaking**.

`dna-empop` keeps no network client and no corpus index (ADR-0077). The export
can be checked against the laboratory's exporter.

### 2. Labelled report conventions (`dna-normalize`)

The profile declares conventions that do not record a call. For example:

```toml
[[nomenclature.conventions]]
name = "hvs1_c_length"
label = "not_recorded"
kind = "run_length"
# A length change of an uninterrupted C run of at least 10 in this window.
window = [16180, 16193]
base = "C"
minimum_run = 10
```

The nomenclature plugin labels each represented call that a convention covers
with that convention (`"convention": {"name": "hvs1_c_length", "label":
"not_recorded"}`). The labelled call stays
in `notation`, the haplotype is unchanged, and no output drops a call. A rule
matches only on the represented call and the reference, never on read
evidence, and asserts no base: SRS-NOM rules apply.

Downstream:

- **The evaluation harness** gives a sequence-class difference whose DNA-only
  calls all carry `not_recorded` the cause `lab`, evidence
  `convention_not_recorded`, scope `deferred`. This is the computed rule of the
  cause follow-up.
- **The EMPOP export** keeps labelled calls in the mitotype. EMPOP's standard
  query disregards hotspot length differences, and leaving the call out would
  misstate the sequence inside the range.
- **A delivered report** can show or hide the label. That choice belongs to
  reporting, not to DNA's calls.

## Alternatives considered

- **A separate `dna empop` command and `empop/v1` contract.** Rejected (user
  decision, 2026-10-17). The `notation` output already composes normalization
  and conformance, and another contract would duplicate its provenance.
- **Fixed profile regions as ranges.** Rejected. Every `N` inside them would
  become a reference base in the search.
- **Deleting the HVS-I length calls.** Rejected. It changes the haplotype
  (ADR-0077 §5) and hides evidence from the reviewers.
- **A FASTA-like export per range.** Not chosen. EMPOP accepts it, and it
  avoids notation entirely, but the laboratory's exporter and reviewers use
  rCRS-coded mitotypes. It can be added later from the consensus segments if
  the rCRS-coded form proves fragile.
- **Rewriting to EMPOP's phylogenetic alignment.** Rejected by ADR-0077.

## Validation plan

- **Unit tests (`dna-empop`):**
  - an `N` splits a range;
  - an excluded variant and an excluded insertion cut their positions;
  - an insertion whose anchor ends a range is written, and one whose anchor
    is cut is not;
  - tokens are distinct and sorted;
  - an IUPAC code passes through;
  - regions keep the profile's order.
- **Unit tests (`dna-normalize`):**
  - a convention labels exactly the covered calls;
  - an unlabelled call is unchanged;
  - the haplotype rebuilt from the labelled notation is identical.
- **Invariant, development panel v5:**
  - every position inside an exported range is decided by the consensus
    (no `N`, no excluded call); the positions on a tie the consensus resolved
    as the reference are counted;
  - the mitotype holds the represented calls inside the ranges, which
    reconstruct the record's haplotype (SRS-NOM-003); the export only
    filters them.
- **Equivalence, Python check (ADR-0066):** given the same calls and
  intervals, DNA's export is byte-equal to `mtdna_raw`'s `EmpopProfile`
  (`ranges_text`, `mitotype`). Run it on every development-panel sample.
- **Truth comparison:** for concordant samples, report:
  - how often DNA's mitotype equals the truth's;
  - how often DNA's ranges lie inside the truth's intervals;
  - the differences, by class.
- **Harness:** the 31 HVS-I length samples move to `lab`
  `convention_not_recorded`, deferred. No haplotype measure changes.

## Rollout and rollback

Both parts are post-calling. The core and its outputs are unchanged.

- Rollout: new plugin `empop` at version 1, and a new nomenclature plugin
  version for the labels.
- Rollback: remove the profile's `[empop]` section or its conventions. No
  stored data depends on either.

## Observability and operations

`notation/v1` carries `cut`, so a split range is explained in the document.
No new logging is needed.

## Security and privacy

The export is a local string inside the sample's own output. DNA makes no
network call, and nothing leaves the machine unless a user submits it. The
data policy is unchanged.

## Compatibility and migration

`notation/v1` is revised in place under the versioning policy, recorded as
**Breaking**. No migration: earlier documents are regenerated.

## Decision

2026-10-17, user:
- the EMPOP export is a section of `notation/v1`;
- the HVS-I C-tract length convention is a labelled policy in `dna-normalize`
  that never deletes a call.

2026-10-17, user, on the Researcher's review of the implementation plan:
- a `run_length` convention covers a call only inside an uninterrupted run of
  at least `minimum_run` bases (10 for HVS-I). With `T16189`, the reads
  resolve the two short runs, and a recorded length there is real;
- a range never covers a position DNA left undecided or excluded. A tie the
  consensus resolved is exported as DNA's decision.

2026-10-17, on the Researcher's review of the increment 2 measurement: an
insertion is written by its anchor, as EMPOP numbers it, instead of only when
both junction positions lie in one range.

## Implementation

Planned in two increments, each its own pull request.

**Increment 1, labelled conventions** (`nomenclature` 3):

- Profile keys:

  ```toml
  [[nomenclature.conventions]]
  name = "hvs1_c_length"
  label = "not_recorded"
  kind = "run_length"
  window = [16180, 16193]
  base = "C"
  minimum_run = 10
  ```

  - `label` takes one value, `not_recorded`.
  - `kind` takes one value, `run_length`. It covers an inserted `base` whose
    anchor lies in `window`, and a deleted `base` in `window`, both judged
    per rendered call. The read's represented haplotype must have an
    uninterrupted run of `base` of at least `minimum_run` at the call,
    counting a deleted base in the run it shortens. Substitutions are never
    covered, so `16194C` and `16189C` stay content.
  - Conventions may not overlap.
- Labels are applied per rendered call, after representation, by
  `dna_normalize::variant_nomenclature`. Rendering itself stays in the
  facade.
- Every notation view carries the label as `convention` on a labelled call:
  `variants/v1`, `sample_evidence/v10` and `notation/v1`, each revised in
  place and recorded as **Breaking**. `SampleAnalysis::call` carries it too.
- Harness: a `sequence`-class sample whose DNA-only calls all carry
  `not_recorded`, and which has no truth-only call, takes the cause `lab`,
  scope `deferred`, evidence `convention_not_recorded`.
- Measured on development panel v5 (main `4e6f8d1`): 31 of the 47 HVS-I-only
  sequence samples differ only by a C-length call (29 by `16193.kC`, 2 by
  `16193DEL`). The others differ in content (for example `16194C`, or a call
  where the truth has `16189C`) and keep their cause.
  - With `minimum_run = 10`, 29 of the 31 move to `lab`, deferred. The two
    `16193DEL` samples have `16185T`, which splits the run (C8), so they
    stay unlabelled.
  - 2 more samples move: their only other difference is a truth IUPAC code
    that the class's IUPAC-compatible measure admits.
  - 31 samples move in all, every one from `undetermined`. Statuses, classes
    and every other `scores.tsv` field are unchanged.
  - 35 samples with a DNA-only HVS-I C-length call keep another cause: 22
    `sequence` and 13 `representation`. Each also differs otherwise, or its
    call sits in a run shorter than 10.

**Increment 2, the EMPOP export** (#99; `empop` 1, `notation/v1` **Breaking**):

- `[empop]` is an empty table: declaring it turns the export on, and it
  requires `[notation]`. It has no keys. The region order is the profile's,
  the ranges come from the evidence, and no tunable is needed.
- The range rule "less every position of an excluded variant" is already met
  by the core's callable segments of a called record (SRS-VAR-019). The
  plugin takes the records' callable segments as they are, without a second
  cut.
- The `cut` reasons come from the variants document only:
  - `excluded_variant`: a record position of an excluded variant;
  - `unresolved`: any other record position that no record calls (`N`);
  - `overlap`: a position that two records call, which never happens in the
    consensus path.
- A tie is handled as the consensus handled it. A tie with no decision is
  already `N`, so it is cut as `unresolved`. A one-against-one tie that the
  consensus resolved as the reference is DNA's decision, in the haplotype
  and in the export alike.
- `dna-empop` takes rendered calls as plain position, ordinal and text. It
  does not depend on `dna-normalize`. `SampleAnalysis::call` gives the same
  query, and the harness records it (`empop_*` in `scores.tsv`).
- Measured on development panel v5 (main `b075635`), 16,401 queries:
  - Equivalence: given DNA's ranges and calls, `mtdna_raw`'s `EmpopProfile`
    renders byte-equal `ranges_text` and `mitotype` for every query.
  - Resolved ties inside a range: 1,129 positions in 263 samples, exported
    as DNA's decision.
  - Undecided inside a range: 16 positions in 13 samples. In each, `dna call`
    re-aligns the consensus record and places its `N` at another reference
    position than the consensus site, so the variants document asserts a
    position the consensus left undecided:
    - 11 sites at a run's end (`303-310` as `CCCCCCCCCN`, `568-574`,
      HVS-I);
    - 4 inside an indel or repeat context (`103`, `516`/`523`, `16026`).
    The export follows the document. This is the consensus-to-`call`
    interface loss of PROP-0006, now visible through `empop_undecided`.
  - Truth comparison: the reviewers' `EmpopProfile` from the truth record,
    inside the positions both declare, for 16,400 exportable truth records.
    - DNA's ranges lie inside the truth's intervals in 15,947 samples.
    - The mitotypes are equal in 14,305, and in 15,162 once hotspot length
      differences (which EMPOP's standard query disregards) are set aside.
    - Most of the rest are truth IUPAC codes (`heteroplasmy`, 992) and
      `sequence` samples. In 'none' samples, the remaining differences are
      `309.kC` and `573.kC` lengths that the lab's comparator tolerates.
