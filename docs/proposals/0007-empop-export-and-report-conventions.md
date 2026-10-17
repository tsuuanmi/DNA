---
id: PROP-0007
type: proposal
status: draft
owners: []
created: 2026-10-17
related-requirements: [SRS-NOM-016]
related-decisions: [ADR-0069, ADR-0077]
implementation: []
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
  the haplotype. On development panel v5, 47 sequence-class samples differ
  from the truth only by such a length call, and the evaluation harness cannot
  yet tell them apart from DNA defects.

## Goals

1. An offline EMPOP query in the `notation` output: interpretation ranges and
   an rCRS-coded mitotype, built from the sample's calls.
2. Every position inside an exported range was decided by DNA. No range covers
   an `N`, a contested site, or an excluded call.
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
  - 47 sequence-class samples differ only by an HVS-I C-tract length call.

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
  - the profile's regions intersected with the union of the records'
    callable segments;
  - less every position of an excluded variant;
  - less both positions of a junction with an excluded insertion;
  - merged where they touch.

  They are written in the profile's region order and, within a region, in
  position order (`16024-16365 73-340 438-576`).
- **Mitotype.** The mitotype is built as follows:
  - every eligible represented call inside the ranges, in the profile notation
    (`263G`, `249DEL`, `309.1C`);
  - distinct, sorted by position and then by insertion ordinal;
  - never a reference-equal token.

  An insertion is written only when both of its junction positions are inside
  one range.
- **IUPAC.** A code passes through unchanged when DNA writes one. None is
  written today.
- **Document.** The plugin adds an `empop` object to `notation/v1`:
  - `ranges`: 1-based inclusive pairs, plus the rendered `ranges_text`;
  - `mitotype`: the space-separated tokens;
  - `cut`: the excluded positions with their reason (`unresolved`,
    `contested`, `excluded_variant`), so a reviewer sees why a range is split.

  The contract change follows the
  [versioning policy](../governance/versioning.md#public-json-contracts):
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
# A length change of the C run in this window, after window representation.
window = [16180, 16193]
kind = "run_length"
```

The nomenclature plugin labels each represented call that a convention covers
with that convention (`convention: "hvs1_c_length"`). The labelled call stays
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
  - an insertion at a range edge is left out;
  - tokens are distinct and sorted;
  - an IUPAC code passes through;
  - regions keep the profile's order.
- **Unit tests (`dna-normalize`):**
  - a convention labels exactly the covered calls;
  - an unlabelled call is unchanged;
  - the haplotype rebuilt from the labelled notation is identical.
- **Invariant, development panel v5:**
  - every position inside an exported range is decided by the consensus
    (no `N`, no contested site, no excluded call);
  - the sequence rebuilt from rCRS and the mitotype over the ranges equals
    the consensus sequence there.
- **Equivalence, Python check (ADR-0066):** given the same calls and
  intervals, DNA's export is byte-equal to `mtdna_raw`'s `EmpopProfile`
  (`ranges_text`, `mitotype`). Run it on every development-panel sample.
- **Truth comparison:** for concordant samples, report:
  - how often DNA's mitotype equals the truth's;
  - how often DNA's ranges lie inside the truth's intervals;
  - the differences, by class.
- **Harness:** the 47 HVS-I length samples move to `lab`
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

The rest is under review.

## Implementation

None yet. The range rule "less every position of an excluded variant" is met
by the core's callable segments of a called record (SRS-VAR-019): the `empop`
plugin takes the records' callable segments as they are, without a second
cut.
