---
implements: [SRS-CONS, SRS-NOM, SRS-VAR]
decided-by: [ADR-0069]
---

# DNA Notation JSON

`dna notation <sample-id> <sample-id>.variants.json --reference <reference.fasta>`
runs the post-calling plugins over one variants document
([[reference/variants|`dna.variants/v1`]],
[[ADR-0069]])
and writes one deterministic `results/<sample-id>.notation.json` document
identified as `dna.notation/v1`. The authoritative schema is
[`schemas/notation-v1.schema.json`](schemas/notation-v1.schema.json), and the
synthetic example is
[`examples/notation-v1.example.json`](examples/notation-v1.example.json).

The target profile named by the configuration must declare notation. The
document's sample must be the requested sample, and its reference must be the
supplied reference. Only eligible variants are represented.

## Fields

- `provenance`:
  - the reference `name`, `topology`, and sequence `sha256`;
  - `configuration_sha256`;
  - the target `profile`;
  - `plugins`: `normalization` and `nomenclature`, then `conformance` and
    `empop` when the profile declares them;
  - `source`: the `schema_version` and the SHA-256 of the variants document's
    bytes.
- `notation`: the per-read view of [[reference/sample-evidence/notation|sample notation]],
  with reads named as in the source document. For the same profile it equals
  the `notation` of the source document.
- `conformance` (only when the profile declares it):
  - `rules`: the profile's notation conventions, in order;
  - `findings`: every represented variant of a read that does not follow a
    rule, as its rendered `calls`, the `read`, and the `rule`.

  Findings are reported, never applied. A finding names a convention that the
  representation does not follow; it is not an error in the haplotype.

## EMPOP query

When the profile declares `[empop]`, `empop` is an offline query for EMPOP
(SRS-NOM-018, PROP-0007). DNA sends it nowhere. EMPOP reads the omission of
sequence inside a declared range as the reference, so a range covers only
positions that DNA decided:

- `ranges`: 1-based inclusive interpretation ranges, in the profile's region
  order and, within a region, in position order. A range covers the positions
  exactly one record calls: its callable reference segments, which already
  leave out every `N` and the positions of every excluded variant
  (SRS-VAR-019). Ranges that touch merge.
- `ranges_text`: the same, as EMPOP takes it (`16024-16365 73-340 438-576`).
- `mitotype`: the represented calls inside the ranges, space-separated,
  distinct, and ordered by position and then insertion ordinal. EMPOP numbers
  an insertion `X.k` at `X`, so an insertion is written when its anchor `X`
  lies in a range. An excluded insertion cuts `X` itself, so a range that
  includes `X` lists every insertion after `X`. A labelled call
  ([[reference/sample-evidence/notation|report convention]]) is written too: the
  query must not misstate the sequence, and EMPOP's standard query disregards
  hotspot length differences anyway.
- `ignored_in_search` (empop version 3): the mitotype's insertions and
  deletions of a profile length tract's base inside the tract's `search`
  window ([[reference/profiles|profiles]]), which EMPOP's standard query
  disregards (16189-16194, 302-316, 567-574 on the human profile): written,
  but ignored in a search. The field documents EMPOP's search behaviour, not
  DNA's confidence in a call: `315.1C` is listed because the 302-316 window
  holds it, and it is still called and checked normally.
- `cut`: positions inside a profile region and a record that no range
  covers, as `start`, `end`, and `reason`:
  - `unresolved`: the record does not call it (`N`);
  - `excluded_variant`: a position of a variant excluded from reporting;
  - `not_analysed`: a reportable region's end the call does not analyse,
    observed only inside reads' unreliable spans (SRS-CONS-012, empop
    version 2);
  - `overlap`: two records call it.

A tie the consensus resolved is DNA's decision and is exported as such. The
query covers only what the variants document states: when `call` places a
consensus `N` at another reference position than the consensus site, the
range follows the document ([[validation/discrepancy-catalogue|catalogue]]).

## Rules

- `insertion_at_run_end`: an insertion of one repeated base that touches a
  reference homopolymer of that base, at least `minimum_run_length` long, must
  follow the homopolymer's 3' base (`315.1C`, not `311.1C`).
- `insertion_matches_run`: bases inserted between two bases of such a
  homopolymer must be its base. Another difference there is written as a
  substitution plus a length change (`314G 315.1C`, not `313.1G`).
