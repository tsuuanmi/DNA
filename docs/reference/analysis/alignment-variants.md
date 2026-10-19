---
decided-by: [ADR-0076]
---

# Analysis Alignment, Aligned Calls, and Variants

Part of the canonical [analysis contract](README.md).

## Alignment summary

The selected alignment reports:

- `orientation` (`forward` or `reverse`);
- `callable_bases` and callable `identity`;
- `unresolved_bases` (unresolved calls) and `masked_bases` (masked calls of
  the read's callability: dephased ones align with their call and profile, the
  others as unresolved) and `gap_opens`;
- one or two 0-based half-open `reference_segments`;
- `callable_reference_segments`: the 0-based half-open parts of the mapped
  segments observed by unmasked calls that read a resolved base (not `N`),
  and by deletions between them, in alignment order;
- `wraps_origin`.

Gapped query/reference rows, operation runs, fixed-point profile score, exact-match/mismatch redundancy, and deletion columns remain internal; each call's column is in [`calls`](#aligned-calls). Placement consumes the retained post-trim `EvidenceProfile` sequence, while the retained primary sequence remains attached to traceback columns for callable/identity metrics and downstream primary-sequence variant extraction. `reference_segments` therefore describe post-trim mapped coverage.

## Aligned calls

`calls` holds every call of the selected alignment, one per call, in call
order (the read's own order, also for a reverse read), so each joins the
[`dna.read_evidence/v1`](../read-evidence.md) call of the same `index` and the
document shows the chain sanger → align → call for its read
([PROP-0011](../../proposals/0011-step-contracts.md)). A call the alignment
does not hold has no entry; the aligned calls may reach beyond `read.trim`.
Each holds:

- `index`: the call's 0-based index in the read's own order;
- `column`: where the alignment places the call, before normalization:
  - `position`: the 1-based reference position of a matched call; for an
    inserted call, the last reference position before it in reference order;
  - `ordinal`: 0 for a matched call; for an inserted call, its 1-based place
    among the calls inserted after `position`, in reference order;

  absent for a call inserted before the first reference position;
- `trust`: the core's trust in the call, as its variant calling and
  consensus judge it (`dna_call::read_trust`, the same trust `dna review`
  reads):
  - `placement`: as in `dna.read_evidence/v1`, the `kind` (`scored`,
    `anchoring`, `unresolved`) and, for a masked call, the modality's
    `reason`;
  - `trusted`: the modality places the call clear of its signal edges, and
    it lies outside the read-end margin of a read end it does not vouch for;
  - `exempt`: a bounded-run exemption lifts `read_end` from the call
    ([ADR-0076](../../decisions/adr/0076-bounded-run-read-end-exemption.md));
  - `withheld_by`: what withholds trust, before then after the call in read
    order (at most one on each side), empty when the call is trusted:
    `{"kind": "read_end"}` (the read's own end, which its modality does not
    vouch for), `{"kind": "untrimmed"}` (a call outside the read's
    informative interval), or `{"kind": "masked", "reason": …}` (a masked
    call that does not anchor the alignment, with the modality's reason
    verbatim). An uninformative call is withheld by itself; an untrusted
    call can have no withholder;
  - `scope`: `full` (the call's base, a run's length and ends, and the
    indels it shows) or `base` (its base only).

Pending: the read's component calls of each phase mixture
([PROP-0010](../../proposals/0010-length-mixture-decomposition.md)) join this
document in a later version, once PROP-0010's first increment is merged.

## Variant calls

Each normalized variant contains only `position`, `reference`, `alternate`,
`kind`, and direct `calls`. The variant position is 1-based on the supplied
reference strand.

Each call is optimized for scientific review rather than implementation
traceability:

- `role`: `supporting` or `flanking`;
- `base`: called base projected onto the reference strand;
- `peaks`: raw analyzed A/C/G/T channel heights sampled together at the unique
  primary-event coordinate and projected to reference orientation. A flanking
  call whose primary base is unresolved (`N`) is omitted from `calls`;
- `quality`: the existing uncalibrated relative score under a concise public
  field name.

Original call indexes, ABIF PLOC coordinates, mapped call positions, selected-peak
positions/sources, penalties, calibration flags, and vendor scores remain internal.

For reverse reads, both `base` and `peaks` are reference-oriented, so the
reviewer can compare the normalized alternate allele directly with the signal
channels.

### SNVs

An SNV contains one or more supporting calls. The top-level variant `position`
already identifies the biological coordinate, so that coordinate is not repeated
inside each call.

### Insertions

An insertion contains one supporting call per inserted base plus available
flanking calls. Supporting-call bases and peaks show the observed inserted
sequence directly.

### Deletions

A deletion has no trace signal at the deleted reference base. It therefore
contains available flanking calls with their real peak/quality evidence and never
fabricates deleted-base dna.

For repeat-associated indels, normalization can move the reported allele
representation away from the observed alignment gap. The normalized
`position/reference/alternate` remain authoritative for the biological variant.
