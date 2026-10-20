---
implements: [SRS-ALN, SRS-CALL, SRS-VAR]
decided-by:
  - ADR-0027
  - ADR-0047
  - ADR-0057
  - ADR-0060
  - ADR-0062
  - ADR-0067
  - ADR-0069
  - ADR-0071
  - ADR-0076
covers:
  - crates/dna-call/src/
  - crates/dna-call/src/model/variant.rs
  - crates/dna-call/src/read_call.rs
  - crates/dna-kernel/src/read_evidence.rs
  - crates/dna-sanger/src/read_processing/evidence.rs
  - src/input/sequence.rs
---

# Variant Calling Method

Part of the canonical [[design/pipeline|DNA pipeline]].

Extracts evidence-backed primary-sequence differences from the selected
alignment and constructs validated anchored REF/ALT alleles. Only differences
in the primary sequence are considered; no allele-frequency, genotype,
heteroplasmy, post-calling haplotype canonicalization, or target nomenclature
inference is performed.

### Substep 7.1 — Difference extraction

Walking the alignment columns:

- a column where the query is `-` is a **deletion** of the reference bases;
- a column where the reference is `-` is an **insertion** of the query bases;
- a column with unequal canonical query and reference bases is an **SNV**.

A column on a masked call ([[design/callability|callability]]) is no SNV candidate
and produces no record. An insertion with a masked inserted call is excluded
with the mask reasons of those calls. Other differences whose allele contains a
non-canonical base, or whose indel length exceeds `max_indel_length`, increment
the excluded-candidate warning count rather than being reported.

### Substep 7.2 — Allele construction and anchoring

The caller converts each selected-alignment difference into a validated
reference-oriented allele representation.

For indels, the current production implementation preserves the canonical gap
placement selected by alignment under `SRS-ALN-012`; it does not independently
left-shift, rotate, or right-shift the event through a repeat. When an aligned
left flank is unavailable, the caller derives the real reference predecessor
where possible. A true leading linear indel uses a real right anchor. Circular
origin-spanning representation preserves the alignment-selected side of the
rCRS seam.

Internally each variant retains its contig, 1-based position, reference/alternate
alleles, kind, and direct call mappings. Compact analysis v9 emits only `position`,
`reference`, `alternate`, `kind`, and `calls`. Every public call contains
only its supporting/flanking `role`, reference-oriented called `base`,
co-located reference-oriented A/C/G/T primary-event channel heights in `peaks`,
and uncalibrated `quality`. Original call index, PLOC, mapped call position,
trace-strand symbols, selected-peak positions/sources, penalties, and vendor
evidence remain internal. Deletions carry real aligned flanks only and never
fabricate deleted-base dna. A flanking call whose primary base is unresolved
(`N`: its strongest channels tie, carry no signal, or all qualify) is omitted
from the public calls, which admit only A/C/G/T bases. The emitted reference allele is validated against
the supplied reference.

The implemented post-calling `variant_normalization` capability defined by
ADR-0060 may move sequence-equivalent events only after reconstructing and
preserving the complete represented haplotype. That capability is intentionally
separate from this caller.

### Substep 7.3 — Configured eligibility

A called candidate is retained only when its 1-based anchor `position` lies
inside at least one configured inclusive region. SNV supporting calls and every
inserted-base supporting call must each have a highest A/C/G/T peak greater than
or equal to `sanger_evidence.minimum_peak_height`, unless its uncalibrated
relative score reaches `sanger_evidence.peak_rescue_quality` and its second
base holds less than `sanger_evidence.peak_rescue_second_share` of the two
strongest peaks, and an
uncalibrated relative score strictly greater than
`sanger_evidence.relative_quality_threshold`. Insertion flanks are not evaluated.
Deletions have no supporting trace base, so their flanks are not subjected to
peak or quality thresholds; their caller anchor must still be in a region.
Read eligibility (ADR-0062, ADR-0067, SRS-VAR-013) then marks a variant
ineligible when any of its mapped calls is one the modality places at an edge
of the read's usable signal (`read_end`), or lies within `read_end_margin`
calls of a physical read end the modality does not vouch for. The Sanger
modality places at an edge every call within its own `signal_edge_margin`
calls of an
uninformative call: a call beyond the trim interval, a masked call aligned as
unresolved, or a read end (SRS-CALL-013). A difference therefore needs that
many informative calls between it and the nearest place where the alignment
has no information, at the read's ends and around internal masked segments
alike. A call the modality vouches for by its base only (a leading call,
SRS-CALL-014, or a call beside a junction whose spacing marks a missing
peak, SRS-CALL-026) supports a substitution, but no insertion or deletion and no
run's end; and where one of its neighbours is masked, `N`, or vouched for
by its base only, and its base is the reference base next to it, it may sit
one position off, so it supports no substitution either (`unresolved_flank`,
SRS-VAR-018). The calls of a run the read bounds are exempt where the
margin's edge is the read losing phase right after a repeat run (ADR-0076):
an insertion or deletion of that run is then not `read_end`, and the
consensus trusts the run's calls where they match the reference or are
inserted, and the run's bounds where they match the reference, while
substituted calls and bounds that read another base keep the margin. A read's physical ends count only when its modality does not
vouch for them: Sanger reads never do, and reviewed consensus sequences always
do, so a consensus difference at its first or last base stays eligible. A masked
evidence call — a supporting call, or a flanking call of a deletion — adds the
reason of its segment: `post_homopolymer` when the segment starts in the window
after a repeat run, otherwise `dephased_signal`, `mixed_signal`, `weak_signal`,
or `irregular_spacing`. Trace order is the sequencing direction for both
strands, so forward and reverse reads are handled alike without primer
knowledge.

A read must also resolve the end of any run of one base whose length a
variant changes (ADR-0071, SRS-VAR-016), or the variant gets `run_boundary`.
The run is the read's maximal stretch of the edited base around the edit: the
substituted base, inserted bases that are all one base, or what a deletion of
one repeated base leaves of its run. The call beyond the run resolves its end
when it is informative and reads another non-`N` base, either unmasked or
dephased: a read that dephases right after a run still shows where the run
ends. An unvouched physical read end resolves nothing; a vouched one does.
Insertions and deletions need both ends resolved and the whole run unmasked;
a substitution needs only the ends it lies on. A lone substituted or inserted
base, an edit of several distinct bases, and the deletion of a whole run
change no run's length. The variant's reference-strand alleles are compared
with the read's trace-strand calls through the read's placement.

Vendor PCON is not used by this filter. For SNVs, a supporting call with more than one co-localized qualifying channel is retained as an observed called difference but is ineligible for clean-SNV reporting with `mixed_supporting_dna`. Insertions and deletions are not subjected to this point-mixed-signal gate; persistent mixed-length evidence is a separate method boundary.

The caller reads each read only as `ReadEvidence`
([[ADR-0069]]).
The region gate, `indel_length_exceeded`, `non_canonical_allele`,
`read_end`, `run_boundary`, `unresolved_flank`, and `placement_conflict` are
core rules. Everything else arrives from the modality as data:
- the Sanger adapter (`read_processing::evidence`) raises the support vetoes
  `peak_below_minimum` and `relative_quality_not_above_threshold`, which apply
  to SNVs and insertions, and `mixed_supporting_dna`, which applies to SNVs
  only;
- it labels masked calls with their segment reason.

An insertion whose inserted calls extend a repeat of their unit in the read,
where the read bounds the run, is judged for support vetoes on its placements
along the repeat (SRS-VAR-020, `Eligibility::repeat_unvetoed`): any window
of as many repeat calls that is unmasked and free of insertion vetoes clears
it. The alignment takes the repeat's 3'-most calls as the inserted ones
(SRS-ALN-012); which identical calls those are is arbitrary, so a weak call
among them must not decide. The bounded run is the one ADR-0076's exemptions
read (`length_only`), and the consensus applies the same rule to a junction's
inserted calls.

A variant's reasons are reported in this order:
1. `placement_conflict`, for every variant of a consensus record placed away
   from the interval its name declares (SRS-VAR-017);
2. `outside_target_region`;
3. the union of its supporting calls' vetoes, in the adapter's vocabulary
   order (the order listed above);
4. `read_end`;
5. `run_boundary`;
6. `unresolved_flank`, when a mapped call is `N` (SRS-VAR-018); a deletion
   neither of whose flanks reads a base is not observed at all, only
   recorded as excluded with this reason;
7. the mask reasons of its evidence calls, in mapping order and without
   duplicates.

Each removed candidate increments `excluded_variant_candidates` once, even when
it fails more than one eligibility condition. The pure variant stage also returns
a concise exclusion diagnostic containing kind, contig, caller position when
available, and all failed rules. The core's one-read path (`dna_call::read_call`)
writes one WARN record per diagnostic without reference/alternate alleles. Sample aggregation logs aggregate
counts of differential-locus observations and variant-associated calls that
retain a basecall-independent profile; those operational counts do not alter the
scientific result.

### Substep 7.4 — Ordering

Reported variants are sorted by `(contig, position, reference, alternate)` and
deduplicated.
