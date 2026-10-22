---
implements: [SRS-ALN, SRS-API, SRS-CALL, SRS-COMPAT, SRS-NOM, SRS-SAMPLE]
decided-by: [ADR-0065, ADR-0075, ADR-0076, ADR-0077]
---

# Evaluation measures and discrepancy classes

Part of the [[validation/evaluation-harness|evaluation harness]]: what `dna-eval` measures and how it classifies discrepancies.

## Measures

The truth side is compared as stored in the database; DNA's side uses its
notation. Every sample is compared under six measures: primary, strict,
IUPAC-compatible, the two haplotype measures, and haplotype LoD.

Run lengths a profile convention leaves unrecorded (`not_recorded`,
SRS-NOM-017) are not compared under any measure, on either side: an
insertion of the convention's base anchored in its window, or a deletion of
it there, is left out when that side's haplotype with its substitutions only
has a run of the base of at least `minimum_run` at the call. The reviewers do
not write such a length (`hvs1_c_length`: 1.1 % of the truth's samples with
`16189C` carry an indel at 16180-16195), so DNA's length and the reference
length the truth keeps are no difference. Counting the run without either
side's insertions and deletions keeps an over-estimated length from carrying
a run across the minimum.

**Primary.** The lab's concordance (`mtdna_raw`, `mtdna_raw:src/core/comparison.py`,
`SampleComparator.is_concordant`), with DNA as tool A (the automated caller)
and the truth as tool B:

- an allele is `(position, reference, alternate)`. Positions are normalised
  (`16193.0` → `16193`, thousands separators removed), and alleles whose
  reference equals the alternate are dropped;
- `N/A` when either side has no allele;
- insertions and deletions at 309, 455, 463, and 573, and at their `.n`
  insertion positions, are left out ("special" positions); substitutions
  there count. This is the lab's rule, and only the primary measure keeps
  it. DNA is target-agnostic, so the strict, IUPAC-compatible, and
  haplotype measures judge a run length at these positions like any other;
- `N` when a position both sides call holds different allele sets
  (`comparison.py` also skips 309.1, 309.2, and 573.1, which hold insertions
  and are special already);
- `N` when DNA alone calls an allele inside a reportable region of the
  target profile: HV1 (16024–16365), HV2 (73–340), or HV3 (438–576) for
  rCRS, named in the profile's order;
- `N` when the truth alone calls an allele inside an HV region and inside
  DNA's callable intervals, or DNA reports none;
- `Y` otherwise. Interval and region tests use an allele's anchor base, so
  `16193.1` lies in `[16024, 16193]`.

A truth allele where DNA has no callable base is therefore not counted, and
a truth IUPAC allele differs from any definite base DNA calls.

**Strict.** The primary rules with the special positions judged, and a truth-only allele in an HV region also
counts inside the truth's analysed intervals: where DNA has no callable base
it is missed (rule `uncovered`). This makes DNA's coverage gaps visible.

**IUPAC-compatible.** The primary rules with the special positions judged, and DNA's base matches a truth IUPAC
code that admits it, like `iupac_bases_compatible`: two alleles with the same
position and reference match when their called bases share a possible base
(`16189C` and `16189Y`), and where DNA calls nothing inside its callable
intervals its base is the reference, which matches a truth IUPAC allele whose
code admits the reference (`16189Y` over a reference `T`).

**Haplotype and haplotype IUPAC-compatible.** The primary and the
IUPAC-compatible rules, with the special positions judged, applied only to alleles that do not describe the same
sequence on both sides, so a haplotype written or placed differently is not a
difference (`16190.1C` and `16193.1C` in the C run 16190–16193, `248DEL` and
`249DEL` over `AA`, `514DEL 515DEL` and `523DEL 524DEL` in the AC repeat,
`16182.1C 16189C 16193DEL` and `16183C 16184A 16189C`):

- the comparable intervals are the HV regions within DNA's callable
  intervals (merged where they touch), or the HV regions when DNA reports
  none, and within the truth's analysed intervals when it records any: a
  stretch the reviewers did not analyse holds no truth to compare. There is
  no circular seam inside them;
- a DNA allele outside the truth's analysed intervals raises no `dna_only`
  (the primary rules, the lab's own comparator, still count it);
- in each, both sides' alleles are applied to the reference: an allele
  replaces as many bases as its reference has (`514`, `CA`, `-` deletes 514
  and 515; equally many alternate bases substitute them one by one),
  insertions follow their anchor base in ordinal order, and IUPAC codes stay
  as written. Only the reference's length is used, not its bases;
- the interval is walked base by base and cut wherever both local sequences
  have the same length, so the sequences are equal exactly when every piece
  is. A piece that differs is a *block*. A block that reaches the end of its
  interval with unequal lengths is *open*: the length difference may
  continue where DNA has no callable base;
- the haplotype measure compares symbols exactly. The haplotype
  IUPAC-compatible measure also lets a symbol that a truth allele wrote match
  a DNA symbol that shares a base with it (`16191Y` against DNA's `T`, or
  against the reference base where DNA calls nothing), as two alleles do
  under the IUPAC-compatible rules. A symbol that comes from the reference on
  the truth side matches only itself, so a DNA code (`Y`, `N`) where the
  truth calls nothing stays a difference;
- alleles whose bases lie in the interval outside every block describe the
  same sequence and leave the comparison; the rules decide the rest, special
  positions included, and `N/A` is decided on all alleles.
  An interval where one side has two alleles at one base or insertion
  ordinal, or an allele that cannot be applied (an insertion without an
  ordinal, an ordinal on a reference base, an inserted `-`, a reference or
  alternate that is empty or not ASCII, or bases past the interval), cannot
  be built and stays under the rules.

Only representation therefore changes a status relative to the base rules
with the special positions judged: a sample is concordant on a haplotype
measure exactly when every allele the rules count describes the same sequence
as the other side's, or the rules find it concordant anyway.

**Haplotype LoD.** The haplotype measure, except at a truth heteroplasmy
the traces cannot show (`lod.rs`). A truth substitution to a two-base IUPAC
code is *observed* when the covering admitted reads show both of its bases
there, each as some read's call or as some read's second channel (a channel
other than its own call) at 0.10 or more of its normalised A/C/G/T profile
(rounded to thousandths, so 0.0995 counts), or when one of its bases shows
as such a second channel; otherwise it is *below the limit of detection*.
Every observation of an admitted read counts, masked ones included, so the
limit errs toward counting a heteroplasmy as observed. Below the limit the allele compares as under the haplotype
IUPAC-compatible measure: DNA's base there, or the reference where DNA calls
nothing but reads the position, matches when it is one of the code's bases,
in the alleles and in the sequences. At or above the limit the code must
match exactly. A minor reference base where DNA calls the alternate shows the
heteroplasmy as well. The limit is 10 % (user, 2026-10-16; see the
[[PROP-0005|approval protocol]]); each
sample's count of truth IUPAC positions below it is reported (`below_lod`).

**DNA's callable intervals.** The comparable span is DNA's callable
intervals, which leave out unresolved consensus bases (`N`; SRS-ALN-015), so
neither DNA's sequence nor the IUPAC-compatible rules' "where DNA calls
nothing" ever stands for a position DNA left unresolved.

For each sample the harness also reports, on the primary measure, the alleles
unique to each side (special positions left out), every rule that applies,
the positions those rules name, whether every such position holds a truth
IUPAC allele (`iupac_only`), and their windows. The windows are the target
profile's nomenclature windows, then its length tracts outside them, each
labelled by its profile name in lower case with any other character than a
letter or digit written `_`, and `other`: for rCRS `hvs_ii` 303–315,
`hvs_iii` 512–524, `hvs_i` 16180–16193, and `hvs3_c_tract` 568–573. It
lists the `uncovered` truth alleles as well.

The harness holds no reference coordinate of its own: the regions and
windows come from the configuration's target profile, read once
(`geometry.rs`). The laboratory's comparison conventions it reproduces, the
special positions and the region names, are in one module with their
sources (`conventions.rs`).

## Discrepancy classes

Each sample takes one class, which names the part of DNA (or the lab) that
owns its difference from the truth (ADR-0077). A sample that failed under a
measure is `failed`, and one not applicable under a measure (a side without an
allele) is `not_comparable`; any other takes the first of these whose measure
is discordant:

| Class | Measure | Difference | Owner |
|---|---|---|---|
| `failed` | any | DNA failed, or the truth has no record | core or the input |
| `not_comparable` | any | DNA or the truth has no allele | the input, or coverage |
| `sequence` | haplotype IUPAC-compatible | the sequences differ | core, or the truth (adjudication) |
| `heteroplasmy` | haplotype | only at a truth IUPAC code | core minor signal (deferred) |
| `unanalysed` | primary | every primary difference is a DNA call outside the truth's analysed intervals | a review where the truth has no analysis |
| `representation` | primary | the same sequence written differently | `dna-normalize`, or a truth minority form |
| `coverage` | strict | a truth allele where DNA has no callable base | core callability, or the lab's assay |
| `none` | — | every measure concordant | — |

`class_regions` lists the windows of the discordant positions of the class's
measure. Whether a representation is DNA's to fix or a minority form in the
truth is not decided per sample.

## Discrepancy causes

Each sample with a discrepancy also has a cause: who accounts for it, so the
difference reaches the one who can act on it. The core is judged on the
working hypothesis that its input signal is good, or clearly poor through the
lab. Each cause says which case applies:
- `lab`: the lab did poorly (a missing or short trace, a weak run, or a truth
  the reads contradict);
- `technology`: the lab did its part, but the read shows dephasing, noise, or
  length heterogeneity that the core cannot resolve and signal processing in
  `dna-sanger` may later;
- `dna`: a core defect;
- `undetermined`: a reviewer decides.

Every discordant position of
the sample's class gives a finding, `cause:evidence@position` (or `@window`)
in `cause_evidence`; the sample's `cause` is the first of `dna`, `lab`,
`technology`, and `undetermined` among them. Each finding also has a scope
(`cause_scope` for the sample): `defect`, or `deferred` for an item left
undone by design. Three kinds are deferred:
- the IUPAC codes DNA does not write;
- a call the reviewers' convention does not record;
- signal processing in `dna-sanger`, which takes every `technology` finding.

A defect comes before a deferred item, so a sample with both a `lab` and a
`technology` finding is the lab's. The causes come from the
admitted reads' evidence at each region position (with the core's own trust
verdict on each call, `SampleAnalysis::trust`), the limit-of-detection
signals, and how the reviewers write each sequence over the whole truth
database; nothing is listed by hand, and the truth is only read.

| Class | Evidence | Cause |
|---|---|---|
| `coverage` (each uncovered truth allele) | no admitted read covers the position (`no_read`) | `lab`: a missing or short trace |
| | a truth IUPAC code of two bases that a read shows, judged as under `heteroplasmy` (`iupac_not_reported`, `below_lod`) | `dna` deferred, or `technology` deferred |
| | a read covers it with a call masked for phase loss: dephased or mixed signal (a mixed segment counts as phase loss), or a segment after a repeat run (`phase_loss`) | `technology`, deferred: signal processing |
| | a read covers it with a call masked for weak signal or irregular spacing (`weak_signal`) | `lab`: the run's quality |
| | clean reads tie and the consensus writes `N` (`contested_reads`) | `undetermined`: for the reviewer |
| | it lies in a stretch of runs the consensus left undecided, no read showing both run ends (`run_unread`) | `technology`, deferred: signal processing |
| | a read's call there is `N` (`unresolved_call`) | `dna`: its callability policy |
| | a clean call the core does not trust, withheld by a call masked for phase loss: dephased or mixed signal, or a segment DNA starts after a repeat run (`beside_phase_loss`) | `technology`, deferred: signal processing |
| | a weak-signal veto denies a call that reads the reference (`veto_out_of_scope`) | `dna`: the veto guards variant support |
| | a call there shows more than one qualifying channel (`mixed_call`) | `undetermined` |
| | every untrusted call is withheld by its read's own end, its trim, or a call masked for weak signal (`read_edge`), or carries a weak-signal veto (`weak_call`); an ADR-0076 exemption counts as trust | `lab`: read quality |
| | otherwise (`untrusted_call`) | `dna`: its callability policy |
| `unanalysed` (each window) | DNA calls where the reviewers did not analyse (`outside_truth_intervals`) | `undetermined`: only a review there verifies them |
| `representation` (each window) | DNA writes the reviewers' majority form of the window sequence (`truth_minority_form`) | `lab`: a minority form in the truth |
| | the truth writes it (`dna_minority_form`) | `dna`: its notation |
| | outside the windows, the majority of truth samples carrying each side's calls (`…_proxy`) | as above |
| `heteroplasmy` (each truth IUPAC code) | a read shows the code's second base at the limit of detection (`iupac_not_reported`) | `dna`, deferred: DNA writes no IUPAC code by design |
| | a read shows it only at 3 % or more (`below_lod`) | `technology`, deferred: signal processing |
| | clean reads on both strands show no second base at 3 % (`no_second_base`) | `lab`: the truth, for a re-check |
| | one strand only (`one_strand`), or a code of more than two bases (`iupac_code`) | `undetermined` |
| `sequence` (each window) | every DNA-only call is one a profile convention labels `not_recorded` (SRS-NOM-017), and the truth has no call of its own (`convention_not_recorded`) | `lab`, deferred: the reviewers' convention |
| `sequence` (each SNV or deletion) | otherwise, a DNA-only call where clean reads on both strands show only the reference (`phantom`) | `dna` |
| | a truth-only call there (`reads_show_reference`) | `lab`: the truth, for a re-check |
| `failed` | no record in the truth, an unreadable trace, or every trace rejected | `lab` |
| | a DNA error | `dna` |

Any other discrepancy is `undetermined` (`adjudicate`, `not_comparable`).
The written forms group the truth's IUPAC-free samples by the sequence of
each repeat window, so DNA's form and the truth's are compared for the same
sequence; outside the windows no sequence groups them, and the count of
samples carrying each side's calls is a proxy. Each run writes this index as
`forms.tsv`, keyed by the truth database's SHA-256. With `--evidence-cache`
the index is cached per truth database, reference, and repeat windows
([[validation/evaluation-harness#Running DNA|forms cache]]), so a run on the
same inputs reads it instead of indexing the whole database. A sample can need more than
one party: the summary counts the samples per class and cause, defects
first and then deferred items, each in cause order, and per cause, scope, and evidence the samples with such a
finding. Insertions are not judged under `sequence` yet: reads showing the
reference at an anchor do not exclude one.
