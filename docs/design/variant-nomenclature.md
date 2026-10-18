# Variant Nomenclature

Variant nomenclature is an optional target-specific representation layer after
variant calling and haplotype-preserving normalization. The target knowledge —
which windows exist, their reference sequence and structure, and which rules
apply in which order — comes from a target profile (ADR-0063,
[profile reference](../reference/profiles.md)); this layer owns only the generic
engine.

It does not call variants, reinterpret Sanger signal, or change the reconstructed
biological sequence.

## Example: the human-mtDNA HVS-II window

The shipped human-mtDNA profile declares an `anchored_homopolymer` window at
rCRS positions 303-315:

~~~text
303       309 310 311   315
 CCCCCCC   T   CCCCC
 |------|  |   |---|
 left run anchor right run
~~~

The validated reference window is:

~~~text
CCCCCCCTCCCCC
~~~

The `T` at rCRS 310 separates two C runs. Sequence-equivalent alignments may
describe movement of this anchor as substitutions, insertions, or deletions.
Nomenclature instead describes the same haplotype by the lengths of the two C
runs.

For example, a normalized pair:

~~~text
309 C>T
310 T>C
~~~

moves the observed `T` one position left without changing the local sequence
relative to deleting one C from the left run and adding one C to the right run.
The nomenclature layer therefore represents the same haplotype as:

~~~text
delete one C from the 309-side run
insert one C after 315
~~~

The exact public `Variant` representation remains anchored REF/ALT. A later
notation layer may mechanically render those selected insertions as forms such
as `309.1C`, `309.2C`, and `315.1C`.

## Algorithm

`variant_nomenclature::apply(reference, profile, input)`:

1. loads the reference with the profile topology, requires the profile's
   reference when it pins one, and validates the input reference identity;
2. reconstructs the complete alternate sequence from normalized variants;
3. verifies it equals the alternate sequence retained by normalization;
4. for each profile window inside the reference, requires the profile window
   sequence and isolates the edits wholly contained in it;
5. reconstructs the local alternate window;
6. tries the window's rules in order and takes the first candidate that
   reconstructs exactly that local haplotype; with no candidate, the normalized
   edits stay;
7. preserves all variants outside the windows;
8. reconstructs the complete represented haplotype and requires byte-for-byte
   equality with the normalization result before returning.

The rule vocabulary and its structure preconditions are in the
[profile reference](../reference/profiles.md). The human-mtDNA profile uses:

| Window | Structure | Rules, in order |
| --- | --- | --- |
| HVS-II 303-315 | C runs around T310 | run lengths at 309 and 315, the rest left-anchored (`309.1C 315.1C`, `312T 315.1C 315.2C 315.3C`); duplicated anchor `311T 315.1C`; anchor deletion `310C 315DEL`; one-base gain ending in C as substitutions + `315.1C` |
| HVS-III 512-524 | `AC` tandem repeat | whole-copy loss or gain at the 3' end plus substitutions (`523DEL 524DEL`, `513A 523DEL 524DEL`, `524.1A 524.2C`) |
| HVS-I 16180-16193 | A-tract, then C runs around T16189 | phylogenetic `16183C 16184A 16189C`; the A-tract and C-tract rule (`16183DEL`, `16188.1C`, `16182C 16183C 16189C`, `16189C 16190T 16193.1C 16193.2C`); anchor deletion `16189C 16193DEL` |

The rules follow the lab's reviewer convention (ADR-0077). Each reproduces the
form most reviewers wrote for every window sequence of at least 5 corpus
samples (`crates/dna-normalize/tests/fixtures/lab-convention.tsv`); the research
prototype of the same rules reproduced the form of about 99.97 % of all corpus
samples. Where the corpus departs
from SWGDAM or EMPOP, the corpus wins: reviewers never write `310C` for a
missing T310. The windows start before the repeats where the corpus writes a
flank base together with them (`513A 523DEL 524DEL`; `16181C 16182C 16183C
16189C`).

A window that an edit crosses keeps its normalized form, and the read's other
windows are still represented. In `sample`, `call`, and `notation`, a rule may
also rely only on bases the read reads, inside its callable segments: a rule
rewrites the whole window's sequence, and over an unread base that sequence is
the reference's, not the read's. A window therefore keeps its normalized form
when its representation would call an unread base (a substituted or deleted
base, or an insertion's anchor base), or when the bases its choice of form
rests on are unread:

- `anchored_run_lengths`: the anchor and the run after it, and the base before
  the anchor when the left run's length changes;
- `anchored_tract`: the last lead base, the first left-run base, the anchor,
  and the run after it;
- the other anchor rules: the anchor and the run after it;
- `motif_shift`: the window's `read_copies` terminal motif copies (2 in the
  shipped profile) when the length changes;
- `canonical_haplotype`: the whole window, which its choice compares with the
  declared haplotype.

Without that, a read that leaves 16189-16193 unread would see its readable
`16182C 16183C` rewritten as `16182DEL 16183DEL 16188.1C 16188.2C`, a form that
holds only if 16189 is T.

## Adding a target

A target whose repeats fit the existing structures and rules needs only a
profile file. A new structure or rule is a code change: it must propose
candidates that the engine can verify against the window haplotype, and its
preconditions must be checked when the profile loads, so an invalid pairing
never reaches the engine.

## Shared mechanics

Sequence-edit conversion, whole-haplotype application, deterministic edit
ordering, and anchored public-variant rendering are target-independent mechanics
owned by the `variant_representation` module of `dna-normalize`.

`variant_normalization` and `variant_nomenclature` both reuse those mechanics
while owning separate policy and error boundaries.

## Workflow composition and notation

The `sample`, `call`, and `notation` commands compose this per read when the profile declares notation
(SRS-NOM-010). Each read's eligible variants are normalized under the profile's
indel placement by `variant_normalization` and then passed through the profile
windows, so that, for example, the HV2F and HV3R descriptions of one poly-C
haplotype converge before they are compared. Composition lives in
`pipeline::represent`, following ADR-0060 §7. A window that an edit
straddles keeps its normalized form, because SRS-NOM-007 forbids partially
rewriting it; the read's other windows are still represented, so one
crossing edit, such as a deletion of 512-513 beside HVS-III, no longer leaves
an HVS-II run-length change in its raw `309.1C 309.2T 310C` form.

Rendering in the profile's notation style (`per_base_decimal`) is mechanical
serialization (ADR-0060 §8) and lives in `report::notation` (SRS-NOM-011); the report groups
identical calls with their supporting reads (SRS-NOM-012). The same chain
(Tracy alignment → right alignment → mtDNA policy → per-base rendering) is the
approach of the legacy `mtdna_raw` pipeline.

Measured outside the repository (ADR-0065) on 160 held-out samples with
reviewer consensus sequences, the notation reproduces the reviewer calls
exactly for 158 samples. The other two differ only in how the HVS-II window
names the same haplotype, which the conformance checks below report. Remaining
disagreements on Sanger traces come from the Sanger evidence, not from this
representation layer.

## Notation command and conformance

`notation` runs the same per-read composition over the eligible variants of a
`dna.variants/v1` document, so a variants document carried to another process
gets the notation that `call` derives in-process (ADR-0069 phase 4).

When the profile declares `conformance`, the `conformance` module then checks
each read's represented variants against the declared conventions
(SRS-NOM-016). It reads homopolymer runs of at least `minimum_run_length` bases
from the reference sequence and reports findings without changing a call:

- `insertion_at_run_end`: an insertion of one repeated base that touches a run
  of that base must follow the run's 3' base;
- `insertion_matches_run`: bases inserted between two bases of a run must be
  the run's base.

The shipped human-mtDNA profile declares both rules with runs of four or more
bases. They follow how EMPOP/ISFG practice writes length variation in the
control-region C stretches. In the HVS-II window, the run-length
representation can write a G inside the right C run as `313.1G`, while
reviewers write `314G 315.1C`. Both reconstruct the same haplotype, and the
finding names that difference.

Measured outside the repository (ADR-0065) on 160 held-out samples with
reviewer consensus sequences:
- `notation` reproduced the in-process notation for every sample;
- the checker reported exactly the two samples whose notation differs from
  the reviewers' in this way, and nothing else.

## Report conventions

A report convention labels calls that the target's reviewers write
differently, and keeps them (SRS-NOM-017). It is not a window rule: a window
rule chooses among sequence-equivalent forms, while a convention covers calls
that change content, so it may label a call but never drop one (ADR-0077 §5).

The labelling runs on each read's rendered calls after representation. For
`run_length`, it rebuilds the read's represented haplotype, keyed by reference
position and insertion ordinal. For a covered call, it measures the
uninterrupted run of the convention's base:
- an inserted base: the run that holds it;
- a deleted base: the run across the junction it leaves, plus the deleted
  base.

A run shorter than `minimum_run` leaves the call unlabelled. A call that two
reads label differently, or that one read leaves unlabelled, is unlabelled.

The shipped human-mtDNA profile declares `hvs1_c_length`. The reviewers do not
record the length of the uninterrupted HVS-I C run that `16189C` creates
(`16184-16193`, C10 and longer), which carries length heterogeneity. With
`T16189`, the two short runs are read through. A length call there is
evidence, so it stays unlabelled. On development panel v5 (main `4e6f8d1`), 31
`sequence` samples differ from the truth only by labelled calls, and the
[evaluation harness](../validation/evaluation-harness.md#discrepancy-causes)
gives them to the lab as deferred. Two `16193DEL` samples whose `16185T`
splits the run (C8) stay unlabelled.

## EMPOP query

`dna-empop`'s `export` builds the offline EMPOP query of `notation`
(SRS-NOM-018). The same query comes from `SampleAnalysis::call`. It formats and
decides nothing (ADR-0077 §5). EMPOP reads an omission inside a range as the
reference, so a range covers only the positions exactly one record calls: the
records' callable segments, which already leave out every `N` and every
excluded variant's positions (SRS-VAR-019). A record position no range covers
is listed as `unresolved`, `excluded_variant`, or `overlap`.

The mitotype follows EMPOP's numbering, which places an insertion `X.k` at
`X`. A range that holds `X` claims the complete list of insertions after `X`,
so an insertion is written when its anchor lies in a range. That list is
complete because an excluded insertion cuts `X` itself: where `X` survives,
every insertion DNA observed after `X` is eligible. An insertion whose anchor
is cut is not written, even when `X+1` is in a range.

A tie the consensus resolved is DNA's decision. If it is wrong, it is wrong in
the haplotype and in the query alike, so the query exports it as decided.

On development panel v5:
- the query is byte-equal to the laboratory's exporter (`mtdna_raw`'s
  `EmpopProfile`) given the same ranges and calls, for every sample;
- inside the positions both declare, DNA's mitotype equals the reviewers' in
  most samples once hotspot length differences, which EMPOP's standard query
  disregards, are set aside. The rest is mainly truth IUPAC codes.

## Scientific basis

This policy follows the forensic mtDNA convention that insertion notation is
anchored to the preceding reference position with decimal suffixes and that
length variation in the HVS-II C tracts is represented consistently around the
309/315 boundaries. It is intentionally narrower than a complete mtDNA
nomenclature implementation.

Relevant background includes the ISFG mtDNA typing recommendations and the
standardization literature for rCRS-based mtDNA notation.

## Non-goals

This implementation does not yet cover:

- Sanger-specific repeat artifact interpretation and primer callable ranges;
- sample reconciliation or consensus beyond listing supporting reads;
- VCF/HGVS formatting;
- NGS-specific behavior.
