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
| HVS-II | C runs around T310 | run lengths at 309/315; duplicated anchor `311T 315.1C`; anchor deletion `310C 315DEL`; one-base gain ending in C as substitutions + EMPOP `315.1C` |
| HVS-III | `AC` tandem repeat | one-motif loss as substitutions + `523DEL 524DEL` |
| HVS-I | C runs around T16189 | phylogenetic `16183C 16184A 16189C`; anchor deletion `16189C 16193DEL` |

These encode the `mtdna_raw` profiles. Its "preserve interpreted
substitutions" rule is the same as keeping the normalized form here, so it has
no separate implementation. Its EMPOP terminal rule also accepts multi-base
gains through a full edit-script derivation; DNA accepts only the validated
one-base gain and otherwise keeps the normalized form.

An edit crossing a window boundary fails closed rather than being partly
rewritten.

## Adding a target

A target whose repeats fit the existing structures and rules needs only a
profile file. A new structure or rule is a code change: it must propose
candidates that the engine can verify against the window haplotype, and its
preconditions must be checked when the profile loads, so an invalid pairing
never reaches the engine.

## Shared mechanics

Sequence-edit conversion, whole-haplotype application, deterministic edit
ordering, and anchored public-variant rendering are target-independent mechanics
owned by the `variant_representation` module of `dna-post`.

`variant_normalization` and `variant_nomenclature` both reuse those mechanics
while owning separate policy and error boundaries.

## Workflow composition and notation

The `sample`, `call`, and `notation` commands compose this per read when the profile declares notation
(SRS-NOM-010). Each read's eligible variants are normalized under the profile's
indel placement by `variant_normalization` and then passed through the profile
windows, so that, for example, the HV2F and HV3R descriptions of one poly-C
haplotype converge before they are compared. Composition lives in
`pipeline::represent`, following ADR-0060 §7. An edit that straddles a
window keeps its normalized form, because SRS-NOM-007 forbids partially
rewriting it.

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
