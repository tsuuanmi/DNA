# Variant Nomenclature

Variant nomenclature is an optional target-specific representation layer after
variant calling and, for the current human-mtDNA path, after haplotype-preserving
normalization.

It does not call variants, reinterpret Sanger signal, or change the reconstructed
biological sequence.

## Current implemented scope

The first implemented rule is the human-mtDNA HVS-II poly-C window at rCRS
positions 303-315:

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

`variant_nomenclature::mtdna::apply_hv2_polyc`:

1. loads and validates the supplied reference identity;
2. reconstructs the complete alternate sequence from normalized variants;
3. verifies it equals the alternate sequence retained by normalization;
4. isolates edits wholly contained in the HVS-II 303-315 window;
5. reconstructs the local alternate window;
6. when that window contains only C bases around exactly one T anchor, derives
   left- and right-run length deltas relative to rCRS;
7. replaces only the local decomposition with run-length edits at the 309 and
   315 boundaries;
8. preserves all variants outside the window;
9. reconstructs the complete represented haplotype and requires byte-for-byte
   equality with the normalization result before returning.

An edit crossing the window boundary fails closed rather than being partly
rewritten.

## Shared mechanics

Sequence-edit conversion, whole-haplotype application, deterministic edit
ordering, and anchored public-variant rendering are target-independent mechanics
owned by the crate-internal `variant_representation` module.

`variant_normalization` and `variant_nomenclature` both reuse those mechanics
while owning separate policy and error boundaries.

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

- HVS-III 513-524 AC-repeat nomenclature;
- HVS-I 16189/16193 poly-C nomenclature;
- Sanger-specific repeat artifact interpretation;
- decimal-string rendering;
- sample reconciliation;
- VCF/HGVS formatting;
- NGS-specific behavior.
