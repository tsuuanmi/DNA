# Variant Normalization Method

Part of the canonical [DNA pipeline architecture](../architecture/variant-lifecycle.md).

Variant normalization is an optional post-calling representation transform. It
does not call variants from sequencing evidence and does not apply
target-specific nomenclature.

## Inputs

The current public capability receives:

- one single-record FASTA reference;
- one `CalledVariantSet` carrying the expected reference identity;
- one explicit `NormalizationPolicy`.

The implemented policy is `RightAligned`, selected by a target profile's
`indel_placement = "right"` in the `sample` workflow (ADR-0063).

## Step 1 — Validate reference identity and source alleles

The supplied FASTA is loaded through the canonical reference loader. Its name and
sequence SHA-256 must match the `CalledVariantSet` reference identity.

Every source variant must:

- use the same contig name;
- use a positive 1-based anchor;
- contain uppercase A/C/G/T/N alleles;
- have a kind consistent with its REF/ALT lengths;
- carry a REF allele that exactly matches the supplied reference at its current
  anchored span.

Origin-spanning source REF alleles are currently rejected by this standalone
normalizer rather than silently rotated across the canonical seam.

## Step 2 — Convert anchored variants to minimal sequence edits

Caller REF/ALT anchoring is removed mechanically by trimming the common prefix
and common suffix from each allele pair.

Conceptually:

~~~text
anchored deletion: CAA -> C
minimal edit: delete AA after C

anchored insertion: C -> CAT
minimal edit: insert AT after C
~~~

The resulting edit uses zero-based half-open reference coordinates plus an
alternate sequence.

## Step 3 — Reconstruct the complete source haplotype

All source edits are applied jointly to the reference.

Edits that overlap incompatibly, duplicate one insertion boundary, exceed the
reference, or otherwise cannot produce one deterministic sequence are rejected.

The reconstructed sequence is retained in `VariantNormalizationResult` as
`alternate_sequence`.

## Step 4 — Apply the right-alignment policy

Only pure insertions and pure deletions are positionally shifted.

### Deletion shift

A deletion can shift one base right when the first deleted reference base equals
the first reference base following the deletion.

~~~text
reference: C A A A A G
delete:      A
                 ↓ repeated one-base shifts
canonical:         A
~~~

### Insertion shift

An insertion at one reference boundary can shift one base right when its first
inserted base equals the crossed reference base. The inserted sequence is
rotated while the boundary moves so the represented haplotype remains
sequence-equivalent.

For tandem repeats this one-base operation naturally composes into whole-motif
movement.

### Global haplotype guard

A candidate one-base movement is accepted only if applying the **entire**
candidate edit set reconstructs the exact source alternate sequence.

This is important for nearby edits: an indel is not independently pushed through
another SNV or indel when that movement would change phase or sequence.

Gap candidates are considered from the right-most event toward the left so the
result is deterministic.

### Canonical seam

No candidate shift is generated beyond the final FASTA base. A circular
reference such as the mtDNA genome therefore does not permit indefinite
representation rotation across the fixed FASTA coordinate seam.

## Step 5 — Render normalized edits

Minimal edits are converted back to anchored public `Variant` values:

- SNV: direct one-base REF/ALT;
- insertion: adjacent real reference base anchors the inserted sequence;
- deletion: adjacent real reference base anchors the deleted sequence.

The final variants are sorted deterministically by
`(contig, position, reference, alternate)`.

## Invariant

For reference `R`, source edits `S`, and normalized edits `N`:

~~~text
apply(R, S) == apply(R, N)
~~~

If this invariant fails, normalization returns a typed error.

## Ownership boundary

This module owns sequence-equivalent representation movement only.

It does not own:

- chromatogram interpretation or caller eligibility;
- Sanger alignment traceback canonicalization;
- mtDNA 309/315, 513-524, 16189/16193 nomenclature-window rules;
- sample consensus/reconciliation;
- VCF/HGVS serialization or interchange normalization.

Those remain separate capabilities under
[ADR-0060](../decisions/adr/0060-separate-variant-canonicalization-nomenclature.md).
