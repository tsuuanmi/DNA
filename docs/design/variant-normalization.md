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

Only pure insertions and pure deletions are positionally shifted; a deletion
trade may rewrite the substitutions a deletion passes.

### Deletion shift

A deletion can shift one base right when the first deleted reference base equals
the first reference base following the deletion.

~~~text
reference: C A A A A G
delete:      A
                 ↓ repeated one-base shifts
canonical:         A
~~~

### Deletion trade

When a deletion cannot shift further, it may still move up to 20 bases 3' if
only substitutions lie in the stretch it crosses and the same haplotype needs
no more edits: the bases in between are written as substitutions, which may
replace substitutions the deletion passes. The 3'-most such placement wins.
A deletion whose one-base shift would run into a substitution of the next
base cannot shift, so it trades: `16297DEL 16298C` over `TT` becomes
`16297C 16298DEL`.
With `16297C`, `16293DEL 16298C` becomes `16293C 16298DEL`, as the lab's
reviewers write it; a plain run shift is the trade without substitutions.

A trade asserts every base of the stretch it crosses, so it has two more
bounds, given as `TradeBounds`:

- **Read bases.** Every base of the stretch must lie in the read's callable
  segments; otherwise the rewritten substitutions would describe bases the
  read does not show.
- **Window boundaries.** The stretch may not cross a boundary of a profile
  nomenclature window, in either direction: inside a window, the window's
  rule owns the representation (`16193DEL 16194T` keeps its deletion at the
  HVS-I window's last base). Normalization takes the window spans as data and
  applies no window rule (SRS-VN-010).

The sample, call, and notation workflows pass both; `normalize_with` bounds
neither. Calling every
IUPAC-free reviewed corpus haplotype (about 101,000 samples) from its
sequence, the trades DNA makes span 1, 5, 6, or 10 bases (the longest is
`513DEL` to `523DEL` inside HVS-III), so the 20-base bound is twice the
longest. The trade makes DNA write 36 more of those samples the reviewers'
way (35 of them `16293C 16298DEL`) and none fewer; on development panel v5,
3 more samples are concordant and none regress.

A candidate whose edits would overlap, such as a deletion shifted onto a
substituted base, describes no haplotype and is rejected like one that
changes the sequence.

### Insertion shift

An insertion at one reference boundary can shift one base right when its first
inserted base equals the crossed reference base. While the insertion does not
yet follow a reference base equal to its first inserted base, it can also cross
a base that a substitution of the same haplotype makes equal to it, so it
reaches the 3' end of a reference run of its base: with `16172C`, an inserted C
moves from 16171 to `16174.1C`, as the lab's reviewers write it. An insertion
already at the end of a run stays there (`16262.1C 16263C`). The inserted
sequence is rotated while the boundary moves so the represented haplotype
remains sequence-equivalent.

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
- the mtDNA nomenclature-window rules (HVS-II 303-315, HVS-III 512-524,
  HVS-I 16180-16193);
- sample consensus/reconciliation;
- VCF/HGVS serialization or interchange normalization.

Those remain separate capabilities under
[ADR-0060](../decisions/adr/0060-separate-variant-canonicalization-nomenclature.md).
