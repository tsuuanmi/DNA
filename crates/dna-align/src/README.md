# dna-align

Owns placement, the core's first step
([PROP-0008](../../../docs/proposals/0008-pipeline-responsibilities.md),
[ADR-0079](../../../docs/decisions/adr/0079-core-crates-by-step.md)): where
does a read lie on the reference? It depends only on `dna-kernel`.

Public surface: `align_best`, [placed](placed.rs) (the `dna.placed_read/v1`
contract: `write` and `read`, which re-checks a placement against the read's
evidence and the reference, [ADR-0087](../../../docs/decisions/adr/0087-placed-read-contract.md),
proposed), the `[alignment]` section (`AlignmentConfig`,
`RawAlignmentConfig`), and the [model](model/README.md): a read-only `Alignment`,
columns included, whose callable segments only its own methods narrow, so
`dna-align` is its only writer. `Alignment::placed` is the second producer
of alignments, for a source that places a read itself (the consensus's
decided reads): it checks the invariants `align_best` guarantees and
computes the metrics from the columns. `placed::read` is the third: it
reads back a placement `align_best` selected and `placed::write` wrote, and
re-checks it against the read's evidence and the reference. The
`test-support` feature adds
`Alignment::fixture` for the guards of downstream crates.

Output contract: the `alignment` and each call's `column` of
[`dna.analysis/v10`](../../../docs/reference/schemas/analysis-v10.schema.json) (`dna analyze`).

Owns bounded, deterministic affine-gap alignment, strand selection, traceback,
and canonical repeat-equivalent gap placement.

Entry point: `align_best` from `orient.rs`, over the read's `ReadEvidence`
([ADR-0069](../../../docs/decisions/adr/0069-plugin-first-modality-core-post-calling.md)).
The module imports no Sanger type.

Key children: `scoring.rs`, [`exact/`](exact/README.md), `gotoh.rs`,
`traceback.rs`, `canonical.rs`, `runs.rs`, and `orient.rs`.

The query is the evidence's informative interval. Calls with an *unresolved*
mask enter it as `N` without a profile; calls with an *anchoring* mask keep
their base and profile. `orient.rs` also derives the masked-base counts and the
callable reference segments: the reference observed by unmasked calls that
read a resolved base. `canonical.rs` places every gap at its 3'-most
sequence-equivalent placement with the shared move of
`dna_kernel::model::placement`, keeping a move when the primary-base score is
unchanged (SRS-ALN-012). Before that, on a read whose ends are vouched (a
consensus record), `canonical.rs` moves a deletion beside an `N` without
evidence 3' across it (SRS-ALN-017), and `runs.rs` re-expresses masked
repeat stretches of the selected alignment as one length edit per run
(ADR-0072). When the selected orientation's two best placements tie only on
calls at one end of the query that the read asserts unreliable (inside
`ReadEvidence::unreliable`, or masked), `orient.rs` leaves those calls
unplaced and aligns the rest again (SRS-ALN-006, ADR-0086).

This module does not extract variants or mutate upstream signal evidence.

The evidence-profile Gotoh implementation is intentionally first-party after reuse evaluation rather
than by default. `exact/` owns proof-based acceleration: direct upper-bound placement for exact profile-optimal reads and score-bounded q-gram pruning for low-edit reads. Candidate seeds can provide a lower bound, but only a completeness proof may exclude reference space or an orientation. Both tiers reuse `memchr` for maintained SIMD substring localization and fall back to full Gotoh whenever proof or cost conditions are not satisfied; an orientation evaluated exactly, by proof or by full Gotoh, sets the threshold that can reject the other. Current `rust-bio` pairwise alignment does not directly satisfy
DNA's required evidence-profile scoring and `open + k * extension` gap contract;
DNA also requires deterministic repeat-equivalent placement, explicit placement
ambiguity handling, and circular-reference span semantics. The owning method
document records the full rationale. A future maintained implementation may
replace this code only if it satisfies the same scientific contract and is
validated independently.

`config.rs` owns the `[alignment]`, with the alignment-cell cap configuration section: its raw record, validated
record, and rules (ADR-0069).

See [alignment requirements](../../../docs/requirements/alignment.md),
[alignment method](../../../docs/design/alignment.md), and
[alignment invariants](../../../docs/architecture/invariants/alignment.md).

See [crates](../../README.md).
