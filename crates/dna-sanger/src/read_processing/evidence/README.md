# Evidence

The Sanger modality adapter: one read's `ReadEvidence` from its processed
products ([ADR-0069](../../../../../docs/decisions/adr/0069-plugin-first-modality-core-post-calling.md)).
What the evidence holds, and in what order the rules apply, is described in
[read_processing](../README.md).

- `mod.rs`: `read_evidence`, which builds the per-call verdicts and applies
  the rules below in order; the Sanger reason labels and `reason`, which
  reads a published label back; each masked call's placement; and the
  read's poor start (SRS-CALL-024).
- `identity.rs`: the calls the read's own signal vouches for by their base:
  the leading calls (SRS-CALL-014), the edge calls (SRS-CALL-017), and the
  start calls (SRS-CALL-018). Each call is judged against the norms
  (median primary peak, width, and spacing) of its own read.
- `reference_only.rs`: the read's reference-only calls (SRS-CALL-022,
  ADR-0080).
- `post_repeat.rs`: the bounded runs before a post-repeat phase loss
  (ADR-0076).
- `tests.rs`: the adapter's tests.

The submodules read the read's own calls, loci, and callability, never a
reference position, alignment, or truth.
