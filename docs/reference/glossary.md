# Glossary

- **Analyzed channel:** one of the processed A/C/G/T fluorescence arrays stored in the canonical ABIF DATA.9-12 tags used by DNA.
- **Call index:** 0-based identity of a locus in the ordered read-call sequence.
- **PLOC:** vendor-defined 0-based trace-sample position for a called locus; the current MVP uses PLOC.2 as the locus authority.
- **Primary call:** DNA's strongest nucleotide interpretation at one locus. It is a read-level interpretation, not guaranteed biological truth.
- **Ambiguity:** DNA's representation of significant multi-channel evidence at a locus; it is not automatically genotype or heteroplasmy.
- **Unresolved:** evidence for which DNA intentionally declines a canonical interpretation, commonly represented by `N`.
- **Source evidence:** validated decoded chromatogram data and required metadata from which downstream observations are derived.
- **Read-level evidence:** observations and differences inferred from one chromatogram.
- **Sample-level evidence:** conclusions aggregated from multiple independently acquired read observations; not part of the current single-trace core.
- **Primary-sequence difference:** a difference between the selected read primary sequence and the supplied reference. It is not, by itself, a genotype or clinical claim.
- **Core confidence floor:** the smallest scientific path that must be well understood and validated first.
- **Current supported baseline:** capabilities the product already implements intentionally and retains while validation continues.
- **Research:** non-normative exploratory documentation under `docs/research/`.
- **Read callability:** the signal-derived, per-read view of where a read is still one ladder: phase-state segments, a per-position mask, and the callable span (`dna.read_callability/v1`).
- **Phase state:** the behaviour of the measured signal over a run of calls — `in_phase`, `dephased` (double peaks explained by the main ladder plus shadows of the primary calls one to three calls away, at least one of them one call away: slippage or a ladder smeared by one call), `mixed` (double peaks the shadow model does not explain), `weak`, or `irregular` (call spacing).
- **Masked call:** a call inside a non-in-phase segment; it keeps its call and evidence but is declared not callable by the read's own signal: it supports no variant, aligns as unresolved unless dephased, and observes a sample locus as `masked`.
- **Callable span:** the 0-based half-open call interval from the first to the last unmasked call; empty when every call is masked. The trim interval is the callable span plus at most `read_end_margin` adjacent dephased calls per side.
- **Callable reference segments:** the reference intervals a read observes with unmasked calls, and with deletions between them; inside its mapped segments, a position outside them is covered by a masked call.
