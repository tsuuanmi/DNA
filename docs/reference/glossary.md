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
