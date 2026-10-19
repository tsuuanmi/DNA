---
implements: [SRS-IN]
decided-by:
  - ADR-0001
  - ADR-0002
  - ADR-0004
  - ADR-0005
  - ADR-0006
  - ADR-0008
  - ADR-0010
  - ADR-0014
  - ADR-0022
  - ADR-0026
  - ADR-0044
  - ADR-0045
  - ADR-0048
  - ADR-0049
  - ADR-0050
  - ADR-0051
  - ADR-0052
  - ADR-0054
  - ADR-0055
  - ADR-0056
  - ADR-0070
---

# Retired Documents Ledger

Documents removed from the tree but kept findable
([[governance/documentation-architecture]] section 8). Each row names what was
decided or learned, what replaced it, and the commit holding the full text:

```bash
git show <commit>:<path>
```

A retired ID keeps resolving: `[[ADR-0004]]` and a bare `ADR-0004` point here.
Rows are appended by `docs_graph.py retire`; never edit a row, and never reuse an
ID.

| id | title | outcome | superseded-by | retired | full-text |
|---|---|---|---|---|---|
| ADR-0044 | (number skipped) | the number was never committed to this repository's history | - | 2026-10-19 | never-committed |
| ADR-0045 | (number skipped) | the number was never committed to this repository's history | - | 2026-10-19 | never-committed |
| ADR-0048 | (number skipped) | the number was never committed to this repository's history | - | 2026-10-19 | never-committed |
| ADR-0049 | (number skipped) | the number was never committed to this repository's history | - | 2026-10-19 | never-committed |
| ADR-0050 | (number skipped) | the number was never committed to this repository's history | - | 2026-10-19 | never-committed |
| ADR-0051 | (number skipped) | the number was never committed to this repository's history | - | 2026-10-19 | never-committed |
| ADR-0052 | (number skipped) | the number was never committed to this repository's history | - | 2026-10-19 | never-committed |
| ADR-0054 | (number skipped) | the number was never committed to this repository's history | - | 2026-10-19 | never-committed |
| ADR-0055 | (number skipped) | the number was never committed to this repository's history | - | 2026-10-19 | never-committed |
| ADR-0056 | (number skipped) | the number was never committed to this repository's history | - | 2026-10-19 | never-committed |
| ADR-0002 | Use a Single Layered Crate for the MVP | One crate until reuse or build isolation justifies a split; replaced by the plugin-family workspace split. | ADR-0070 | 2026-10-19 | `65c221e560e1:docs/decisions/adr/0002-single-crate-layering.md` |
| ADR-0004 | Align Directly to Bundled rCRS in the MVP | Direct semi-global alignment to a bundled rCRS with a 50,000-base reference limit (kept as SRS-IN-005); replaced by circular alignment with origin wrapping. | ADR-0010 | 2026-10-19 | `65c221e560e1:docs/decisions/adr/0004-rcrs-direct-alignment.md` |
| ADR-0005 | Use Versioned JSON and Standards-Compliant VCF | Planned versioned JSON plus VCF output; the VCF half was dropped for JSON-only auditable output. | ADR-0008 | 2026-10-19 | `65c221e560e1:docs/decisions/adr/0005-versioned-output-contracts.md` |
| ADR-0006 | Mirror Rust Source Documentation One-to-One | Mirrored Rust source documentation one-to-one under docs/src; replaced by colocated source READMEs and rustdoc. | ADR-0022 | 2026-10-19 | `65c221e560e1:docs/decisions/adr/0006-source-documentation-mirroring.md` |
| ADR-0014 | Project a Compact Analysis Result Summary | Compact analysis result summary (v5); replaced by reviewer-facing signal evidence. | ADR-0026 | 2026-10-19 | `65c221e560e1:docs/decisions/adr/0014-compact-result-summary.md` |
| governance/lifecycle | Documentation Lifecycle | Folded into the documentation architecture standard (classes, status vocabulary, promotion) and the DNA policy. | - | 2026-10-19 | `65c221e560e1:docs/governance/lifecycle.md` |
| governance/ownership | Documentation and Code Ownership | Folded into the DNA policy area table and the standard's change-impact rule. | - | 2026-10-19 | `65c221e560e1:docs/governance/ownership.md` |
| ADR-0001 | Deliver an End-to-End MVP First | Deliver an end-to-end Sanger MVP first (one AB1, one rCRS FASTA, re-calling, QC, alignment, variants); delivered and since extended by the sample, consensus, and plugin architecture. | - | 2026-10-19 | `65c221e560e1:docs/decisions/adr/0001-mvp-vertical-slice.md` |
<!-- ledger:end -->
