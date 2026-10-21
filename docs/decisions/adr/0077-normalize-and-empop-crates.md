---
id: ADR-0077
status: accepted
date: 2026-10-17
related: [ADR-0047, ADR-0060, ADR-0069, ADR-0070]
---

# ADR-0077: Post-calling split into normalize and EMPOP crates

<!-- docs:begin lifecycle -->
- **Status:** Accepted
- **Date:** 2026-10-17
- **Related:** [ADR-0047](0047-canonical-right-aligned-mtdna-gaps.md), [ADR-0060](0060-separate-variant-canonicalization-nomenclature.md), [ADR-0069](0069-plugin-first-modality-core-post-calling.md), [ADR-0070](0070-workspace-split-by-plugin-family.md)
<!-- docs:end lifecycle -->

<!-- docs:begin backlinks -->
- **Governs:** [architecture/overview](../../architecture/overview.md), [architecture/variant-lifecycle](../../architecture/variant-lifecycle.md), [design/variant-nomenclature](../../design/variant-nomenclature.md), [reference/rust-api](../../reference/rust-api.md), [reference/rust-api-sample](../../reference/rust-api-sample.md), [validation/discrepancy-catalogue](../../validation/discrepancy-catalogue.md), [validation/evaluation-harness](../../validation/evaluation-harness.md), [validation/evaluation-measures](../../validation/evaluation-measures.md), [validation/evaluation-triage](../../validation/evaluation-triage.md), [validation/evaluation-workload](../../validation/evaluation-workload.md)
<!-- docs:end backlinks -->

## Context

`dna-post` holds every post-calling plugin: representation, normalization,
nomenclature, and conformance (ADR-0070 §1). Two different notation targets
now apply to it, and they disagree.

- **The lab's reviewer convention.** Delivered reports and reviewer truth are
  written this way. Within the repeat windows, a length change is a C
  insertion or deletion at the 3′ end of the C run it belongs to. A moved
  interrupting base is written as a substitution. AC-unit changes go 3′.
  Measured on about 111,000 reviewed corpus samples:
  - `309.1C 315.1C` (45,695) and `315.1C` (37,221);
  - `523DEL 524DEL` (61,013) against `514DEL 515DEL` (19);
  - `310C`, which SWGDAM recommends when T310 is absent: 0.
- **Forensic database conventions (EMPOP).** These use phylogenetic alignment
  (Bandelt and Parson 2008; Huber, Parson, and Dür 2018) and SWGDAM's ordered
  rules. Some HVS-I states depend on the nearest phylogenetic motif.

On development panel v5, 90 samples differ from the truth only in
representation. DNA's notation is the lab's, so the normalizer must follow
the lab convention. An EMPOP form is a separate output with its own rules, and
it must never be guessed from the lab form.

## Decision

1. `dna-post` is split into two post-calling crates:
   - **`dna-normalize`** holds haplotype-preserving representation,
     sequence-equivalent normalization, and target nomenclature. It owns the
     `normalization` and `nomenclature` plugins. Its target is the lab's
     reviewer convention. A rule follows the convention the reviewed corpus
     shows, with its corpus count as evidence, and the corpus wins where it
     departs from SWGDAM or EMPOP.
   - **`dna-empop`** holds conformance to forensic notation conventions and,
     later, an offline EMPOP export of a sample's calls and ranges. It owns
     the `conformance` plugin. It reports and never rewrites
     (ADR-0069). It is deterministic and offline: it neither queries EMPOP nor
     sends sample data anywhere (data policy), and it keeps no corpus
     frequency index.
2. Both crates depend only on `dna-kernel`, and only the `dna` facade composes
   them (ADR-0070 §2). The module validator and CI enforce the new map.
3. The split moves code only. Plugin identifiers, versions, contracts,
   public re-exports (`dna::variant_normalization`,
   `dna::variant_nomenclature`), and every output stay the same.
4. Rule changes come in their own changes, each measured on the development
   panel. Every output must rebuild the identical haplotype; the haplotype
   measures must not change; the primary measure must gain. Reviewer truth
   is never normalized: the harness compares DNA's notation with the truth
   as written. The planned rule changes are:
   - window representation per window rather than per read (SRS-NOM-010);
   - the HVS-III window starting at 512;
   - an HVS-I C-tract rule that keeps 16189C and places length changes at
     16193;
   - broader HVS-II rules.

5. **Responsibility boundary** (amended 2026-10-17).
   - The core caller (`dna-core` with the modality plugin) decides **what**
     the sample's sequence is: placement, callability, consensus, and the
     called variants. It is judged by the haplotype measures. A
     representation it chooses internally, such as a re-expressed run
     (ADR-0072), only makes evidence comparable; it is not the delivered
     notation.
   - `dna-normalize` decides **how** a decided sequence is written, in the
     lab's convention. Every output rebuilds the identical haplotype, and no
     rule asserts a base the read does not read.
   - `dna-empop` checks and exports the **format** and never rewrites.
   - A convention that changes content is not normalization. An example is
     the reviewers not recording the HVS-I C-tract length. Such a convention
     needs its own labelled report policy, which classifies and never deletes.
   - Every difference from the truth is classified by its owner in the
     [[validation/discrepancy-catalogue|discrepancy catalogue]]. A
     class not handled yet stays there as a TODO. Reviewer truth can be wrong:
     a difference that the read evidence attributes to the truth goes to the
     reviewers for a re-check and does not change DNA.

## Alternatives considered

- **Keep `dna-post` and add two crates beside it.** Rejected. Normalization
  code would then live in two crates.
- **One normalizer targeting EMPOP.** Rejected. Notation would then differ
  from the delivered reports and the reviewer truth in about 61,000 corpus
  samples (`523DEL 524DEL`).
- **An online EMPOP client in DNA.** Rejected. It would send sample
  haplotypes to an external service and make results depend on the network.

## Consequences

- The workspace has six member crates. `dna-normalize` replaces `dna-post`'s
  representation, normalization, and nomenclature modules; `dna-empop` takes
  conformance.
- Documentation names the two crates. ADR-0070 §1's `dna-post` entry is
  amended by this decision.
