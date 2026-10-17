# ADR-0077: Post-calling split into normalize and EMPOP crates

- **Status:** Accepted
- **Date:** 2026-10-17
- **Related decisions:** [ADR-0047](0047-canonical-right-aligned-mtdna-gaps.md),
  [ADR-0060](0060-separate-variant-canonicalization-nomenclature.md),
  [ADR-0069](0069-plugin-first-modality-core-post-calling.md),
  [ADR-0070](0070-workspace-split-by-plugin-family.md) (amends §1)

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
