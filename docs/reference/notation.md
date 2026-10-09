# DNA Notation JSON

`signal notation <sample-id> <sample-id>.variants.json --reference <reference.fasta>`
runs the post-calling plugins over one variants document
([`dna.variants/v1`](variants.md),
[ADR-0069](../decisions/adr/0069-plugin-first-modality-core-post-calling.md))
and writes one deterministic `results/<sample-id>.notation.json` document
identified as `dna.notation/v1`. The authoritative schema is
[`schemas/notation-v1.schema.json`](schemas/notation-v1.schema.json), and the
synthetic example is
[`examples/notation-v1.example.json`](examples/notation-v1.example.json).

The target profile named by the configuration must declare notation. The
document's sample must be the requested sample, and its reference must be the
supplied reference. Only eligible variants are represented.

## Fields

- `provenance`:
  - the reference `name`, `topology`, and sequence `sha256`;
  - `configuration_sha256`;
  - the target `profile`;
  - `plugins`: `normalization` and `nomenclature`, then `conformance` when the
    profile declares it;
  - `source`: the `schema_version` and the SHA-256 of the variants document's
    bytes.
- `notation`: the per-read view of [sample notation](sample-evidence/notation.md),
  with reads named as in the source document. For the same profile it equals
  the `notation` of the source document.
- `conformance` (only when the profile declares it):
  - `rules`: the profile's notation conventions, in order;
  - `findings`: every represented variant of a read that does not follow a
    rule, as its rendered `calls`, the `read`, and the `rule`.

  Findings are reported, never applied. A finding names a convention that the
  representation does not follow; it is not an error in the haplotype.

## Rules

- `insertion_at_run_end`: an insertion of one repeated base that touches a
  reference homopolymer of that base, at least `minimum_run_length` long, must
  follow the homopolymer's 3' base (`315.1C`, not `311.1C`).
- `insertion_matches_run`: bases inserted between two bases of such a
  homopolymer must be its base. Another difference there is written as a
  substitution plus a length change (`314G 315.1C`, not `313.1G`).
