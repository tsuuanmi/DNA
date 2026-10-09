# Target Profiles

A target profile holds the knowledge about one sequencing target; the
[configuration](configuration.md) holds method parameters and names one profile
with its root `profile` key. Requirements are [`SRS-PRF-*`](../requirements/profiles.md);
the rationale is [ADR-0063](../decisions/adr/0063-target-profiles.md). The
shipped profile is [`config/profiles/human-mtdna-rcrs.toml`](../../config/profiles/human-mtdna-rcrs.toml).

Profiles are strict TOML of at most 1 MiB. Unknown or missing keys and every
violated constraint fail at load. Results record the profile `id` and the
SHA-256 of the file bytes in `provenance.profile`.

## Schema version 1

| Section | Key | Required | Validation |
|---|---|---|---|
| root | `schema_version` | yes | exactly `1` |
| | `id` | yes | 1-64 of `a-z 0-9 . _ -`, starting with a letter or digit |
| `reference` | `topology` | yes | `linear` or `circular` |
| | `sequence_sha256` | no | 64 lowercase hex digits; another reference fails closed |
| `variant_calling` | `regions` | yes | non-empty inclusive 1-based ranges within `1..=50000`, used as a union |
| `normalization` | `indel_placement` | with `notation` | `right` (3'/right-most, never across the FASTA seam) |
| `notation` | `style` | with `normalization` | `per_base_decimal` (`73G`, `249DEL`, `309.1C`) |
| `nomenclature` | `windows` | no | non-empty array of tables, in reference order without overlap |
| `conformance` | `rules` | with `notation` | non-empty, without repetition: `insertion_at_run_end`, `insertion_matches_run` ([notation result](notation.md#rules)) |
| | `minimum_run_length` | with `rules` | at least `2`: the shortest reference homopolymer the rules apply to |

`normalization` and `notation` form the notation chain of `sample`, `call`, and
`notation`, and are declared together or not at all. `conformance` declares the
notation conventions that `notation` checks and reports; it never changes a
call. Nomenclature windows are used by that chain and by the
`variant_nomenclature::apply` API.

Each `[[nomenclature.windows]]` table:

| Key | Validation |
|---|---|
| `name` | non-empty, unique |
| `start` | 1-based position of the first window base |
| `sequence` | uppercase A/C/G/T; a reference that does not carry it at `start` fails when it is loaded |
| `structure` | `{ kind = "anchored_homopolymer", repeat_base, anchor, anchor_base }` (anchor is a 1-based position inside the window holding `anchor_base`, which differs from `repeat_base`) or `{ kind = "tandem_repeat", motif }` (shorter than the window and ending it) |
| `rules` | non-empty, without repetition, tried in order |
| `canonical_haplotype` | exactly when `canonical_haplotype` is a rule: `[{ position, base }, ...]`, unique in-window substitutions that change the reference base |

| Rule | Structure | Representation |
|---|---|---|
| `anchored_run_lengths` | anchored homopolymer of repeat bases around one anchor | anchor movement as run-length change at the run ends |
| `anchor_duplication` | anchored homopolymer, repeat base after the anchor | duplicated anchor as the base after it plus a terminal repeat insertion |
| `anchor_deletion` | anchored homopolymer, repeat base after the anchor | deleted anchor as a repeat base plus deletion of the last run base |
| `terminal_repeat_insertion` | anchored homopolymer | one-base gain ending in the repeat base as substitutions plus a terminal insertion, when no longer than the normalized form |
| `canonical_haplotype` | any | the declared substitutions, when they reproduce the window haplotype |
| `motif_shift` | tandem repeat | loss of one motif copy as substitutions plus deletion of the terminal copy |

A rule only proposes; the first candidate that reproduces the window haplotype
exactly is used, and a window no rule matches keeps its normalized form
([design](../design/variant-nomenclature.md)).
