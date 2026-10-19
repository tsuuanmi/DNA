---
implements: [SRS-NOM, SRS-PRF, SRS-VN]
decided-by: [ADR-0063]
covers: [config/profiles/human-mtdna-rcrs.toml, crates/dna-kernel/src/profile/, src/input/mod.rs]
---

# Target Profiles

A target profile holds the knowledge about one sequencing target; the
[[reference/configuration|configuration]] holds method parameters and names one profile
with its root `profile` key. Requirements are [[requirements/profiles|`SRS-PRF-*`]];
the rationale is [[ADR-0063]]. The
shipped profile is [`config/profiles/human-mtdna-rcrs.toml`](../../config/profiles/human-mtdna-rcrs.toml).

Profiles are strict TOML of at most 1 MiB. Unknown or missing keys and every
violated constraint fail at load. Results record the profile `id` and the
SHA-256 of the file bytes in `provenance.profile`.

## Schema version 1

| Section | Key | Required | Validation |
|---|---|---|---|
| root | `schema_version` | yes | exactly `2` |
| | `id` | yes | 1-64 of `a-z 0-9 . _ -`, starting with a letter or digit |
| `reference` | `topology` | yes | `linear` or `circular` |
| | `sequence_sha256` | no | 64 lowercase hex digits; another reference fails closed |
| `variant_calling` | `regions` | yes | non-empty inclusive 1-based ranges within `1..=50000`, used as a union |
| `normalization` | `indel_placement` | with `notation` | `right` (3'/right-most, never across the FASTA seam) |
| | `max_trade` | with `notation` | non-negative integer: farthest a deletion moves 3' past substitutions in one trade (SRS-VN-008); `0` disables trades |
| `notation` | `style` | with `normalization` | `per_base_decimal` (`73G`, `249DEL`, `309.1C`) |
| `nomenclature` | `windows` | no | non-empty array of tables, in reference order without overlap |
| | `conventions` | no, with `notation` | array of tables, in reference order without overlap |
| `conformance` | `rules` | with `notation` | non-empty, without repetition: `insertion_at_run_end`, `insertion_matches_run` ([notation result](notation.md#rules)) |
| | `minimum_run_length` | with `rules` | at least `2`: the shortest reference homopolymer the rules apply to |
| `empop` | (none) | no, with `notation` | an empty table: `notation` exports an EMPOP query ([notation result](notation.md#empop-query)) |

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
| `structure` | `{ kind = "anchored_homopolymer", repeat_base, anchor, anchor_base }` (anchor is a 1-based position inside the window holding `anchor_base`, which differs from `repeat_base`) or `{ kind = "tandem_repeat", motif, read_copies }` (`motif` shorter than the window and ending it; `read_copies`, positive, the terminal motif copies a read must show for a length change) |
| `rules` | non-empty, without repetition, tried in order |
| `canonical_haplotype` | exactly when `canonical_haplotype` is a rule: `[{ position, base }, ...]`, unique in-window substitutions that change the reference base |

| Rule | Structure | Representation |
|---|---|---|
| `anchored_run_lengths` | anchored homopolymer of repeat bases around one anchor | the lab convention: the first anchor base splits the window; the left run's length changes at its 3' end, and the rest is written left-anchored against the right run (substitutions, then repeat-base insertions or deletions at the window's 3' end) |
| `anchored_tract` | anchored homopolymer led by a run of another base, then repeat bases around one anchor | the lab convention for a lead run before anchored runs: with the anchor before a right run of its reference length, the lead run's length changes at its 3' end, the left run is written left-anchored with its length change at its 3' end, and the right run as substitutions; otherwise the whole window is written left-anchored |
| `anchor_duplication` | anchored homopolymer, repeat base after the anchor | duplicated anchor as the base after it plus a terminal repeat insertion |
| `anchor_deletion` | anchored homopolymer, repeat base after the anchor | deleted anchor as a repeat base plus deletion of the last run base |
| `terminal_repeat_insertion` | anchored homopolymer | one-base gain ending in the repeat base as substitutions plus a terminal insertion, when no longer than the normalized form |
| `canonical_haplotype` | any | the declared substitutions, when they reproduce the window haplotype |
| `motif_shift` | tandem repeat | loss or gain of whole motif copies at the window's 3' end, plus position-wise substitutions against the remaining reference; an unchanged length as substitutions |

A rule only proposes; the first candidate that reproduces the window haplotype
exactly is used, and a window no rule matches keeps its normalized form
([[design/variant-nomenclature|design]]).

Each `[[nomenclature.conventions]]` table declares a report convention: rendered
calls that the target's reviewers write differently. A covered call is
labelled with the convention in every notation view and kept; no call is
dropped or rewritten (SRS-NOM-017).

| Key | Validation |
|---|---|
| `name` | non-empty, unique |
| `label` | `not_recorded`: the reviewers do not record the call |
| `kind` | `run_length` |
| `window` | `[start, end]`, inclusive 1-based positions with `1 <= start <= end` |
| `base` | A/C/G/T: the run's base |
| `minimum_run` | at least `2` |

`run_length` covers an inserted `base` whose anchor lies in `window`, and a
deleted `base` in `window`, when the read's represented haplotype has an
uninterrupted run of `base` of at least `minimum_run` bases at the call,
counting a deleted base in the run it shortens. A substitution is never
covered. The shipped human-mtDNA profile declares `hvs1_c_length`
(`16180-16193`, `C`, `10`): the reviewers do not record the length of the
uninterrupted HVS-I C run that `16189C` creates.
