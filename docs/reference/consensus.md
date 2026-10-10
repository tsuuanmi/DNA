# DNA Consensus JSON

`dna consensus <sample-id> <trace>... --reference <reference.fasta>` writes two
files:
- `results/<sample-id>.consensus.json`, identified as `dna.consensus/v1`;
- `results/<sample-id>.consensus.fasta`.

The authoritative schema is
[`schemas/consensus-v1.schema.json`](schemas/consensus-v1.schema.json); the
synthetic example is
[`examples/consensus-v1.example.json`](examples/consensus-v1.example.json).
The method is described in [design](../design/consensus.md).

## Fields

- **`provenance`**, as in [variants](variants.md). Its `plugins` are
  `sanger`, `core`, and `consensus`.
- **`reads[]`:** the admitted reads in the sample's read order (content
  SHA-256), each with `name` (the trace stem) and `sha256`.
- **`rejected_reads[]`:** reads set aside before the consensus, with `name`,
  `sha256`, and `rejection` as in the
  [sample read registry](sample-evidence/reads-coverage.md).
- **`segments[]`:** the consensus sequence.
  - `name`: the FASTA record identifier, `<sample-id>_<first>-<last>` in
    1-based reference positions.
  - `reference`: the covered interval, 0-based half-open.
  - `sequence`: upper-case `A`, `C`, `G`, `T`, and `N` where undecided.
    Inserted bases are spliced in, and deleted positions are left out.
- **`sites[]`:** decided intervals whose call differs from the reference,
  that were not called outright, or whose clean observations disagreed.
  - `start`, `end`: the interval in 1-based inclusive reference positions.
    It covers one position, or a stretch of whole reference runs.
  - `reference`: the reference sequence over the interval.
  - `call`: the decided sequence, with insertions at the interval's edges
    included. It is empty when every position is deleted, and `null` when
    nothing was decided. `N` marks an undecided run or position.
  - `state`:
    - `called`: agreement or a strict majority;
    - `contested`: a tie, with the reference taken for one read against one;
    - `unresolved`: no clean observation.
  - `supporting_reads`, `opposing_reads`, `uninformative_reads`: reads whose
    clean sequence matches the call, whose clean sequence differs, and that
    cover the interval without a clean sequence.
- **`summary`:** counts over the segments:
  - `called_positions`;
  - `contested_sites`;
  - `unresolved_positions` (the `N` count).

## FASTA

The FASTA has one record per segment, with 70-column lines. It is valid input
for `dna call <sample-id> results/<sample-id>.consensus.fasta --reference
<reference.fasta>`, which writes the consensus's `dna.variants/v1` document.

## Workflow

To turn a sample's traces into variants:

1. `dna consensus <id> <trace>... --reference <fasta>` writes the consensus
   and its FASTA.
2. `dna call <id> results/<id>.consensus.fasta --reference <fasta>` writes
   `results/<id>.variants.json` (`dna.variants/v1`) with its notation.
3. Optionally, `dna notation <id> results/<id>.variants.json --reference
   <fasta>` adds the conformance findings.

`dna sample` remains the read-level evidence used to review a decision: per
read placement, callability, loci, support and opposition. Its notation lists
every read's calls and makes no decision. On two disjoint 320-sample held-out
subsets, the consensus path matched more samples exactly and was more precise
than that union ([ADR-0073](../decisions/adr/0073-sample-consensus.md)).

## Interpretation boundary

- The consensus is an adjudicated sequence. It is not a genotype, and it
  carries no heteroplasmy call or fraction.
- Read-level evidence stays in `dna sample`.
