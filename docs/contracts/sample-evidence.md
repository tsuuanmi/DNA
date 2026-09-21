# DNA Sample Evidence JSON

`signal sample <sample-id> <trace.ab1>... --reference <reference.fasta>`
writes one deterministic `results/<sample-id>.sample.json` document identified as
`dna.sample_evidence/v8`. The authoritative schema is
[`schemas/sample-evidence-v8.schema.json`](schemas/sample-evidence-v8.schema.json)
and the example is
[`examples/sample-evidence-v8.example.json`](examples/sample-evidence-v8.example.json).

The sample identifier and read names are reviewer-facing provenance. They never
constrain scientific placement, orientation, overlap discovery, or variant
reconciliation.

## Contract sections

- [Reads and coverage](sample-evidence/reads-coverage.md)
- [Overlaps and differential loci](sample-evidence/overlaps-loci.md)
- [Variants and differential-locus signal evidence](sample-evidence/variants-signal.md)

## Why the evidence layers are separate

The four arrays intentionally preserve different evidence layers:

```text
selected per-read alignments
      ↓
coverage[]              local mapped-read denominator/orientation topology
      ↓
overlaps[]              which mapped read pairs are eligible for later reconciliation
      ↓
locus_differences[]     what each read observed at differential reference coordinates
      ↓
normalization/filtering
      ↓
variants[]              normalized alleles + eligibility + trace evidence
```

`locus_differences[]` is useful for disagreement and coverage reasoning,
including reference-vs-alternate or unresolved evidence. `variants[]` is the
reviewer-facing normalized biological call layer and is where per-read peak
evidence belongs.

None of these arrays is a consensus result.

## Contract boundary

v8 remains compact and difference-focused. It does not serialize per-base
evidence for loci where every covering read agrees with the reference. Pairwise
overlap records summarize only admission-relevant counts rather than dense
per-coordinate comparisons. The scientific pipeline still processes each read
independently before sample aggregation.

The current implementation emits v8 only. Earlier sample-evidence contracts are
not emitted as aliases or compatibility output.

## Non-goals

The v8 contract contains no consensus sequence, sample-level adjudicated variant
verdict, majority-vote result, genotype, heteroplasmy estimate, haplogroup
interpretation, F/R pair object, primer/HV placement rule, or filename-derived
placement. `overlaps[]` is an evidence/admission graph, not a pair-first merge
structure.
