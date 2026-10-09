# Known Limitations

These are current interpretation/validation boundaries, not future promises.

- A single chromatogram does not establish genotype, quantitative heteroplasmy, phase, contamination, pathogenicity, or clinical significance.
- Rolling SNR and relative quality are observational/relative measures, not calibrated Phred error probabilities.
- Read callability detects where a read's signal stops being one ladder and masks it; it does not recover the sequence behind a phase shift, and a slipped population above one half flips the primary inside a run with no in-run signature, so only cross-strand disagreement at sample scope can reveal it. Its `dephased` label means that shadows of neighbouring calls explain the double peaks; it does not separate slippage from a ladder smeared by one call, and dinucleotide slippage, whose shadows lie two calls away, is labelled `mixed`. A masked call supports no variant even when its primary call is right, and a variant needs `read_end_margin` informative calls between it and an uninformative end, so a difference right where a read dephases is reported only from a strand that reads it in phase. Until rejected reads exist, a read with too few callable calls fails the `sample` operation.
- Mixed or secondary signal remains evidence; it is not automatically biological heteroplasmy.
- Unresolved evidence remains unresolved rather than being forced into a definitive biological claim.
- External-tool agreement is differential evidence, not a compatibility requirement or ground-truth substitute.
- The standalone post-calling mtDNA normalizer currently rejects a source variant whose anchored REF allele itself spans the FASTA/rCRS end-to-start seam; it does not silently rotate or reinterpret that event.
- Production readiness requires approved real-trace validation in addition to engineering test success.

Canonical behavioral boundaries are defined by [requirements](../requirements/README.md), [design](../design/README.md), and [architecture invariants](../architecture/invariants/README.md).
