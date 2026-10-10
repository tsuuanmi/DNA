# Known Limitations

These are current interpretation/validation boundaries, not future promises.

- A single chromatogram does not establish genotype, quantitative heteroplasmy, phase, contamination, pathogenicity, or clinical significance.
- Rolling SNR and relative quality are observational/relative measures, not calibrated Phred error probabilities.
- Read callability detects where a read's signal stops being one ladder and masks it; it does not recover the sequence behind a phase shift, and a slipped population above one half flips the primary inside a run with no in-run signature, so only cross-strand disagreement at sample scope can reveal it. Its `dephased` label means that shadows of neighbouring calls explain the double peaks; it does not separate slippage from a ladder smeared by one call, and dinucleotide slippage, whose shadows lie two calls away, is labelled `mixed`. A masked call supports no variant even when its primary call is right, and a variant needs `read_end_margin` informative calls between it and an uninformative end, so a difference right where a read dephases is reported only from a strand that reads it in phase. When several length populations of similar size meet at a poly-C run, no strand reads the run length, and length insertions such as 309.1C are not reported ([phase-recovery research](../research/phase-recovery/README.md)). A run-length edit whose run end a read does not resolve is excluded (`run_boundary`, ADR-0071), A masked repeat stretch whose calls spell the reference's runs is re-expressed as one length edit per run (ADR-0072), so a read that reads the run lengths before dephasing reports them; a read that loses the anchor base among several length populations still reports none.
- Mixed or secondary signal remains evidence; it is not automatically biological heteroplasmy.
- Unresolved evidence remains unresolved rather than being forced into a definitive biological claim.
- The sample consensus (`dna consensus`, ADR-0073) writes `N` where no read observes a position or run cleanly. In particular, a poly-C run that every read meets near its end (`read_end`) or beside an unresolved call stays `N`, so most remaining consensus misses are HVS-II length calls. One read against one read takes the reference and is marked `contested`. A point heteroplasmy is never written as an IUPAC code.
- External-tool agreement is differential evidence, not a compatibility requirement or ground-truth substitute.
- The standalone post-calling mtDNA normalizer currently rejects a source variant whose anchored REF allele itself spans the FASTA/rCRS end-to-start seam; it does not silently rotate or reinterpret that event.
- Production readiness requires approved real-trace validation in addition to engineering test success.

Canonical behavioral boundaries are defined by [requirements](../requirements/README.md), [design](../design/README.md), and [architecture invariants](../architecture/invariants/README.md).
