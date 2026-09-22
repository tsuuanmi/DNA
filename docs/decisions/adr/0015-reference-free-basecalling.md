# ADR-0015: Add One Reference-Free Basecall JSON Contract

- **Status:** Accepted
- **Date:** 2026-08-23

## Context

DNA's ABIF decode, signal-derived calling, rolling signal analysis, and
quality trimming do not require a reference, but they were reachable only through
`signal analyze`. Apollo exposed a broad basecall command with JSON, TSV, FASTA,
and FASTQ outputs. Copying that surface would create several publication
contracts and would misrepresent DNA's uncalibrated relative quality as FASTQ
quality.

## Decision

Add `signal basecall <trace.ab1>` as a second one-file command. It reuses one
shared Rust read-processing path and atomically creates exactly one
`results/<trace-stem>.basecalls.json` document identified as
`dna.basecalls/v1`.

The result contains input/configuration provenance, complete primary, ambiguity,
and retained sequences, trim bounds, merged observational noisy regions, and
unresolved-call warning counts and vendor-disagreement counts when optional
vendor calls are available. It omits reference, alignment, variants, vendor
evidence, per-call tables, full peaks, and rolling windows. It makes no Phred,
genotype, heteroplasmy, or clinical claim.

The complete strict configuration schema remains version 4. Basecall uses its
basecalling, signal-processing, and quality-control sections while retaining the
whole-file checksum as deterministic provenance. DNA adds no alternate config,
output override, stdout mode, format switch, compatibility alias, FASTA, TSV,
FASTQ, VCF, or BCF output.

## Consequences

Users can obtain auditable reference-free sequence results without duplicating
scientific algorithms. Analysis and basecall results coexist because their
derived suffixes differ, and each retains atomic no-overwrite publication.
Complete sequence strings make the basecall result more identifying than compact
analysis v5, so it follows the AB1 data policy.

A future interchange exporter or calibrated FASTQ contract requires a separate
decision and validation evidence; it must derive from typed results rather than
reimplement base calling.
