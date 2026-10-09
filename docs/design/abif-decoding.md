# ABIF Decoding Method

Part of the canonical [DNA pipeline](pipeline.md).

Within the source architecture, ABIF is a **format layer** in the Sanger modality crate (`dna_sanger::abif`), which the facade's Sanger input adapter (`input::sanger`) calls. Applied Biosystems Sanger sequencing sample files commonly use the `.ab1` extension, but decoding is determined by the ABIF container signature and required sequencing tags rather than the filename suffix.

## Implementation sourcing

ABIF is an external format, so ADR-0059 requires ecosystem reuse evaluation.
The current first-party decoder is nevertheless retained deliberately.

The container layout and sequencing tags are evaluated against the Applied
Biosystems Genetic Analysis Data File Format specification. An available Rust
candidate, `bio_files::ab1`, was also evaluated; its upstream source describes
the implementation as directly adapted from Biopython and currently carries an
unresolved `data_offset` handling TODO. That is not sufficient evidence to
replace a parser on DNA's untrusted-input boundary.

Any replacement must demonstrate, at minimum, equivalent behavior for bounded
offset/count arithmetic, inline versus external payloads, reserved/padded
allocations, duplicate tags, exact DATA/FWO_/PLOC/PBAS/PCON layout semantics,
malformed/truncated files, and representative production traces. Differential
agreement alone does not replace DNA's scientific and security validation.

Parses the ABIF container and validates every directory entry, offset, element
size, and element count before access. It extracts:

- the four `DATA.9`–`DATA.12` channels as signed 16-bit samples, reordered into
  canonical A/C/G/T order using the `FWO_.1` channel-order string;
- the `PLOC.2` basecall positions (strictly increasing, within the sample
  range), projected to canonical Sanger `locus_positions`;
- optional vendor base strings and quality values. Their decoded cardinality may
  differ from PLOC and is retained as trace-integrity evidence rather than
  changing the PLOC-defined call series.

The decoded `Chromatogram` is the canonical Sanger evidence boundary. It records the source file name and SHA-256, the canonical four channel arrays, canonical locus positions, and optional vendor evidence. ABIF directory entries and tag names—including `PLOC`—do not cross that boundary. ABIF
version, channel order, and sample count are validated during decode but are not
duplicated as retained metadata.
