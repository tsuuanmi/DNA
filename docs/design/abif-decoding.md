# ABIF Decoding Method

Part of the canonical [DNA pipeline](pipeline.md).

Parses the ABIF container and validates every directory entry, offset, element
size, and element count before access. It extracts:

- the four `DATA.9`–`DATA.12` channels as signed 16-bit samples, reordered into
  canonical A/C/G/T order using the `FWO_.1` channel-order string;
- the `PLOC.2` basecall positions (strictly increasing, within the sample
  range);
- optional vendor base strings and quality values. Their decoded cardinality may
  differ from PLOC and is retained as trace-integrity evidence rather than
  changing the PLOC-defined call series.

The decoded chromatogram records the source file name and SHA-256, the canonical
four channel arrays, the basecall positions, and optional vendor evidence. ABIF
version, channel order, and sample count are validated during decode but are not
duplicated as retained metadata.
