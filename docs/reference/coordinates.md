# Coordinate Contract

DNA uses several coordinate domains. They are deliberately distinct.

| Domain | Base | Interval style | Example |
|---|---:|---|---|
| original call index | 0-based | internal implementation coordinate |
| ABIF PLOC / trace sample | 0-based | internal implementation coordinate |
| call/sample windows | 0-based | half-open | `[start, end)` |
| trim interval | 0-based | half-open | `[trim.start, trim.end)` |
| reference segments | 0-based | half-open | `[start, end)` |
| reported biological variant position | 1-based | scalar | `position = 73` |

## Rules

1. Coordinate domains are not interchangeable simply because they are represented by integers.
2. Reverse-strand alignment changes biological orientation but does not change the identity of the original trace call.
3. Inserted supporting calls have no biological reference position at the inserted base.
4. Deleted reference bases have no supporting trace-base call and must not be assigned fabricated signal evidence.
5. Variant normalization may move the reported anchor in a repeat while observed call mappings remain tied to the alignment evidence.
6. Circular projection uses the reference modulo length; wrapped alignments may be represented by two 0-based half-open reference segments.

See also [system invariants](../architecture/invariants/README.md) and the versioned output contracts.

Reviewer-facing `dna.analysis/v7` and `dna.sample_evidence/v8` variant evidence intentionally does not serialize original call index or PLOC. Those coordinates remain internal for mapping and tests; public variant evidence uses normalized biological position plus reference-oriented base/peaks/quality.
