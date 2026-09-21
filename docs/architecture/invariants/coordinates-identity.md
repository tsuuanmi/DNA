# Coordinates and identity Invariants

These invariants are part of the canonical [system invariant set](../invariants.md).

- **INV-COORD-001:** Trace sample positions, original call indexes, PLOC values, trim bounds, signal-window bounds, and reference segments are 0-based unless a contract explicitly says otherwise.
- **INV-COORD-002:** Reported biological variant positions are 1-based.
- **INV-COORD-003:** Half-open intervals use `[start, end)`.
- **INV-COORD-004:** Call index, trace-sample position, and biological reference position are different coordinate domains and must not be substituted implicitly.
- **INV-ID-001:** Reverse-strand processing preserves the original trace call identity and PLOC.
- **INV-ID-002:** Variant normalization may change reported allele placement without rewriting the observed trace-call mappings.
