# Read Placement Invariants

Part of the canonical [read/sample invariant family](sample-boundaries.md).

- **INV-READ-001:** Every trace is scientifically processed and placed independently before any cross-read reconciliation.
- **INV-READ-002:** Read orientation and covered reference segments are derived from alignment evidence; filename, amplicon/HV label, primer label, and declared F/R direction do not constrain default placement.
- **INV-READ-003:** Cross-read reconciliation uses normalized reference-coordinate/variant space rather than canonical F/R pairs as exclusive merge keys.
- **INV-READ-004:** A missing canonical partner does not invalidate an otherwise admitted read; optional assay metadata remains provenance or post-mapping QC unless a separately specified method explicitly says otherwise.
