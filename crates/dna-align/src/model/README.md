# Alignment model

Owns `Alignment`: a read's placement, orientation, covered and callable
reference segments, and metrics, all read-only, with a read-only view of its
columns (`AlignmentColumn`). Only `align_best` and `Alignment::placed` build
one, and only its own methods narrow its callable segments:
`without_callable` leaves excluded positions out (SRS-VAR-019) and
`uncallable` drops them all for a placement in doubt (SRS-VAR-017).

See [dna-align](../README.md).
