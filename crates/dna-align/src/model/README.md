# Alignment model

Owns `Alignment`: a read's placement, orientation, covered and callable
reference segments, and metrics, with a read-only view of its columns
(`AlignmentColumn`). Only `align_best` and `Alignment::placed` build one.

See [dna-align](../README.md).
