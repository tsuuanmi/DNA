# Output and operations Invariants

These invariants are part of the canonical [system invariant set](README.md).

- **INV-OUT-001:** A failed core analysis publishes no scientific result.
- **INV-OUT-002:** Core result publication is atomic and does not overwrite an existing result.
- **INV-OUT-003:** Operational logs are separate from deterministic scientific result contracts.
- **INV-OUT-004:** A versioned public schema is not mutated retroactively; incompatible contract changes require a new schema version.
- **INV-OUT-005:** Validation measurement exports are separate ignored local artifacts. They reuse authoritative scientific evidence but cannot alter public result schemas, become production compatibility outputs, or apply research thresholds during export.
- **INV-OUT-006:** Validation event diagnostics expose existing call and signal provenance only. Source chromatogram sample coordinates remain unchanged, while all diagnostic A/C/G/T arrays are projected to the selected reference orientation. Diagnostics cannot reselect events or mutate scientific evidence.
