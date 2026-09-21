# Pipeline Invariants

These invariants are part of the canonical [system invariant set](../invariants.md).

- **INV-PIPE-001:** Scientific stages consume validated output from earlier stages and do not silently re-parse or reinterpret external inputs.
- **INV-PIPE-002:** Reference-aware stages do not alter upstream signal-derived base calls.
- **INV-PIPE-003:** Expected external-input failures are explicit typed failures, not panics or silent fallbacks.
- **INV-PIPE-004:** Identical scientific inputs, configuration, algorithm versions, and supported execution environment produce deterministic scientific output.
